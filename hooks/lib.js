// Shared helpers for the lucid-coder hooks. Every reader fails open: a bad
// store or flag never removes the core ruleset.
const crypto = require('crypto');
const fs = require('fs');
const os = require('os');
const path = require('path');

const RULES_BUDGET = 4000; // learned-rules share of the injected text (PLAN.md 5.5)
const OUTPUT_CAP = 9500; // additionalContext is capped at 10,000 chars
const RULE_MAX = 300;
const SCOPE_ORDER = ['file-glob', 'repo', 'language', 'global']; // narrowest first

const FALLBACK_CORE = [
  '# lucid-coder',
  'Write scientific R/Python for the next reader: one-pass understanding and a clean rerun.',
  'Add a check only if it stops silent data loss, a silent mislabel or a silent overwrite.',
].join('\n');

function stateDir() {
  return process.env.CLAUDE_CONFIG_DIR || path.join(os.homedir(), '.claude');
}

function flagPath() {
  return path.join(stateDir(), '.lucid-coder-active');
}

function pointerPath() {
  return path.join(stateDir(), '.lucid-coder-data');
}

// Skills' Bash calls may lack CLAUDE_PLUGIN_DATA, so the hooks leave a pointer
// to the dir the plugin really uses.
function dataDir() {
  if (process.env.CLAUDE_PLUGIN_DATA) return process.env.CLAUDE_PLUGIN_DATA;
  try {
    const pointed = fs.readFileSync(pointerPath(), 'utf8').trim();
    if (pointed) return pointed;
  } catch (e) {}
  return path.join(stateDir(), 'lucid-coder');
}

function rememberDataDir() {
  if (!process.env.CLAUDE_PLUGIN_DATA) return;
  try {
    fs.mkdirSync(stateDir(), { recursive: true });
    fs.writeFileSync(pointerPath(), process.env.CLAUDE_PLUGIN_DATA);
  } catch (e) {}
}

// 'git@host:owner/repo.git' and 'https://host/owner/repo' both become 'host/owner/repo'.
function normaliseRemote(url) {
  const u = url.trim().replace(/\/+$/, '').replace(/\.git$/, '');
  const full = u.match(/^[a-z+]+:\/\/(?:[^@/]+@)?([^/:]+)(?::\d+)?\/(.+)$/i);
  if (full) return full[1].toLowerCase() + '/' + full[2];
  const scp = u.match(/^(?:[^@/]+@)?([^/:]+):(.+)$/);
  if (scp) return scp[1].toLowerCase() + '/' + scp[2].replace(/^\/+/, '');
  return u;
}

// The git dir behind <dir>/.git, or null. Like git, a .git without HEAD
// (e.g. a stray empty /tmp/.git) is not a repo.
function gitDirOf(dir) {
  const dotGit = path.join(dir, '.git');
  let gitdir = dotGit;
  try {
    if (fs.statSync(dotGit).isFile()) {
      const m = fs.readFileSync(dotGit, 'utf8').match(/^gitdir:\s*(.+)$/m);
      gitdir = m ? path.resolve(dir, m[1].trim()) : '';
    }
  } catch (e) {
    return null;
  }
  return gitdir && fs.existsSync(path.join(gitdir, 'HEAD')) ? gitdir : null;
}

// Reads .git directly (no git spawn): id is the origin remote, else the root path.
function repoInfo(cwd) {
  for (let dir = path.resolve(cwd); ; dir = path.dirname(dir)) {
    const gitdir = gitDirOf(dir);
    if (gitdir) {
      // Worktrees keep config in the common dir.
      let configDir = gitdir;
      try { configDir = path.resolve(gitdir, fs.readFileSync(path.join(gitdir, 'commondir'), 'utf8').trim()); } catch (e) {}
      let id = dir;
      try {
        const config = fs.readFileSync(path.join(configDir, 'config'), 'utf8');
        const section = config.match(/^\s*\[remote\s+"origin"\]([^[]*)/m);
        const url = section && section[1].match(/^\s*url\s*=\s*(.+)$/m);
        if (url) id = normaliseRemote(url[1]);
      } catch (e) {}
      return { root: dir, id };
    }
    if (path.dirname(dir) === dir) return null;
  }
}

// Absent flag means on: lucid-coder is always-on until the user turns it off.
function readMode() {
  try {
    return fs.readFileSync(flagPath(), 'utf8').trim() === 'off' ? 'off' : 'on';
  } catch (e) {
    return 'on';
  }
}

function setMode(mode) {
  if (mode === 'off') {
    fs.mkdirSync(stateDir(), { recursive: true });
    fs.writeFileSync(flagPath(), 'off');
    return;
  }
  try { fs.unlinkSync(flagPath()); } catch (e) {}
}

function readCore() {
  try {
    return fs.readFileSync(path.join(__dirname, '..', 'rules', 'core.md'), 'utf8').trim();
  } catch (e) {
    return FALLBACK_CORE;
  }
}

// Language tiers: rules/<lang>.md, injected after core for the languages
// the repo uses. A missing tier file is just skipped.
function readTiers(langs) {
  const tiers = [];
  for (const lang of langs) {
    try { tiers.push(fs.readFileSync(path.join(__dirname, '..', 'rules', `${lang}.md`), 'utf8').trim()); } catch (e) {}
  }
  return tiers;
}

const SKIP_DIRS = new Set(['.git', 'node_modules', 'renv', '.venv', '__pycache__', 'results']);
const LANG_EXT = { '.r': 'r', '.rmd': 'r', '.qmd': 'r', '.py': 'python', '.ipynb': 'python' };
const DETECT_DEPTH = 3;

// Which of r/python the tree under cwd uses, as a Set. Sync scan to depth 3,
// stops once both are found; any fs error leaves what was found so far.
function detectLanguages(cwd) {
  const found = new Set();
  const walk = (dir, depth) => {
    if (found.size === 2 || depth > DETECT_DEPTH) return;
    let entries = [];
    try { entries = fs.readdirSync(dir, { withFileTypes: true }); } catch (e) { return; }
    for (const entry of entries) {
      if (found.size === 2) return;
      if (entry.isDirectory()) {
        if (!SKIP_DIRS.has(entry.name)) walk(path.join(dir, entry.name), depth + 1);
        continue;
      }
      const lang = LANG_EXT[path.extname(entry.name).toLowerCase()];
      if (lang) found.add(lang);
    }
  };
  walk(cwd, 0);
  return found;
}

const LANG_VALUES = { r: ['r'], python: ['python'], both: ['r', 'python'] };

function langPath() {
  return path.join(stateDir(), '.lucid-coder-lang');
}

// The override is per repo: the flag file is a JSON map from repo id (or the
// real path outside a repo) to r, python or both. Anything unreadable, including
// an old plain-text flag, means auto-detect.
function langKey(cwd) {
  let dir = path.resolve(cwd);
  try { dir = fs.realpathSync(dir); } catch (e) {}
  const repo = repoInfo(dir);
  return repo ? repo.id : dir;
}

function readLangMap() {
  try {
    const map = JSON.parse(fs.readFileSync(langPath(), 'utf8'));
    return map && typeof map === 'object' && !Array.isArray(map) ? map : {};
  } catch (e) {
    return {};
  }
}

function readLang(cwd = process.cwd()) {
  try {
    const value = readLangMap()[langKey(cwd)];
    return LANG_VALUES[value] ? value : 'auto';
  } catch (e) {
    return 'auto';
  }
}

function setLang(value, cwd = process.cwd()) {
  const map = readLangMap();
  const key = langKey(cwd);
  if (LANG_VALUES[value]) map[key] = value;
  else delete map[key];
  fs.mkdirSync(stateDir(), { recursive: true });
  fs.writeFileSync(langPath(), JSON.stringify(map));
}

// Tier languages for this session: the override flag wins over detection.
function languages(cwd) {
  const lang = readLang(cwd);
  if (lang !== 'auto') return LANG_VALUES[lang];
  return ['r', 'python'].filter(l => detectLanguages(cwd).has(l));
}

const UNTRUSTED = 'lucid-coder.rules.json in this repo is not trusted; run /lucid-coder:rules trust after reviewing it.';

function sha256(buf) {
  return crypto.createHash('sha256').update(buf).digest('hex');
}

// Returns the active rules that apply here. A repo file whose hash is not
// trusted sets rules.notice instead of loading, so render can say so.
// langs, when given, drops language-scoped rules for languages not in use.
function readActiveRules(langs) {
  let store = {};
  try { store = JSON.parse(fs.readFileSync(path.join(dataDir(), 'rules.json'), 'utf8')) || {}; } catch (e) {}
  let repo = null;
  try { repo = repoInfo(process.cwd()); } catch (e) {}
  const applies = s => !s ||
    !(s.kind === 'repo' && !(repo && s.repo === repo.id)) &&
    !(s.kind === 'language' && langs && !langs.includes(s.language));
  const rules = (Array.isArray(store.rules) ? store.rules : []).filter(r =>
    r && r.status === 'active' && applies(r.scope));
  if (!repo) return rules;
  try {
    const text = fs.readFileSync(path.join(repo.root, 'lucid-coder.rules.json'));
    const hash = sha256(text);
    const trusted = (Array.isArray(store.trusted) ? store.trusted : [])
      .some(t => t && t.repo === repo.id && t.sha256 === hash);
    if (!trusted) {
      rules.notice = UNTRUSTED;
      return rules;
    }
    for (const r of JSON.parse(text).rules || []) {
      if (r && r.status === 'active') rules.push({ ...r, scope: { kind: 'repo', repo: repo.id } });
    }
  } catch (e) {
    // no repo file, or a broken one: skip it
  }
  return rules;
}

// One line, no leading '#', so rule text cannot open a fake heading.
function escapeRule(text) {
  return String(text || '').replace(/\s+/g, ' ').trim().replace(/^#+\s*/, '').slice(0, RULE_MAX);
}

function scopeRank(rule) {
  const i = SCOPE_ORDER.indexOf(rule.scope && rule.scope.kind);
  return i === -1 ? SCOPE_ORDER.length : i;
}

function scopeLabel(scope) {
  const s = scope || {};
  if (s.kind === 'language') return `language: ${s.language}`;
  if (s.kind === 'file-glob') return `file-glob: ${s.glob}`;
  return s.kind;
}

function render(core, rules, tiers = []) {
  let out = 'LUCID-CODER ACTIVE\n\n' + core;
  for (const tier of tiers) out += '\n\n' + tier;
  if (rules.length) out += learnedBlock(out.length, rules);
  if (rules.notice) out += '\n\n' + rules.notice;
  return out;
}

function learnedBlock(used, rules) {
  const lines = rules
    .slice()
    .sort((a, b) => scopeRank(a) - scopeRank(b))
    .map(r => `- [${escapeRule(r.topic)}, ${escapeRule(scopeLabel(r.scope))}] ${escapeRule(r.statement)}`);
  // The 200-char margin also covers the header, the hidden-count line and the trust notice.
  const budget = Math.min(RULES_BUDGET, OUTPUT_CAP - used - 200);
  let block = '';
  let shown = 0;
  for (const line of lines) {
    if (block.length + line.length + 1 > budget) break;
    block += line + '\n';
    shown++;
  }
  let out = '\n\n## Learned rules\n' + block;
  if (shown < lines.length) {
    out += `${lines.length - shown} rules not shown; run /lucid-coder:rules.\n`;
  }
  return out;
}

// Native Claude Code: SubagentStart drops raw stdout, so it needs the JSON form.
function writeOutput(event, text) {
  try {
    if (event === 'SubagentStart') {
      process.stdout.write(JSON.stringify(
        { hookSpecificOutput: { hookEventName: event, additionalContext: text } }));
      return;
    }
    process.stdout.write(text);
  } catch (e) {
    // stdout closed at hook exit must not surface as a hook failure
  }
}

// Calls cb once with the parsed payload, or {} on bad input, stdin error or timeout.
function readStdin(cb) {
  let input = '';
  let done = false;
  const finish = () => {
    if (done) return;
    done = true;
    let payload = {};
    try { payload = JSON.parse(input.replace(/^﻿/, '')) || {}; } catch (e) {}
    cb(payload);
  };
  process.stdin.on('data', chunk => { input += chunk; });
  process.stdin.on('end', finish);
  process.stdin.on('error', finish);
  setTimeout(finish, 1000).unref();
}

module.exports = {
  FALLBACK_CORE, readMode, setMode, readCore, readActiveRules,
  readTiers, detectLanguages, readLang, setLang, languages,
  escapeRule, render, writeOutput, readStdin,
  dataDir, rememberDataDir, repoInfo, sha256, RULE_MAX,
};
