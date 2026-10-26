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

const UNTRUSTED = 'lucid-coder.rules.json in this repo is not trusted; run /lucid-coder:rules trust after reviewing it.';

function sha256(buf) {
  return crypto.createHash('sha256').update(buf).digest('hex');
}

// Returns the active rules that apply here. A repo file whose hash is not
// trusted sets rules.notice instead of loading, so render can say so.
function readActiveRules() {
  let store = {};
  try { store = JSON.parse(fs.readFileSync(path.join(dataDir(), 'rules.json'), 'utf8')) || {}; } catch (e) {}
  let repo = null;
  try { repo = repoInfo(process.cwd()); } catch (e) {}
  const rules = (Array.isArray(store.rules) ? store.rules : []).filter(r =>
    r && r.status === 'active' &&
    !(r.scope && r.scope.kind === 'repo' && !(repo && r.scope.repo === repo.id)));
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

function render(core, rules) {
  let out = 'LUCID-CODER ACTIVE\n\n' + core;
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
  escapeRule, render, writeOutput, readStdin,
  dataDir, rememberDataDir, repoInfo, sha256, RULE_MAX,
};
