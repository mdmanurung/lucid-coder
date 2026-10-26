#!/usr/bin/env node
// The one writer of the learned-rules store (PLAN.md 5.2). Exit 0 ok,
// 1 usage or validation error, 2 conflict with an active rule.
const fs = require('fs');
const path = require('path');
const { parseArgs } = require('util');
const { dataDir, repoInfo, sha256, RULE_MAX } = require('../hooks/lib');

const USAGE = `Usage: node scripts/rules.js <command>
  list [--status proposed|active|rejected|superseded|retired] [--json]
  propose --topic T --scope global|language|repo|file-glob [--language r|python] [--glob G]
          --statement S --why W [--good G] [--bad B] [--source user-correction|user-preference|repo-review]
  approve <id>
  supersede <new-id> <old-id>
  reject <id>
  retire <id>
  trust`;

const STATUSES = ['proposed', 'active', 'rejected', 'superseded', 'retired'];
const SOURCES = ['user-correction', 'user-preference', 'repo-review'];
const TOPIC = /^[a-z]+(\.[a-z0-9_]+)+$/;
const SECRET = /Bearer\s|AKIA[0-9A-Z]{16}|ghp_|sk-[A-Za-z0-9]{20,}|-----BEGIN|(secret|password|token)\s*=/i;
const ABS_PATH = /\/home\/|\/Users\/|\b[A-Za-z]:\\/;
const REPO_FILE = 'lucid-coder.rules.json';

function die(msg, code = 1) {
  console.error(msg);
  process.exit(code);
}

const today = () => new Date().toISOString().slice(0, 10);
const storePath = () => path.join(dataDir(), 'rules.json');
const scopeKey = s => [s.kind, s.language || s.repo || s.glob || ''].join(':');

function scopeText(s) {
  return s.kind === 'global' ? 'global' : `${s.kind}: ${s.language || s.repo || s.glob}`;
}

// Returns the first problem with a rule, or null.
function ruleProblem(r) {
  if (!r || typeof r !== 'object') return 'a rule is not an object';
  if (!/^r-\d{4,}$/.test(r.id || '')) return `bad id ${r.id}`;
  if (!TOPIC.test(r.topic || '')) return `${r.id}: topic "${r.topic}" must look like r.pipe or user.my_topic`;
  const s = r.scope || {};
  if (s.kind === 'language' && !['r', 'python'].includes(s.language)) return `${r.id}: language scope needs --language r or python`;
  if (s.kind === 'file-glob' && !s.glob) return `${r.id}: file-glob scope needs --glob`;
  if (s.kind === 'repo' && !s.repo) return `${r.id}: repo scope needs a repo id`;
  if (!['global', 'language', 'repo', 'file-glob'].includes(s.kind)) return `${r.id}: scope must be global, language, repo or file-glob`;
  if (!r.statement) return `${r.id}: --statement is required`;
  if (!r.why) return `${r.id}: --why is required`;
  for (const field of ['statement', 'why', 'good', 'bad']) {
    const v = r[field];
    if (v == null) continue;
    if (typeof v !== 'string') return `${r.id}: ${field} must be text`;
    if (/[\r\n]/.test(v)) return `${r.id}: ${field} must be one line; remove the newline`;
    if (field !== 'why' && v.length > RULE_MAX) return `${r.id}: ${field} is ${v.length} chars; keep it to ${RULE_MAX}`;
    if (SECRET.test(v)) return `${r.id}: ${field} looks like it holds a secret; remove it`;
    if (field !== 'why' && ABS_PATH.test(v)) return `${r.id}: ${field} holds an absolute path; use a repo-relative path`;
  }
  if (/^\s*#/.test(r.statement) || /^\s*#/.test(r.why)) return `${r.id}: statement and why cannot start with #`;
  if (!STATUSES.includes(r.status)) return `${r.id}: unknown status ${r.status}`;
  if (!SOURCES.includes(r.source)) return `${r.id}: source must be one of ${SOURCES.join(', ')}`;
  return null;
}

function load() {
  let text;
  try {
    text = fs.readFileSync(storePath(), 'utf8');
  } catch (e) {
    if (e.code === 'ENOENT') return { version: 1, rules: [], trusted: [] };
    die(`Cannot read ${storePath()}: ${e.message}`);
  }
  let store;
  try { store = JSON.parse(text); } catch (e) {}
  if (!store || !Array.isArray(store.rules) || !Array.isArray(store.trusted || [])) {
    die(`${storePath()} is not a valid store; it was left unchanged. Fix it by hand or restore rules.json.bak.`);
  }
  store.trusted = store.trusted || [];
  return store;
}

function save(store) {
  for (const r of store.rules) {
    const problem = ruleProblem(r);
    if (problem) die(`Not saved: ${problem}`);
  }
  const seen = new Set();
  for (const r of store.rules.filter(x => x.status === 'active')) {
    const key = r.topic + '|' + scopeKey(r.scope);
    if (seen.has(key)) die(`Not saved: two active rules for ${r.topic} (${scopeText(r.scope)})`);
    seen.add(key);
  }
  const file = storePath();
  fs.mkdirSync(path.dirname(file), { recursive: true });
  if (fs.existsSync(file)) fs.copyFileSync(file, file + '.bak');
  const tmp = `${file}.${process.pid}.tmp`;
  fs.writeFileSync(tmp, JSON.stringify(store, null, 2) + '\n');
  fs.renameSync(tmp, file);
}

function find(store, id, status) {
  const r = store.rules.find(x => x.id === id);
  if (!r) die(`No rule ${id}; run: node scripts/rules.js list`);
  if (status && r.status !== status) die(`${id} is ${r.status}, not ${status}`);
  return r;
}

function show(r) {
  return `${r.id} [${r.status}] ${r.topic} (${scopeText(r.scope)}): ${r.statement}`;
}

function currentRepo() {
  const repo = repoInfo(process.cwd());
  if (!repo) die('Not inside a git repository; run this from the repo the rule is for.');
  return repo;
}

const commands = {
  list(store, args, opts) {
    if (opts.status && !STATUSES.includes(opts.status)) die(`--status must be one of ${STATUSES.join(', ')}`);
    const rules = store.rules.filter(r => !opts.status || r.status === opts.status);
    if (opts.json) return console.log(JSON.stringify(rules, null, 2));
    console.log(rules.length ? rules.map(show).join('\n') : 'No rules.');
  },

  propose(store, args, opts) {
    const scope = { kind: opts.scope };
    if (opts.scope === 'language') scope.language = opts.language;
    if (opts.scope === 'file-glob') scope.glob = opts.glob;
    if (opts.scope === 'repo') scope.repo = currentRepo().id;
    const last = Math.max(0, ...store.rules.map(r => parseInt(String(r.id).slice(2), 10) || 0));
    const rule = {
      id: 'r-' + String(last + 1).padStart(4, '0'),
      topic: opts.topic,
      scope,
      statement: opts.statement,
      why: opts.why,
      good: opts.good || null,
      bad: opts.bad || null,
      status: 'proposed',
      supersedes: null,
      source: opts.source || 'user-correction',
      created: today(),
      approved: null,
    };
    const problem = ruleProblem(rule);
    if (problem) die(`Not proposed: ${problem.replace(rule.id + ': ', '')}`);
    store.rules.push(rule);
    save(store);
    console.log(JSON.stringify(rule, null, 2));
  },

  approve(store, [id]) {
    const rule = find(store, id, 'proposed');
    const clash = store.rules.find(r => r.status === 'active' && r.topic === rule.topic &&
      scopeKey(r.scope) === scopeKey(rule.scope));
    if (clash) {
      die(`${clash.id} is already active for ${rule.topic} (${scopeText(rule.scope)}). ` +
        `To replace it run: node scripts/rules.js supersede ${id} ${clash.id}`, 2);
    }
    Object.assign(rule, { status: 'active', approved: today() });
    save(store);
    console.log('Approved ' + show(rule));
  },

  supersede(store, [newId, oldId]) {
    const next = find(store, newId, 'proposed');
    const old = find(store, oldId, 'active');
    if (next.topic !== old.topic || scopeKey(next.scope) !== scopeKey(old.scope)) {
      die(`${newId} and ${oldId} differ in topic or scope; supersede needs both to match`);
    }
    Object.assign(next, { status: 'active', approved: today(), supersedes: oldId });
    old.status = 'superseded';
    save(store);
    console.log(`Approved ${show(next)}\nSuperseded ${oldId}`);
  },

  reject(store, [id]) {
    find(store, id, 'proposed').status = 'rejected';
    save(store);
    console.log('Rejected ' + id);
  },

  retire(store, [id]) {
    find(store, id, 'active').status = 'retired';
    save(store);
    console.log('Retired ' + id);
  },

  trust(store) {
    const repo = currentRepo();
    const file = path.join(repo.root, REPO_FILE);
    if (!fs.existsSync(file)) die(`No ${REPO_FILE} at the root of this repo (${repo.root}).`);
    const hash = sha256(fs.readFileSync(file));
    store.trusted = store.trusted.filter(t => t.repo !== repo.id);
    store.trusted.push({ repo: repo.id, sha256: hash, trusted: today() });
    save(store);
    console.log(`Trusted ${REPO_FILE} for ${repo.id} (sha256 ${hash}). Any edit to it needs trust again.`);
  },
};

function main() {
  let parsed;
  try {
    parsed = parseArgs({
      allowPositionals: true,
      options: Object.fromEntries(
        ['status', 'topic', 'scope', 'language', 'glob', 'statement', 'why', 'good', 'bad', 'source']
          .map(k => [k, { type: 'string' }]).concat([['json', { type: 'boolean' }]])),
    });
  } catch (e) {
    die(`${e.message}\n${USAGE}`);
  }
  const [command, ...args] = parsed.positionals;
  if (!commands[command]) die(USAGE);
  const needs = { approve: 1, reject: 1, retire: 1, supersede: 2 }[command] || 0;
  if (args.length !== needs) die(USAGE);
  commands[command](load(), args, parsed.values);
}

main();
