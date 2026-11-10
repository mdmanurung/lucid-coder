#!/usr/bin/env node
// The one writer of the learned-rules store (PLAN.md 5.2). Exit 0 ok,
// 1 usage or validation error, 2 conflict with an active rule.
const fs = require('fs');
const path = require('path');
const { parseArgs } = require('util');
const { dataDir, repoInfo, sha256, ruleProblem, STATUSES } = require('../hooks/lib');
const { harvest, markSeen } = require('../hooks/memory');

const USAGE = `Usage: node scripts/rules.js <command>
  list [--status proposed|active|rejected|superseded|retired] [--json]
  propose --topic T --scope global|language|repo|file-glob [--language r|python] [--glob G] [--repo ID]
          --statement S --why W [--good G] [--bad B]
          [--source user-correction|user-preference|repo-review|user-memory:<project>/<file>]
  approve <id>
  supersede <new-id> <old-id>
  reject <id>
  retire <id>
  trust --show | trust --sha256 <hash>
  harvest [--since-last] [--json] | harvest --mark-seen`;

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
    if (opts.scope === 'repo') scope.repo = opts.repo || currentRepo().id;
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

  // --show prints the file and its hash from one read; --sha256 trusts only
  // that content, so an edit between review and trust is refused.
  trust(store, args, opts) {
    const repo = currentRepo();
    const file = path.join(repo.root, REPO_FILE);
    if (!fs.existsSync(file)) die(`No ${REPO_FILE} at the root of this repo (${repo.root}).`);
    const text = fs.readFileSync(file);
    const hash = sha256(text);
    if (opts.show) return console.log(`${text.toString('utf8').trimEnd()}\nsha256 ${hash}`);
    if (opts.sha256 != null && opts.sha256.toLowerCase() !== hash) {
      die(`Not trusted: ${REPO_FILE} now has sha256 ${hash}, not the reviewed ${opts.sha256}. Review it again.`);
    }
    store.trusted = store.trusted.filter(t => t.repo !== repo.id);
    store.trusted.push({ repo: repo.id, sha256: hash, trusted: today() });
    save(store);
    console.log(`Trusted ${REPO_FILE} for ${repo.id} (sha256 ${hash}). Any edit to it needs trust again.`);
  },
};

// Reads Claude Code memory, never writes it; --mark-seen writes only hashes and a time.
function harvestCommand(opts) {
  if (opts['mark-seen']) return console.log(`Marked ${markSeen()} feedback memories as seen.`);
  const found = harvest({ sinceLast: opts['since-last'] });
  if (opts.json) return console.log(JSON.stringify(found, null, 2));
  if (!found.length) return console.log('No feedback memories to review.');
  for (const c of found) {
    console.log(`${c.source} [${c.projects_seen.length > 1 ? 'global' : c.repo || 'no repo'}] ${c.description}`);
  }
}

function main() {
  let parsed;
  try {
    parsed = parseArgs({
      allowPositionals: true,
      options: Object.fromEntries(
        ['status', 'topic', 'scope', 'language', 'glob', 'repo', 'statement', 'why', 'good', 'bad', 'source', 'sha256']
          .map(k => [k, { type: 'string' }]).concat(['json', 'since-last', 'mark-seen', 'show'].map(k => [k, { type: 'boolean' }]))),
    });
  } catch (e) {
    die(`${e.message}\n${USAGE}`);
  }
  const [command, ...args] = parsed.positionals;
  if (command === 'harvest') return harvestCommand(parsed.values);
  if (!commands[command]) die(USAGE);
  const needs = { approve: 1, reject: 1, retire: 1, supersede: 2 }[command] || 0;
  if (args.length !== needs) die(USAGE);
  commands[command](load(), args, parsed.values);
}

main();
