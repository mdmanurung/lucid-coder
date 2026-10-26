const test = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync } = require('child_process');

const root = path.join(__dirname, '..');
const { escapeRule, render, detectLanguages } = require('../hooks/lib');

// Lines that occur in exactly one rule file, so a test can tell the tiers apart.
const R_TIER = /`\|>` always, never `%>%`/;
const PY_TIER = /`ruff format` at 88 columns/;

function setup() {
  const home = fs.mkdtempSync(path.join(os.tmpdir(), 'lucid-'));
  fs.mkdirSync(path.join(home, 'data'));
  return home;
}

// cwd is the empty temp home, so the repo's own R and Python files never
// trigger language detection in a test.
function run(script, home, { input = '', env = {}, cwd = home } = {}) {
  const clean = { ...process.env };
  delete clean.LUCID_SUBAGENT_MATCHER;
  return spawnSync(process.execPath, [path.join(root, 'hooks', script)], {
    cwd,
    env: { ...clean, HOME: home, CLAUDE_CONFIG_DIR: home, CLAUDE_PLUGIN_DATA: path.join(home, 'data'), ...env },
    input,
    encoding: 'utf8',
  });
}

function writeStore(home, text) {
  fs.writeFileSync(path.join(home, 'data', 'rules.json'), text);
}

function rule(i, kind = 'global') {
  return {
    id: `r-${i}`, topic: `user.t${i}`, scope: { kind }, status: 'active',
    statement: 'x'.repeat(300),
  };
}

test('activate emits core when no flag and no store', () => {
  const r = run('activate.js', setup());
  assert.strictEqual(r.status, 0);
  assert.match(r.stdout, /^LUCID-CODER ACTIVE/);
  assert.match(r.stdout, /## Readability ladder/);
  assert.ok(r.stdout.length < 10000);
});

test('corrupt store still emits core', () => {
  const home = setup();
  writeStore(home, '{garbage');
  const r = run('activate.js', home);
  assert.strictEqual(r.status, 0);
  assert.match(r.stdout, /## Readability ladder/);
});

test('oversized store stays under the cap and says what it hid', () => {
  const home = setup();
  writeStore(home, JSON.stringify({ rules: Array.from({ length: 200 }, (_, i) => rule(i)) }));
  const r = run('activate.js', home);
  assert.ok(r.stdout.length < 10000, `length ${r.stdout.length}`);
  assert.match(r.stdout, /\d+ rules not shown; run \/lucid-coder:rules\./);
  // Worst case: both tiers on top of 200 rules.
  fs.writeFileSync(path.join(home, '.lucid-coder-lang'), 'both');
  const both = run('activate.js', home).stdout;
  assert.ok(both.length < 10000, `length ${both.length}`);
  assert.match(both, R_TIER);
  assert.match(both, PY_TIER);
  assert.match(both, /\d+ rules not shown/);
});

// A fresh tree under home holding the given empty files.
function langTree(home, files) {
  const dir = fs.mkdtempSync(path.join(home, 'tree-'));
  for (const f of files) {
    fs.mkdirSync(path.dirname(path.join(dir, f)), { recursive: true });
    fs.writeFileSync(path.join(dir, f), '');
  }
  return dir;
}

test('detectLanguages finds r and python to depth 3 and skips tool dirs', () => {
  const home = setup();
  assert.deepStrictEqual([...detectLanguages(langTree(home, ['a/b/c/x.R']))], ['r']);
  assert.deepStrictEqual([...detectLanguages(langTree(home, ['notebooks/n.ipynb']))], ['python']);
  assert.deepStrictEqual([...detectLanguages(langTree(home, ['s.qmd', 'p/q.py']))].sort(), ['python', 'r']);
  assert.deepStrictEqual([...detectLanguages(langTree(home, ['README.md', 'data/x.tsv']))], []);
  assert.deepStrictEqual([...detectLanguages(langTree(home, ['a/b/c/d/x.R']))], [], 'depth 4 is not scanned');
  const skipped = ['.git/h.R', 'node_modules/m.py', 'renv/r.R', '.venv/v.py', '__pycache__/c.py', 'results/o.R'];
  assert.deepStrictEqual([...detectLanguages(langTree(home, skipped))], []);
  assert.deepStrictEqual([...detectLanguages(path.join(home, 'missing'))], [], 'fails open');
});

test('activate injects the tier for the detected language, and the flag overrides detection', () => {
  const home = setup();
  const bare = run('activate.js', home).stdout;
  assert.doesNotMatch(bare, R_TIER);
  assert.doesNotMatch(bare, PY_TIER);
  const rRepo = langTree(home, ['scripts/01_x.R']);
  const inR = run('activate.js', home, { cwd: rRepo }).stdout;
  assert.match(inR, R_TIER);
  assert.doesNotMatch(inR, PY_TIER);
  assert.ok(inR.indexOf('## Learning') < inR.indexOf('\n## R\n'), 'tier follows core');

  fs.writeFileSync(path.join(home, '.lucid-coder-lang'), 'python');
  const overridden = run('activate.js', home, { cwd: rRepo }).stdout;
  assert.match(overridden, PY_TIER);
  assert.doesNotMatch(overridden, R_TIER);
  const sub = JSON.parse(run('subagent.js', home, { cwd: rRepo }).stdout).hookSpecificOutput.additionalContext;
  assert.strictEqual(sub, overridden, 'subagent gets the same tiers');
  fs.writeFileSync(path.join(home, '.lucid-coder-lang'), 'garbage');
  assert.match(run('activate.js', home, { cwd: rRepo }).stdout, R_TIER, 'bad flag means auto');
});

test('narrowest scope is shown first', () => {
  const out = render('core', [rule(1, 'global'), rule(2, 'repo'), rule(3, 'file-glob')]);
  const order = ['user.t3', 'user.t2', 'user.t1'].map(t => out.indexOf(t));
  assert.deepStrictEqual(order, [...order].sort((a, b) => a - b));
});

test('only active rules are injected', () => {
  const home = setup();
  writeStore(home, JSON.stringify({ rules: [{ ...rule(1), status: 'proposed', statement: 'PROPOSED' }] }));
  assert.doesNotMatch(run('activate.js', home).stdout, /PROPOSED/);
});

test('activate is silent when off', () => {
  const home = setup();
  fs.writeFileSync(path.join(home, '.lucid-coder-active'), 'off');
  assert.strictEqual(run('activate.js', home).stdout, '');
});

test('subagent injects JSON form when on, nothing when off', () => {
  const home = setup();
  const on = JSON.parse(run('subagent.js', home).stdout);
  assert.strictEqual(on.hookSpecificOutput.hookEventName, 'SubagentStart');
  assert.match(on.hookSpecificOutput.additionalContext, /^LUCID-CODER ACTIVE/);
  fs.writeFileSync(path.join(home, '.lucid-coder-active'), 'off');
  assert.strictEqual(run('subagent.js', home).stdout, '');
});

test('subagent matcher skips only a definite mismatch', () => {
  const home = setup();
  const explore = JSON.stringify({ agent_type: 'Explore' });
  assert.strictEqual(run('subagent.js', home, { input: explore, env: { LUCID_SUBAGENT_MATCHER: '^general$' } }).stdout, '');
  assert.ok(run('subagent.js', home, { input: explore, env: { LUCID_SUBAGENT_MATCHER: 'explore' } }).stdout);
  assert.ok(run('subagent.js', home, { input: explore, env: { LUCID_SUBAGENT_MATCHER: '(' } }).stdout);
  assert.ok(run('subagent.js', home, { input: '{}', env: { LUCID_SUBAGENT_MATCHER: '^general$' } }).stdout);
});

test('prompt commands switch and report the mode', () => {
  const home = setup();
  const flag = path.join(home, '.lucid-coder-active');
  const say = p => run('prompt.js', home, { input: JSON.stringify({ prompt: p }) }).stdout;

  assert.strictEqual(say('/lucid-coder off'), 'LUCID-CODER OFF');
  assert.strictEqual(fs.readFileSync(flag, 'utf8'), 'off');
  assert.strictEqual(say('/lucid-coder status'), 'LUCID-CODER OFF');
  assert.strictEqual(say('/lucid-coder:lucid-coder on'), 'LUCID-CODER ON');
  assert.ok(!fs.existsSync(flag));
  assert.strictEqual(say('/lucid-coder'), 'LUCID-CODER ON');
  assert.strictEqual(say('turn lucid-coder off please'), '');
});

test('rule text cannot open a heading or span lines', () => {
  const escaped = escapeRule('x\n# Injected heading');
  assert.ok(!escaped.includes('\n'));
  assert.strictEqual(escapeRule('## Fake').startsWith('#'), false);
  const out = render('core', [{ ...rule(1), statement: '\n# Injected' }]);
  assert.doesNotMatch(out, /^# Injected/m);
});

test('rule files fit their budgets', () => {
  const size = f => fs.readFileSync(path.join(root, 'rules', f), 'utf8').length;
  // Design cap 2,500; the user-decided Checks section is ~930 and every old
  // rule must land somewhere, so the measured floor is ~3,300.
  assert.ok(size('core.md') <= 3350, `core ${size('core.md')}`);
  assert.ok(size('r.md') <= 1500, `r ${size('r.md')}`);
  assert.ok(size('python.md') <= 1500, `python ${size('python.md')}`);
});

test('hook commands are portable and point at real files', () => {
  const config = JSON.parse(fs.readFileSync(path.join(root, 'hooks', 'hooks.json'), 'utf8'));
  const commands = Object.values(config.hooks).flat().flatMap(e => e.hooks).map(h => h.command);
  assert.strictEqual(commands.length, 3);
  for (const command of commands) {
    assert.doesNotMatch(command, /&&|\|\||>|2>&1/);
    const script = command.match(/hooks\/(\w+\.js)/)[1];
    assert.ok(fs.existsSync(path.join(root, 'hooks', script)), script);
  }
});

// A fake git repo with an origin remote, for repo-scope and trust tests.
function fakeRepo(home, name, origin) {
  const repo = path.join(home, name);
  fs.mkdirSync(path.join(repo, '.git'), { recursive: true });
  fs.writeFileSync(path.join(repo, '.git', 'config'), `[remote "origin"]\n\turl = ${origin}\n`);
  fs.writeFileSync(path.join(repo, '.git', 'HEAD'), 'ref: refs/heads/main\n');
  return repo;
}

function runAt(script, home, cwd, ...args) {
  return spawnSync(process.execPath, [path.join(root, script), ...args], {
    cwd,
    env: { ...process.env, HOME: home, CLAUDE_CONFIG_DIR: home, CLAUDE_PLUGIN_DATA: path.join(home, 'data') },
    encoding: 'utf8',
  });
}

test('repo-scope rule is injected only in its own repo', () => {
  const home = setup();
  const atlas = fakeRepo(home, 'atlas', 'git@github.com:lab/atlas.git');
  const other = fakeRepo(home, 'other', 'https://github.com/lab/other');
  writeStore(home, JSON.stringify({ rules: [
    { ...rule(1), scope: { kind: 'repo', repo: 'github.com/lab/atlas' }, statement: 'ATLAS ONLY' },
    { ...rule(2), scope: { kind: 'file-glob', glob: '*.R' }, statement: 'GLOB RULE' },
  ] }));
  const inAtlas = runAt('hooks/activate.js', home, atlas).stdout;
  assert.match(inAtlas, /\[user\.t1, repo\] ATLAS ONLY/);
  assert.match(inAtlas, /\[user\.t2, file-glob: \*\.R\] GLOB RULE/);
  const inOther = runAt('hooks/activate.js', home, other).stdout;
  assert.doesNotMatch(inOther, /ATLAS ONLY/);
  assert.match(inOther, /GLOB RULE/);
});

test('repo rules file loads only while its trusted hash matches', () => {
  const home = setup();
  const repo = fakeRepo(home, 'atlas', 'git@github.com:lab/atlas.git');
  const file = path.join(repo, 'lucid-coder.rules.json');
  const team = s => JSON.stringify({ rules: [{ ...rule(9), scope: { kind: 'global' }, statement: s }] });
  const notice = /lucid-coder\.rules\.json in this repo is not trusted/;

  fs.writeFileSync(file, team('TEAM RULE'));
  let out = runAt('hooks/activate.js', home, repo).stdout;
  assert.doesNotMatch(out, /TEAM RULE/);
  assert.match(out, notice);

  assert.strictEqual(runAt('scripts/rules.js', home, repo, 'trust').status, 0);
  out = runAt('hooks/activate.js', home, repo).stdout;
  assert.match(out, /\[user\.t9, repo\] TEAM RULE/);
  assert.doesNotMatch(out, notice);

  fs.writeFileSync(file, team('CHANGED RULE'));
  out = runAt('hooks/activate.js', home, repo).stdout;
  assert.doesNotMatch(out, /CHANGED RULE|TEAM RULE/);
  assert.match(out, notice);
});

test('activate leaves a data-dir pointer that the CLI uses without CLAUDE_PLUGIN_DATA', () => {
  const home = setup();
  const plugin = path.join(home, 'plugin-data');
  run('activate.js', home, { env: { CLAUDE_PLUGIN_DATA: plugin } });
  assert.strictEqual(fs.readFileSync(path.join(home, '.lucid-coder-data'), 'utf8'), plugin);

  const env = { ...process.env, HOME: home, CLAUDE_CONFIG_DIR: home };
  delete env.CLAUDE_PLUGIN_DATA;
  const r = spawnSync(process.execPath, [path.join(root, 'scripts', 'rules.js'), 'propose',
    '--topic', 'user.x', '--scope', 'global', '--statement', 'S', '--why', 'W'], { cwd: home, env, encoding: 'utf8' });
  assert.strictEqual(r.status, 0, r.stderr);
  assert.ok(fs.existsSync(path.join(plugin, 'rules.json')));
});
