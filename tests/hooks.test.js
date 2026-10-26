const test = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync } = require('child_process');

const root = path.join(__dirname, '..');
const { escapeRule, render } = require('../hooks/lib');

function setup() {
  const home = fs.mkdtempSync(path.join(os.tmpdir(), 'lucid-'));
  fs.mkdirSync(path.join(home, 'data'));
  return home;
}

function run(script, home, { input = '', env = {} } = {}) {
  const clean = { ...process.env };
  delete clean.LUCID_SUBAGENT_MATCHER;
  return spawnSync(process.execPath, [path.join(root, 'hooks', script)], {
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

test('core ruleset fits its budget', () => {
  assert.ok(fs.readFileSync(path.join(root, 'rules', 'core.md'), 'utf8').length <= 4000);
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
