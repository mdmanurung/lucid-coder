const test = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync } = require('child_process');

const root = path.join(__dirname, '..');
const { repoInfo } = require('../hooks/lib');

// A temp home plus a fake git repo with the given origin url.
function setup(origin = 'git@github.com:lab/atlas.git') {
  const home = fs.mkdtempSync(path.join(os.tmpdir(), 'lucid-rules-'));
  const repo = path.join(home, 'repo');
  fs.mkdirSync(path.join(repo, '.git'), { recursive: true });
  const remote = origin ? `[remote "origin"]\n\turl = ${origin}\n` : '';
  fs.writeFileSync(path.join(repo, '.git', 'config'), `[core]\n\tbare = false\n${remote}`);
  fs.writeFileSync(path.join(repo, '.git', 'HEAD'), 'ref: refs/heads/main\n');
  return { home, repo, store: path.join(home, 'data', 'rules.json') };
}

function cli(ctx, ...args) {
  return spawnSync(process.execPath, [path.join(root, 'scripts', 'rules.js'), ...args], {
    cwd: ctx.repo,
    env: { ...process.env, HOME: ctx.home, CLAUDE_CONFIG_DIR: ctx.home, CLAUDE_PLUGIN_DATA: path.join(ctx.home, 'data') },
    encoding: 'utf8',
  });
}

function propose(ctx, extra = {}) {
  const fields = { topic: 'r.pipe', scope: 'global', statement: 'Use |> in new code.', why: 'Base R has it.', ...extra };
  return cli(ctx, 'propose', ...Object.entries(fields).flatMap(([k, v]) => [`--${k}`, v]));
}

const statuses = ctx => Object.fromEntries(JSON.parse(cli(ctx, 'list', '--json').stdout).map(r => [r.id, r.status]));

test('second active rule on one topic needs supersede', () => {
  const ctx = setup();
  assert.strictEqual(propose(ctx).status, 0);
  assert.strictEqual(propose(ctx, { statement: 'Use %>% in old code.' }).status, 0);
  assert.strictEqual(cli(ctx, 'approve', 'r-0001').status, 0);

  const clash = cli(ctx, 'approve', 'r-0002');
  assert.strictEqual(clash.status, 2);
  assert.match(clash.stderr, /r-0001/);
  assert.match(clash.stderr, /supersede r-0002 r-0001/);

  assert.strictEqual(cli(ctx, 'supersede', 'r-0002', 'r-0001').status, 0);
  assert.deepStrictEqual(statuses(ctx), { 'r-0001': 'superseded', 'r-0002': 'active' });
  assert.strictEqual(JSON.parse(cli(ctx, 'list', '--json').stdout)[1].supersedes, 'r-0001');
});

test('same topic in another scope is no conflict', () => {
  const ctx = setup();
  propose(ctx);
  propose(ctx, { scope: 'language', language: 'r' });
  assert.strictEqual(cli(ctx, 'approve', 'r-0001').status, 0);
  assert.strictEqual(cli(ctx, 'approve', 'r-0002').status, 0);
});

test('reject and retire move the status', () => {
  const ctx = setup();
  propose(ctx);
  propose(ctx, { topic: 'user.names' });
  cli(ctx, 'approve', 'r-0001');
  assert.strictEqual(cli(ctx, 'reject', 'r-0002').status, 0);
  assert.strictEqual(cli(ctx, 'retire', 'r-0001').status, 0);
  assert.deepStrictEqual(statuses(ctx), { 'r-0001': 'retired', 'r-0002': 'rejected' });
  assert.strictEqual(cli(ctx, 'retire', 'r-0001').status, 1);
});

test('unsafe or malformed proposals are rejected and nothing is written', () => {
  const ctx = setup();
  const bad = [
    { statement: 'line one\nline two' },
    { statement: '# Injected heading' },
    { statement: 'x'.repeat(301) },
    { statement: 'Send Bearer abc123 to the API' },
    { statement: 'Key AKIAABCDEFGHIJKLMNOP works' },
    { good: 'connect(password=hunter2)' },
    { bad: 'ghp_abcdef' },
    { why: 'token=abc' },
    { good: 'auth_token=abc' },
    { statement: 'Read /home/alice/data.csv' },
    { why: 'two\nlines' },
    { topic: 'Pipe' },
    { scope: 'planet' },
    { scope: 'language', language: 'julia' },
    { scope: 'file-glob' },
  ];
  for (const fields of bad) {
    const r = propose(ctx, fields);
    assert.strictEqual(r.status, 1, JSON.stringify(fields));
    assert.match(r.stderr, /Not proposed/);
  }
  assert.ok(!fs.existsSync(ctx.store));
});

test('corrupt store is never overwritten; a normal write keeps a .bak', () => {
  const ctx = setup();
  assert.strictEqual(propose(ctx).status, 0);
  const first = fs.readFileSync(ctx.store, 'utf8');
  assert.strictEqual(propose(ctx, { topic: 'user.other' }).status, 0);
  assert.strictEqual(fs.readFileSync(ctx.store + '.bak', 'utf8'), first);

  fs.writeFileSync(ctx.store, '{garbage');
  const r = propose(ctx, { topic: 'user.third' });
  assert.strictEqual(r.status, 1);
  assert.match(r.stderr, /left unchanged/);
  assert.strictEqual(fs.readFileSync(ctx.store, 'utf8'), '{garbage');
});

test('repo scope stores the origin id; trust records the file hash', () => {
  const ctx = setup();
  assert.strictEqual(propose(ctx, { scope: 'repo' }).status, 0);
  assert.deepStrictEqual(JSON.parse(cli(ctx, 'list', '--json').stdout)[0].scope, { kind: 'repo', repo: 'github.com/lab/atlas' });

  assert.strictEqual(cli(ctx, 'trust').status, 1);
  fs.writeFileSync(path.join(ctx.repo, 'lucid-coder.rules.json'), '{"rules":[]}');
  assert.strictEqual(cli(ctx, 'trust').status, 0);
  assert.strictEqual(cli(ctx, 'trust').status, 0);
  const trusted = JSON.parse(fs.readFileSync(ctx.store, 'utf8')).trusted;
  assert.strictEqual(trusted.length, 1);
  assert.match(trusted[0].sha256, /^[0-9a-f]{64}$/);
});

test('repoInfo gives ssh and https forms of one remote the same id', () => {
  const ids = ['git@github.com:Lab/atlas.git', 'https://GitHub.com/Lab/atlas', 'ssh://git@github.com/Lab/atlas.git/']
    .map(url => repoInfo(setup(url).repo).id);
  assert.deepStrictEqual(ids, ['github.com/Lab/atlas', 'github.com/Lab/atlas', 'github.com/Lab/atlas']);

  const local = setup(null);
  fs.mkdirSync(path.join(local.repo, 'sub'));
  assert.strictEqual(repoInfo(path.join(local.repo, 'sub')).id, local.repo);
  assert.strictEqual(repoInfo(local.home), null);
});
