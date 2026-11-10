const test = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync } = require('child_process');

const root = path.join(__dirname, '..');

// A temp home plus a fake git repo holding a team rules file.
function setup() {
  const home = fs.mkdtempSync(path.join(os.tmpdir(), 'lucid-trust-'));
  const repo = path.join(home, 'repo');
  fs.mkdirSync(path.join(repo, '.git'), { recursive: true });
  fs.writeFileSync(path.join(repo, '.git', 'config'), '[remote "origin"]\n\turl = git@github.com:lab/atlas.git\n');
  fs.writeFileSync(path.join(repo, '.git', 'HEAD'), 'ref: refs/heads/main\n');
  const file = path.join(repo, 'lucid-coder.rules.json');
  fs.writeFileSync(file, '{"rules": []}\n');
  return { home, repo, file, store: path.join(home, 'data', 'rules.json') };
}

function cli(ctx, ...args) {
  return spawnSync(process.execPath, [path.join(root, 'scripts', 'rules.js'), ...args], {
    cwd: ctx.repo,
    env: { ...process.env, HOME: ctx.home, CLAUDE_CONFIG_DIR: ctx.home, CLAUDE_PLUGIN_DATA: path.join(ctx.home, 'data') },
    encoding: 'utf8',
  });
}

const trusted = ctx => (fs.existsSync(ctx.store) ? JSON.parse(fs.readFileSync(ctx.store, 'utf8')).trusted : []);

test('trust --show prints the file and its hash and writes nothing', () => {
  const ctx = setup();
  const r = cli(ctx, 'trust', '--show');
  assert.strictEqual(r.status, 0, r.stderr);
  assert.match(r.stdout, /^\{"rules": \[\]\}\nsha256 [0-9a-f]{64}\n$/);
  assert.ok(!fs.existsSync(ctx.store));
});

test('trust --sha256 trusts only the reviewed content', () => {
  const ctx = setup();
  const hash = cli(ctx, 'trust', '--show').stdout.match(/sha256 ([0-9a-f]{64})/)[1];

  // The file changes between review and trust: refused, nothing recorded.
  fs.writeFileSync(ctx.file, '{"rules": [{"injected": true}]}\n');
  const stale = cli(ctx, 'trust', '--sha256', hash);
  assert.strictEqual(stale.status, 1);
  assert.match(stale.stderr, /Not trusted: .* not the reviewed/);
  assert.deepStrictEqual(trusted(ctx), []);

  const fresh = cli(ctx, 'trust', '--show').stdout.match(/sha256 ([0-9a-f]{64})/)[1];
  assert.strictEqual(cli(ctx, 'trust', '--sha256', fresh.toUpperCase()).status, 0);
  assert.deepStrictEqual(trusted(ctx).map(t => t.sha256), [fresh]);
});
