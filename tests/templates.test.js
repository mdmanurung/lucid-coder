const test = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync } = require('child_process');

const root = path.join(__dirname, '..');
const templates = path.join(root, 'templates');
const examples = path.join(templates, 'examples');
// Host tools (ruff, air, uv) are installed under ~/.local/bin, which a bare spawn may not see.
// PYTHONDONTWRITEBYTECODE keeps the example run from leaving __pycache__ in templates/.
const env = { ...process.env, PYTHONDONTWRITEBYTECODE: '1', PATH: `${path.join(os.homedir(), '.local', 'bin')}:${process.env.PATH}` };
const R_MODULE = 'module load statistical/R/4.5.2/gcc.8.5.0';

function sh(cmd, args, cwd = templates) {
  return spawnSync(cmd, args, { cwd, env, encoding: 'utf8' });
}

function rscript(args, cwd = templates) {
  return sh('bash', ['-lc', `${R_MODULE} && Rscript ${args}`], cwd);
}

function has(cmd) {
  return sh(cmd, ['--version']).status === 0;
}

// Scripts refuse an existing run dir, so hand them a child of the temp dir.
function newRunDir() {
  return path.join(fs.mkdtempSync(path.join(os.tmpdir(), 'lucid-run-')), 'run');
}

function assertRunOutputs(runDir) {
  for (const f of ['contrasts.tsv', 'top-gene_by-condition.pdf', 'top-gene_by-condition.png',
    'top-gene_by-condition.csv', 'session_info.txt']) {
    assert.ok(fs.statSync(path.join(runDir, f)).size > 0, f);
  }
}

const EXAMPLE_R = 'examples/01_fit-contrasts.R';
const EXAMPLE_PY = 'examples/01_fit-contrasts.py';
const hasR = rscript('-e "invisible(1)"').status === 0;
const rTimeout = 120000;

test('template files exist', () => {
  for (const f of ['air.toml', '.lintr', 'ruff.toml', 'palettes.R', 'palettes.py', 'theme.R', 'theme.py',
    EXAMPLE_R, EXAMPLE_PY]) {
    assert.ok(fs.existsSync(path.join(templates, f)), f);
  }
});

test('ruff format and check pass on templates', (t) => {
  if (!has('ruff')) return t.skip('ruff not installed');
  const fmt = sh('ruff', ['format', '--no-cache', '--check', '--config', 'ruff.toml', '.']);
  assert.strictEqual(fmt.status, 0, fmt.stdout + fmt.stderr);
  const lint = sh('ruff', ['check', '--no-cache', '--config', 'ruff.toml', '.']);
  assert.strictEqual(lint.status, 0, lint.stdout + lint.stderr);
});

test('air format check passes on the R templates', (t) => {
  if (!has('air')) return t.skip('air not installed');
  const r = sh('air', ['format', '--check', EXAMPLE_R, 'palettes.R', 'theme.R']);
  assert.strictEqual(r.status, 0, r.stdout + r.stderr);
});

test('lintr finds nothing in the R example', { timeout: rTimeout }, (t) => {
  if (!hasR) return t.skip('R module not available');
  const r = rscript(`-e 'l <- lintr::lint("${EXAMPLE_R}"); print(l); quit(status = length(l))'`);
  assert.strictEqual(r.status, 0, r.stdout + r.stderr);
});

test('palettes.R self-check passes', { timeout: rTimeout }, (t) => {
  if (!hasR) return t.skip('R module not available');
  const r = rscript('palettes.R');
  assert.strictEqual(r.status, 0, r.stderr);
});

test('R example runs into a new run dir', { timeout: rTimeout }, (t) => {
  if (!hasR) return t.skip('R module not available');
  const runDir = newRunDir();
  const r = rscript(`${EXAMPLE_R} ${runDir}`);
  assert.strictEqual(r.status, 0, r.stderr);
  assertRunOutputs(runDir);
  assert.notStrictEqual(rscript(`${EXAMPLE_R} ${runDir}`).status, 0, 'existing run dir must be refused');
});

const uvTimeout = 600000; // first call downloads a Python and the wheels
const uvPython = ['run', '--python', '3.12', '--with', 'matplotlib,pandas,scipy', 'python'];

test('palettes.py self-check passes', { timeout: uvTimeout }, (t) => {
  if (!has('uv')) return t.skip('uv not installed');
  const r = sh('uv', [...uvPython, 'palettes.py']);
  assert.strictEqual(r.status, 0, r.stderr);
});

test('Python example runs into a new run dir', { timeout: uvTimeout }, (t) => {
  if (!has('uv')) return t.skip('uv not installed');
  const runDir = newRunDir();
  const r = sh('uv', [...uvPython, EXAMPLE_PY, runDir]);
  assert.strictEqual(r.status, 0, r.stderr);
  assertRunOutputs(runDir);
  assert.notStrictEqual(sh('uv', [...uvPython, EXAMPLE_PY, runDir]).status, 0, 'existing run dir must be refused');
});

// PLAN.md phase 2 done-when: after align copies the configs into a repo, the
// formatters pass on the example scripts with no --config flag.
test('copied configs are picked up by ruff and air with no flags', (t) => {
  const repo = fs.mkdtempSync(path.join(os.tmpdir(), 'lucid-align-'));
  fs.mkdirSync(path.join(repo, '.git'));
  fs.writeFileSync(path.join(repo, '.git', 'HEAD'), 'ref: refs/heads/main\n');
  for (const f of fs.readdirSync(templates)) {
    if (fs.statSync(path.join(templates, f)).isFile()) fs.copyFileSync(path.join(templates, f), path.join(repo, f));
  }
  fs.mkdirSync(path.join(repo, 'scripts'));
  for (const f of fs.readdirSync(examples)) {
    fs.copyFileSync(path.join(examples, f), path.join(repo, 'scripts', f));
  }
  let checked = 0;
  if (has('ruff')) {
    const r = sh('ruff', ['format', '--no-cache', '--check', 'scripts'], repo);
    assert.strictEqual(r.status, 0, r.stdout + r.stderr);
    checked++;
  }
  if (has('air')) {
    const r = sh('air', ['format', '--check', 'scripts'], repo);
    assert.strictEqual(r.status, 0, r.stdout + r.stderr);
    checked++;
  }
  if (checked === 0) t.skip('neither ruff nor air installed');
});
