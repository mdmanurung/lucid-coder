// Gate self-tests: fixtures are deterministic, and each task's check.R passes the reference
// script and fails the naive one.
const test = require('node:test');
const assert = require('node:assert');
const crypto = require('crypto');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync } = require('child_process');

const root = path.join(__dirname, '..', '..');
const tasksDir = path.join(root, 'eval', 'tasks');
const R_MODULE = 'module load statistical/R/4.5.2/gcc.8.5.0';
const TIMEOUT = 180000;

function rscript(args, cwd) {
  return spawnSync('bash', ['-lc', `${R_MODULE} && Rscript ${args}`], { cwd, encoding: 'utf8' });
}

const hasR = rscript('-e "invisible(1)"', root).status === 0;

function fileHashes(dir) {
  const hashes = {};
  for (const entry of fs.readdirSync(dir, { recursive: true })) {
    const file = path.join(dir, entry);
    if (fs.statSync(file).isFile()) {
      hashes[entry] = crypto.createHash('sha256').update(fs.readFileSync(file)).digest('hex');
    }
  }
  return hashes;
}

test('make-fixtures.R rewrites identical files', { timeout: TIMEOUT }, (t) => {
  if (!hasR) return t.skip('R module not available');
  const before = fileHashes(tasksDir);
  const run = rscript('eval/fixtures/make-fixtures.R', root);
  assert.strictEqual(run.status, 0, run.stderr);
  // Proves the rerun wrote here, so equal hashes mean "regenerated identically".
  assert.match(run.stderr, /task 01 trap/);
  assert.deepStrictEqual(fileHashes(tasksDir), before);
  for (const [file] of Object.entries(before)) {
    assert.ok(fs.statSync(path.join(tasksDir, file)).size < 1e6, `${file} is over 1 MB`);
  }
});

// Mirror run.sh: inputs copied into a fresh git repo, the script in scripts/, run once.
function gate(taskId, scriptPath) {
  const taskDir = path.join(tasksDir, taskId);
  const repo = fs.mkdtempSync(path.join(os.tmpdir(), `gate-${taskId}-`));
  fs.cpSync(path.join(taskDir, 'inputs'), repo, { recursive: true });
  fs.mkdirSync(path.join(repo, 'scripts'), { recursive: true });
  const script = path.join('scripts', path.basename(scriptPath));
  fs.copyFileSync(scriptPath, path.join(repo, script));
  spawnSync('git', ['init', '-q'], { cwd: repo });
  const run = rscript(`${script} results/run1`, repo);
  assert.strictEqual(run.status, 0, `${scriptPath} failed:\n${run.stderr}`);
  const check = rscript(`${path.join(taskDir, 'check.R')} ${repo} ${taskDir}`, repo);
  assert.strictEqual(check.status, 0, check.stderr);
  fs.rmSync(repo, { recursive: true, force: true });
  return JSON.parse(check.stdout.trim());
}

const cases = [
  ['01_pseudobulk-de', '01_compare-conditions.R', /GENE17 called/],
  ['02_metadata-join', '02_summarise-samples.R', /wrong or missing condition\/batch/],
  ['03_figure-palette', '03_plot-proportions.R', /proportions_by-condition\.csv/],
  ['04_refactor-legacy', '04_group-means.R', /values differ|rows, expected/],
];

for (const [taskId, script, naiveDetail] of cases) {
  const referenceDir = path.join(tasksDir, taskId, 'reference');

  test(`${taskId}: reference passes the gate`, { timeout: TIMEOUT }, (t) => {
    if (!hasR) return t.skip('R module not available');
    const verdict = gate(taskId, path.join(referenceDir, script));
    assert.strictEqual(verdict.gate, 'pass', verdict.detail);
  });

  test(`${taskId}: naive answer fails the gate`, { timeout: TIMEOUT }, (t) => {
    if (!hasR) return t.skip('R module not available');
    const verdict = gate(taskId, path.join(referenceDir, 'naive', script));
    assert.strictEqual(verdict.gate, 'fail');
    assert.match(verdict.detail, naiveDetail);
  });
}
