// Python twins: each task's check.R passes the reference script and fails the naive one,
// and score.js scores the reference script through ruff and uv.
const test = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync } = require('child_process');

const { scoreRun } = require('../score');

const root = path.join(__dirname, '..', '..');
const tasksDir = path.join(root, 'eval', 'tasks');
const R_MODULE = 'module load statistical/R/4.5.2/gcc.8.5.0';
const UV = path.join(os.homedir(), '.local', 'bin', 'uv');
const UV_ARGS = ['run', '--python', '3.12', '--with', 'pandas,scipy,matplotlib', 'python'];
// The first uv call may fetch the interpreter and packages.
const TIMEOUT = 600000;

const rscript = (args, cwd) =>
  spawnSync('bash', ['-lc', `${R_MODULE} && Rscript ${args}`], { cwd, encoding: 'utf8' });
const hasR = rscript('-e "invisible(1)"', root).status === 0;
const hasUv = spawnSync(UV, ['--version']).status === 0;
const skipReason = !hasR ? 'R module not available' : !hasUv ? 'uv not available' : null;

// Mirror run.sh: inputs copied into a fresh git repo, the script in scripts/, run once.
function repoWith(taskId, scriptPath) {
  const repo = fs.mkdtempSync(path.join(os.tmpdir(), `py-${taskId}-`));
  fs.cpSync(path.join(tasksDir, taskId, 'inputs'), repo, { recursive: true });
  fs.mkdirSync(path.join(repo, 'scripts'));
  const script = path.join('scripts', path.basename(scriptPath));
  fs.copyFileSync(scriptPath, path.join(repo, script));
  spawnSync('git', ['init', '-q'], { cwd: repo });
  return { repo, script };
}

function gate(taskId, scriptPath) {
  const taskDir = path.join(tasksDir, taskId);
  const { repo, script } = repoWith(taskId, scriptPath);
  const run = spawnSync(UV, [...UV_ARGS, script, 'results/run1'], { cwd: repo, encoding: 'utf8', timeout: TIMEOUT });
  assert.strictEqual(run.status, 0, `${scriptPath} failed:\n${run.stderr}`);
  const check = rscript(`${path.join(taskDir, 'check.R')} ${repo} ${taskDir}`, repo);
  assert.strictEqual(check.status, 0, check.stderr);
  fs.rmSync(repo, { recursive: true, force: true });
  return JSON.parse(check.stdout.trim());
}

const cases = [
  ['01py_pseudobulk-de', '01_compare-conditions.py', /GENE17 called/],
  ['02py_metadata-join', '02_summarise-samples.py', /wrong or missing condition\/batch/],
  ['03py_figure-palette', '03_plot-proportions.py', /proportions_by-condition\.csv/],
];

for (const [taskId, script, naiveDetail] of cases) {
  const referenceDir = path.join(tasksDir, taskId, 'reference');

  test(`${taskId}: reference passes the gate`, { timeout: TIMEOUT }, (t) => {
    if (skipReason) return t.skip(skipReason);
    const verdict = gate(taskId, path.join(referenceDir, script));
    assert.strictEqual(verdict.gate, 'pass', verdict.detail);
  });

  test(`${taskId}: naive answer fails the gate`, { timeout: TIMEOUT }, (t) => {
    if (skipReason) return t.skip(skipReason);
    const verdict = gate(taskId, path.join(referenceDir, 'naive', script));
    assert.strictEqual(verdict.gate, 'fail');
    assert.match(verdict.detail, naiveDetail);
  });

  test(`${taskId}: score.js scores the reference script`, { timeout: TIMEOUT }, (t) => {
    if (skipReason) return t.skip(skipReason);
    const runDir = fs.mkdtempSync(path.join(os.tmpdir(), `py-score-${taskId}-`));
    const taskDir = path.join(tasksDir, taskId);
    const { repo } = repoWith(taskId, path.join(referenceDir, script));
    fs.cpSync(repo, path.join(runDir, 'repo'), { recursive: true });
    fs.writeFileSync(path.join(runDir, 'meta.json'), JSON.stringify({ task: taskId, arm: 'A', rep: 1, exit_code: 0, task_dir: taskDir }));
    const row = scoreRun(runDir);
    assert.deepStrictEqual(
      [row.script_found, row.rerun_ok, row.air_ok, row.lint_count, row.commented_code_lines, row.header_complete],
      [true, true, true, 0, 0, true]);
    assert.ok(row.check_count >= 1, `check_count ${row.check_count}`);
    assert.deepStrictEqual([row.uptake_pipe, row.uptake_theme, row.native_pipe, row.rule_warn1], [null, null, null, null]);
    assert.strictEqual(row.uptake_no_title, taskId.startsWith('03') ? true : null);
    fs.rmSync(repo, { recursive: true, force: true });
    fs.rmSync(runDir, { recursive: true, force: true });
  });
}
