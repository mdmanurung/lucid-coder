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
// The R that run.sh and score.js use: R4_51 holds R, papermill and the ir kernel.
const R_BIN = path.join(process.env.R_ROOT || '/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51', 'bin');
const TIMEOUT = 180000;
const SKILL_SCRIPT = 'skills/pseudobulk-de/pseudobulk_de.R';

function shell(command, cwd) {
  return spawnSync('bash', ['-c', command], { cwd, encoding: 'utf8', env: { ...process.env, PATH: `${R_BIN}:${process.env.PATH}` } });
}
const rscript = (args, cwd) => shell(`Rscript ${args}`, cwd);

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
  if (!hasR) return t.skip('R4_51 R not available');
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
// expectExit: the script's exit status; the task 05 reference refuses the inputs by design.
// A skill task (scriptPath is a dir holding skills/) copies that tree and runs the skill contract.
function gate(taskId, scriptPath, expectExit = 0) {
  const taskDir = path.join(tasksDir, taskId);
  const repo = fs.mkdtempSync(path.join(os.tmpdir(), `gate-${taskId}-`));
  fs.cpSync(path.join(taskDir, 'inputs'), repo, { recursive: true });
  let command;
  if (fs.statSync(scriptPath).isDirectory()) {
    fs.cpSync(path.join(scriptPath, 'skills'), path.join(repo, 'skills'), { recursive: true });
    command = `${SKILL_SCRIPT} --input data --output results/run1`;
  } else {
    fs.mkdirSync(path.join(repo, 'scripts'), { recursive: true });
    const script = path.join('scripts', path.basename(scriptPath));
    fs.copyFileSync(scriptPath, path.join(repo, script));
    command = `${script} results/run1`;
  }
  spawnSync('git', ['init', '-q'], { cwd: repo });
  const run = rscript(command, repo);
  assert.strictEqual(run.status, expectExit, `${scriptPath} exited ${run.status}:\n${run.stderr}`);
  // A skill renders report.md from its notebook, as the 01sk prompt runs it.
  if (fs.existsSync(path.join(repo, 'skills/pseudobulk-de/report.ipynb'))) {
    const render = shell('papermill skills/pseudobulk-de/report.ipynb results/run1/report.ipynb -p out_dir results/run1 && jupyter nbconvert --to markdown --no-input results/run1/report.ipynb', repo);
    assert.strictEqual(render.status, 0, render.stderr);
  }
  const check = rscript(`${path.join(taskDir, 'check.R')} ${repo} ${taskDir}`, repo);
  assert.strictEqual(check.status, 0, check.stderr);
  fs.rmSync(repo, { recursive: true, force: true });
  return { ...JSON.parse(check.stdout.trim()), stderr: run.stderr };
}

const cases = [
  ['01_pseudobulk-de', '01_compare-conditions.R', /GENE17 called/],
  ['02_metadata-join', '02_summarise-samples.R', /wrong or missing condition\/batch/],
  ['03_figure-palette', '03_plot-proportions.R', /proportions_by-condition\.csv/],
  ['04_refactor-legacy', '04_group-means.R', /values differ|rows, expected/],
  ['05_recode-trap', '05_summarise-cohorts.R', /wrong age group:.*S07 \(BR1 -> older\)/],
  ['06_grep-trap', '06_filter-controls.R', /has 0 rows, expected 280/],
  // '.': the skill tree under reference/ (and reference/naive/), not one script.
  ['01sk_pseudobulk-de', '.', /GENE17 called.*SKILL\.md lacks trigger_keywords, Why This Exists/],
];
// Reference scripts that stop on purpose: exit status and the stderr the stop() must carry.
const referenceStops = { '05_recode-trap': [1, /unmapped labels BR 1 \(S07\), br2 \(S19\)/] };

for (const [taskId, script, naiveDetail] of cases) {
  const referenceDir = path.join(tasksDir, taskId, 'reference');
  const [referenceExit, referenceStderr] = referenceStops[taskId] || [0, null];

  test(`${taskId}: reference passes the gate`, { timeout: TIMEOUT }, (t) => {
    if (!hasR) return t.skip('R4_51 R not available');
    const verdict = gate(taskId, path.join(referenceDir, script), referenceExit);
    assert.strictEqual(verdict.gate, 'pass', verdict.detail);
    if (referenceStderr) assert.match(verdict.stderr, referenceStderr);
  });

  test(`${taskId}: naive answer fails the gate`, { timeout: TIMEOUT }, (t) => {
    if (!hasR) return t.skip('R4_51 R not available');
    const verdict = gate(taskId, path.join(referenceDir, 'naive', script));
    assert.strictEqual(verdict.gate, 'fail');
    assert.match(verdict.detail, naiveDetail);
  });
}
