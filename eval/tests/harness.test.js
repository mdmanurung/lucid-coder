const test = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync } = require('child_process');

const evalDir = path.join(__dirname, '..');
const { COLUMNS, scoreRun } = require('../score');
const { makePacket } = require('../make-packet');
const { kappa, parseJudgeReply, readTsv } = require('../report');

const TS = '20260101T000000Z';
const KEY = 'sk-ant-test-DO-NOT-PRINT-0123456789';

const GOOD_R = `# Purpose: count rows per group in the mock data.
# Inputs: the mock table x.csv in data.
# Outputs: <run_dir>/counts.csv
# Run: Rscript scripts/01_good.R results/run1
# Notes: mock script for the harness tests.
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
out_dir <- args[1]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x <- read.csv(file.path("data", "x.csv"))
counts <- x %>% split(x$group) %>% vapply(nrow, integer(1))
plot_style <- quote(ggplot() + theme_bw(base_size = 8))
write.csv(
  data.frame(group = names(counts), n = counts),
  file.path(out_dir, "counts.csv"),
  row.names = FALSE
)
`;

const BAD_R = `x<-read.csv("data/x.csv")
# y <- x[x$a > 1, ]
y = x |> subset(T)
p <- ggplot() + ggtitle("x")
stop("boom")
`;

const tmp = () => fs.mkdtempSync(path.join(os.tmpdir(), 'lucid-eval-'));

function mockTask(root) {
  const dir = path.join(root, 'tasks', '99_mock');
  fs.mkdirSync(path.join(dir, 'inputs', 'data'), { recursive: true });
  fs.writeFileSync(path.join(dir, 'check.R'), 'cat(\'{"gate":"pass","detail":"mock gate"}\\n\')\n');
  fs.writeFileSync(path.join(dir, 'prompt.md'), 'Write scripts/01_x.R.\n');
  fs.writeFileSync(path.join(dir, 'inputs', 'data', 'x.csv'), 'group,a\ng1,1\ng1,2\ng2,3\n');
  return dir;
}

function mockRun(tsDir, name, { arm, script, claude, taskDir }) {
  const runDir = path.join(tsDir, name);
  fs.mkdirSync(path.join(runDir, 'repo', 'scripts'), { recursive: true });
  fs.cpSync(path.join(taskDir, 'inputs'), path.join(runDir, 'repo'), { recursive: true });
  fs.writeFileSync(path.join(runDir, 'meta.json'), JSON.stringify({
    task: '99_mock', arm, rep: 1, model: 'm', budget: 2, exit_code: 0, task_dir: taskDir,
    arm_channel: 'append-system-prompt-file',
  }));
  if (claude !== undefined) fs.writeFileSync(path.join(runDir, 'claude.json'), claude);
  if (script) fs.writeFileSync(path.join(runDir, 'repo', 'scripts', '01_x.R'), script);
  return runDir;
}

// One scored mock grid shared by the score, packet and report tests.
let scored;
function scoredGrid() {
  if (scored) return scored;
  const root = tmp();
  const taskDir = mockTask(root);
  const tsDir = path.join(root, 'results', TS);
  const envelope = JSON.stringify({ type: 'result', total_cost_usd: 0.42, duration_ms: 61000, num_turns: 7 });
  mockRun(tsDir, '99_mock_A_r1', { arm: 'A', script: BAD_R, claude: '{}', taskDir });
  mockRun(tsDir, '99_mock_C_r1', { arm: 'C', script: GOOD_R, claude: envelope, taskDir });
  mockRun(tsDir, '99_mock_D_r1', { arm: 'D', script: null, claude: 'not json', taskDir });
  const r = spawnSync(process.execPath, [path.join(evalDir, 'score.js'), tsDir], { encoding: 'utf8' });
  assert.strictEqual(r.status, 0, r.stderr);
  const rows = Object.fromEntries(readTsv(path.join(root, 'results', 'results.tsv')).map(x => [x.arm, x]));
  scored = { root, tsDir, taskDir, rows };
  return scored;
}

test('score.js: a good script scores well on every column', () => {
  const g = scoredGrid().rows.C;
  assert.deepStrictEqual(
    [g.run_id, g.cost_usd, g.duration_s, g.turns, g.script_found, g.rerun_ok, g.gate, g.gate_detail],
    [`${TS}/99_mock_C_r1`, '0.42', '61', '7', 'true', 'true', 'pass', 'mock gate']);
  assert.deepStrictEqual(
    [g.air_ok, g.lint_count, g.header_complete, g.check_count, g.commented_code_lines],
    // lint_count 0: its two %>% lints are pipe_consistency, which score.js drops.
    ['true', '0', 'true', '1', '0']);
  assert.deepStrictEqual([g.uptake_pipe, g.uptake_no_title, g.uptake_theme], ['true', 'true', 'true']);
  // It uses %>% and out_dir, and never sets options(warn = 1).
  assert.deepStrictEqual([g.native_pipe, g.rule_out_dir, g.rule_warn1], ['false', 'true', 'false']);
});

test('score.js: a bad script fails style, rerun and uptake', () => {
  const b = scoredGrid().rows.A;
  assert.deepStrictEqual([b.cost_usd, b.duration_s, b.turns], ['NA', 'NA', 'NA'], 'missing fields are NA');
  assert.deepStrictEqual([b.rerun_ok, b.air_ok, b.header_complete], ['false', 'false', 'false']);
  assert.ok(Number(b.lint_count) >= 3, `lint_count ${b.lint_count}`);
  assert.ok(Number(b.commented_code_lines) >= 1);
  assert.strictEqual(b.check_count, '1');
  assert.deepStrictEqual([b.uptake_pipe, b.uptake_no_title, b.uptake_theme], ['false', 'false', 'false']);
  // It pipes with |> only, so the house pipe rule passes even in a bad script.
  assert.deepStrictEqual([b.native_pipe, b.rule_out_dir, b.rule_warn1], ['true', 'false', 'false']);
});

test('score.js: no script means script_found false and dependent columns NA', () => {
  const d = scoredGrid().rows.D;
  assert.strictEqual(d.script_found, 'false');
  for (const c of ['rerun_ok', 'air_ok', 'lint_count', 'header_complete', 'check_count',
    'commented_code_lines', 'uptake_pipe', 'uptake_no_title', 'uptake_theme', 'native_pipe',
    'rule_out_dir', 'rule_warn1', 'cost_usd']) {
    assert.strictEqual(d[c], 'NA', c);
  }
});

test('score.js: header once, fixed columns, rescoring adds no rows', () => {
  const { root, tsDir } = scoredGrid();
  const file = path.join(root, 'results', 'results.tsv');
  spawnSync(process.execPath, [path.join(evalDir, 'score.js'), tsDir]);
  const lines = fs.readFileSync(file, 'utf8').trim().split('\n');
  assert.strictEqual(lines[0], COLUMNS.join('\t'));
  assert.strictEqual(lines.length, 4);
  assert.ok(lines.every(l => l.split('\t').length === COLUMNS.length));
});

test('score.js: a check.R that crashes gives gate NA with the error', () => {
  const root = tmp();
  const taskDir = mockTask(root);
  fs.writeFileSync(path.join(taskDir, 'check.R'), 'stop("harness broke")\n');
  const row = scoreRun(mockRun(path.join(root, TS), 'r', { arm: 'A', script: GOOD_R, claude: '{}', taskDir }));
  assert.strictEqual(row.gate, null);
  assert.match(row.gate_detail, /harness broke/);
});

function runSh(args, env) {
  const clean = { ...process.env };
  delete clean.ANTHROPIC_API_KEY;
  return spawnSync('bash', [path.join(evalDir, 'run.sh'), ...args], { encoding: 'utf8', env: { ...clean, ...env } });
}

test('run.sh --dry-run prints 32 commands with per-arm flags and creates nothing', () => {
  const home = tmp();
  const r = runSh(['--tasks', '01,02,03,04', '--arms', 'A,B,C,D', '--reps', '2', '--dry-run'],
    { HOME: home, EVAL_RESULTS_DIR: path.join(home, 'results') });
  assert.strictEqual(r.status, 0, r.stderr);
  assert.ok(!fs.existsSync(path.join(home, 'results')), 'dry-run creates nothing');
  const lines = r.stdout.trim().split('\n');
  assert.strictEqual(lines.length, 32);
  for (const line of lines) {
    const arm = line.match(/_([ABCD])_r[12]\/repo /)[1];
    assert.match(line, /claude --safe-mode --no-session-persistence -p "\$\(cat .*prompt\.md\)" --model claude-opus-5-5 --max-budget-usd 2 --output-format json --permission-mode acceptEdits/);
    assert.doesNotMatch(line, /--bare|CLAUDE_CONFIG_DIR/);
    assert.match(line, /--allowedTools Read Write Edit Bash\\\(Rscript:\\\*\\\)/);
    assert.doesNotMatch(line, /--plugin-dir/);
    if (arm === 'A') {
      assert.doesNotMatch(line, /append-system-prompt-file|arm-prompt/);
    } else {
      assert.match(line, new RegExp(`arm-prompt\\.js ${arm} > \\S+_${arm}_r[12]/arm-prompt\\.md && `));
      assert.match(line, new RegExp(`--append-system-prompt-file \\S+_${arm}_r[12]/arm-prompt\\.md`));
    }
  }
  assert.strictEqual(lines.filter(l => /_A_r/.test(l)).length, 8);
});

test('run.sh real mode with a stub claude: layout, login auth, arm file, meta, refusal', () => {
  const root = tmp();
  mockTask(root);
  const home = tmp();
  const bin = path.join(root, 'bin');
  fs.mkdirSync(bin);
  // The stub runs in <run_dir>/repo, so ../ is the run dir.
  fs.writeFileSync(path.join(bin, 'claude'), [
    '#!/usr/bin/env bash',
    'printf "%s\\n" "$@" > ../args.txt',
    'echo "${ANTHROPIC_API_KEY:-unset}" > ../key.txt',
    'echo "{\\"total_cost_usd\\": 0.1, \\"duration_ms\\": 1000, \\"num_turns\\": 2}"',
  ].join('\n'), { mode: 0o755 });
  // An exported key must not reach claude: runs use the login, not the API.
  const env = { HOME: home, PATH: `${bin}:${process.env.PATH}`, EVAL_TASKS_DIR: path.join(root, 'tasks'),
    EVAL_RESULTS_DIR: path.join(root, 'results'), ANTHROPIC_API_KEY: KEY };
  const r = runSh(['--tasks', '99', '--arms', 'A,C', '--reps', '1'], env);
  assert.strictEqual(r.status, 0, r.stderr);
  const [ts] = fs.readdirSync(path.join(root, 'results'));
  const runA = path.join(root, 'results', ts, '99_mock_A_r1');
  const runC = path.join(root, 'results', ts, '99_mock_C_r1');
  assert.ok(fs.existsSync(path.join(runA, 'repo', '.git', 'HEAD')));
  assert.ok(fs.existsSync(path.join(runA, 'repo', 'data', 'x.csv')));
  assert.strictEqual(fs.readFileSync(path.join(runA, 'key.txt'), 'utf8').trim(), 'unset');
  assert.ok(!fs.existsSync(path.join(runA, 'arm-prompt.md')));
  const argsA = fs.readFileSync(path.join(runA, 'args.txt'), 'utf8');
  assert.ok(argsA.startsWith('--safe-mode\n--no-session-persistence\n-p\n'));
  assert.doesNotMatch(argsA, /append-system-prompt-file/);
  const argsC = fs.readFileSync(path.join(runC, 'args.txt'), 'utf8');
  assert.ok(argsC.includes(`--append-system-prompt-file\n${path.join(runC, 'arm-prompt.md')}\n`));
  assert.match(fs.readFileSync(path.join(runC, 'arm-prompt.md'), 'utf8'), /^LUCID-CODER ACTIVE/);
  const meta = JSON.parse(fs.readFileSync(path.join(runC, 'meta.json'), 'utf8'));
  assert.deepStrictEqual([meta.task, meta.arm, meta.rep, meta.exit_code, meta.arm_channel],
    ['99_mock', 'C', 1, 0, 'append-system-prompt-file']);
  assert.strictEqual(JSON.parse(fs.readFileSync(path.join(runA, 'meta.json'), 'utf8')).arm_channel, 'none');
  // A task listed twice maps to the same run dir: the second must be refused.
  const twice = runSh(['--tasks', '99,99', '--arms', 'A', '--reps', '1'], { ...env, EVAL_RESULTS_DIR: path.join(root, 'results2') });
  assert.strictEqual(twice.status, 1);
  assert.match(twice.stderr, /99_mock_A_r1 already exists/);
});

test('arm-prompt.js: D carries the three rules, C does not, A is empty', () => {
  const arm = a => spawnSync(process.execPath, [path.join(evalDir, 'arm-prompt.js'), a], { encoding: 'utf8' }).stdout;
  const statements = JSON.parse(fs.readFileSync(path.join(evalDir, 'arms', 'rules-d.json'), 'utf8')).rules.map(r => r.statement);
  assert.strictEqual(statements.length, 3);
  const c = arm('C');
  const d = arm('D');
  assert.match(c, /^LUCID-CODER ACTIVE/);
  for (const s of statements) {
    assert.ok(d.includes(s), s);
    assert.ok(!c.includes(s), s);
  }
  assert.strictEqual(arm('A'), '');
});

test('make-packet.js: blinded copies, empty sheet, key outside, refuses a reshuffle', () => {
  const { tsDir } = scoredGrid();
  assert.strictEqual(makePacket(tsDir, 'packet'), 2, 'the run with no script is left out');
  const files = fs.readdirSync(path.join(tsDir, 'packet')).sort();
  assert.deepStrictEqual(files, ['P01.R', 'P02.R', 'sheet.tsv']);
  const key = readTsv(path.join(tsDir, 'packet-key.tsv'));
  assert.deepStrictEqual(key.map(k => k.run_id).sort(), [`${TS}/99_mock_A_r1`, `${TS}/99_mock_C_r1`]);
  for (const f of ['P01.R', 'P02.R']) {
    const text = fs.readFileSync(path.join(tsDir, 'packet', f), 'utf8');
    assert.ok(text === GOOD_R || text === BAD_R);
  }
  assert.throws(() => makePacket(tsDir, 'packet'), /exists/);
});

test('judge.sh with a stub claude writes judges.tsv from fenced JSON replies', () => {
  const root = tmp();
  const taskDir = mockTask(root);
  const tsDir = path.join(root, TS);
  mockRun(tsDir, 'a', { arm: 'A', script: BAD_R, claude: '{}', taskDir });
  mockRun(tsDir, 'c', { arm: 'C', script: GOOD_R, claude: '{}', taskDir });
  const stub = path.join(root, 'claude-stub');
  fs.writeFileSync(stub, [
    '#!/usr/bin/env bash',
    'model=""; prev=""; for a in "$@"; do [[ "$prev" == --model ]] && model="$a"; prev="$a"; done',
    'q1=4; [[ "$model" == *sonnet* ]] && q1=2',
    'printf \'{"type":"result","total_cost_usd":0.01,"result":"```json\\\\n{\\\\"q1_one_read\\\\": %s, \\\\"q2_find_decisions\\\\": 3, \\\\"q3_change_checks\\\\": \\\\"no\\\\", \\\\"q3_which\\\\": \\\\"\\\\", \\\\"notes\\\\": \\\\"clear\\\\"}\\\\n```"}\' "$q1"',
  ].join('\n'), { mode: 0o755 });
  const r = spawnSync('bash', [path.join(evalDir, 'judge.sh'), tsDir], {
    encoding: 'utf8', env: { ...process.env, ANTHROPIC_API_KEY: KEY, JUDGE_CLAUDE: stub },
  });
  assert.strictEqual(r.status, 0, r.stderr);
  const rows = readTsv(path.join(tsDir, 'judges.tsv'));
  assert.strictEqual(rows.length, 4);
  assert.deepStrictEqual([...new Set(rows.map(x => x.neutral_id))].sort(), ['J01', 'J02']);
  for (const row of rows) {
    assert.strictEqual(row.q1_one_read, row.judge === 'claude-sonnet-5-5' ? '2' : '4');
    assert.deepStrictEqual([row.q2_find_decisions, row.q3_change_checks, row.notes], ['3', 'no', 'clear']);
  }
  assert.strictEqual(readTsv(path.join(tsDir, 'judge-key.tsv')).length, 2);
});

test('parseJudgeReply: bare JSON passes, bad replies give an error', () => {
  const wrap = result => JSON.stringify({ type: 'result', result });
  assert.deepStrictEqual(
    parseJudgeReply(wrap('{"q1_one_read":5,"q2_find_decisions":1,"q3_change_checks":"Yes","q3_which":"add a row check","notes":""}')),
    { q1_one_read: 5, q2_find_decisions: 1, q3_change_checks: 'yes', q3_which: 'add a row check', notes: '' });
  assert.ok(parseJudgeReply('oops').error);
  assert.ok(parseJudgeReply(wrap('I think it is a 4.')).error);
  assert.ok(parseJudgeReply(wrap('{"q1_one_read":7,"q2_find_decisions":1,"q3_change_checks":"no"}')).error);
});

test('kappa: textbook and hand-derived values', () => {
  // Wikipedia "Cohen's kappa" worked example: 50 items, yes/yes 20, yes/no 5,
  // no/yes 10, no/no 15. p_o = 0.7, p_e = 0.5*0.6 + 0.5*0.4 = 0.5, kappa = 0.4.
  const a = [...Array(25).fill('yes'), ...Array(25).fill('no')];
  const b = [...Array(20).fill('yes'), ...Array(5).fill('no'), ...Array(10).fill('yes'), ...Array(15).fill('no')];
  assert.ok(Math.abs(kappa(a, b, ['no', 'yes'], 'none') - 0.4) < 1e-12);
  // With two categories every weighting equals the unweighted kappa.
  assert.ok(Math.abs(kappa(a, b, ['no', 'yes'], 'quadratic') - 0.4) < 1e-12);
  // Ordinal, k = 3: a = 1,1,2,3 and b = 1,2,2,3. Quadratic weights d^2/4:
  // observed 0.25/4 = 0.0625, expected 0.3125, kappa = 1 - 0.2 = 0.8.
  // Linear weights d/2: observed 0.125, expected 0.4375, kappa = 5/7.
  assert.ok(Math.abs(kappa([1, 1, 2, 3], [1, 2, 2, 3], [1, 2, 3], 'quadratic') - 0.8) < 1e-12);
  assert.ok(Math.abs(kappa([1, 1, 2, 3], [1, 2, 2, 3], [1, 2, 3], 'linear') - 5 / 7) < 1e-12);
  assert.strictEqual(kappa([1, 2, 3], [1, 2, 3], [1, 2, 3]), 1);
  assert.strictEqual(kappa([2, 2], [2, 2], [1, 2, 3]), null, 'one category used: undefined');
});

test('report.js renders without and with human sheets, never a combined score', () => {
  const { root, tsDir } = scoredGrid();
  if (!fs.existsSync(path.join(tsDir, 'packet'))) makePacket(tsDir, 'packet');
  const report = args => {
    const r = spawnSync(process.execPath, [path.join(evalDir, 'report.js'), tsDir, ...args], { encoding: 'utf8' });
    assert.strictEqual(r.status, 0, r.stderr);
    return fs.readFileSync(path.join(root, 'results', 'report.md'), 'utf8');
  };
  const bare = report([]);
  for (const h of ['## L1', '## L2', '## L3', '## L4', '## LLM judges (proxy, not human evidence)']) {
    assert.ok(bare.includes(h), h);
  }
  assert.match(bare, /Arms delivered by --append-system-prompt-file in --safe-mode sessions, because plugin hooks do not run in an isolated session \(step 0, 2026-10-26\)/);
  assert.match(bare, /Pending: no filled reader sheets/);
  assert.match(bare, /n per task × arm cell = 1/);
  assert.match(bare, /\| gate pass \| 1\/1 \| 1\/1 \| 1\/1 \|/);
  assert.match(bare, /\| rerun ok \| 0\/1 \| 1\/1 \| 0\/0 \|/);
  assert.match(bare, /no combined score/);

  const sheet = readTsv(path.join(tsDir, 'packet', 'sheet.tsv'));
  const fill = (name, q1) => {
    const file = path.join(root, name);
    fs.writeFileSync(file, 'neutral_id\tq1_one_read\tq2_find_decisions\tq3_change_checks\tq3_which\tnotes\n' +
      sheet.map((s, i) => `${s.neutral_id}\t${q1[i]}\t3\t${i ? 'yes' : 'no'}\t\t`).join('\n') + '\n');
    return file;
  };
  const withHuman = report(['--human', fill('r1.tsv', [2, 5]), '--human', fill('r2.tsv', [3, 5])]);
  assert.doesNotMatch(withHuman, /Pending/);
  assert.match(withHuman, /Agreement, reader 1 \(r1\.tsv\) vs reader 2 \(r2\.tsv\)/);
  // q1 pairs (2,3) and (5,5) on 1-5, quadratic: 1 - (1/32) / (7/32) = 6/7.
  assert.match(withHuman, /q1 0\.86 \(n = 2\)/);
});
