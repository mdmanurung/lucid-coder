#!/usr/bin/env node
// Scores every run under eval/results/<timestamp>/ and appends one row per run
// to eval/results/results.tsv (the dir that holds <timestamp>). Runs already in
// the file are skipped, so rescoring never duplicates rows.
// Usage: node eval/score.js eval/results/<timestamp>
const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync } = require('child_process');
const { codeMeasures } = require('./measures');

const COLUMNS = ['run_id', 'task', 'arm', 'rep', 'exit_code', 'cost_usd', 'duration_s', 'turns',
  'script_found', 'rerun_ok', 'gate', 'gate_detail', 'air_ok', 'lint_count', 'header_complete',
  'check_count', 'commented_code_lines', 'uptake_pipe', 'uptake_no_title', 'uptake_theme',
  // Added after the pilot (D1 revised, new arm D rules); pilot rows hold NA here.
  'native_pipe', 'rule_out_dir', 'rule_warn1',
  // Added after rerun4; older rows lack it and readers treat the missing cell as NA.
  'join_guard_count',
  // Rubric G items (eval/rubric.md, eval/measures.js); added after grid 5.
  'code_lines', 'joins_unkeyed', 'warning_suppress', 'longest_pipe', 'max_nesting',
  'vague_names', 'wrappers_once', 'session_info', 'late_loads', 'base_df_calls', 'upper_constants',
  // Added with the 3+ uses function rule (2026-10-27).
  'functions_under3'];
const EVAL_DIR = __dirname;
const TEMPLATES = path.join(EVAL_DIR, '..', 'templates');
// R4_51 conda env (R 4.5.1): tidyverse, limma, edgeR, IRkernel, papermill, pandoc. Override with R_ROOT.
const R_ROOT = process.env.R_ROOT || '/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51';
const RERUN_MS = 120000;
const PY_RERUN_MS = 600000;
const GATE_MS = 600000;
// rules/core.md asks for Notes only when a decision is invisible, so it is not required.
const HEADER_WORDS = ['Purpose', 'Inputs', 'Outputs', 'Run'];
// Two lints are dropped so no arm is penalised for following its own setup:
// pipe_consistency_linter fires on every %>%, which arm D's r.pipe rule asks for;
// commented_code_linter fires on header paths ("# Inputs: data/x.tsv"), which grow
// as headers get more complete. Only header lines (before the first code line) are exempt.
const LINT_R = 'f <- commandArgs(TRUE)[1]; l <- lintr::lint(f); txt <- readLines(f); ' +
  'first_code <- c(which(!grepl("^\\\\s*(#|$)", txt)), length(txt) + 1)[1]; ' +
  'k <- vapply(l, function(x) x$linter, ""); ln <- vapply(l, function(x) as.integer(x$line_number), 1L); ' +
  'keep <- k != "pipe_consistency_linter" & !(k == "commented_code_linter" & ln < first_code); ' +
  'cat(sum(keep), sum(keep & k == "commented_code_linter"), "\\n")';

const toolEnv = () => ({
  ...process.env,
  PATH: `${R_ROOT}/bin:${process.env.PATH}:${path.join(os.homedir(), '.local', 'bin')}`,
  LD_LIBRARY_PATH: `${R_ROOT}/lib64/R/lib:${process.env.LD_LIBRARY_PATH || ''}`,
});

const run = (cmd, args, opts) => spawnSync(cmd, args, { encoding: 'utf8', env: toolEnv(), ...opts });

function readJson(file) {
  try { return JSON.parse(fs.readFileSync(file, 'utf8')); } catch (e) { return null; }
}

const oneLine = s => String(s).replace(/[\t\r\n]+/g, ' ').trim();

// The agent's script: the first scripts/*.R or *.py other than the legacy input of task 04.
function findScript(repo) {
  let files = [];
  try { files = fs.readdirSync(path.join(repo, 'scripts')); } catch (e) {}
  const found = files.filter(f => /\.(R|py)$/i.test(f) && !/^00_legacy\./.test(f)).sort()[0];
  if (found) return path.join('scripts', found);
  // ClawBio-style tasks (01sk): the script lives in skills/<name>/.
  let skills = [];
  try { skills = fs.readdirSync(path.join(repo, 'skills')).sort(); } catch (e) {}
  for (const skill of skills) {
    let inSkill = [];
    try { inSkill = fs.readdirSync(path.join(repo, 'skills', skill)); } catch (e) {}
    const script = inSkill.filter(f => /\.(R|py)$/i.test(f)).sort()[0];
    if (script) return path.join('skills', skill, script);
  }
  return null;
}

const isSkill = script => script.split(path.sep)[0] === 'skills';

const isPython = script => /\.py$/i.test(script);
// uv may have to fetch the interpreter and packages on its first call.
const PY_CMD = ['uv', ['run', '--python', '3.12', '--with', 'pandas,scipy,matplotlib', 'python']];

// `claude --output-format json` prints one result object; tolerate an array of events.
function claudeStats(file) {
  let out = readJson(file);
  if (Array.isArray(out)) out = out.filter(e => e && e.type === 'result').pop();
  const num = v => (typeof v === 'number' ? v : null);
  return {
    cost_usd: num(out && out.total_cost_usd),
    duration_s: num(out && out.duration_ms) === null ? null : out.duration_ms / 1000,
    turns: num(out && out.num_turns),
  };
}

// Reruns the script on a copy of the repo without results/, so leftovers cannot help.
function rerun(repo, script) {
  const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'lucid-rerun-'));
  try {
    fs.cpSync(repo, tmp, { recursive: true, filter: src => src !== path.join(repo, 'results') });
    const [cmd, args, timeout] = isPython(script) ? [...PY_CMD, PY_RERUN_MS] : ['Rscript', [], RERUN_MS];
    // A skill takes the ClawBio contract flags; a task script takes the run dir only.
    const out = path.join('results', 'rerun');
    const scriptArgs = isSkill(script) ? [script, '--input', 'data', '--output', out] : [script, out];
    const r = run(cmd, [...args, ...scriptArgs], { cwd: tmp, timeout });
    return r.status === 0;
  } finally {
    fs.rmSync(tmp, { recursive: true, force: true });
  }
}

function gate(repo, taskDir) {
  const checkR = path.join(taskDir || '', 'check.R');
  if (!taskDir || !fs.existsSync(checkR)) return { gate: null, gate_detail: `no check.R in ${taskDir}` };
  const r = run('Rscript', [checkR, repo, taskDir], { cwd: taskDir, timeout: GATE_MS });
  const parsed = (r.stdout || '').split('\n').reverse().map(l => {
    try { return JSON.parse(l); } catch (e) { return null; }
  }).find(j => j && j.gate);
  if (r.status !== 0 || !parsed) {
    return { gate: null, gate_detail: `check.R error (exit ${r.status}): ${String(r.stderr || '').slice(-300)}` };
  }
  return { gate: parsed.gate, gate_detail: parsed.detail || '' };
}

// ruff check exits 1 when it finds something, so the status says nothing; only a
// missing binary (r.error) gives NA.
function ruffCount(repo, args) {
  const r = run('ruff', ['check', '--no-cache', '--output-format', 'json', ...args], { cwd: repo });
  if (r.error) return null;
  try { return JSON.parse(r.stdout).length; } catch (e) { return null; }
}

// run.sh copies air.toml and ruff.toml in before the run, so the agent formats
// to the same widths; copying again here keeps scoring on the template configs.
function style(repo, script) {
  if (isPython(script)) {
    fs.copyFileSync(path.join(TEMPLATES, 'ruff.toml'), path.join(repo, 'ruff.toml'));
    const fmt = run('ruff', ['format', '--check', '--no-cache', script], { cwd: repo });
    return {
      air_ok: fmt.error ? null : fmt.status === 0,
      lint_count: ruffCount(repo, [script]),
      commented_code_lines: ruffCount(repo, ['--select', 'ERA001', script]),
    };
  }
  fs.copyFileSync(path.join(TEMPLATES, 'air.toml'), path.join(repo, 'air.toml'));
  fs.copyFileSync(path.join(TEMPLATES, '.lintr'), path.join(repo, '.lintr'));
  const air = run('air', ['format', '--check', script], { cwd: repo });
  const lint = run('Rscript', ['-e', LINT_R, script], { cwd: repo });
  const counts = lint.status === 0 ? (lint.stdout.trim().match(/^(\d+) (\d+)$/) || []) : [];
  return {
    air_ok: air.error ? null : air.status === 0,
    lint_count: counts[1] === undefined ? null : Number(counts[1]),
    commented_code_lines: counts[2] === undefined ? null : Number(counts[2]),
  };
}

function textMeasures(code, script) {
  const head = code.split('\n').slice(0, 20).join('\n');
  // The plot rules only apply to a script that plots; elsewhere they are NA, not a free pass.
  // `import matplotlib`, not bare `matplotlib`: the prompt's uv command names it in every header.
  const plots = /ggplot\(|(import|from)\s+matplotlib|pyplot/.test(code);
  // Comment lines are dropped so a guard mentioned in a comment does not count.
  const codeLines = code.split('\n').filter(l => !/^\s*#/.test(l)).join('\n');
  const joinGuard = isPython(script) ? /\bvalidate\s*=(?!=)/g : /\brelationship\s*=(?!=)/g;
  const common = {
    ...codeMeasures(code, isPython(script)),
    join_guard_count: (codeLines.match(joinGuard) || []).length,
    header_complete: HEADER_WORDS.every(w => new RegExp(`\\b${w}\\b`, 'i').test(head)),
    rule_out_dir: /\bout_dir\s*(<-|=)/.test(code),
  };
  // The pipe, theme and warn rules are R-only; a Python script gets NA there.
  if (isPython(script)) {
    return {
      ...common,
      check_count: (code.match(/\braise\b|\bsys\.exit\(|\bassert\b/g) || []).length,
      uptake_pipe: null,
      uptake_no_title: plots ? !/\.set_title\(|plt\.title\(|suptitle\(/.test(code) : null,
      uptake_theme: null,
      native_pipe: null,
      rule_warn1: null,
    };
  }
  return {
    ...common,
    check_count: (code.match(/\bstop\(|\bstopifnot\(|\bvalidate_that\(|assert/g) || []).length,
    uptake_pipe: code.includes('%>%') && !code.includes('|>'),
    uptake_no_title: plots ? !/ggtitle\(|labs\([^)]*\btitle\s*=/.test(code) : null,
    uptake_theme: plots ? /theme_bw\(\s*base_size\s*=\s*8\s*\)/.test(code) : null,
    // House rule: |> and never %>%. A script with no pipe at all is NA, not a pass.
    native_pipe: code.includes('|>') || code.includes('%>%') ? code.includes('|>') && !code.includes('%>%') : null,
    rule_warn1: /options\(\s*warn\s*=\s*1\s*\)/.test(code),
  };
}

function scoreRun(runDir) {
  const meta = readJson(path.join(runDir, 'meta.json')) || {};
  const repo = path.join(runDir, 'repo');
  const taskDir = meta.task_dir || path.join(EVAL_DIR, 'tasks', String(meta.task));
  const row = {
    run_id: `${path.basename(path.dirname(runDir))}/${path.basename(runDir)}`,
    task: meta.task, arm: meta.arm, rep: meta.rep, exit_code: meta.exit_code,
    ...claudeStats(path.join(runDir, 'claude.json')),
    ...gate(repo, taskDir),
  };
  const script = findScript(repo);
  row.script_found = Boolean(script);
  if (!script) return row;
  row.rerun_ok = rerun(repo, script);
  Object.assign(row, style(repo, script), textMeasures(fs.readFileSync(path.join(repo, script), 'utf8'), script));
  return row;
}

function toTsv(row) {
  return COLUMNS.map(c => {
    const v = row[c];
    if (v === null || v === undefined) return 'NA';
    return oneLine(v);
  }).join('\t');
}

function main() {
  const tsDir = process.argv[2];
  if (!tsDir || !fs.existsSync(tsDir)) {
    console.error('usage: node eval/score.js eval/results/<timestamp>');
    process.exit(1);
  }
  const out = path.join(path.dirname(path.resolve(tsDir)), 'results.tsv');
  const done = new Set();
  if (fs.existsSync(out)) {
    for (const line of fs.readFileSync(out, 'utf8').split('\n').slice(1)) done.add(line.split('\t')[0]);
    // A file from before a column was added gets the new header; old rows stay short.
    const [head, ...rest] = fs.readFileSync(out, 'utf8').split('\n');
    if (head !== COLUMNS.join('\t') && COLUMNS.join('\t').startsWith(head + '\t')) {
      fs.writeFileSync(out, [COLUMNS.join('\t'), ...rest].join('\n'));
    }
  } else {
    fs.writeFileSync(out, COLUMNS.join('\t') + '\n');
  }
  const runs = fs.readdirSync(tsDir).map(d => path.join(path.resolve(tsDir), d))
    .filter(d => fs.existsSync(path.join(d, 'meta.json'))).sort();
  let added = 0;
  for (const runDir of runs) {
    const id = `${path.basename(path.dirname(runDir))}/${path.basename(runDir)}`;
    if (done.has(id)) continue;
    fs.appendFileSync(out, toTsv(scoreRun(runDir)) + '\n');
    added++;
  }
  console.log(`score.js: ${added} rows added to ${out} (${runs.length - added} already scored)`);
}

if (require.main === module) main();

module.exports = { COLUMNS, findScript, scoreRun, textMeasures, toTsv, toolEnv };
