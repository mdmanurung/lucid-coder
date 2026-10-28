#!/usr/bin/env node
// One-page report by evidence level (PLAN.md 8). Measures stay separate:
// there is never a combined score.
// Usage: node eval/report.js eval/results/<timestamp> [--human sheet.tsv]... [--step0 notes.txt]
//        node eval/report.js --judges eval/results/<timestamp>   (parse judge-raw/ into judges.tsv)
// Writes report.md next to results.tsv.
const fs = require('fs');
const path = require('path');
const { parseArgs } = require('util');

const ARM_NAMES = { A: 'baseline', B: 'ponytail', C: 'lucid-coder', D: 'lucid-coder + 3 rules' };
const CHANNEL_NOTE = 'Arms delivered by --append-system-prompt-file in --safe-mode sessions, because plugin hooks do not run in an isolated session (step 0, 2026-10-26).';
const JUDGE_COLUMNS = ['neutral_id', 'judge', 'q1_one_read', 'q2_find_decisions', 'q3_change_checks', 'q3_which', 'notes'];

function readTsv(file) {
  if (!fs.existsSync(file)) return [];
  const [head, ...lines] = fs.readFileSync(file, 'utf8').split('\n').filter(l => l.trim());
  const cols = head.split('\t');
  return lines.map(l => {
    const cells = l.split('\t');
    return Object.fromEntries(cols.map((c, i) => [c, (cells[i] || '').trim()]));
  });
}

// Weighted Cohen's kappa for two raters over the same items.
// weight: 'quadratic' or 'linear' disagreement weights for ordinal scales, 'none' for nominal.
// Returns null when it is undefined (no items, or both raters used one category).
function kappa(a, b, categories, weight = 'quadratic') {
  const k = categories.length;
  const n = a.length;
  if (!n || n !== b.length) return null;
  const idx = v => categories.indexOf(v);
  const w = (i, j) => {
    if (weight === 'none') return i === j ? 0 : 1;
    const d = Math.abs(i - j) / (k - 1);
    return weight === 'linear' ? d : d * d;
  };
  const pa = Array(k).fill(0);
  const pb = Array(k).fill(0);
  let observed = 0;
  for (let t = 0; t < n; t++) {
    const i = idx(a[t]);
    const j = idx(b[t]);
    if (i < 0 || j < 0) throw new Error(`value outside categories: ${a[t]}, ${b[t]}`);
    pa[i] += 1 / n;
    pb[j] += 1 / n;
    observed += w(i, j) / n;
  }
  let expected = 0;
  for (let i = 0; i < k; i++) for (let j = 0; j < k; j++) expected += w(i, j) * pa[i] * pb[j];
  return expected === 0 ? null : 1 - observed / expected;
}

// Outer layer: the `claude --output-format json` envelope; inner: the judge's own JSON,
// possibly fenced. Returns the answer fields or { error }.
function parseJudgeReply(text) {
  let outer;
  try { outer = JSON.parse(text); } catch (e) { return { error: 'reply is not JSON' }; }
  if (Array.isArray(outer)) outer = outer.filter(e => e && e.type === 'result').pop();
  const body = outer && typeof outer.result === 'string' ? outer.result : '';
  const start = body.indexOf('{');
  const end = body.lastIndexOf('}');
  let ans;
  try { ans = JSON.parse(body.slice(start, end + 1)); } catch (e) { return { error: 'no JSON object in result' }; }
  const score = v => (Number.isInteger(v) && v >= 1 && v <= 5 ? v : null);
  const q3 = String(ans.q3_change_checks || '').toLowerCase();
  if (!score(ans.q1_one_read) || !score(ans.q2_find_decisions) || !['yes', 'no'].includes(q3)) {
    return { error: 'answer fields missing or out of range' };
  }
  return {
    q1_one_read: ans.q1_one_read, q2_find_decisions: ans.q2_find_decisions,
    q3_change_checks: q3, q3_which: ans.q3_which || '', notes: ans.notes || '',
  };
}

const clean = v => String(v == null ? '' : v).replace(/[\t\r\n]+/g, ' ').trim();

// judge-raw/<id>.<model>.json -> judges.tsv; an unparseable reply becomes an NA row.
function writeJudges(tsDir) {
  const raw = path.join(tsDir, 'judge-raw');
  const lines = [JUDGE_COLUMNS.join('\t')];
  for (const f of fs.readdirSync(raw).filter(f => f.endsWith('.json')).sort()) {
    const [, id, judge] = f.match(/^([^.]+)\.(.+)\.json$/) || [];
    if (!id) continue;
    const ans = parseJudgeReply(fs.readFileSync(path.join(raw, f), 'utf8'));
    const row = ans.error
      ? { neutral_id: id, judge, q1_one_read: 'NA', q2_find_decisions: 'NA', q3_change_checks: 'NA', q3_which: '', notes: `unparseable: ${ans.error}` }
      : { neutral_id: id, judge, ...ans };
    lines.push(JUDGE_COLUMNS.map(c => clean(row[c])).join('\t'));
  }
  fs.writeFileSync(path.join(tsDir, 'judges.tsv'), lines.join('\n') + '\n');
  return lines.length - 1;
}

const isNA = v => v === undefined || v === '' || v === 'NA';

function median(values) {
  const v = values.filter(x => !isNA(x)).map(Number).sort((x, y) => x - y);
  if (!v.length) return 'NA';
  const m = v.length % 2 ? v[(v.length - 1) / 2] : (v[v.length / 2 - 1] + v[v.length / 2]) / 2;
  return String(Math.round(m * 100) / 100);
}

function rate(values, hit = 'true') {
  const v = values.filter(x => !isNA(x));
  return `${v.filter(x => x === hit).length}/${v.length}`;
}

const fmtKappa = k => (k === null ? 'undefined' : k.toFixed(2));

function table(header, rows) {
  return [`| ${header.join(' | ')} |`, `|${header.map(() => '---').join('|')}|`]
    .concat(rows.map(r => `| ${r.join(' | ')} |`)).join('\n');
}

const L3_MEASURES = [
  ['gate pass', r => rate(r.map(x => x.gate), 'pass')],
  ['script found', r => rate(r.map(x => x.script_found))],
  ['rerun ok', r => rate(r.map(x => x.rerun_ok))],
  ['format check pass', r => rate(r.map(x => x.air_ok))],
  ['lint count, median', r => median(r.map(x => x.lint_count))],
  ['header complete', r => rate(r.map(x => x.header_complete))],
  ['check count, median', r => median(r.map(x => x.check_count))],
  ['join guards, median', r => median(r.map(x => x.join_guard_count))],
  ['commented-code lines, median', r => median(r.map(x => x.commented_code_lines))],
  ['code lines, median (rubric R1)', r => median(r.map(x => x.code_lines))],
  ['joins without keys, median (S1)', r => median(r.map(x => x.joins_unkeyed))],
  ['blanket warning suppression, median (S5)', r => median(r.map(x => x.warning_suppress))],
  ['longest pipe, median steps (R2)', r => median(r.map(x => x.longest_pipe))],
  ['deepest nesting, median (R3)', r => median(r.map(x => x.max_nesting))],
  ['uninformative names, median (R4)', r => median(r.map(x => x.vague_names))],
  ['wrappers used once, median (R10)', r => median(r.map(x => x.wrappers_once))],
  ['named functions used < 3 times, median (R11)', r => median(r.map(x => x.functions_under3))],
  ['session info written (P2)', r => rate(r.map(x => x.session_info))],
  ['late library/import, median (P4)', r => median(r.map(x => x.late_loads))],
  ['base-R data-frame calls, median (H5)', r => median(r.map(x => x.base_df_calls))],
  ['UPPER_SNAKE R parameters, median (H6, house style)', r => median(r.map(x => x.upper_constants))],
  ['uptake: %>% only', r => rate(r.map(x => x.uptake_pipe))],
  ['uptake: no plot title', r => rate(r.map(x => x.uptake_no_title))],
  ['uptake: theme_pubr(base_size = 8)', r => rate(r.map(x => x.uptake_theme))],
  ['house: |> only (pipe used)', r => rate(r.map(x => x.native_pipe))],
  ['D rule: out_dir variable', r => rate(r.map(x => x.rule_out_dir))],
  ['D rule: options(warn = 1)', r => rate(r.map(x => x.rule_warn1))],
];

// Ratings (judges or readers) joined to arms through a key file; one table per rater.
function ratingSection(ratings, keyFile, armOf, arms, raterCol) {
  const runOf = Object.fromEntries(readTsv(keyFile).map(k => [k.neutral_id, k.run_id]));
  const raters = [...new Set(ratings.map(r => r[raterCol]))];
  const out = [];
  for (const rater of raters) {
    const mine = ratings.filter(r => r[raterCol] === rater);
    const byArm = a => mine.filter(r => armOf[runOf[r.neutral_id]] === a);
    out.push(`**${rater}**`, '', table(['measure', ...arms.map(a => `${a} ${ARM_NAMES[a] || ''}`.trim())], [
      ['n scripts', ...arms.map(a => String(byArm(a).length))],
      ['q1 one read, median (1-5)', ...arms.map(a => median(byArm(a).map(r => r.q1_one_read)))],
      ['q2 find decisions, median (1-5)', ...arms.map(a => median(byArm(a).map(r => r.q2_find_decisions)))],
      ['q3 would change a check (yes)', ...arms.map(a => rate(byArm(a).map(r => r.q3_change_checks), 'yes'))],
    ]), '');
  }
  if (raters.length < 2) {
    out.push('Agreement: needs two raters.');
    return out;
  }
  const [r1, r2] = raters;
  const second = Object.fromEntries(ratings.filter(r => r[raterCol] === r2).map(r => [r.neutral_id, r]));
  const pairs = q => ratings.filter(r => r[raterCol] === r1 && second[r.neutral_id] &&
    !isNA(r[q]) && !isNA(second[r.neutral_id][q])).map(r => [r[q], second[r.neutral_id][q]]);
  const agree = (q, cats, weight) => {
    const p = pairs(q);
    return `${fmtKappa(kappa(p.map(x => x[0]), p.map(x => x[1]), cats, weight))} (n = ${p.length})`;
  };
  const scale = ['1', '2', '3', '4', '5'];
  out.push(`Agreement, ${r1} vs ${r2} (weighted Cohen's kappa; quadratic weights for q1 and q2, unweighted for q3): ` +
    `q1 ${agree('q1_one_read', scale, 'quadratic')}, q2 ${agree('q2_find_decisions', scale, 'quadratic')}, ` +
    `q3 ${agree('q3_change_checks', ['no', 'yes'], 'none')}.`);
  if (raters.length > 2) out.push(`Only the first two raters enter the kappa; ${raters.length} rated.`);
  return out;
}

function buildReport({ tsDir, rows, human = [], step0 = '' }) {
  const ts = path.basename(tsDir);
  const runDirs = fs.readdirSync(tsDir).filter(d => fs.existsSync(path.join(tsDir, d, 'meta.json')));
  const metas = runDirs.map(d => JSON.parse(fs.readFileSync(path.join(tsDir, d, 'meta.json'), 'utf8')));
  const arms = [...new Set(rows.map(r => r.arm))].sort();
  const tasks = [...new Set(rows.map(r => r.task))].sort();
  const armOf = Object.fromEntries(rows.map(r => [r.run_id, r.arm]));
  const cell = (t, a) => rows.filter(r => r.task === t && r.arm === a);
  const cellNs = tasks.flatMap(t => arms.map(a => cell(t, a).length));
  const nNote = `n per task × arm cell = ${cellNs.length ? [...new Set(cellNs)].join(' or ') : 0}; ` +
    `n per arm = ${arms.map(a => `${a}: ${rows.filter(r => r.arm === a).length}`).join(', ')}.`;
  const armHeader = arms.map(a => `${a} ${ARM_NAMES[a] || ''}`.trim());
  const channels = [...new Set(metas.map(m => m.arm_channel || 'unrecorded'))].join(', ');
  const cost = rows.reduce((s, r) => s + (isNA(r.cost_usd) ? 0 : Number(r.cost_usd)), 0);

  const out = [
    `# lucid-coder eval: ${ts}`,
    '',
    `Model ${[...new Set(metas.map(m => m.model))].join(', ')}, budget $${[...new Set(metas.map(m => m.budget))].join(', ')} per run. ` +
    'Each measure is reported on its own; there is no combined score, and a failed gate is never offset by style.',
    '',
    '## L1: the harness ran',
    '',
    `${runDirs.length} runs found, ${rows.length} scored. ` +
    `claude exited non-zero in ${metas.filter(m => m.exit_code !== 0).length}. ` +
    `No usable claude.json (cost NA) in ${rows.filter(r => isNA(r.cost_usd)).length}. ` +
    `No script written in ${rows.filter(r => r.script_found === 'false').length}. ` +
    `Total cost $${cost.toFixed(2)}.`,
    '',
    '## L2: the text reached the session',
    '',
    CHANNEL_NOTE,
    `Channel recorded in meta.json: ${channels}.`,
    step0 ? `\nStep 0 notes:\n\n${step0.trim()}` : '',
    '',
    '## L3: the output changed',
    '',
    nNote + ' Rate cells show k/n with NA excluded.',
    '',
    table(['measure', ...armHeader], [
      ['n runs', ...arms.map(a => String(rows.filter(r => r.arm === a).length))],
      ...L3_MEASURES.map(([label, f]) => [label, ...arms.map(a => f(rows.filter(r => r.arm === a)))]),
    ]),
    '',
    `Gate pass by task (${nNote.split(';')[0]}):`,
    '',
    table(['task', ...armHeader], tasks.map(t => [t, ...arms.map(a => rate(cell(t, a).map(r => r.gate), 'pass'))])),
    '',
  ];
  const failed = rows.filter(r => r.gate !== 'pass');
  out.push(failed.length ? 'Gates not passed:\n\n' + failed.map(r => `- ${r.run_id}: ${r.gate} (${r.gate_detail || 'no detail'})`).join('\n') : 'Every gate passed.', '');

  out.push('## L4: humans find it more readable', '');
  const humanRows = human.flatMap((sheet, i) => sheet.rows
    .filter(r => !isNA(r.q1_one_read) || !isNA(r.q2_find_decisions) || !isNA(r.q3_change_checks))
    .map(r => ({ ...r, reader: `reader ${i + 1} (${path.basename(sheet.file)})` })));
  if (!humanRows.length) {
    out.push('Pending: no filled reader sheets yet.', '');
  } else {
    out.push(`Blinded readers; ${nNote}`, '');
    out.push(...ratingSection(humanRows, path.join(tsDir, 'packet-key.tsv'), armOf, arms, 'reader'), '');
  }

  out.push('## LLM judges (proxy, not human evidence)', '');
  const judges = readTsv(path.join(tsDir, 'judges.tsv'));
  if (!judges.length) {
    out.push('Not run.', '');
  } else {
    out.push(`Two blinded model judges; ${nNote}`, '');
    out.push(...ratingSection(judges, path.join(tsDir, 'judge-key.tsv'), armOf, arms, 'judge'), '');
  }
  return out.join('\n').replace(/\n{3,}/g, '\n\n');
}

function main() {
  const { values, positionals } = parseArgs({
    allowPositionals: true,
    options: { judges: { type: 'boolean' }, human: { type: 'string', multiple: true }, step0: { type: 'string' } },
  });
  const tsDir = positionals[0] && path.resolve(positionals[0]);
  if (!tsDir || !fs.existsSync(tsDir)) {
    console.error('usage: node eval/report.js eval/results/<timestamp> [--human sheet.tsv]... [--step0 notes.txt]\n' +
      '       node eval/report.js --judges eval/results/<timestamp>');
    process.exit(1);
  }
  if (values.judges) {
    console.log(`report.js: ${writeJudges(tsDir)} judge rows written to ${path.join(tsDir, 'judges.tsv')}`);
    return;
  }
  const resultsDir = path.dirname(tsDir);
  const rows = readTsv(path.join(resultsDir, 'results.tsv'))
    .filter(r => r.run_id.startsWith(path.basename(tsDir) + '/'));
  const human = (values.human || []).map(file => ({ file, rows: readTsv(file) }));
  const step0 = values.step0 ? fs.readFileSync(values.step0, 'utf8') : '';
  const out = path.join(resultsDir, 'report.md');
  fs.writeFileSync(out, buildReport({ tsDir, rows, human, step0 }));
  console.log(`report.js: wrote ${out}`);
}

if (require.main === module) main();

module.exports = { kappa, parseJudgeReply, writeJudges, buildReport, readTsv };
