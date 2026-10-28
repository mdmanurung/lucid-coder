#!/usr/bin/env node
// Builds a side-by-side review page: one panel per analysis step, the four arms as
// columns, shuffled per panel so position never names an arm. Optional human-corpus
// panels (four repos per step) follow. The reader picks a column and says why.
// The arm key goes to <out>/side-by-side-key.tsv; eval/results/ is gitignored.
// Usage: node eval/side-by-side.js <results-dir> <out-dir> [--rep 1] [--arms A,D] [--tasks 02,04] [--whole] [--human panels.json]
// --whole compares whole scripts instead of banner sections.
const crypto = require('crypto');
const fs = require('fs');
const path = require('path');
const { findScript } = require('./score');

// Banner text -> section; the first match wins, so "filter, normalise, fit" is fit.
const SECTIONS = [
  ['fit', 'Fit the model', /model|fit|test|edger|limma|normali[sz]/i],
  ['figure', 'Figures in the script', /figure|plot/i],
  ['write', 'Write outputs', /write|output|reproduc|session|json/i],
  ['report', 'Report text in the script', /report/i],
  ['load', 'Load inputs', /load|read|input/i],
  ['prep', 'Filter and aggregate', /filter|aggregat|pseudobulk|sum|donor|cell/i],
  ['setup', 'Header, parameters, arguments', /./],
];
const BANNER = /^#.*(-{4,}|={4,})\s*$/;

function splitScript(code) {
  const parts = {};
  let key = 'setup';
  for (const line of code.split('\n')) {
    if (BANNER.test(line)) key = SECTIONS.find(([, , re]) => re.test(line.replace(/[-=#]/g, ' ')))[0];
    (parts[key] ??= []).push(line);
  }
  return Object.fromEntries(Object.entries(parts).map(([k, v]) => [k, v.join('\n').trim()]));
}

// A notebook shown as text: markdown cells as quoted blocks, code cells as they are.
function notebookText(file) {
  const nb = JSON.parse(fs.readFileSync(file, 'utf8'));
  return nb.cells.map(cell => {
    const src = [].concat(cell.source).join('');
    const tags = cell.metadata?.tags?.length ? ` [${cell.metadata.tags.join(', ')}]` : '';
    return cell.cell_type === 'markdown'
      ? `## ── markdown ──\n${src.split('\n').map(l => '#> ' + l).join('\n')}`
      : `## ── code${tags} ──\n${src}`;
  }).join('\n\n');
}

function markdownSection(file, heading) {
  const text = fs.readFileSync(file, 'utf8');
  const m = text.match(new RegExp(`^## +${heading}[^\\n]*\\n([\\s\\S]*?)(?=^## |(?![\\s\\S]))`, 'm'));
  return m ? `## ${heading}\n${m[1].trim()}` : '';
}

const shuffle = xs => xs.map(x => [crypto.randomInt(1e9), x]).sort((a, b) => a[0] - b[0]).map(([, x]) => x);
const LABELS = ['W', 'X', 'Y', 'Z'];

// One task's runs, one per arm, as panels. Several tasks in one folder get task-prefixed ids.
function armPanels(resultsDir, rep, arms, onlyTasks, whole) {
  const runs = fs.readdirSync(resultsDir)
    .filter(d => new RegExp(`_[${arms.join('')}]_r${rep}$`).test(d))
    .filter(d => !onlyTasks || onlyTasks.some(t => d.startsWith(t + '_')))
    .sort();
  const tasks = [...new Set(runs.map(run => run.replace(/_[A-D]_r\d+$/, '')))];
  const key = [];
  const panels = [];
  for (const task of tasks) {
    const taskRuns = runs.filter(run => run.startsWith(task + '_'));
    if (taskRuns.length !== arms.length) throw new Error(`${task}: expected arms ${arms} at rep ${rep}, found ${taskRuns.length}`);
    const pieces = Object.fromEntries(taskRuns.map(run => {
      const repo = path.join(resultsDir, run, 'repo');
      const script = path.join(repo, findScript(repo));
      const dir = path.dirname(script);
      const nb = fs.readdirSync(dir).find(f => f.endsWith('.ipynb'));
      const skillMd = path.join(dir, 'SKILL.md');
      const code = fs.readFileSync(script, 'utf8');
      return [run, {
        ...splitScript(code),
        whole: code.trim(),
        notebook: nb ? notebookText(path.join(dir, nb)) : '',
        workflow: fs.existsSync(skillMd) ? markdownSection(skillMd, 'Workflow') : '',
      }];
    }));
    // A script without banners puts everything in setup, so its other columns would be empty:
    // compare that task as whole scripts instead.
    const bannerless = taskRuns.some(run => Object.keys(pieces[run]).filter(k => k !== 'setup' && SECTIONS.some(([id]) => id === k)).length === 0);
    // Whole-script mode still shows the notebook and SKILL.md workflow as their own panels.
    const extras = [['notebook', 'Report notebook', 'r'], ['workflow', 'SKILL.md: Workflow section', 'markdown']];
    const steps = (whole || bannerless) ? [['whole', 'Whole script', 'r'], ...extras] : [
      // Script order, not the keyword priority of SECTIONS.
      ...['setup', 'load', 'prep', 'fit', 'figure', 'write', 'report'].map(id => [id, SECTIONS.find(s => s[0] === id)[1], 'r']),
      ['notebook', 'Report notebook', 'r'],
      ['workflow', 'SKILL.md: Workflow section', 'markdown'],
    ];
    const prefix = tasks.length > 1 ? task.replace(/_.*/, '') : 'arm';
    const language = taskRuns.some(run => /py_/.test(run)) ? 'python' : 'r';
    // A step only some versions have would compare code against an empty column; skip it.
    for (const [id, title, lang] of steps.filter(([id]) => taskRuns.every(run => pieces[run][id]))) {
      const columns = shuffle(taskRuns).map((run, i) => {
        key.push([`${prefix}_${id}`, LABELS[i], run, run.match(/_([A-D])_r\d+$/)[1]]);
        return { label: LABELS[i], code: pieces[run][id] || '(nothing in this step)' };
      });
      const taskTitle = tasks.length > 1 ? `${task.replace(/_/, ' ').replace(/-/g, ' ')}: ` : '';
      panels.push({ id: `${prefix}_${id}`, group: 'arms', title: taskTitle + title, language: lang === 'r' ? language : lang, columns });
    }
  }
  return { panels, key };
}

function humanPanels(file) {
  return JSON.parse(fs.readFileSync(file, 'utf8')).panels.map(p => ({
    id: `human_${p.id}`,
    group: 'human',
    title: p.title,
    language: p.language,
    columns: shuffle(p.snippets).map((s, i) => ({
      label: LABELS[i],
      code: s.code,
      // Shown only after the reader saves a pick.
      source: `${s.repo}: ${s.path}${s.cell != null ? ` cell ${s.cell}` : `:${s.start_line}-${s.end_line}`}`,
    })),
  }));
}

function main() {
  const argv = process.argv.slice(2);
  const opt = name => (argv.includes(name) ? argv[argv.indexOf(name) + 1] : null);
  const [resultsDir, outDir] = argv;
  if (!resultsDir || !outDir) {
    console.error('usage: node eval/side-by-side.js <results-dir> <out-dir> [--rep 1] [--arms A,D] [--tasks 02,04] [--whole] [--notes notes.json] [--human panels.json]');
    process.exit(1);
  }
  // --arms A,D gives a two-column head-to-head; columns are still shuffled per panel.
  const { panels, key } = armPanels(resultsDir, opt('--rep') || '1', (opt('--arms') || 'A,B,C,D').split(','), opt('--tasks')?.split(','), argv.includes('--whole'));
  const human = opt('--human') ? humanPanels(opt('--human')) : [];
  // --notes file.json: {panelId: text} shown under a panel title, e.g. a scope statement.
  const notes = opt('--notes') ? JSON.parse(fs.readFileSync(opt('--notes'), 'utf8')) : {};
  for (const p of panels) if (notes[p.id]) p.note = notes[p.id];
  fs.mkdirSync(outDir, { recursive: true });
  const template = fs.readFileSync(path.join(__dirname, 'side-by-side.html'), 'utf8');
  const data = JSON.stringify({ source: path.basename(resultsDir), panels: [...panels, ...human] }).replace(/</g, '\\u003c');
  fs.writeFileSync(path.join(outDir, 'side-by-side.html'), template.replace('/*DATA*/null', () => data));
  fs.writeFileSync(path.join(outDir, 'side-by-side-key.tsv'),
    ['panel\tcolumn\trun_id\tarm', ...key.map(r => r.join('\t'))].join('\n') + '\n');
  console.log(`side-by-side.js: ${panels.length} arm panels, ${human.length} human panels -> ${outDir}`);
}

if (require.main === module) main();

module.exports = { splitScript, notebookText, markdownSection };
