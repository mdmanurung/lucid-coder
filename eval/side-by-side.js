#!/usr/bin/env node
// Builds a side-by-side review page: one panel per analysis step, the four arms as
// columns, shuffled per panel so position never names an arm. Optional human-corpus
// panels (four repos per step) follow. The reader picks a column and says why.
// The arm key goes to <out>/side-by-side-key.tsv; eval/results/ is gitignored.
// Usage: node eval/side-by-side.js <results-dir> <out-dir> [--rep 1] [--human panels.json]
const crypto = require('crypto');
const fs = require('fs');
const path = require('path');

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

function armPanels(resultsDir, rep) {
  const runs = fs.readdirSync(resultsDir).filter(d => new RegExp(`_[A-D]_r${rep}$`).test(d)).sort();
  if (runs.length !== 4) throw new Error(`expected 4 arms at rep ${rep} in ${resultsDir}, found ${runs.length}`);
  const skill = run => path.join(resultsDir, run, 'repo', 'skills', 'pseudobulk-de');
  const pieces = Object.fromEntries(runs.map(run => {
    const dir = skill(run);
    const script = fs.readdirSync(dir).find(f => /\.(R|py)$/.test(f));
    const nb = path.join(dir, 'report.ipynb');
    return [run, {
      ...splitScript(fs.readFileSync(path.join(dir, script), 'utf8')),
      notebook: fs.existsSync(nb) ? notebookText(nb) : '',
      workflow: markdownSection(path.join(dir, 'SKILL.md'), 'Workflow'),
    }];
  }));
  const steps = [
    // Script order, not the keyword priority of SECTIONS.
    ...['setup', 'load', 'prep', 'fit', 'figure', 'write', 'report'].map(id => [id, SECTIONS.find(s => s[0] === id)[1], 'r']),
    ['notebook', 'Report notebook (report.ipynb)', 'r'],
    ['workflow', 'SKILL.md: Workflow section', 'markdown'],
  ];
  const key = [];
  const panels = steps
    .filter(([id]) => runs.some(run => pieces[run][id]))
    .map(([id, title, language]) => {
      const columns = shuffle(runs).map((run, i) => {
        key.push([`arm_${id}`, LABELS[i], run, run.match(/_([A-D])_r\d+$/)[1]]);
        return { label: LABELS[i], code: pieces[run][id] || '(nothing in this step)' };
      });
      return { id: `arm_${id}`, group: 'arms', title, language, columns };
    });
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
    console.error('usage: node eval/side-by-side.js <results-dir> <out-dir> [--rep 1] [--human panels.json]');
    process.exit(1);
  }
  const { panels, key } = armPanels(resultsDir, opt('--rep') || '1');
  const human = opt('--human') ? humanPanels(opt('--human')) : [];
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
