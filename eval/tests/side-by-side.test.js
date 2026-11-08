const test = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync } = require('child_process');
const { splitScript, notebookText, armPanels } = require('../side-by-side');

test('banners map to sections; fit wins over filter; text before the first banner is setup', () => {
  const parts = splitScript([
    '# Purpose: x',
    '# ---- Arguments ----',
    'args <- 1',
    '# ---- filter, normalise, fit ----',
    'fit <- 2',
    '# Write outputs ----',
    'write(3)',
  ].join('\n'));
  assert.deepStrictEqual(Object.keys(parts), ['setup', 'fit', 'write']);
  assert.match(parts.setup, /Purpose[\s\S]*args <- 1/);
});

test('page embeds code with $ sequences intact and writes a per-panel key', () => {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), 'sbs-'));
  for (const arm of ['A', 'B', 'C', 'D']) {
    const dir = path.join(root, `01sk_x_${arm}_r1`, 'repo', 'skills', 'pseudobulk-de');
    fs.mkdirSync(dir, { recursive: true });
    fs.writeFileSync(path.join(dir, 'x.R'), `# Load data ----\nx <- df$\`a b\` # ${arm}\n`);
    fs.writeFileSync(path.join(dir, 'SKILL.md'), '## Workflow\n1. run\n## Safety\n');
    fs.writeFileSync(path.join(dir, 'report.ipynb'), JSON.stringify({ cells: [
      { cell_type: 'code', metadata: { tags: ['parameters'] }, source: ['out_dir <- "r"'] },
      { cell_type: 'markdown', metadata: {}, source: '# Report' },
    ] }));
  }
  const out = path.join(root, 'out');
  const r = spawnSync('node', [path.join(__dirname, '..', 'side-by-side.js'), root, out], { encoding: 'utf8' });
  assert.strictEqual(r.status, 0, r.stderr);
  const html = fs.readFileSync(path.join(out, 'side-by-side.html'), 'utf8');
  const data = JSON.parse(html.match(/const DATA = (.*);\n/)[1]);
  const load = data.panels.find(p => p.id === 'arm_load');
  assert.ok(load.columns.every(c => c.code.includes('df$`a b`')));
  assert.match(data.panels.find(p => p.id === 'arm_notebook').columns[0].code, /code \[parameters\][\s\S]*#> # Report/);
  const key = fs.readFileSync(path.join(out, 'side-by-side-key.tsv'), 'utf8').trim().split('\n');
  assert.strictEqual(key.length, 1 + 4 * data.panels.length);
});

test('--arms A,D builds two-column panels per task from scripts/', () => {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), 'sbs-pair-'));
  for (const task of ['02_join', '03_plot']) {
    for (const arm of ['A', 'B', 'C', 'D']) {
      const dir = path.join(root, `${task}_${arm}_r1`, 'repo', 'scripts');
      fs.mkdirSync(dir, { recursive: true });
      fs.writeFileSync(path.join(dir, `${task}.R`), `# Load data ----\nx <- 1 # ${arm}\n`);
    }
  }
  const out = path.join(root, 'out');
  const r = spawnSync('node', [path.join(__dirname, '..', 'side-by-side.js'), root, out, '--arms', 'A,D'], { encoding: 'utf8' });
  assert.strictEqual(r.status, 0, r.stderr);
  const data = JSON.parse(fs.readFileSync(path.join(out, 'side-by-side.html'), 'utf8').match(/const DATA = (.*);\n/)[1]);
  assert.deepStrictEqual(data.panels.map(p => p.id), ['02_load', '03_load']);
  assert.ok(data.panels.every(p => p.columns.length === 2 && p.columns.every(c => /# [AD]$/.test(c.code))));
});

// Three tasks, each with a load and a write panel (even), or load only (odd).
function pairRoot(prefix, tasks) {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), prefix));
  for (const [task, banners] of tasks) {
    for (const arm of ['A', 'D']) {
      const dir = path.join(root, `${task}_${arm}_r1`, 'repo', 'scripts');
      fs.mkdirSync(dir, { recursive: true });
      fs.writeFileSync(path.join(dir, `${task}.R`), banners.map(b => `# ${b} ----\nx <- 1\n`).join(''));
    }
  }
  return root;
}
const aLeftCount = (key, prefix) => key.filter(([id, label, , arm]) => id.startsWith(prefix) && label === 'W' && arm === 'A').length;

test('two arms: arm A leads in half the panels of every task and of the pool', () => {
  const root = pairRoot('sbs-balance-', [
    ['01_a', ['Load data', 'Write outputs']], ['02_b', ['Load data', 'Write outputs']],
    ['03_c', ['Load data']], ['04_d', ['Load data']], ['05_e', ['Load data']],
  ]);
  // Independent per-panel shuffles pass one run with p ~ 0.1, all 20 with p ~ 1e-20.
  for (let i = 0; i < 20; i++) {
    const { panels, key } = armPanels(root, '1', ['A', 'D'], null, false, i);
    assert.strictEqual(panels.length, 7);
    for (const task of ['01', '02']) assert.strictEqual(aLeftCount(key, task), 1, `A leads task ${task} in ${aLeftCount(key, task)} of 2`);
    const pooled = aLeftCount(key, '');
    assert.ok(pooled === 3 || pooled === 4, `A on the left in ${pooled} of 7`);
  }
});

test('the same seed rebuilds the same key and --seed is written to the meta file', () => {
  const root = pairRoot('sbs-seed-', [['01_a', ['Load data', 'Write outputs']], ['02_b', ['Load data', 'Write outputs']], ['03_c', ['Load data']]]);
  const keys = [11, 11, 12, 13, 14].map(seed => JSON.stringify(armPanels(root, '1', ['A', 'D'], null, false, seed).key));
  assert.strictEqual(keys[0], keys[1]);
  assert.ok(new Set(keys).size > 1, 'different seeds give different keys');
  const out = path.join(root, 'out');
  const r = spawnSync('node', [path.join(__dirname, '..', 'side-by-side.js'), root, out, '--arms', 'A,D', '--seed', '7'], { encoding: 'utf8' });
  assert.strictEqual(r.status, 0, r.stderr);
  assert.strictEqual(JSON.parse(fs.readFileSync(path.join(out, 'side-by-side-meta.json'), 'utf8')).seed, 7);
});

test('a script without banners turns that task into one whole-script panel', () => {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), 'sbs-whole-'));
  for (const arm of ['A', 'D']) {
    const dir = path.join(root, `02_join_${arm}_r1`, 'repo', 'scripts');
    fs.mkdirSync(dir, { recursive: true });
    const code = arm === 'A' ? 'library(dplyr)\nx <- read()\n' : '# Load data ----\nx <- read()\n# Write outputs ----\nwrite(x)\n';
    fs.writeFileSync(path.join(dir, '02_join.R'), code);
  }
  const out = path.join(root, 'out');
  const r = spawnSync('node', [path.join(__dirname, '..', 'side-by-side.js'), root, out, '--arms', 'A,D'], { encoding: 'utf8' });
  assert.strictEqual(r.status, 0, r.stderr);
  const data = JSON.parse(fs.readFileSync(path.join(out, 'side-by-side.html'), 'utf8').match(/const DATA = (.*);\n/)[1]);
  assert.deepStrictEqual(data.panels.map(p => p.id), ['arm_whole']);
  assert.ok(data.panels[0].columns.some(c => c.code.includes('# Write outputs ----')));
});
