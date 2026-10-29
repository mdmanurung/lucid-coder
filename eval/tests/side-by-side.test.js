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

test('two arms: each arm is on the left in half the panels, odd count off by one', () => {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), 'sbs-balance-'));
  for (const task of ['01_a', '02_b', '03_c', '04_d', '05_e', '06_f', '07_g']) {
    for (const arm of ['A', 'D']) {
      const dir = path.join(root, `${task}_${arm}_r1`, 'repo', 'scripts');
      fs.mkdirSync(dir, { recursive: true });
      fs.writeFileSync(path.join(dir, `${task}.R`), `# Load data ----\nx <- 1\n`);
    }
  }
  // Independent per-panel shuffles pass one run with p = 70/128, all 20 with p ~ 6e-6.
  for (let i = 0; i < 20; i++) {
    const { panels, key } = armPanels(root, '1', ['A', 'D'], null, false);
    assert.strictEqual(panels.length, 7);
    const aLeft = key.filter(([, label, , arm]) => label === 'W' && arm === 'A').length;
    assert.ok(aLeft === 3 || aLeft === 4, `A on the left in ${aLeft} of 7`);
  }
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
