const test = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync } = require('child_process');
const { splitScript, notebookText } = require('../side-by-side');

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
