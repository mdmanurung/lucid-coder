const test = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { collect } = require('../collect-comments');

test('collects own-line and trailing #$ comments with original line numbers and scores', () => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'lucid-collect-'));
  fs.writeFileSync(path.join(dir, 'P01.R'), [
    '#$ q1 one read (1-5): 4',
    '#$ q2 find decisions (1-5):',
    '#$ q3 change checks (yes/no): yes',
    '#$ q3 which: drop the column check',
    '#$ ---- your #$ comments below ----',
    '#$ q2=3',
    'library(dplyr)',
    '#$ why this cutoff?',
    'cutoff <- 0.05',
    'x <- 1 #$ vague name',
    '# an ordinary comment',
    '#$ end note',
  ].join('\n'));
  fs.writeFileSync(path.join(dir, 'P02.R'), 'a <- 1\n');
  const { rows, scores } = collect(dir);
  assert.deepStrictEqual(rows.map(r => [r.id, r.line, r.text, r.code]), [
    ['P01', 1, 'q2=3', 'library(dplyr)'],
    ['P01', 2, 'why this cutoff?', 'cutoff <- 0.05'],
    ['P01', 3, 'vague name', 'x <- 1'],
    ['P01', 4, 'end note', ''],
  ]);
  assert.deepStrictEqual(scores.P01, { q1: '4', q2: '3', q3: 'yes', q3_which: 'drop the column check' });
  assert.deepStrictEqual(scores.P02, {});
});

test('make-packet --inline header round-trips through the collector', () => {
  const { INLINE_HEADER } = require('../make-packet');
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'lucid-inline-'));
  const filled = INLINE_HEADER.replace(/(\(1-5; can you say[^:]*:) $/m, '$1 5').replace(/(\(yes\/no[^:]*:) $/m, '$1 no');
  fs.writeFileSync(path.join(dir, 'P01.R'), filled + 'a <- 1 #$ fine\n');
  const { rows, scores } = collect(dir);
  assert.deepStrictEqual(scores.P01, { q1: '5', q3: 'no' });
  assert.deepStrictEqual(rows.map(r => [r.line, r.text, r.code]), [[1, 'fine', 'a <- 1']]);
});

test('skill packets: report and SKILL.md comments, #$P procedure tags and q4 fields', () => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'lucid-skill-'));
  fs.writeFileSync(path.join(dir, 'P01.R'), [
    '#$ q4 procedure sound (yes/no; design, replication unit, method, filters, contrast): no',
    '#$ q4 which (one sentence): tests cells, not donors',
    'fit <- lm(y ~ g) #$P cell-level model',
  ].join('\n'));
  fs.writeFileSync(path.join(dir, 'P01.report.md'), '# Report\n#$P why drop genes below 10 counts?\nFiltered 12 genes.\n');
  fs.writeFileSync(path.join(dir, 'P01.SKILL.md'), '---\nname: x\n---\n');
  const { rows, scores } = collect(dir);
  assert.deepStrictEqual(scores.P01, { q4: 'no', q4_which: 'tests cells, not donors' });
  assert.deepStrictEqual(rows.map(r => [r.part, r.kind, r.line, r.text, r.code]), [
    ['R', 'procedure', 1, 'cell-level model', 'fit <- lm(y ~ g)'],
    ['report.md', 'procedure', 2, 'why drop genes below 10 counts?', 'Filtered 12 genes.'],
  ]);
});
