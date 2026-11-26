const test = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { compare } = require('../compare-readers');

test('picks map to arms per panel; JSON and TSV answers; agreement over panels both answered', () => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'readers-'));
  const file = (name, text) => (fs.writeFileSync(path.join(dir, name), text), path.join(dir, name));
  const key = file('key.tsv', 'panel\tcolumn\trun_id\tarm\np1\tW\tt_A_r1\tA\np1\tX\tt_D_r1\tD\np2\tW\tt_D_r1\tD\np2\tX\tt_A_r1\tA\np3\tW\tt_A_r1\tA\np3\tX\tt_D_r1\tD\n');
  const first = file('first.tsv', 'panel\tchoice\np1\tX\np2\tW\np3\tnone\n');
  const second = file('second.json', JSON.stringify([{ panel: 'p1', choice: 'X' }, { panel: 'p2', choice: 'X' }]));
  const { rows, agree, both } = compare(key, first, second);
  assert.deepStrictEqual(rows.map(r => [r.reader1, r.reader2, r.agree]), [['D', 'D', true], ['D', 'A', false], ['none', '', null]]);
  assert.deepStrictEqual([agree, both], [1, 2]);
});

test('a missing panel or choice column is a clear error, not an empty comparison', () => {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'readers-cols-'));
  const file = (name, text) => (fs.writeFileSync(path.join(dir, name), text), path.join(dir, name));
  const key = file('key.tsv', 'panel\tcolumn\trun_id\tarm\np1\tW\tt_A_r1\tA\n');
  const good = file('good.tsv', 'panel\tchoice\np1\tW\n');
  assert.throws(() => compare(key, file('nopanel.tsv', 'step\tchoice\np1\tW\n'), good), /nopanel\.tsv lacks column panel/);
  assert.throws(() => compare(key, good, file('nochoice.json', JSON.stringify([{ panel: 'p1' }]))), /nochoice\.json lacks column choice/);
  assert.throws(() => compare(file('badkey.tsv', 'panel\tcol\tarm\np1\tW\tA\n'), good, good), /badkey\.tsv lacks column column/);
});
