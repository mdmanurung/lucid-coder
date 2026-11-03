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
