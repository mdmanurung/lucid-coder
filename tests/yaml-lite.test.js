const test = require('node:test');
const assert = require('node:assert');
const { parse, stringify } = require('../scripts/yaml-lite');

test('parses the supported subset', () => {
  const doc = parse([
    '# header',
    'version: 1',
    'name: cd14  # trailing comment',
    'q: "a: b # not a comment"',
    "s: 'it''s'",
    'date: 2026-10-29',
    'flags: [true, false, null, 0.05, "x, y"]',
    'empty: []',
    'nested:',
    '  deeper:',
    '    key: ask',
    'filters:',
    '  - step: filterByExpr',
    '    cutoff: 10',
    '  - step: other',
    'outputs:',
    '- a.tsv',
    '- b.pdf',
    'none:',
  ].join('\n'));
  assert.deepStrictEqual(doc, {
    version: 1, name: 'cd14', q: 'a: b # not a comment', s: "it's", date: '2026-10-29',
    flags: [true, false, null, 0.05, 'x, y'], empty: [], nested: { deeper: { key: 'ask' } },
    filters: [{ step: 'filterByExpr', cutoff: 10 }, { step: 'other' }],
    outputs: ['a.tsv', 'b.pdf'], none: null,
  });
});

test('rejects anything outside the subset with a line number', () => {
  const cases = {
    'a: 1\nb: &x 2': /^Error: line 2: unsupported/,
    'a: *x': /^Error: line 1: unsupported/,
    'a: |\n  text': /^Error: line 1: unsupported/,
    'a: >+\n  text': /^Error: line 1: unsupported/,
    'a: >\nb: 1': /^Error: line 1: empty folded/,
    'a: {b: 1}': /^Error: line 1: unsupported/,
    'a:\n\tb: 1': /^Error: line 2: tabs/,
    'a: 1\na: 2': /^Error: line 2: duplicate key a/,
    'a: b: c': /^Error: line 1: a plain value cannot hold/,
    'a: one\n  two': /^Error: line 2: unexpected indent/,
    'a: [x, [y]]': /^Error: line 1: nested flow/,
    'a: "open': /^Error: line 1: unclosed double quote/,
    'a: [x, y': /^Error: line 1: unclosed flow list/,
    '---\na: 1': /^Error: line 1: document markers/,
    '  a: 1': /^Error: line 1: the top level/,
    'just text': /^Error: line 1: expected "key: value"/,
    '__proto__: 1': /^Error: line 1: key __proto__/,
    'a:\n  - - x': /^Error: line 2: unsupported/,
    'a: [x,,y]': /^Error: line 1: empty item/,
    'a: [,x]': /^Error: line 1: empty item/,
    'a: [x, y # c]': /^Error: line 1: a comment inside/,
    'a: >\n    x\n  y': /^Error: line 3: folded line is less/,
    'a: >\n  x\n    y': /^Error: line 3: more-indented/,
  };
  for (const [text, err] of Object.entries(cases)) {
    assert.throws(() => parse(text), err, text);
  }
});

test('stringify reads back, quoting what would not', () => {
  const value = { approved: { date: '2026-10-29', by: 'Ann Lee', sha256: '0'.repeat(63) + '5', n: '12', t: 'a: b' } };
  const text = stringify(value);
  assert.match(text, /^approved:\n {2}date: 2026-10-29\n {2}by: Ann Lee\n/);
  assert.deepStrictEqual(parse(text), value);
});

test('stringify rejects what it cannot write and quotes flow items', () => {
  for (const by of ['Ann\nLee', 'Ann\rLee', 'a\u2028b', 'a\x01b']) {
    assert.throws(() => stringify({ by }), /line break or control/, JSON.stringify(by));
  }
  for (const bad of [NaN, undefined, { a: 1 }]) assert.throws(() => stringify({ l: [bad] }));
  const value = { l: ['a,b', 'c]', '[d', '{e}', 'f'], k: '- x' };
  assert.deepStrictEqual(parse(stringify(value)), value);
});

test('comment strip is linear and folded keeps content indent', () => {
  const t = Date.now();
  parse('a: x' + ' '.repeat(200000) + 'y');
  assert.ok(Date.now() - t < 1000);
  assert.deepStrictEqual(parse('a: >-\n    x\n    y'), { a: 'x y' });
});
