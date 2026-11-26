#!/usr/bin/env node
// Maps two readers' side-by-side picks to arms through the page key and counts agreement.
// Answers: the page's "Copy all answers" JSON, or a TSV with panel and choice columns.
// Agreement counts panels both readers answered; "none" (no preference) is its own pick.
// Usage: node eval/compare-readers.js <side-by-side-key.tsv> <reader1> <reader2>
const fs = require('fs');

// A required column that is absent would read as undefined and silently drop every pick.
function requireColumns(head, needed, file) {
  const missing = needed.filter(c => !head.includes(c));
  if (missing.length) throw new Error(`${file} lacks column ${missing.join(', ')}; found ${head.join(', ') || 'none'}`);
}

function readPicks(file) {
  const text = fs.readFileSync(file, 'utf8').trim();
  if (text.startsWith('[')) {
    const rows = JSON.parse(text);
    rows.forEach(r => requireColumns(Object.keys(r), ['panel', 'choice'], file));
    return Object.fromEntries(rows.map(r => [r.panel, r.choice]));
  }
  const [head, ...lines] = text.split('\n').map(l => l.split('\t'));
  requireColumns(head.map(h => h.trim()), ['panel', 'choice'], file);
  return Object.fromEntries(lines.map(c => [c[head.indexOf('panel')], c[head.indexOf('choice')]?.trim()]));
}

function compare(keyFile, file1, file2) {
  const arm = {};
  const [keyHead, ...keyRows] = fs.readFileSync(keyFile, 'utf8').trim().split('\n').map(l => l.split('\t'));
  requireColumns(keyHead, ['panel', 'column', 'arm'], keyFile);
  const at = name => keyHead.indexOf(name);
  for (const row of keyRows) arm[`${row[at('panel')]}\t${row[at('column')]}`] = row[at('arm')];
  const [p1, p2] = [readPicks(file1), readPicks(file2)];
  const toArm = (panel, choice) => (choice === 'none' ? 'none' : arm[`${panel}\t${choice}`] ?? '');
  const panels = [...new Set(Object.keys(arm).map(k => k.split('\t')[0]))];
  const rows = panels.map(panel => {
    const [a1, a2] = [toArm(panel, p1[panel]), toArm(panel, p2[panel])];
    return { panel, reader1: a1, reader2: a2, agree: a1 && a2 ? a1 === a2 : null };
  });
  const both = rows.filter(r => r.agree !== null);
  return { rows, agree: both.filter(r => r.agree).length, both: both.length };
}

if (require.main === module) {
  const [keyFile, file1, file2] = process.argv.slice(2);
  if (!file2) { console.error('usage: node eval/compare-readers.js <key.tsv> <reader1> <reader2>'); process.exit(1); }
  let result;
  try {
    result = compare(keyFile, file1, file2);
  } catch (e) {
    console.error(`compare-readers.js: ${e.message}`);
    process.exit(1);
  }
  const { rows, agree, both } = result;
  console.log(['panel\treader1_arm\treader2_arm\tagree', ...rows.map(r => [r.panel, r.reader1, r.reader2, r.agree ?? ''].join('\t'))].join('\n'));
  console.log(`agree ${agree}/${both} panels both answered`);
}

module.exports = { compare };
