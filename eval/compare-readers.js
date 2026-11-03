#!/usr/bin/env node
// Maps two readers' side-by-side picks to arms through the page key and counts agreement.
// Answers: the page's "Copy all answers" JSON, or a TSV with panel and choice columns.
// Agreement counts panels both readers answered; "none" (no preference) is its own pick.
// Usage: node eval/compare-readers.js <side-by-side-key.tsv> <reader1> <reader2>
const fs = require('fs');

function readPicks(file) {
  const text = fs.readFileSync(file, 'utf8').trim();
  if (text.startsWith('[')) return Object.fromEntries(JSON.parse(text).map(r => [r.panel, r.choice]));
  const [head, ...lines] = text.split('\n').map(l => l.split('\t'));
  return Object.fromEntries(lines.map(c => [c[head.indexOf('panel')], c[head.indexOf('choice')]?.trim()]));
}

function compare(keyFile, file1, file2) {
  const arm = {};
  for (const [panel, column, , a] of fs.readFileSync(keyFile, 'utf8').trim().split('\n').slice(1).map(l => l.split('\t'))) arm[`${panel}\t${column}`] = a;
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
  const { rows, agree, both } = compare(keyFile, file1, file2);
  console.log(['panel\treader1_arm\treader2_arm\tagree', ...rows.map(r => [r.panel, r.reader1, r.reader2, r.agree ?? ''].join('\t'))].join('\n'));
  console.log(`agree ${agree}/${both} panels both answered`);
}

module.exports = { compare };
