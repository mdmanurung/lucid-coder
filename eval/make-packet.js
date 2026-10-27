#!/usr/bin/env node
// Blinds each run's script for readers. Writes <ts>/<name>/<id>.R or .py in a fresh
// random order, <ts>/<name>/sheet.tsv for the reader to fill, and the private
// key <ts>/<name>-key.tsv (neutral_id -> run_id) outside the packet dir.
// Ids start with the first letter of <name> (P for packet, J for judge), so
// human and judge ids never match.
// Usage: node eval/make-packet.js eval/results/<timestamp> [packet|judge]
const crypto = require('crypto');
const fs = require('fs');
const path = require('path');
const { findScript } = require('./score');

// Rating fields written at the top of each script with --inline; the reader
// types after each colon and adds own #$ comments anywhere. collect-comments.js reads both.
const INLINE_HEADER = [
  "#$ ===== Reviewer: type after each colon, then add #$ comments anywhere below =====",
  "#$ q1 one read (1-5; can you say what this script does after one read?): ",
  "#$ q2 find decisions (1-5; can you find each threshold, filter, model choice and path?): ",
  "#$ q3 change checks (yes/no; would you add, remove or change a check?): ",
  "#$ q3 which (one sentence): ",
  "#$ q4 procedure sound (yes/no; design, replication unit, method, filters, contrast): ",
  "#$ q4 which (one sentence): ",
  "#$ ---- your comments: #$ for code, #$P for procedure; own line above the code, or after it ----",
].join("\n") + "\n";

const SHEET_COLUMNS = ['neutral_id', 'q1_one_read', 'q2_find_decisions', 'q3_change_checks', 'q3_which', 'notes'];

function shuffle(items) {
  const a = items.slice();
  for (let i = a.length - 1; i > 0; i--) {
    const j = crypto.randomInt(i + 1);
    [a[i], a[j]] = [a[j], a[i]];
  }
  return a;
}

function makePacket(tsDir, name = 'packet', { inline = false } = {}) {
  const ts = path.resolve(tsDir);
  const outDir = path.join(ts, name);
  if (fs.existsSync(outDir)) throw new Error(`${outDir} exists; a new shuffle would orphan filled sheets`);
  const runs = fs.readdirSync(ts).sort()
    .map(d => ({ runDir: path.join(ts, d), runId: `${path.basename(ts)}/${d}` }))
    .filter(r => fs.existsSync(path.join(r.runDir, 'meta.json')))
    .map(r => ({ ...r, script: findScript(path.join(r.runDir, 'repo')) }))
    .filter(r => r.script);
  const prefix = name[0].toUpperCase();
  fs.mkdirSync(outDir);
  const key = ['neutral_id\trun_id'];
  const sheet = [SHEET_COLUMNS.join('\t')];
  shuffle(runs).forEach((r, i) => {
    const id = prefix + String(i + 1).padStart(2, '0');
    const text = fs.readFileSync(path.join(r.runDir, 'repo', r.script), 'utf8');
    fs.writeFileSync(path.join(outDir, id + path.extname(r.script)), (inline ? INLINE_HEADER : '') + text);
    // A ClawBio-style skill also ships its SKILL.md and the report its run wrote,
    // so the reader can judge the procedure; comments go in them as #$ / #$P too.
    const skillDir = path.dirname(r.script);
    if (skillDir.split(path.sep)[0] === 'skills') {
      const extra = { 'SKILL.md': path.join(r.runDir, 'repo', skillDir, 'SKILL.md'),
        'report.md': path.join(r.runDir, 'repo', 'results', 'run1', 'report.md') };
      for (const [suffix, src] of Object.entries(extra)) {
        if (fs.existsSync(src)) fs.copyFileSync(src, path.join(outDir, `${id}.${suffix}`));
      }
    }
    key.push(`${id}\t${r.runId}`);
    sheet.push(id + '\t'.repeat(SHEET_COLUMNS.length - 1));
  });
  fs.writeFileSync(path.join(outDir, 'sheet.tsv'), sheet.join('\n') + '\n');
  fs.writeFileSync(path.join(ts, `${name}-key.tsv`), key.join('\n') + '\n');
  return runs.length;
}

if (require.main === module) {
  const args = process.argv.slice(2);
  const inline = args.includes('--inline');
  const [tsDir, name = 'packet'] = args.filter(a => a !== '--inline');
  if (!tsDir || !/^[a-z]+$/.test(name)) {
    console.error('usage: node eval/make-packet.js eval/results/<timestamp> [packet|judge] [--inline]');
    process.exit(1);
  }
  try {
    console.log(`make-packet.js: ${makePacket(tsDir, name, { inline })} scripts blinded into ${path.join(tsDir, name)}`);
  } catch (e) {
    console.error(`make-packet.js: ${e.message}`);
    process.exit(1);
  }
}

module.exports = { makePacket, SHEET_COLUMNS, INLINE_HEADER };
