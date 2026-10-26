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

const SHEET_COLUMNS = ['neutral_id', 'q1_one_read', 'q2_find_decisions', 'q3_change_checks', 'q3_which', 'notes'];

function shuffle(items) {
  const a = items.slice();
  for (let i = a.length - 1; i > 0; i--) {
    const j = crypto.randomInt(i + 1);
    [a[i], a[j]] = [a[j], a[i]];
  }
  return a;
}

function makePacket(tsDir, name = 'packet') {
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
    fs.copyFileSync(path.join(r.runDir, 'repo', r.script), path.join(outDir, id + path.extname(r.script)));
    key.push(`${id}\t${r.runId}`);
    sheet.push(id + '\t'.repeat(SHEET_COLUMNS.length - 1));
  });
  fs.writeFileSync(path.join(outDir, 'sheet.tsv'), sheet.join('\n') + '\n');
  fs.writeFileSync(path.join(ts, `${name}-key.tsv`), key.join('\n') + '\n');
  return runs.length;
}

if (require.main === module) {
  const [tsDir, name = 'packet'] = process.argv.slice(2);
  if (!tsDir || !/^[a-z]+$/.test(name)) {
    console.error('usage: node eval/make-packet.js eval/results/<timestamp> [packet|judge]');
    process.exit(1);
  }
  try {
    console.log(`make-packet.js: ${makePacket(tsDir, name)} scripts blinded into ${path.join(tsDir, name)}`);
  } catch (e) {
    console.error(`make-packet.js: ${e.message}`);
    process.exit(1);
  }
}

module.exports = { makePacket, SHEET_COLUMNS };
