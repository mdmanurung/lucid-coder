#!/usr/bin/env node
// Collects a reader's in-line `#$` comments from a blinded packet into
// comments.tsv, one row per comment: neutral id, line, comment, and the code it
// refers to (the same line for a trailing comment, else the next code line).
// `#$ q1=4 q2=3` anywhere in a file records scores; `q3=yes|no` too.
// With the key next to the packet, rows also get run_id and arm; read the key
// only after the reader has finished.
// Usage: node eval/collect-comments.js eval/results/<packet-ts>/packet [--key]
const fs = require('fs');
const path = require('path');

const MARK = /#\$(P?)\s?(.*)$/;
const SCORE = /\b(q[123])\s*=\s*(\d|yes|no)\b/gi;
// The rating fields make-packet.js --inline writes at the top of each script.
const FIELD = /^\s*#\$\s*(q1|q2|q3 which|q3|q4 which|q4)\b[^:]*:\s*(.*)$/i;
const FIELD_KEYS = { q1: "q1", q2: "q2", q3: "q3", "q3 which": "q3_which", q4: "q4", "q4 which": "q4_which" };
const RULER = /^\s*#\$\s*(-{3,}|={3,})/;

function collect(packetDir) {
  const rows = [];
  const scores = {};
  // P01.R, plus P01.SKILL.md and P01.report.md for a skill task.
  for (const file of fs.readdirSync(packetDir).filter(f => /^[A-Z]\d+\.(R|py|SKILL\.md|report\.md)$/i.test(f)).sort()) {
    const id = file.split(".")[0];
    const part = file.slice(id.length + 1);
    const lines = fs.readFileSync(path.join(packetDir, file), 'utf8').split('\n');
    let original = 0; // line number in the script as delivered, before the reader's own lines
    const pending = [];
    scores[id] ??= {};
    lines.forEach(line => {
      const field = line.match(FIELD);
      if (field) {
        if (field[2].trim()) scores[id][FIELD_KEYS[field[1].toLowerCase()]] = field[2].trim();
        return;
      }
      if (RULER.test(line)) return;
      const m = line.match(MARK);
      const own = m && /^\s*#\$/.test(line);
      if (own) {
        pending.push([m[1] ? "procedure" : "code", m[2].trim()]);
        return;
      }
      original++;
      const code = line.replace(MARK, '').trimEnd();
      for (const [kind, text] of pending.splice(0)) rows.push({ id, part, kind, line: original, text, code: code.trim() });
      if (m) rows.push({ id, part, kind: m[1] ? "procedure" : "code", line: original, text: m[2].trim(), code: code.trim() });
    });
    for (const [kind, text] of pending) rows.push({ id, part, kind, line: original, text, code: '' });
  }
  // Inline q1=4 style scores fill any field left blank in the header.
  for (const r of rows) {
    for (const s of r.text.matchAll(SCORE)) scores[r.id][s[1].toLowerCase()] ??= s[2].toLowerCase();
  }
  return { rows, scores };
}

function readKey(packetDir) {
  const keyFile = path.join(path.dirname(packetDir), `${path.basename(packetDir)}-key.tsv`);
  const key = {};
  for (const line of fs.readFileSync(keyFile, 'utf8').trim().split('\n').slice(1)) {
    const [id, runId] = line.split('\t');
    key[id] = { run_id: runId, arm: (runId.match(/_([A-D])_r\d+$/) || [])[1] || '' };
  }
  return key;
}

const oneLine = v => String(v ?? '').replace(/[\t\n]/g, ' ');

function main() {
  const [packetDir, flag] = process.argv.slice(2);
  if (!packetDir || !fs.existsSync(packetDir)) {
    console.error('usage: node eval/collect-comments.js eval/results/<packet-ts>/packet [--key]');
    process.exit(1);
  }
  const { rows, scores } = collect(packetDir);
  const key = flag === '--key' ? readKey(packetDir) : null;
  const extra = id => (key ? [key[id]?.arm, key[id]?.run_id] : []);
  const keyCols = key ? ['arm', 'run_id'] : [];
  const suffix = key ? '-keyed' : '';
  const write = (name, cols, lines) =>
    fs.writeFileSync(path.join(packetDir, name + suffix + '.tsv'), [cols.join('\t'), ...lines.map(l => l.map(oneLine).join('\t'))].join('\n') + '\n');
  write('comments', ['neutral_id', 'file', 'kind', 'line', 'comment', 'code', ...keyCols],
    rows.map(r => [r.id, r.part, r.kind, r.line, r.text, r.code, ...extra(r.id)]));
  write('scores', ['neutral_id', 'q1', 'q2', 'q3', 'q3_which', 'q4', 'q4_which', ...keyCols],
    Object.entries(scores).map(([id, s]) => [id, s.q1, s.q2, s.q3, s.q3_which, s.q4, s.q4_which, ...extra(id)]));
  const file = path.join(packetDir, 'comments' + suffix + '.tsv');
  console.log(`collect-comments.js: ${rows.length} comments from ${new Set(rows.map(r => r.id)).size} scripts -> ${file}`);
}

if (require.main === module) main();

module.exports = { collect };
