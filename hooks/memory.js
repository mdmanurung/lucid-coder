// Read-only scan of Claude Code auto-memory for `feedback` entries, plus the
// session-start nudge. Only a scan time and content hashes are stored, never
// memory text (PLAN.md 5.4). Nothing here writes to Claude Code memory.
const fs = require('fs');
const os = require('os');
const path = require('path');
const { dataDir, repoInfo, sha256 } = require('./lib');

const configDir = () => process.env.CLAUDE_CONFIG_DIR || path.join(os.homedir(), '.claude');
const seenPath = () => path.join(dataDir(), 'memory-seen.json');
const encode = name => name.replace(/[^A-Za-z0-9]/g, '-');

// Same flag file prompt.js writes with `nudge off`; absent means on.
function nudgeOn() {
  try {
    return fs.readFileSync(path.join(configDir(), '.lucid-coder-nudge'), 'utf8').trim() !== 'off';
  } catch (e) {
    return true;
  }
}

// Claude Code names a project dir by its path with every non-alphanumeric
// character as '-', which is ambiguous; walk the disk for the path that encodes to it.
function decodeProject(enc, dir = '/') {
  let names = [];
  try { names = fs.readdirSync(dir); } catch (e) { return null; }
  for (const name of names) {
    const e = '-' + encode(name);
    if (enc === e) return path.join(dir, name);
    if (enc.startsWith(e + '-')) {
      const found = decodeProject(enc.slice(e.length), path.join(dir, name));
      if (found) return found;
    }
  }
  return null;
}

// type is matched at line start so `node_type:` never counts.
function parseMemory(text) {
  const m = text.match(/^---\r?\n([\s\S]*?)\r?\n---\r?\n?([\s\S]*)$/);
  if (!m) return null;
  const type = (m[1].match(/^[ \t]*type:[ \t]*["']?([\w-]+)/m) || [])[1] || null;
  const d = m[1].match(/^([ \t]*)description:[ \t]*(.*)$/m);
  let description = d ? d[2].trim() : '';
  if (d && /^([>|][+-]?)?$/.test(description)) {
    // Block scalar or empty value: join the lines indented under the key.
    const rest = m[1].slice(d.index + d[0].length).split(/\r?\n/).slice(1);
    const lines = [];
    for (const l of rest) {
      if (l.trim() && l.search(/\S/) <= d[1].length) break;
      if (l.trim()) lines.push(l.trim());
    }
    description = lines.join(' ');
  }
  if (!description) description = m[2].trim().split(/\r?\n/)[0].slice(0, 200);
  if (/^".*"$/.test(description)) {
    try { description = JSON.parse(description); } catch (e) { description = description.slice(1, -1); }
  } else if (/^'.*'$/.test(description)) {
    description = description.slice(1, -1).replace(/''/g, "'");
  }
  return { type, description, body: m[2].trim() };
}

// Every memory file as { file, dir }, from directory listings only (no reads).
function memoryFiles() {
  const root = path.join(configDir(), 'projects');
  let dirs = [];
  try { dirs = fs.readdirSync(root).sort(); } catch (e) { return []; }
  const out = [];
  for (const dir of dirs) {
    const memDir = path.join(root, dir, 'memory');
    let files = [];
    try { files = fs.readdirSync(memDir).filter(f => f.endsWith('.md') && f !== 'MEMORY.md').sort(); } catch (e) { continue; }
    for (const f of files) out.push({ file: path.join(memDir, f), dir });
  }
  return out;
}

// The feedback memories among files as { file, dir, text, ...parsed }.
function feedbackMemories(files = memoryFiles()) {
  const out = [];
  for (const { file, dir } of files) {
    let text;
    try { text = fs.readFileSync(file, 'utf8'); } catch (e) { continue; }
    const parsed = parseMemory(text);
    if (parsed && parsed.type === 'feedback') out.push({ file, dir, text, ...parsed });
  }
  return out;
}

function readSeen() {
  try {
    const s = JSON.parse(fs.readFileSync(seenPath(), 'utf8'));
    return { scanned: s.scanned || null, seen: new Set(Array.isArray(s.seen) ? s.seen : []) };
  } catch (e) {
    return { scanned: null, seen: new Set() };
  }
}

// One candidate per description (case, whitespace and trailing punctuation ignored); projects_seen >= 2 marks a global candidate.
// Semantic duplicates are left to the skill.
function harvest({ sinceLast = false } = {}) {
  const { seen } = readSeen();
  const projects = {};
  const project = dir => {
    if (!projects[dir]) {
      const where = decodeProject(dir);
      const repo = where && repoInfo(where);
      projects[dir] = { name: where ? path.basename(where) : dir, path: where, repo: repo ? repo.id : null };
    }
    return projects[dir];
  };
  const groups = new Map();
  for (const m of feedbackMemories()) {
    if (sinceLast && seen.has(sha256(m.text))) continue;
    const key = m.description.toLowerCase().replace(/\s+/g, ' ').replace(/[\s.!?;:,]+$/, '') || m.file;
    const p = project(m.dir);
    const pid = p.path || m.dir;
    const group = groups.get(key);
    if (group) {
      group.memory_paths.push(m.file);
      if (!group.dirs.has(pid)) {
        group.dirs.add(pid);
        group.projects_seen.push(p.name);
      }
      continue;
    }
    groups.set(key, {
      memory_path: m.file,
      memory_paths: [m.file],
      project: p.name,
      project_path: p.path,
      repo: p.repo,
      source: `user-memory:${p.name}/${path.basename(m.file)}`.replace(/\s+/g, '_'),
      description: m.description,
      body_excerpt: m.body.slice(0, 300),
      projects_seen: [p.name],
      dirs: new Set([pid]),
    });
  }
  return [...groups.values()].map(({ dirs, ...g }) => g);
}

function markSeen() {
  const hashes = feedbackMemories().map(m => sha256(m.text));
  const file = seenPath();
  fs.mkdirSync(path.dirname(file), { recursive: true });
  fs.writeFileSync(file, JSON.stringify({ scanned: new Date().toISOString(), seen: hashes }) + '\n');
  return hashes.length;
}

// One line for SessionStart, or ''. Never throws. Reading every memory costs
// seconds on a network home, so it lists dirs, stats files and reads at most
// NUDGE_CAP files newer than the last scan; before any scan it reads none.
const NUDGE_CAP = 20;
function memoryNudge() {
  try {
    if (!nudgeOn()) return '';
    const files = memoryFiles();
    if (!files.length) return '';
    const { scanned, seen } = readSeen();
    const since = scanned ? Date.parse(scanned) : NaN;
    if (Number.isNaN(since)) return 'Preferences in Claude Code memory not yet reviewed; /lucid-coder:learn --from-memory';
    const newer = files.filter(({ file }) => {
      try { return fs.statSync(file).mtimeMs > since; } catch (e) { return false; }
    });
    if (newer.length > NUDGE_CAP) return `${NUDGE_CAP}+ new preferences in memory; /lucid-coder:learn --from-memory to review`;
    const n = feedbackMemories(newer).filter(m => !seen.has(sha256(m.text))).length;
    if (!n) return '';
    return `${n} new ${n === 1 ? 'preference' : 'preferences'} in memory; /lucid-coder:learn --from-memory to review`;
  } catch (e) {
    return '';
  }
}

module.exports = { harvest, markSeen, memoryNudge, parseMemory, decodeProject };
