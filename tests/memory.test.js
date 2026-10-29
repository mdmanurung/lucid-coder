const test = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync } = require('child_process');

const root = path.join(__dirname, '..');

const memo = (type, description, body = 'Body text.', key = 'type') =>
  `---\nname: x\ndescription: "${description}"\nmetadata:\n  ${key}: ${type}\n---\n\n${body}\n`;

// A temp config dir with memories for two projects: a git repo (atlas) and a plain dir (notes).
function setup() {
  const home = fs.realpathSync(fs.mkdtempSync(path.join(os.tmpdir(), 'lucid-mem-')));
  const atlas = path.join(home, 'work', 'my-atlas');
  const notes = path.join(home, 'work', 'notes.v2');
  fs.mkdirSync(path.join(atlas, '.git'), { recursive: true });
  fs.writeFileSync(path.join(atlas, '.git', 'HEAD'), 'ref: refs/heads/main\n');
  fs.writeFileSync(path.join(atlas, '.git', 'config'), '[remote "origin"]\n\turl = git@github.com:lab/atlas.git\n');
  fs.mkdirSync(notes, { recursive: true });
  const memDir = p => {
    const d = path.join(home, 'projects', p.replace(/[^A-Za-z0-9]/g, '-'), 'memory');
    fs.mkdirSync(d, { recursive: true });
    return d;
  };
  const a = memDir(atlas);
  const n = memDir(notes);
  fs.writeFileSync(path.join(a, 'MEMORY.md'), '- [Fonts](fonts.md) - index line\n');
  fs.writeFileSync(path.join(a, 'fonts.md'), memo('feedback', 'Plot fonts >= 8 pt', 'Use base_size = 8.'));
  fs.writeFileSync(path.join(a, 'ihw.md'), memo('feedback', 'IHW as default p-adjustment'));
  fs.writeFileSync(path.join(a, 'proj.md'), memo('project', 'Deadline is Friday'));
  fs.writeFileSync(path.join(a, 'decoy.md'), `---\nname: d\ndescription: decoy\nnode_type: feedback\ntype: reference\n---\nx\n`);
  fs.writeFileSync(path.join(n, 'fonts-too.md'), memo('feedback', 'Plot fonts >= 8 pt.'));
  return { home, atlas, notes, a };
}

function cli(ctx, ...args) {
  return spawnSync(process.execPath, [path.join(root, 'scripts', 'rules.js'), ...args], {
    cwd: ctx.home,
    env: { ...process.env, HOME: ctx.home, CLAUDE_CONFIG_DIR: ctx.home, CLAUDE_PLUGIN_DATA: path.join(ctx.home, 'data') },
    encoding: 'utf8',
  });
}

const harvest = (ctx, ...args) => JSON.parse(cli(ctx, 'harvest', '--json', ...args).stdout);

// memoryNudge reads the env on each call; a child process keeps it from leaking.
function nudge(ctx) {
  const code = `process.stdout.write(require(${JSON.stringify(path.join(root, 'hooks', 'memory.js'))}).memoryNudge())`;
  return spawnSync(process.execPath, ['-e', code], {
    env: { ...process.env, HOME: ctx.home, CLAUDE_CONFIG_DIR: ctx.home, CLAUDE_PLUGIN_DATA: path.join(ctx.home, 'data') },
    encoding: 'utf8',
  }).stdout;
}

test('harvest keeps feedback only, groups duplicates across projects and resolves repos', () => {
  const ctx = setup();
  const found = harvest(ctx);
  assert.deepStrictEqual(found.map(c => c.description).sort(), ['IHW as default p-adjustment', 'Plot fonts >= 8 pt']);
  const fonts = found.find(c => c.description === 'Plot fonts >= 8 pt');
  assert.deepStrictEqual(fonts.projects_seen.sort(), ['my-atlas', 'notes.v2']);
  assert.strictEqual(fonts.memory_paths.length, 2);
  assert.strictEqual(fonts.body_excerpt, 'Use base_size = 8.');
  const ihw = found.find(c => c.description === 'IHW as default p-adjustment');
  assert.deepStrictEqual(ihw.projects_seen, ['my-atlas']);
  assert.strictEqual(ihw.repo, 'github.com/lab/atlas');
  assert.strictEqual(ihw.source, 'user-memory:my-atlas/ihw.md');
});

test('mark-seen hides old memories from --since-last and the nudge', () => {
  const ctx = setup();
  // Never scanned: the nudge lists dirs only, so unreadable files still nudge.
  const memFiles = fs.readdirSync(ctx.a).map(f => path.join(ctx.a, f));
  memFiles.forEach(f => fs.chmodSync(f, 0o000));
  assert.strictEqual(nudge(ctx), 'Preferences in Claude Code memory not yet reviewed; /lucid-coder:learn --from-memory');
  memFiles.forEach(f => fs.chmodSync(f, 0o644));
  assert.match(cli(ctx, 'harvest', '--mark-seen').stdout, /Marked 3/);
  const stored = fs.readFileSync(path.join(ctx.home, 'data', 'memory-seen.json'), 'utf8');
  assert.doesNotMatch(stored, /fonts|IHW|Plot/);
  assert.deepStrictEqual(harvest(ctx, '--since-last'), []);
  assert.strictEqual(nudge(ctx), '');

  const cohort = path.join(ctx.a, 'cohort.md');
  fs.writeFileSync(cohort, memo('feedback', 'Cohort as random intercept'));
  const later = new Date(Date.now() + 2000); // past any coarse filesystem mtime
  fs.utimesSync(cohort, later, later);
  assert.deepStrictEqual(harvest(ctx, '--since-last').map(c => c.description), ['Cohort as random intercept']);
  assert.strictEqual(nudge(ctx), '1 new preference in memory; /lucid-coder:learn --from-memory to review');

  for (let i = 0; i < 20; i++) {
    const f = path.join(ctx.a, `new${i}.md`);
    fs.writeFileSync(f, memo('feedback', `Preference ${i}`));
    fs.utimesSync(f, later, later);
  }
  assert.strictEqual(nudge(ctx), '20+ new preferences in memory; /lucid-coder:learn --from-memory to review');

  fs.writeFileSync(path.join(ctx.home, '.lucid-coder-nudge'), 'off');
  assert.strictEqual(nudge(ctx), '');
});

test('a user-memory source and an explicit --repo are accepted by propose', () => {
  const ctx = setup();
  const r = cli(ctx, 'propose', '--topic', 'user.ihw', '--scope', 'repo', '--repo', 'github.com/lab/atlas',
    '--statement', 'Use IHW for p-adjustment.', '--why', 'IHW as default p-adjustment',
    '--source', 'user-memory:my-atlas/ihw.md');
  assert.strictEqual(r.status, 0, r.stderr);
  assert.deepStrictEqual(JSON.parse(r.stdout).scope, { kind: 'repo', repo: 'github.com/lab/atlas' });
  const bad = cli(ctx, 'propose', '--topic', 'user.x', '--scope', 'global', '--statement', 's', '--why', 'w',
    '--source', 'user-memory:a/b/c.md');
  assert.strictEqual(bad.status, 1);
});

test('nudge clears after mark-seen even when harvest found nothing new', () => {
  const ctx = setup();
  cli(ctx, 'harvest', '--mark-seen');
  const f = path.join(ctx.a, 'late.md');
  fs.writeFileSync(f, memo('feedback', 'Late preference'));
  const later = new Date(Date.now() + 2000);
  fs.utimesSync(f, later, later);
  assert.match(nudge(ctx), /^1 new preference/);
  // harvest alone is read-only, so the nudge stays until --mark-seen runs.
  harvest(ctx, '--since-last');
  assert.match(nudge(ctx), /^1 new preference/);
  cli(ctx, 'harvest', '--mark-seen');
  assert.strictEqual(nudge(ctx), '');
});

test('distinct descriptions that differ only in punctuation stay separate', () => {
  const ctx = setup();
  fs.writeFileSync(path.join(ctx.a, 'p1.md'), memo('feedback', 'Use p<0.05 cutoff'));
  fs.writeFileSync(path.join(ctx.a, 'p2.md'), memo('feedback', 'Use p>0.05 cutoff'));
  const found = harvest(ctx).filter(c => /^Use p/.test(c.description));
  assert.strictEqual(found.length, 2);
  assert.ok(found.every(c => c.projects_seen.length === 1));
});

test('projects with the same folder name count as two', () => {
  const ctx = setup();
  const other = path.join(ctx.home, 'other', 'my-atlas');
  fs.mkdirSync(other, { recursive: true });
  const d = path.join(ctx.home, 'projects', other.replace(/[^A-Za-z0-9]/g, '-'), 'memory');
  fs.mkdirSync(d, { recursive: true });
  fs.writeFileSync(path.join(d, 'ihw.md'), memo('feedback', 'IHW as default p-adjustment'));
  const ihw = harvest(ctx).find(c => c.description === 'IHW as default p-adjustment');
  assert.deepStrictEqual(ihw.projects_seen, ['my-atlas', 'my-atlas']);
});

test('sources from folder names with spaces pass propose', () => {
  const ctx = setup();
  const dir = path.join(ctx.home, 'work', 'my notes');
  fs.mkdirSync(dir, { recursive: true });
  const d = path.join(ctx.home, 'projects', dir.replace(/[^A-Za-z0-9]/g, '-'), 'memory');
  fs.mkdirSync(d, { recursive: true });
  fs.writeFileSync(path.join(d, 'x.md'), memo('feedback', 'Spaced folder rule'));
  const c = harvest(ctx).find(c => c.description === 'Spaced folder rule');
  assert.strictEqual(c.source, 'user-memory:my_notes/x.md');
  const r = cli(ctx, 'propose', '--topic', 'user.spaced', '--scope', 'global', '--statement', 's', '--why', 'w', '--source', c.source);
  assert.strictEqual(r.status, 0, r.stderr);
  const bad = cli(ctx, 'propose', '--topic', 'user.x', '--scope', 'global', '--statement', 's', '--why', 'w', '--source', 'user-memory:../x.md');
  assert.strictEqual(bad.status, 1);
});

test('parseMemory handles empty, block-scalar and multi-line descriptions', () => {
  const { parseMemory } = require('../hooks/memory');
  const fm = (d, body = 'Body line.') => `---\nname: x\n${d}\ntype: feedback\n---\n${body}\n`;
  assert.strictEqual(parseMemory(fm('description:')).description, 'Body line.');
  assert.strictEqual(parseMemory(fm('description: >-\n  first\n  second')).description, 'first second');
  assert.strictEqual(parseMemory(fm('description: |\n  one\n  two')).description, 'one two');
  assert.strictEqual(parseMemory(fm('description: plain')).description, 'plain');
});
