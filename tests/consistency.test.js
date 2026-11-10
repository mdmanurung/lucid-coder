// Repo consistency: versions, skill frontmatter and cited paths agree with disk.
const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');

const ROOT = path.resolve(__dirname, '..');
const read = p => fs.readFileSync(path.join(ROOT, p), 'utf8');
const exists = p => fs.existsSync(path.join(ROOT, p));

// Every SKILL.md under skills/, at any depth, as a root-relative path.
function skillFiles(dir = 'skills') {
  return fs.readdirSync(path.join(ROOT, dir), { withFileTypes: true }).flatMap(e => {
    const p = path.join(dir, e.name);
    if (e.isDirectory()) return skillFiles(p);
    return e.name === 'SKILL.md' ? [p] : [];
  });
}
const SKILLS = skillFiles();
const SKILL_DIRS = new Set(fs.readdirSync(path.join(ROOT, 'skills'), { withFileTypes: true })
  .filter(e => e.isDirectory()).map(e => e.name));

test('package.json version matches plugin.json', () => {
  const pkg = JSON.parse(read('package.json'));
  const plugin = JSON.parse(read('.claude-plugin/plugin.json'));
  assert.equal(pkg.version, plugin.version);
});

test('each SKILL.md has a name equal to its directory and a description', () => {
  for (const file of SKILLS) {
    const fm = read(file).match(/^---\n([\s\S]*?)\n---/);
    assert.ok(fm, `${file}: no frontmatter`);
    const name = fm[1].match(/^name:\s*(.+)$/m);
    assert.equal(name && name[1].trim(), path.basename(path.dirname(file)), `${file}: name`);
    // description: inline text, or a block scalar (>-, |) followed by indented lines
    const desc = fm[1].match(/^description:[ \t]*(.*)\n?((?:[ \t]+.*\n?)*)/m);
    const text = desc ? (desc[1].replace(/^[>|][-+]?$/, '') + desc[2]).trim() : '';
    assert.ok(text.length > 0, `${file}: empty description`);
  }
});

test('every references/ path cited in a SKILL.md exists', () => {
  const missing = SKILLS.flatMap(file => {
    const cited = read(file).match(/(?<![\w/.~-])references\/[\w./-]*[\w-]\.\w+/g) || [];
    return cited.filter(ref => !exists(path.join(path.dirname(file), ref)))
      .map(ref => `${file}: ${ref}`);
  });
  assert.deepEqual(missing, []);
});

test('every /lucid-coder:<skill> citation resolves to a skill directory', () => {
  const sources = ['README.md', ...SKILLS,
    ...fs.readdirSync(path.join(ROOT, 'rules')).filter(f => f.endsWith('.md')).map(f => `rules/${f}`)];
  const missing = sources.flatMap(file =>
    [...read(file).matchAll(/\/lucid-coder:([a-z0-9-]+)/g)]
      .filter(m => !SKILL_DIRS.has(m[1]))
      .map(m => `${file}: /lucid-coder:${m[1]}`));
  assert.deepEqual(missing, []);
});

test('every markdown path cited in AGENTS.md exists', () => {
  const text = read('AGENTS.md');
  const cited = [
    ...[...text.matchAll(/`([^`\s]+\.md)`/g)].map(m => m[1]),
    ...[...text.matchAll(/\]\(([^)\s#]+\.md)(?:#[^)]*)?\)/g)].map(m => m[1]),
  ].filter(p => !/^[a-z]+:\/\//.test(p));
  assert.deepEqual(cited.filter(p => !exists(p)), []);
});
