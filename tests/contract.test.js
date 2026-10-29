const test = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync } = require('child_process');

const root = path.join(__dirname, '..');
const TEMPLATE = fs.readFileSync(path.join(root, 'templates', 'contract.yaml'), 'utf8');

// The template, copied to a temp contracts/ dir under its analysis name.
function contract(text = TEMPLATE) {
  const dir = path.join(fs.mkdtempSync(path.join(os.tmpdir(), 'lucid-contract-')), 'contracts');
  fs.mkdirSync(dir);
  const file = path.join(dir, 'cd14_disease_vs_healthy.yaml');
  fs.writeFileSync(file, text);
  return file;
}

function cli(...args) {
  return spawnSync(process.execPath, [path.join(root, 'scripts', 'contract.js'), ...args], { encoding: 'utf8' });
}

const answered = TEMPLATE.replace('lfc_min: ask', 'lfc_min: 1').replace('cutoff: ask', 'cutoff: 10')
  .replace('min_cells: ask', 'min_cells: 10');

test('the template validates and lists its ask fields', () => {
  const file = contract();
  const v = cli('validate', file);
  assert.strictEqual(v.status, 0, v.stderr);
  assert.match(v.stdout, /3 fields still ask/);
  assert.deepStrictEqual(cli('asks', file).stdout.trim().split('\n'),
    ['thresholds.lfc_min', 'filters.0.cutoff', 'params.min_cells']);
});

test('validate reports missing keys, bad types, unknown keys and a wrong file name', () => {
  const file = contract(TEMPLATE.replace('seed: 1', 'seed: one').replace(/^test: .*$/m, 'tset: x'));
  const r = cli('validate', file);
  assert.strictEqual(r.status, 1);
  assert.match(r.stderr, /seed must be an integer/);
  assert.match(r.stderr, /missing test/);
  assert.match(r.stderr, /unknown key tset/);
  fs.renameSync(file, path.join(path.dirname(file), 'other.yaml'));
  assert.match(cli('validate', path.join(path.dirname(file), 'other.yaml')).stderr, /name it contracts\/cd14_disease_vs_healthy\.yaml/);
  assert.strictEqual(cli('validate', contract('design: ask\n')).status, 1, 'ask does not excuse missing keys');
  assert.strictEqual(cli('validate', contract(TEMPLATE.replace('formula: "~ condition + (1 | donor)"', 'formula: ask'))).status, 0);
});

test('approve refuses while any field is ask, then approval goes stale on any edit', () => {
  const pending = contract();
  const refused = cli('approve', pending, '--by', 'Ann');
  assert.strictEqual(refused.status, 1);
  assert.match(refused.stderr, /thresholds\.lfc_min/);
  assert.strictEqual(cli('status', pending).status, 2);

  const file = contract(answered);
  assert.strictEqual(cli('approve', file, '--by', 'Ann Lee').status, 0);
  const s = cli('status', file);
  assert.strictEqual(s.status, 0);
  assert.match(s.stdout, /approved by Ann Lee on \d{4}-\d{2}-\d{2}/);
  assert.strictEqual(cli('approve', file, '--by', 'Ann Lee').status, 0, 're-approval replaces the block');
  assert.strictEqual(fs.readFileSync(file, 'utf8').match(/^approved:/gm).length, 1);

  fs.writeFileSync(file, fs.readFileSync(file, 'utf8').replace('# experimental unit', '# the unit'));
  const stale = cli('status', file);
  assert.strictEqual(stale.status, 3);
  assert.match(stale.stdout, /approval stale/);
});

test('mycelium merges contract columns with the card checks and lists what it cannot check', () => {
  const r = cli('mycelium', contract(), '--table', 'data/samples.tsv', '--counts', 'data/cells.h5ad');
  assert.strictEqual(r.status, 0, r.stderr);
  assert.deepStrictEqual(JSON.parse(r.stdout), {
    schema: 'mycelium-extra.data_contract.v1',
    table: 'data/samples.tsv',
    checks: [
      { id: 'columns', kind: 'schema', columns: ['sample', 'donor', 'cell_type', 'condition'] },
      { id: 'design_complete', kind: 'cohort', no_missing: ['condition', 'donor'] },
      { id: 'reference_condition', kind: 'cohort', column: 'condition', levels: ['healthy'] },
      { id: 'one_row_per_sample', kind: 'unique_key', columns: ['sample'] },
      { id: 'sample_donor_mapping_valid', kind: 'unit_of_replication', unit: 'sample', group: 'donor', min_units_per_group: 1 },
      { id: 'count_integrity', kind: 'h5ad_layer', table: 'data/cells.h5ad', h5_path: 'layers/counts',
        properties: { non_negative: true, integer_valued: true } },
      { id: 'minimum_replication', kind: 'unit_of_replication', unit: 'donor', group: 'condition', units_may_span_groups: true },
      { id: 'paired_structure_preserved', kind: 'pairing', unit: 'donor', within: 'condition' },
    ],
  });
  assert.deepStrictEqual(r.stderr.trim().split('\n'), [
    'unsupported (not checked by mycelium):',
    'design_matrix_full_rank: no data_contract.v1 kind checks it',
    'batch_confounding_checked: no batch entry in the contract columns',
    'sensitivity_to_low_cell_counts: no data_contract.v1 kind checks it',
    'no_complete_confounding_or_rank_deficiency: no data_contract.v1 kind checks it',
    'paired_donor_counts_reported_per_celltype_after_filtering: no data_contract.v1 kind checks it',
  ]);
  const unpaired = cli('mycelium', contract(TEMPLATE.replace('blocking: [donor]', 'blocking: []').replace('+ (1 | donor)', '')), '--table', 't.tsv');
  assert.match(unpaired.stderr, /paired_structure_preserved: design is not paired/);
  assert.match(unpaired.stderr, /count_integrity: needs --counts/);
  assert.strictEqual(cli('mycelium', contract()).status, 1, '--table is required');
  const open = cli('mycelium', contract(TEMPLATE.replace('donor_id: donor', 'donor_id: ask')), '--table', 't.tsv');
  assert.strictEqual(open.status, 1);
  assert.match(open.stderr, /^columns\.donor_id$/m);
});

test('the pseudobulk-de card parses, folded question included, and validates', () => {
  const card = require('../scripts/yaml-lite').parse(
    fs.readFileSync(path.join(root, 'skills', 'pseudobulk-de', 'card.yaml'), 'utf8'));
  assert.strictEqual(card.question, 'Compare expression between conditions within annotated cell populations.\n');
  assert.deepStrictEqual(Object.values(card.verification), [false, false, false, false, false]);
  const r = cli('card', 'pseudobulk-de');
  assert.strictEqual(r.status, 0, r.stderr);
  assert.match(r.stdout, /^scrna\.de\.pseudobulk-donor 0\.1\.0 \(candidate\)/);
  assert.match(r.stdout, /not yet verified: documentation_reviewed, source_verified, executed, simulation_calibrated, scientifically_approved\nvariants: dreamlet/);
  assert.strictEqual(cli('card', 'no-such-skill').status, 1);
});

test('validate errors when a card requirement has no columns entry', () => {
  const r = cli('validate', contract(TEMPLATE.replace(/^ {2}donor_id: donor\n/m, '')));
  assert.strictEqual(r.status, 1);
  assert.match(r.stderr, /columns\.donor_id is missing; the pseudobulk-de card requires it/);
  assert.strictEqual(cli('validate', contract(TEMPLATE.replace('donor_id: donor', 'donor_id: ask'))).status, 0, 'ask counts as mapped');
  assert.strictEqual(cli('validate', contract(TEMPLATE.replace(/^skill: .*\n/m, '').replace(/^columns:\n(  .*\n)+/m, ''))).status, 0,
    'no skill, no card, no columns needed');
});

test('validate requires estimand and interpretation', () => {
  const r = cli('validate', contract(TEMPLATE.replace(/^estimand: .*\n/m, '').replace(/^interpretation: .*\n/m, '')));
  assert.strictEqual(r.status, 1);
  assert.match(r.stderr, /missing estimand/);
  assert.match(r.stderr, /missing interpretation/);
  assert.strictEqual(cli('validate', contract(TEMPLATE.replace(/^estimand: .*$/m, 'estimand: ask'))).status, 0);
});

test('the dreamlet variant parses and gates prints its report gates', () => {
  const v = require('../scripts/yaml-lite').parse(
    fs.readFileSync(path.join(root, 'skills', 'pseudobulk-de', 'variants', 'dreamlet.yaml'), 'utf8'));
  assert.strictEqual(v.parent, 'scrna.de.pseudobulk-donor');
  assert.strictEqual(v.parameters['min.cells'].package_default, 5);
  assert.deepStrictEqual(v.verification.documentation_reviewed, { date: '2026-10-29', by: 'user review', note: 'external review document; execution not run' });
  const r = cli('gates', contract());
  assert.strictEqual(r.status, 0, r.stderr);
  const gates = r.stdout.split('\n').filter(l => l.startsWith('- [ ] '));
  assert.strictEqual(gates.length, v.block_before_report.length);
  assert.ok(gates.every(g => g.endsWith('(scrna.de.pseudobulk-donor.dreamlet)')));
});

// A throwaway plugin root, so the variant can add a check its card lacks.
test('mycelium adds the variant block_before_fit checks to the card checks', () => {
  const plugin = fs.mkdtempSync(path.join(os.tmpdir(), 'lucid-plugin-'));
  for (const f of ['scripts/contract.js', 'scripts/yaml-lite.js', 'hooks/lib.js']) {
    fs.mkdirSync(path.dirname(path.join(plugin, f)), { recursive: true });
    fs.copyFileSync(path.join(root, f), path.join(plugin, f));
  }
  const skill = path.join(plugin, 'skills', 'fx');
  fs.mkdirSync(path.join(skill, 'variants'), { recursive: true });
  fs.mkdirSync(path.join(skill, 'scripts'));
  fs.writeFileSync(path.join(skill, 'scripts', 'fit.R'), '');
  const card = fs.readFileSync(path.join(root, 'skills', 'pseudobulk-de', 'card.yaml'), 'utf8')
    .replace(/^checks:\n( {2}- .*\n)+/m, 'checks: [sample_donor_mapping_valid]\n');
  fs.writeFileSync(path.join(skill, 'card.yaml'), card);
  const variant = fs.readFileSync(path.join(root, 'skills', 'pseudobulk-de', 'variants', 'dreamlet.yaml'), 'utf8')
    .replace('script: scripts/dreamlet.R', 'script: scripts/fit.R')
    .replace(/^block_before_fit:\n( {2}- .*\n)+/m, 'block_before_fit: [minimum_replication, made_up_check]\n');
  fs.writeFileSync(path.join(skill, 'variants', 'dreamlet.yaml'), variant);
  const r = spawnSync(process.execPath, [path.join(plugin, 'scripts', 'contract.js'), 'mycelium',
    contract(TEMPLATE.replace('skill: pseudobulk-de', 'skill: fx')), '--table', 't.tsv'], { encoding: 'utf8' });
  assert.strictEqual(r.status, 0, r.stderr);
  const ids = JSON.parse(r.stdout).checks.map(c => c.id);
  assert.deepStrictEqual(ids.slice(-2), ['sample_donor_mapping_valid', 'minimum_replication']);
  assert.match(r.stderr, /made_up_check: no data_contract\.v1 kind checks it/);
});

test('review fixes: unit is no column, paired from formula, bad skill, null and folded asks, own commands', () => {
  const cols = JSON.parse(cli('mycelium', contract(TEMPLATE.replace('unit: donor', 'unit: patient')), '--table', 't.tsv').stdout).checks[0].columns;
  assert.ok(!cols.includes('patient'), 'unit is a concept, not a column');
  const formula = cli('mycelium', contract(TEMPLATE.replace('blocking: [donor]', 'blocking: []')), '--table', 't.tsv');
  assert.doesNotMatch(formula.stderr, /paired_structure_preserved/, '(1|donor) in the formula is paired');
  const bad = cli('validate', contract(TEMPLATE.replace('skill: pseudobulk-de', 'skill: no-such')));
  assert.strictEqual(bad.status, 1);
  assert.match(bad.stderr, /skill no-such not found/);
  assert.match(cli('validate', contract(TEMPLATE.replace('skill: pseudobulk-de', 'skill: ../x'))).stderr, /not a skill name/);
  assert.match(cli('validate', contract(TEMPLATE.replace(/^engine: .*$/m, 'engine: ../x'))).stderr, /not an engine name/);
  const noVariant = cli('gates', contract(TEMPLATE.replace(/^engine: .*$/m, 'engine: nope')));
  assert.strictEqual(noVariant.status, 0);
  assert.match(noVariant.stderr, /note: no variant nope/);
  assert.deepStrictEqual(cli('asks', contract(answered.replace('seed: 1', 'seed:'))).stdout.trim().split('\n'), ['seed']);
  assert.deepStrictEqual(cli('asks', contract(answered.replace('test: ', 'test: >\n  ask\n#'))).stdout.trim().split('\n'), ['test']);
  assert.strictEqual(cli('constructor', contract()).status, 1);
});
