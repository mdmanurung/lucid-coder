#!/usr/bin/env node
// Analysis contracts (contracts/<analysis>.yaml): the design agreed before any code.
// A contract's skill: links to that skill's method card, skills/<skill>/card.yaml,
// and its engine: to the variant skills/<skill>/variants/<engine>.yaml.
// Exit 0 ok, 1 usage or validation error, 2 not approved, 3 approval stale.
const fs = require('fs');
const path = require('path');
const { parseArgs } = require('util');
const { sha256 } = require('../hooks/lib');
const yaml = require('./yaml-lite');

const USAGE = `Usage: node scripts/contract.js <command> <file>
  validate <file>                 required keys and types; lists fields still "ask"
  asks <file>                     fields still "ask", one dotted path per line
  status <file>                   approved (exit 0), not approved (2), approval stale (3)
  approve <file> --by NAME        writes the approved block; refuses while any field is "ask"
  mycelium <file> --table PATH [--sample-col COL] [--counts H5AD]
                                  prints a mycelium-extra.data_contract.v1 JSON to stdout;
                                  card checks with no matching kind are listed on stderr
  card <skill>                    validates and summarises skills/<skill>/card.yaml
  gates <file>                    prints the card and variant report gates as a checklist`;

const ASK = 'ask';
const str = v => typeof v === 'string' && v.trim() !== '';
const map = v => v !== null && typeof v === 'object' && !Array.isArray(v);
const strList = v => Array.isArray(v) && v.every(str);
const strMap = v => map(v) && Object.values(v).every(x => x === ASK || str(x));
const num = v => typeof v === 'number';

// [type check, what it should be]; optional keys are marked with '?'.
const SPEC = {
  version: [v => v === 1, 'the number 1'],
  analysis: [v => /^[A-Za-z0-9_.-]+$/.test(v), 'a name of letters, digits, _ . -'],
  'skill?': [str, 'text'],
  question: [str, 'text'],
  unit: [str, 'text'],
  'observation_unit?': [str, 'text'],
  'engine?': [str, 'text'],
  estimand: [str, 'text'],
  interpretation: [str, 'text, e.g. positive logFC = higher in disease'],
  replicate: [str, 'text'],
  design: [map, 'a map with formula'],
  'design.formula': [str, 'text'],
  'design.reference_levels?': [strMap, 'a map of factor: level'],
  'design.covariates?': [strList, 'a list of column names'],
  'design.blocking?': [strList, 'a list of column names'],
  contrasts: [v => strMap(v) && Object.keys(v).length > 0, 'a map of name: contrast'],
  test: [str, 'text'],
  multiple_testing: [map, 'a map with method and family'],
  'multiple_testing.method': [str, 'text'],
  'multiple_testing.family': [str, 'text'],
  thresholds: [v => map(v) && Object.values(v).every(x => x === ASK || num(x)), 'a map of name: number'],
  'filters?': [v => Array.isArray(v) && v.every(f => map(f) && str(f.step)), 'a list of maps, each with step'],
  seed: [Number.isInteger, 'an integer'],
  outputs: [v => strList(v) && v.length > 0, 'a list of file names'],
  'columns?': [strMap, 'a map of card requirement: sample-table column'],
  'params?': [map, 'a map of argument: value'],
  'approved?': [v => map(v) && str(v.date) && str(v.by) && /^[0-9a-f]{64}$/.test(v.sha256), 'a map of date, by, sha256'],
};

const listOf = (...keys) => v => Array.isArray(v) && v.every(x => map(x) && keys.every(k => str(x[k])));
const STAGES = ['documentation_reviewed', 'source_verified', 'executed', 'simulation_calibrated', 'scientifically_approved'];
const verification = [v => map(v) && Object.keys(v).sort().join() === [...STAGES].sort().join() &&
  Object.values(v).every(x => x === false || (map(x) && str(String(x.date)) && str(x.by) && str(x.note))),
  `a map of ${STAGES.join(', ')}, each false or {date, by, note}`];
const CARD = {
  id: [str, 'text'],
  version: [v => /^\d+\.\d+\.\d+$/.test(v), 'a version like 0.1.0'],
  status: [str, 'text'],
  question: [str, 'text'],
  modalities: [strList, 'a list'],
  designs: [strList, 'a list'],
  requires: [strList, 'a list'],
  statistical_unit: [str, 'text'],
  observation_unit: [str, 'text'],
  methods: [map, 'a map with preferred'],
  'methods.preferred': [strList, 'a list'],
  'methods.alternatives?': [strList, 'a list'],
  checks: [strList, 'a list'],
  'evidence?': [listOf('repository', 'path', 'commit', 'type'), 'a list of maps of repository, path, commit, type'],
  'methodology?': [listOf('doi'), 'a list of maps with doi'],
  'known_failures?': [strList, 'a list'],
  'block_before_report?': [strList, 'a list'],
  verification,
};
const VARIANT = {
  id: [str, 'text'],
  parent: [str, 'the card id'],
  engine: [str, 'the skill\'s --engine value'],
  script: [str, 'a path under the skill directory'],
  formula_pattern: [str, 'text'],
  estimand: [str, 'text'],
  'statistical_unit?': [str, 'text'],
  observation_unit: [str, 'text'],
  parameters: [v => map(v) && Object.values(v).every(p => map(p) && 'package_default' in p && str(p.source)),
    'a map of name: {package_default, source}'],
  known_failures: [strList, 'a list'],
  block_before_report: [strList, 'a list'],
  block_before_fit: [strList, 'a list of check names'],
  verification,
  provenance: [map, 'a map with upstream, commit, docs'],
  'provenance.upstream': [str, 'a url'],
  'provenance.docs': [strList, 'a list of urls'],
};

// Card check -> data_contract.v1 check, using only the kinds data_contract_check.py
// accepts. c is the contract's columns map. A string return is why no check fits.
const CARD_CHECKS = {
  sample_donor_mapping_valid: c => ({ kind: 'unit_of_replication', unit: c.sample_id, group: c.donor_id, min_units_per_group: 1 }),
  minimum_replication: (c, x) => ({ kind: 'unit_of_replication', unit: c.donor_id, group: c.condition, units_may_span_groups: x.paired }),
  paired_structure_preserved: (c, x) => x.paired ? { kind: 'pairing', unit: c.donor_id, within: c.condition }
    : 'design is not paired (the donor_id column is not in the design)',
  batch_confounding_checked: c => c.batch ? { kind: 'batch_confounding', batch: c.batch, contrast: c.condition }
    : 'no batch entry in the contract columns',
  count_integrity: (c, x) => x.counts
    ? { kind: 'h5ad_layer', table: x.counts, h5_path: c.raw_integer_counts, properties: { non_negative: true, integer_valued: true } }
    : 'needs --counts <h5ad>',
};
// Variant gate names from the dreamlet review that one card check covers.
const CHECK_ALIASES = {
  input_counts_are_raw_integer_and_nonnegative: 'count_integrity',
  donor_condition_metadata_constant_per_sample: 'sample_donor_mapping_valid',
  enough_paired_donors_for_inference_not_just_min_samples: 'minimum_replication',
  each_donor_condition_has_at_most_one_biological_sample_or_explicit_replicate_policy: 'paired_structure_preserved',
  condition_and_batch_not_perfectly_confounded: 'batch_confounding_checked',
  sample_ids_unique_by_biological_sample: 'one_row_per_sample',
};
const LAYERS = new Set(['raw_integer_counts']); // card requirements that name an h5ad path, not a column

function die(msg, code = 1) {
  console.error(msg);
  process.exit(code);
}

function read(file) {
  let text;
  try { text = fs.readFileSync(file, 'utf8'); } catch (e) { die(`Cannot read ${file}: ${e.message}`); }
  try { return { text, doc: yaml.parse(text) }; } catch (e) { die(`${file} ${e.message}`); }
}

// Dotted paths of every leaf that is still "ask".
function asks(value, prefix = '') {
  if (typeof value === 'string' && value.trim() === ASK) return [prefix];
  if (value === null) return prefix ? [prefix] : []; // an empty "key:" is unanswered
  if (typeof value !== 'object') return [];
  return Object.entries(value).flatMap(([k, v]) => asks(v, prefix ? `${prefix}.${k}` : k));
}

function problems(doc, specs) {
  if (!map(doc)) return ['the file is not a map of key: value'];
  const out = [];
  const top = new Set(Object.keys(specs).filter(k => !k.includes('.')).map(k => k.replace('?', '')));
  for (const k of Object.keys(doc)) if (!top.has(k)) out.push(`unknown key ${k}`);
  for (const [spec, [ok, want]] of Object.entries(specs)) {
    const key = spec.replace('?', '');
    const [head, sub] = key.split('.');
    const parent = sub ? doc[head] : doc;
    if (sub && !map(parent)) continue; // the parent's own check reports it
    const v = parent[sub || head];
    if (v === ASK) continue;
    if (v === undefined || v === null) {
      if (!spec.endsWith('?')) out.push(`missing ${key}`);
      continue;
    }
    if (!ok(v)) out.push(`${key} must be ${want}`);
  }
  return out;
}

// The skill's method card, or null when the skill has none. Dies on a bad card.
function readCard(skill) {
  if (!/^[a-z0-9-]+$/.test(skill)) die(`skill ${skill} is not a skill name`);
  const file = path.join(__dirname, '..', 'skills', skill, 'card.yaml');
  if (!fs.existsSync(file)) return null;
  const card = read(file).doc;
  const bad = problems(card, CARD);
  if (bad.length) die(`${file} is not a valid card:\n` + bad.map(p => `- ${p}`).join('\n'));
  return card;
}

// The variant for skill + engine, or null when there is none. Dies on a bad variant.
function readVariant(skill, engine, card) {
  if (!/^[a-z0-9_-]+$/.test(engine)) die(`engine ${engine} is not an engine name`);
  const dir = path.join(__dirname, '..', 'skills', skill);
  const file = path.join(dir, 'variants', `${engine}.yaml`);
  if (!fs.existsSync(file)) return null;
  const variant = read(file).doc;
  const bad = problems(variant, VARIANT);
  if (!bad.length && card && variant.parent !== card.id) bad.push(`parent ${variant.parent} is not the card id ${card.id}`);
  if (!bad.length && variant.engine !== engine) bad.push(`engine ${variant.engine} does not match the file name`);
  if (!bad.length && !fs.existsSync(path.join(dir, variant.script))) bad.push(`script ${variant.script} does not exist`);
  if (bad.length) die(`${file} is not a valid variant:\n` + bad.map(p => `- ${p}`).join('\n'));
  return variant;
}

// The card and variant the contract names; either may be null.
function method(doc) {
  const card = map(doc) && str(doc.skill) && doc.skill !== ASK ? readCard(doc.skill) : null;
  const variant = card && str(doc.engine) && doc.engine !== ASK ? readVariant(doc.skill, doc.engine, card) : null;
  if (card && !variant && str(doc.engine) && doc.engine !== ASK) console.error(`note: no variant ${doc.engine} for ${doc.skill}; only the card checks apply`);
  return { card, variant };
}

function contractProblems(doc, file) {
  const out = problems(doc, SPEC);
  if (!map(doc)) return out;
  const base = path.basename(file).replace(/\.ya?ml$/, '');
  if (str(doc.analysis) && doc.analysis !== ASK && doc.analysis !== base) {
    out.push(`analysis is ${doc.analysis} but the file is ${path.basename(file)}; name it contracts/${doc.analysis}.yaml`);
  }
  const named = (v, re) => str(v) && v !== ASK && (re.test(v) ? true : null);
  const skillOk = named(doc.skill, /^[a-z0-9-]+$/);
  if (skillOk === null) out.push(`skill ${doc.skill} is not a skill name (a-z, 0-9, -)`);
  if (named(doc.engine, /^[a-z0-9_-]+$/) === null) out.push(`engine ${doc.engine} is not an engine name (a-z, 0-9, _, -)`);
  if (skillOk && !fs.existsSync(path.join(__dirname, '..', 'skills', doc.skill))) out.push(`skill ${doc.skill} not found`);
  const card = skillOk && !out.some(p => p.startsWith('skill ')) ? readCard(doc.skill) : null;
  if (card && doc.columns !== ASK) {
    for (const need of card.requires) {
      if (!map(doc.columns) || doc.columns[need] == null) {
        out.push(`columns.${need} is missing; the ${doc.skill} card requires it (name its column, or ask)`);
      }
    }
  }
  return out;
}

function validOrDie(doc, file, lead) {
  const bad = contractProblems(doc, file);
  if (bad.length) die(`${lead}${file} is not a valid contract:\n` + bad.map(p => `- ${p}`).join('\n'));
}

// The text the approval covers: everything but a top-level approved block.
function body(text) {
  const lines = text.split(/\r?\n/);
  const out = [];
  for (let i = 0; i < lines.length; i++) {
    if (!/^approved:\s*(#.*)?$/.test(lines[i])) {
      out.push(lines[i]);
      continue;
    }
    while (i + 1 < lines.length && /^(\s|$)/.test(lines[i + 1])) i++;
  }
  return out.join('\n').trimEnd() + '\n';
}

function status(text, doc) {
  if (!map(doc) || !map(doc.approved)) return ['not approved', 2];
  if (doc.approved.sha256 !== sha256(body(text))) return ['approval stale: the contract changed after approval', 3];
  return [`approved by ${doc.approved.by} on ${doc.approved.date}`, 0];
}

// Column names a formula uses: "~ a + b:c + (1|d)" gives a, b, c, d.
// Ceiling: a function call such as log(x) also yields "log"; parse calls if a contract needs them.
function formulaColumns(formula) {
  return (formula.replace(/^[^~]*~/, '').match(/[A-Za-z_.][A-Za-z0-9_.]*/g) || []);
}

const commands = {
  validate(file) {
    const { doc } = read(file);
    validOrDie(doc, file, '');
    const open = asks(doc);
    console.log(`${file} is valid.` + (open.length ? ` ${open.length} fields still ask: ${open.join(', ')}` : ''));
  },

  asks(file) {
    const open = asks(read(file).doc);
    if (open.length) console.log(open.join('\n'));
  },

  status(file) {
    const { text, doc } = read(file);
    const [msg, code] = status(text, doc);
    console.log(msg);
    process.exit(code);
  },

  approve(file, opts) {
    if (!str(opts.by)) die('approve needs --by NAME');
    const { text, doc } = read(file);
    validOrDie(doc, file, 'Not approved; ');
    const open = asks(doc);
    if (open.length) die(`Not approved; these fields are still ask:\n${open.join('\n')}`);
    const kept = body(text);
    const approved = { date: new Date().toISOString().slice(0, 10), by: opts.by, sha256: sha256(kept) };
    const tmp = `${file}.${process.pid}.tmp`;
    fs.writeFileSync(tmp, kept + '\n' + yaml.stringify({ approved }));
    fs.renameSync(tmp, file);
    console.log(`Approved ${file} (sha256 ${approved.sha256}). Any edit clears the approval.`);
  },

  mycelium(file, opts) {
    if (!str(opts.table)) die('mycelium needs --table PATH (the sample table the analysis uses)');
    const { doc } = read(file);
    validOrDie(doc, file, '');
    const d = doc.design;
    const c = doc.columns || {};
    const open = asks({ design: d, columns: c });
    if (open.length) die(`Answer these first; the data contract needs them:\n${open.join('\n')}`);
    const sample = opts['sample-col'] || c.sample_id;
    const design = [...new Set([...formulaColumns(d.formula), ...(d.blocking || []), ...(d.covariates || [])])];
    const mapped = Object.entries(c).filter(([k]) => !LAYERS.has(k)).map(([, v]) => v);
    const columns = [...new Set([sample, ...mapped, ...design, ...Object.keys(d.reference_levels || {})].filter(Boolean))];
    const checks = [{ id: 'columns', kind: 'schema', columns }];
    if (design.length) checks.push({ id: 'design_complete', kind: 'cohort', no_missing: design });
    for (const [factor, level] of Object.entries(d.reference_levels || {})) {
      checks.push({ id: `reference_${factor}`, kind: 'cohort', column: factor, levels: [level] });
    }
    if (sample) checks.push({ id: 'one_row_per_sample', kind: 'unique_key', columns: [sample] });
    const { card, variant } = method(doc);
    const unsupported = [];
    const ctx = { paired: design.includes(c.donor_id), counts: opts.counts };
    const names = new Set([...(card ? card.checks : []), ...(variant ? variant.block_before_fit : [])]
      .map(n => CHECK_ALIASES[n] || n));
    for (const name of names) {
      if (checks.some(k => k.id === name)) continue;
      const made = CARD_CHECKS[name] ? CARD_CHECKS[name](c, ctx) : 'no data_contract.v1 kind checks it';
      if (typeof made === 'string') unsupported.push(`${name}: ${made}`);
      else checks.push({ id: name, ...made });
    }
    if (unsupported.length) console.error(`unsupported (not checked by mycelium):\n${unsupported.join('\n')}`);
    console.log(JSON.stringify({ schema: 'mycelium-extra.data_contract.v1', table: opts.table, checks }, null, 2));
  },

  card(skill) {
    const card = readCard(skill);
    if (!card) die(`No card for ${skill}; expected skills/${skill}/card.yaml`);
    const unchecked = Object.entries(card.verification).filter(([, v]) => !v).map(([k]) => k);
    let variants = [];
    try { variants = fs.readdirSync(path.join(__dirname, '..', 'skills', skill, 'variants')).map(f => f.replace(/\.yaml$/, '')); } catch (e) {}
    for (const engine of variants) readVariant(skill, engine, card);
    console.log([
      `${card.id} ${card.version} (${card.status})`,
      `question: ${card.question.trim()}`,
      `requires: ${card.requires.join(', ')}`,
      `unit: ${card.statistical_unit}`,
      `methods: ${card.methods.preferred.join(', ')}`,
      `checks: ${card.checks.join(', ')}`,
      `not yet verified: ${unchecked.join(', ') || 'none'}`,
      `variants: ${variants.join(', ') || 'none'}`,
    ].join('\n'));
  },

  gates(file) {
    const { doc } = read(file);
    validOrDie(doc, file, '');
    const { card, variant } = method(doc);
    const lines = [card, variant].filter(Boolean)
      .flatMap(m => (m.block_before_report || []).map(g => `- [ ] ${g} (${m.id})`));
    console.log(lines.length ? `Before reporting ${doc.analysis}:\n${lines.join('\n')}` : `No report gates for ${doc.analysis}.`);
  },
};

function main() {
  let parsed;
  try {
    parsed = parseArgs({
      allowPositionals: true,
      options: Object.fromEntries(['by', 'table', 'sample-col', 'counts'].map(k => [k, { type: 'string' }])),
    });
  } catch (e) {
    die(`${e.message}\n${USAGE}`);
  }
  const [command, file, ...extra] = parsed.positionals;
  if (!Object.hasOwn(commands, command) || !file || extra.length) die(USAGE);
  commands[command](file, parsed.values);
}

main();
