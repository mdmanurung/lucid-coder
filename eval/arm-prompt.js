#!/usr/bin/env node
// Prints the text arm <A|B|C|D> appends to the system prompt. Step 0
// (2026-10-26) showed plugin hooks do not run in an isolated session, so
// run.sh passes this through --append-system-prompt-file instead of --plugin-dir.
const fs = require('fs');
const os = require('os');
const path = require('path');
const { render, readCore, readTiers } = require('../hooks/lib');

const PONYTAIL = path.join(os.homedir(), '.claude/plugins/marketplaces/ponytail/hooks/ponytail-instructions.js');
const RULES_D = path.join(__dirname, 'arms', 'rules-d.json');

// Straight from the file, not readActiveRules, so no data dir or cwd is involved.
const armDRules = () => JSON.parse(fs.readFileSync(RULES_D, 'utf8')).rules.filter(r => r.status === 'active');

// Second argument: comma-separated tier languages. The eval tasks are R, so r is the default.
const tiers = () => readTiers((process.argv[3] || 'r').split(','));

const texts = {
  A: () => '',
  B: () => require(PONYTAIL).getPonytailInstructions('full'),
  C: () => render(readCore(), [], tiers()),
  D: () => render(readCore(), armDRules(), tiers()),
};

const arm = process.argv[2];
if (!texts[arm]) {
  console.error('usage: node eval/arm-prompt.js <A|B|C|D> [r|python|r,python]');
  process.exit(1);
}
process.stdout.write(texts[arm]());
