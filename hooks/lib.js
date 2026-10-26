// Shared helpers for the lucid-coder hooks. Every reader fails open: a bad
// store or flag never removes the core ruleset.
const fs = require('fs');
const os = require('os');
const path = require('path');

const RULES_BUDGET = 4000; // learned-rules share of the injected text (PLAN.md 5.5)
const OUTPUT_CAP = 9500; // additionalContext is capped at 10,000 chars
const RULE_MAX = 300;
const SCOPE_ORDER = ['file-glob', 'repo', 'language', 'global']; // narrowest first

const FALLBACK_CORE = [
  '# lucid-coder',
  'Write scientific R/Python for the next reader: one-pass understanding and a clean rerun.',
  'Add a check only if it stops silent data loss, a silent mislabel or a silent overwrite.',
].join('\n');

function stateDir() {
  return process.env.CLAUDE_CONFIG_DIR || path.join(os.homedir(), '.claude');
}

function flagPath() {
  return path.join(stateDir(), '.lucid-coder-active');
}

function dataDir() {
  return process.env.CLAUDE_PLUGIN_DATA || path.join(stateDir(), 'lucid-coder');
}

// Absent flag means on: lucid-coder is always-on until the user turns it off.
function readMode() {
  try {
    return fs.readFileSync(flagPath(), 'utf8').trim() === 'off' ? 'off' : 'on';
  } catch (e) {
    return 'on';
  }
}

function setMode(mode) {
  if (mode === 'off') {
    fs.mkdirSync(stateDir(), { recursive: true });
    fs.writeFileSync(flagPath(), 'off');
    return;
  }
  try { fs.unlinkSync(flagPath()); } catch (e) {}
}

function readCore() {
  try {
    return fs.readFileSync(path.join(__dirname, '..', 'rules', 'core.md'), 'utf8').trim();
  } catch (e) {
    return FALLBACK_CORE;
  }
}

function readActiveRules() {
  try {
    const store = JSON.parse(fs.readFileSync(path.join(dataDir(), 'rules.json'), 'utf8'));
    return (store.rules || []).filter(r => r && r.status === 'active');
  } catch (e) {
    return [];
  }
}

// One line, no leading '#', so rule text cannot open a fake heading.
function escapeRule(text) {
  return String(text || '').replace(/\s+/g, ' ').trim().replace(/^#+\s*/, '').slice(0, RULE_MAX);
}

function scopeRank(rule) {
  const i = SCOPE_ORDER.indexOf(rule.scope && rule.scope.kind);
  return i === -1 ? SCOPE_ORDER.length : i;
}

function render(core, rules) {
  let out = 'LUCID-CODER ACTIVE\n\n' + core;
  if (!rules.length) return out;

  const lines = rules
    .slice()
    .sort((a, b) => scopeRank(a) - scopeRank(b))
    .map(r => `- [${escapeRule(r.topic)}, ${escapeRule(r.scope && r.scope.kind)}] ${escapeRule(r.statement)}`);
  const budget = Math.min(RULES_BUDGET, OUTPUT_CAP - out.length - 200);
  let block = '';
  let shown = 0;
  for (const line of lines) {
    if (block.length + line.length + 1 > budget) break;
    block += line + '\n';
    shown++;
  }
  out += '\n\n## Learned rules\n' + block;
  if (shown < lines.length) {
    out += `${lines.length - shown} rules not shown; run /lucid-coder:rules.\n`;
  }
  return out;
}

// Native Claude Code: SubagentStart drops raw stdout, so it needs the JSON form.
function writeOutput(event, text) {
  try {
    if (event === 'SubagentStart') {
      process.stdout.write(JSON.stringify(
        { hookSpecificOutput: { hookEventName: event, additionalContext: text } }));
      return;
    }
    process.stdout.write(text);
  } catch (e) {
    // stdout closed at hook exit must not surface as a hook failure
  }
}

// Calls cb once with the parsed payload, or {} on bad input, stdin error or timeout.
function readStdin(cb) {
  let input = '';
  let done = false;
  const finish = () => {
    if (done) return;
    done = true;
    let payload = {};
    try { payload = JSON.parse(input.replace(/^﻿/, '')) || {}; } catch (e) {}
    cb(payload);
  };
  process.stdin.on('data', chunk => { input += chunk; });
  process.stdin.on('end', finish);
  process.stdin.on('error', finish);
  setTimeout(finish, 1000).unref();
}

module.exports = {
  FALLBACK_CORE, readMode, setMode, readCore, readActiveRules,
  escapeRule, render, writeOutput, readStdin,
};
