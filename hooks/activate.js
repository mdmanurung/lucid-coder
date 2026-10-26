// SessionStart: inject the core ruleset plus active learned rules.
const { readMode, readCore, readActiveRules, render, writeOutput, FALLBACK_CORE } = require('./lib');

try {
  if (readMode() !== 'off') {
    writeOutput('SessionStart', render(readCore(), readActiveRules()));
  }
} catch (e) {
  writeOutput('SessionStart', FALLBACK_CORE);
}
