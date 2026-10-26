// SessionStart: inject the core ruleset plus active learned rules.
const { readMode, readCore, readActiveRules, render, writeOutput, rememberDataDir, FALLBACK_CORE } = require('./lib');

rememberDataDir();
try {
  if (readMode() !== 'off') {
    writeOutput('SessionStart', render(readCore(), readActiveRules()));
  }
} catch (e) {
  writeOutput('SessionStart', FALLBACK_CORE);
}
