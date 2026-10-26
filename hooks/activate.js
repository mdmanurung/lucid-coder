// SessionStart: inject the core ruleset plus active learned rules.
const { readMode, readCore, readActiveRules, readTiers, languages, render, writeOutput, rememberDataDir, FALLBACK_CORE } = require('./lib');

rememberDataDir();
try {
  if (readMode() !== 'off') {
    writeOutput('SessionStart', render(readCore(), readActiveRules(), readTiers(languages(process.cwd()))));
  }
} catch (e) {
  writeOutput('SessionStart', FALLBACK_CORE);
}
