// SessionStart: inject the core ruleset plus active learned rules.
const { readMode, readCore, readActiveRules, readTiers, languages, render, writeOutput, rememberDataDir, FALLBACK_CORE } = require('./lib');

rememberDataDir();
try {
  if (readMode() !== 'off') {
    const langs = languages(process.cwd());
    writeOutput('SessionStart', render(readCore(), readActiveRules(langs), readTiers(langs)));
  }
} catch (e) {
  writeOutput('SessionStart', FALLBACK_CORE);
}
