// SessionStart: inject the core ruleset plus active learned rules.
const { readMode, readCore, readActiveRules, readTiers, languages, readEnv, render, writeOutput, rememberDataDir, FALLBACK_CORE } = require('./lib');

rememberDataDir();
try {
  if (readMode() !== 'off') {
    const langs = languages(process.cwd());
    const nudge = require('./memory').memoryNudge();
    const text = render(readCore(), readActiveRules(langs), readTiers(langs), readEnv(process.cwd()));
    writeOutput('SessionStart', nudge ? `${text}\n\n${nudge}` : text);
  }
} catch (e) {
  writeOutput('SessionStart', FALLBACK_CORE);
}
