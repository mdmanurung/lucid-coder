// SessionStart: inject the core ruleset plus active learned rules.
const { readMode, readCore, readActiveRules, readTiers, languages, readEnv, render, writeOutput, rememberDataDir, FALLBACK_CORE } = require('./lib');

// The core text is written before the memory scan, so a slow scan can only
// cost the nudge line.
rememberDataDir();
let written = false;
try {
  if (readMode() !== 'off') {
    const langs = languages(process.cwd());
    writeOutput('SessionStart', render(readCore(), readActiveRules(langs), readTiers(langs), readEnv(process.cwd())));
    written = true;
    const nudge = require('./memory').memoryNudge();
    if (nudge) writeOutput('SessionStart', `\n\n${nudge}`);
  }
} catch (e) {
  if (!written) writeOutput('SessionStart', FALLBACK_CORE);
}
