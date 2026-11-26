// SessionStart: inject the core ruleset, or with `learned` the active learned rules.
// hooks.json runs both parts as separate hooks, because Claude Code caps each
// hook's output on its own: the core text can no longer crowd out learned rules.
const { readMode, readCore, readActiveRules, readTiers, languages, readEnv, renderCore, renderLearned, writeOutput, rememberDataDir, FALLBACK_CORE } = require('./lib');

const learnedPart = process.argv[2] === 'learned';

if (learnedPart) {
  try {
    if (readMode() !== 'off') {
      const text = renderLearned(readActiveRules(languages(process.cwd())));
      if (text) writeOutput('SessionStart', text);
    }
  } catch (e) {
    // fail open: a broken store never blocks the session
  }
} else {
  // The core text is written before the memory scan, so a slow scan can only
  // cost the nudge line.
  rememberDataDir();
  let written = false;
  try {
    if (readMode() !== 'off') {
      const langs = languages(process.cwd());
      writeOutput('SessionStart', renderCore(readCore(), readTiers(langs), readEnv(process.cwd())));
      written = true;
      const nudge = require('./memory').memoryNudge();
      if (nudge) writeOutput('SessionStart', `\n\n${nudge}`);
    }
  } catch (e) {
    if (!written) writeOutput('SessionStart', FALLBACK_CORE);
  }
}
