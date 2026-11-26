// SubagentStart: same text as SessionStart, which does not reach subagents.
// With `learned`, the learned rules instead of the core (a separate hook, as in activate.js).
// LUCID_SUBAGENT_MATCHER (a case-insensitive regex on agent_type) limits which
// subagents get it; anything ambiguous injects, so scoping never drops the persona.
const { readMode, readCore, readActiveRules, readTiers, languages, readEnv, renderCore, renderLearned, writeOutput, readStdin, FALLBACK_CORE } = require('./lib');

const learnedPart = process.argv[2] === 'learned';

function inject() {
  try {
    const langs = languages(process.cwd());
    if (learnedPart) {
      const text = renderLearned(readActiveRules(langs));
      if (text) writeOutput('SubagentStart', text);
      return;
    }
    writeOutput('SubagentStart', renderCore(readCore(), readTiers(langs), readEnv(process.cwd())));
  } catch (e) {
    if (!learnedPart) writeOutput('SubagentStart', FALLBACK_CORE);
  }
}

let matcher = null;
try {
  if (process.env.LUCID_SUBAGENT_MATCHER) matcher = new RegExp(process.env.LUCID_SUBAGENT_MATCHER, 'i');
} catch (e) {
  matcher = null;
}

if (readMode() === 'off') process.exit(0);

// No matcher: stay synchronous and never touch stdin.
if (!matcher) {
  inject();
  process.exit(0);
}

readStdin(payload => {
  const agentType = String(payload.agent_type || '').trim();
  if (!agentType || matcher.test(agentType)) inject();
  process.exit(0);
});
