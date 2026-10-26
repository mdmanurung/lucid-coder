// SubagentStart: same text as SessionStart, which does not reach subagents.
// LUCID_SUBAGENT_MATCHER (a case-insensitive regex on agent_type) limits which
// subagents get it; anything ambiguous injects, so scoping never drops the persona.
const { readMode, readCore, readActiveRules, render, writeOutput, readStdin, FALLBACK_CORE } = require('./lib');

function inject() {
  try {
    writeOutput('SubagentStart', render(readCore(), readActiveRules()));
  } catch (e) {
    writeOutput('SubagentStart', FALLBACK_CORE);
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
