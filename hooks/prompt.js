// UserPromptSubmit: handle /lucid-coder on|off|status. Other prompts get no output.
const { readMode, setMode, writeOutput, readStdin } = require('./lib');

const COMMAND = /^[/@$]lucid-coder(?::lucid-coder)?(?:\s+(on|off|status))?$/;

readStdin(payload => {
  try {
    const match = String(payload.prompt || '').trim().toLowerCase().match(COMMAND);
    if (match) {
      const arg = match[1] || 'status';
      if (arg !== 'status') setMode(arg);
      writeOutput('UserPromptSubmit', readMode() === 'off' ? 'LUCID-CODER OFF' : 'LUCID-CODER ON');
    }
  } catch (e) {
    // fail open: a broken flag write must not block the prompt
  }
  process.exit(0);
});
