// UserPromptSubmit: handle /lucid-coder on|off|status and nudge on|off.
// Other prompts get at most one nudge line. Only the on/off flags are stored;
// prompt text never is (PLAN.md 5.4).
const fs = require('fs');
const os = require('os');
const path = require('path');
const { readMode, setMode, writeOutput, readStdin } = require('./lib');

const COMMAND = /^[/@$]lucid-coder(?::lucid-coder)?(?:\s+(on|off|status))?$/;
const NUDGE_COMMAND = /^[/@$]lucid-coder(?::lucid-coder)?\s+nudge\s+(on|off)$/;
const CORRECTION = /\b(always|never|don'?t use|do not use|stop using|prefer|instead of|in this lab|from now on)\b/i;
const EDIT_TOOLS = new Set(['Edit', 'Write', 'MultiEdit', 'NotebookEdit']);
const TAIL_BYTES = 256 * 1024;
const NUDGE = 'If this is a lasting style correction, offer /lucid-coder:learn once.';

function nudgePath() {
  return path.join(process.env.CLAUDE_CONFIG_DIR || path.join(os.homedir(), '.claude'), '.lucid-coder-nudge');
}

// Absent flag means on (decision D5).
function nudgeOn() {
  try {
    return fs.readFileSync(nudgePath(), 'utf8').trim() !== 'off';
  } catch (e) {
    return true;
  }
}

function setNudge(arg) {
  if (arg === 'off') {
    fs.mkdirSync(path.dirname(nudgePath()), { recursive: true });
    fs.writeFileSync(nudgePath(), 'off');
    return;
  }
  try { fs.unlinkSync(nudgePath()); } catch (e) {}
}

function readTail(file) {
  const fd = fs.openSync(file, 'r');
  try {
    const size = fs.fstatSync(fd).size;
    const length = Math.min(size, TAIL_BYTES);
    const buf = Buffer.alloc(length);
    fs.readSync(fd, buf, 0, length, size - length);
    return buf.toString('utf8');
  } finally {
    fs.closeSync(fd);
  }
}

// A typed prompt, not a tool result and not injected skill text (isMeta).
function isRealPrompt(entry) {
  if (entry.type !== 'user' || entry.isMeta || entry.isSidechain) return false;
  const content = entry.message && entry.message.content;
  if (typeof content === 'string') return true;
  if (!Array.isArray(content)) return false;
  return !content.some(b => b && b.type === 'tool_result') && content.some(b => b && b.type === 'text');
}

function editsTool(entry) {
  const content = entry.type === 'assistant' && !entry.isSidechain && entry.message && entry.message.content;
  return Array.isArray(content) && content.some(b => b && b.type === 'tool_use' && EDIT_TOOLS.has(b.name));
}

// True if the assistant turn before the current prompt edited code.
// The current prompt may already be in the transcript, so prompts are only a
// boundary once an assistant line has been seen.
function previousTurnEdited(transcriptPath) {
  try {
    const lines = readTail(transcriptPath).split('\n');
    let seenAssistant = false;
    for (let i = lines.length - 1; i >= 0; i--) {
      let entry;
      try { entry = JSON.parse(lines[i]); } catch (e) { continue; } // cut first line, blanks
      if (!entry) continue;
      if (entry.type === 'assistant') {
        if (editsTool(entry)) return true;
        seenAssistant = true;
        continue;
      }
      if (seenAssistant && isRealPrompt(entry)) return false;
    }
    return false;
  } catch (e) {
    return false;
  }
}

readStdin(payload => {
  try {
    const prompt = String(payload.prompt || '').trim();
    const nudge = prompt.toLowerCase().match(NUDGE_COMMAND);
    if (nudge) {
      setNudge(nudge[1]);
      writeOutput('UserPromptSubmit', nudgeOn() ? 'LUCID-CODER NUDGE ON' : 'LUCID-CODER NUDGE OFF');
      process.exit(0);
    }
    const match = prompt.toLowerCase().match(COMMAND);
    if (match) {
      const arg = match[1] || 'status';
      if (arg !== 'status') setMode(arg);
      writeOutput('UserPromptSubmit', readMode() === 'off' ? 'LUCID-CODER OFF' : 'LUCID-CODER ON');
      process.exit(0);
    }
    if (readMode() !== 'off' && nudgeOn() && CORRECTION.test(prompt) &&
        payload.transcript_path && previousTurnEdited(payload.transcript_path)) {
      writeOutput('UserPromptSubmit', NUDGE);
    }
  } catch (e) {
    // fail open: a broken flag write or transcript must not block the prompt
  }
  process.exit(0);
});
