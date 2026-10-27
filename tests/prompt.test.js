const test = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const os = require('os');
const path = require('path');
const { spawnSync } = require('child_process');

const script = path.join(__dirname, '..', 'hooks', 'prompt.js');
const NUDGE = 'If this is a lasting style correction, offer /lucid-coder:learn once.';
const CORRECTION = 'Always use |> instead of %>%.';

function setup() {
  return fs.mkdtempSync(path.join(os.tmpdir(), 'lucid-prompt-'));
}

function run(home, payload) {
  return spawnSync(process.execPath, [script], {
    env: { ...process.env, HOME: home, CLAUDE_CONFIG_DIR: home },
    input: JSON.stringify(payload),
    encoding: 'utf8',
  }).stdout;
}

// Lines shaped like real Claude Code transcripts (~/.claude/projects/*/*.jsonl).
const userText = text => ({ type: 'user', message: { role: 'user', content: text } });
const userBlocks = text => ({ type: 'user', message: { role: 'user', content: [{ type: 'text', text }] } });
const toolResult = () => ({
  type: 'user',
  message: { role: 'user', content: [{ tool_use_id: 'toolu_1', type: 'tool_result', content: 'ok' }] },
});
const skillText = () => ({ ...userBlocks('Base directory for this skill: /x'), isMeta: true });
const assistantTool = name => ({
  type: 'assistant',
  message: { role: 'assistant', content: [{ type: 'tool_use', id: 'toolu_1', name, input: {} }] },
});
const assistantText = text => ({ type: 'assistant', message: { role: 'assistant', content: [{ type: 'text', text }] } });

function transcript(home, entries) {
  const file = path.join(home, 'session.jsonl');
  fs.writeFileSync(file, entries.map(e => JSON.stringify(e)).join('\n') + '\n');
  return file;
}

const editedTurn = [
  userText('old prompt'),
  assistantTool('Read'), toolResult(),
  userText('fix the script'), skillText(),
  assistantTool('Edit'), toolResult(),
  assistantText('Done.'),
];

test('nudge fires after an edit turn', () => {
  const home = setup();
  const file = transcript(home, editedTurn);
  assert.strictEqual(run(home, { prompt: CORRECTION, transcript_path: file }), NUDGE);
  // the current prompt may already be in the transcript
  const withCurrent = transcript(home, [...editedTurn, userText(CORRECTION)]);
  assert.strictEqual(run(home, { prompt: CORRECTION, transcript_path: withCurrent }), NUDGE);
  // injected skill text after the edit is not a prompt boundary
  const skillAfterEdit = transcript(home, [userText('x'), assistantTool('Edit'), toolResult(), skillText(), assistantText('ok')]);
  assert.strictEqual(run(home, { prompt: CORRECTION, transcript_path: skillAfterEdit }), NUDGE);
});

test('no nudge when the last turn did not edit', () => {
  const home = setup();
  const file = transcript(home, [
    userText('fix it'), assistantTool('Write'), toolResult(),
    userBlocks('what does this do?'), assistantTool('Read'), toolResult(), assistantText('It reads.'),
  ]);
  assert.strictEqual(run(home, { prompt: CORRECTION, transcript_path: file }), '');
});

test('no nudge for a plain prompt, nudge off, mode off or missing transcript', () => {
  const home = setup();
  const file = transcript(home, editedTurn);
  assert.strictEqual(run(home, { prompt: 'Thanks, now plot it.', transcript_path: file }), '');
  assert.strictEqual(run(home, { prompt: CORRECTION }), '');
  assert.strictEqual(run(home, { prompt: CORRECTION, transcript_path: path.join(home, 'nope.jsonl') }), '');

  fs.writeFileSync(path.join(home, '.lucid-coder-nudge'), 'off');
  assert.strictEqual(run(home, { prompt: CORRECTION, transcript_path: file }), '');
  fs.unlinkSync(path.join(home, '.lucid-coder-nudge'));

  fs.writeFileSync(path.join(home, '.lucid-coder-active'), 'off');
  assert.strictEqual(run(home, { prompt: CORRECTION, transcript_path: file }), '');
});

test('nudge command writes and removes its flag', () => {
  const home = setup();
  const flag = path.join(home, '.lucid-coder-nudge');
  assert.strictEqual(run(home, { prompt: '/lucid-coder nudge off' }), 'LUCID-CODER NUDGE OFF');
  assert.strictEqual(fs.readFileSync(flag, 'utf8'), 'off');
  assert.strictEqual(run(home, { prompt: '/lucid-coder:lucid-coder nudge on' }), 'LUCID-CODER NUDGE ON');
  assert.ok(!fs.existsSync(flag));
  assert.ok(!fs.existsSync(path.join(home, '.lucid-coder-active')), 'mode flag untouched');
});

test('lang command writes and removes its flag', () => {
  const home = setup();
  const flag = path.join(home, '.lucid-coder-lang');
  assert.strictEqual(run(home, { prompt: '/lucid-coder lang python' }), 'LUCID-CODER LANG PYTHON (this repo, from the next session start or /clear)');
  assert.deepStrictEqual(Object.values(JSON.parse(fs.readFileSync(flag, 'utf8'))), ['python']);
  assert.strictEqual(run(home, { prompt: '/lucid-coder:lucid-coder lang both' }), 'LUCID-CODER LANG BOTH (this repo, from the next session start or /clear)');
  assert.strictEqual(run(home, { prompt: '/lucid-coder lang auto' }), 'LUCID-CODER LANG AUTO (this repo, from the next session start or /clear)');
  assert.deepStrictEqual(JSON.parse(fs.readFileSync(flag, 'utf8')), {});
  assert.strictEqual(run(home, { prompt: '/lucid-coder lang julia' }), '');
});

test('mode commands still work', () => {
  const home = setup();
  assert.strictEqual(run(home, { prompt: '/lucid-coder off' }), 'LUCID-CODER OFF (from the next session start or /clear)');
  assert.strictEqual(run(home, { prompt: '/lucid-coder status' }), 'LUCID-CODER OFF');
  assert.strictEqual(run(home, { prompt: '/lucid-coder on' }), 'LUCID-CODER ON (from the next session start or /clear)');
  assert.strictEqual(run(home, { prompt: '/lucid-coder' }), 'LUCID-CODER ON');
});
