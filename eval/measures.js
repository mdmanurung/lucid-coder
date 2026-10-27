// Regex measures for the G items in eval/rubric.md. Each works on the script
// text after strip(): string contents and comments are blanked, so a pipe or
// a `merge(` inside a string or comment is not counted. Line structure is kept.
// ponytail: regexes, not a parser; a call split by an unbalanced bracket inside
// a regex literal can be miscounted. Upgrade to a tokenizer if a count disagrees
// with a reader on more than a rare script.

const ITERATES = /\b(map\w*|walk\w*|[slvm]?apply|Map|Reduce|group_modify|group_map|reframe|summari[sz]e|mutate|across|transform|agg|groupby)\s*\(|~|\\\(|\blambda\b|\bfor\b/;
const VAGUE = new Set(['tmp', 'temp', 'df', 'df2', 'data', 'x', 'res2', 'foo', 'out2']);

// Blank comments and the contents of string literals (quotes kept as "").
function strip(code, python) {
  let out = '';
  let i = 0;
  while (i < code.length) {
    const c = code[i];
    if (c === '#') {
      while (i < code.length && code[i] !== '\n') i++;
      continue;
    }
    const triple = python && (code.startsWith('"""', i) || code.startsWith("'''", i));
    if (triple || c === '"' || c === "'" || (!python && c === '`')) {
      const close = triple ? code.slice(i, i + 3) : c;
      out += close + close;
      i += close.length;
      while (i < code.length && !code.startsWith(close, i)) {
        if (code[i] === '\\') i++;
        else if (code[i] === '\n') out += '\n';
        i++;
      }
      i += close.length;
      continue;
    }
    out += c;
    i++;
  }
  return out;
}

// Argument text of every call to a name matching re, with balanced parentheses.
function callArgs(text, re) {
  const args = [];
  for (const m of text.matchAll(re)) {
    let depth = 0;
    let j = m.index + m[0].length - 1; // at the "("
    const start = j + 1;
    for (; j < text.length; j++) {
      if (text[j] === '(') depth++;
      else if (text[j] === ')' && --depth === 0) break;
    }
    args.push(text.slice(start, j));
  }
  return args;
}

const codeLineList = text => text.split('\n').filter(l => l.trim() !== '');

function codeMeasures(code, python) {
  const text = strip(code, python);
  const lines = codeLineList(text);
  const count = re => (text.match(re) || []).length;

  // S1: joins that do not name their keys.
  const joins = python
    ? callArgs(text, /\.merge\(/g).filter(a => !/\b(on|left_on|left_index)\s*=/.test(a))
    : [
        ...callArgs(text, /\b(left|right|inner|full|semi|anti)_join\(/g).filter(a => !/\bby\s*=|\bjoin_by\(/.test(a)),
        ...callArgs(text, /\bmerge\(/g).filter(a => !/\bby(\.x)?\s*=/.test(a)),
      ];

  // S5: blanket warning suppression.
  const suppress = python
    ? callArgs(text, /\bwarnings\.(filterwarnings|simplefilter)\(/g).filter(a => /""/.test(a) && !/\bcategory\s*=/.test(a)).length
    : count(/\bsuppressWarnings\(/g) + count(/options\(\s*warn\s*=\s*-1/g);

  // R2: the most steps in one pipe (R) or method chain (Python).
  // A chain continues on a line after one ending in a pipe or ggplot +
  // (R), or on a line starting with "." (Python).
  let longest = 0;
  let run = 0;
  let open = false;
  for (const line of lines) {
    const dot = python && /^\s*\./.test(line);
    const steps = python ? (dot ? 1 : 0) : (line.match(/\|>|%>%/g) || []).length;
    run = (python ? dot : open) ? run + steps : steps;
    open = !python && /(\|>|%>%|\+)\s*$/.test(line);
    longest = Math.max(longest, run);
  }
  if (longest) longest += 1; // n pipe operators join n + 1 steps

  // R3: deepest nesting, by braces in R and by indentation in Python.
  let depth = 0;
  let deepest = 0;
  for (const line of lines) {
    if (python) {
      deepest = Math.max(deepest, Math.floor(line.match(/^ */)[0].length / 4));
      continue;
    }
    for (const ch of line) {
      if (ch === '{') deepest = Math.max(deepest, ++depth);
      else if (ch === '}') depth = Math.max(0, depth - 1);
    }
  }

  // R4: uninformative names on the left of an assignment.
  const assignRe = python ? /^\s*([A-Za-z_]\w*)\s*=(?!=)/ : /^\s*([A-Za-z_.][\w.]*)\s*(<-|=)(?!=)/;
  // Only at bracket depth 0, so a named argument on its own line is not an assignment.
  const names = [];
  let brackets = 0;
  for (const line of lines) {
    if (brackets === 0) {
      const m = line.match(assignRe);
      if (m) names.push(m[1]);
    }
    for (const ch of line) {
      if (ch === "(" || ch === "[") brackets++;
      else if (ch === ")" || ch === "]") brackets = Math.max(0, brackets - 1);
    }
  }
  const vague = names.filter(n => VAGUE.has(n.toLowerCase()) || /^[A-Za-z]$/.test(n)).length;

  // R10: functions defined once and referenced once (a wrapper used once).
  const defRe = python ? /^\s*def\s+([A-Za-z_]\w*)\s*\(/gm : /^\s*([A-Za-z_.][\w.]*)\s*(<-|=)\s*function\s*\(/gm;
  const defined = [...text.matchAll(defRe)].map(m => m[1]);
  // A function passed to map() or apply() is reused, so only a lone direct call counts.
  const defLine = python ? /^\s*def\s/ : /<-\s*function\s*\(/;
  const usedOnce = defined.filter(n => {
    const esc = n.replace(/\./g, '\\.');
    const refs = (text.match(new RegExp(`(?<![\\w.])${esc}(?![\\w.])`, 'g')) || []).length - 1;
    const calls = (text.match(new RegExp(`(?<![\\w.])${esc}\\s*\\(`, 'g')) || []).length - (python ? 1 : 0);
    if (refs !== 1 || calls !== 1) return false;
    // A lone call inside an iterator or lambda still runs once per element.
    const callLine = lines.find(l => new RegExp(`(?<![\\w.])${esc}\\s*\\(`).test(l) && !defLine.test(l));
    return !ITERATES.test(callLine || "");
  }).length;

  // R11: named functions used fewer than three times in the script (user rule:
  // name a function only at 3+ uses; a helper file shared by scripts is separate).
  const underThree = defined.filter(n => {
    const esc = n.replace(/\./g, "\\.");
    return (text.match(new RegExp(`(?<![\\w.])${esc}(?![\\w.])`, "g")) || []).length - 1 < 3;
  }).length;

  // P2: versions recorded.
  const sessionInfo = python
    ? /\bprint_versions\(|\bsession_info\.show\(/.test(text)
    : /\bsessionInfo\(|\bsession_info\(/.test(text);

  // P4: library() or import below the first other code line.
  const loadRe = python ? /^\s*(import|from)\s+\w/ : /^\s*(library|require|suppressPackageStartupMessages)\(/;
  // Setup lines (options(), a docstring, a closing bracket) may sit above the loads.
  const setup = l => /^\s*[})]\s*$/.test(l) || /^\s*options\(/.test(l) || (python && /^\s*""/.test(l));
  const firstOther = lines.findIndex(l => !loadRe.test(l) && !setup(l));
  const lateLoads = firstOther === -1 ? 0 : lines.slice(firstOther + 1).filter(l => loadRe.test(l)).length;

  const out = {
    code_lines: lines.length,
    joins_unkeyed: joins.length,
    warning_suppress: suppress,
    longest_pipe: longest,
    max_nesting: deepest,
    vague_names: vague,
    wrappers_once: usedOnce,
    functions_under3: underThree,
    session_info: sessionInfo,
    late_loads: lateLoads,
    // H5 and H6 are R house rules; Python gets NA.
    base_df_calls: python ? null
      : count(/\bdata\.table\(|\bfread\(|\[\s*,[^\]\n]*:=|\bmerge\(|\b(apply|sapply|lapply|mapply)\(/g) +
        count(/\bfor\s*\(\s*\w+\s+in\s+(seq_len\(\s*)?(1\s*:\s*)?nrow\(/g),
    upper_constants: python ? null : names.filter(n => /^[A-Z][A-Z0-9_]+$/.test(n)).length,
  };
  return out;
}

module.exports = { codeMeasures, strip, callArgs };
