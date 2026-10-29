// A strict YAML subset, dependency-free: block maps, block lists of scalars or
// of maps, flow lists of scalars, quoted and plain scalars, folded scalars
// (key: > or >-, as a map value), # comments.
// Anything else throws an Error whose message starts with "line N:".
const NUMBER = /^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$/;
const KEY = /^([A-Za-z0-9_][\w.\-]*)\s*:(?:\s+(.*))?$/;

function fail(n, msg) {
  throw new Error(`line ${n}: ${msg}`);
}

// Text after a quoted value or flow list may only be a comment.
function rest(text, n) {
  const t = text.trim();
  if (t && !t.startsWith('#')) fail(n, `unexpected text after value: ${t}`);
}

// Parses one quoted scalar at the start of s; returns [value, remainder].
function quoted(s, n) {
  if (s[0] === '"') {
    const m = s.match(/^"(?:[^"\\]|\\.)*"/);
    if (!m) fail(n, 'unclosed double quote');
    try { return [JSON.parse(m[0]), s.slice(m[0].length)]; } catch (e) { fail(n, 'bad escape in double quotes'); }
  }
  const m = s.match(/^'(?:[^']|'')*'/);
  if (!m) fail(n, 'unclosed single quote');
  return [m[0].slice(1, -1).replace(/''/g, "'"), s.slice(m[0].length)];
}

function plain(t, n) {
  if (/^([&*!{}|>%@`,\]]|[-?](\s|$))/.test(t)) fail(n, `unsupported YAML syntax: ${t}`);
  if (/:(\s|$)/.test(t)) fail(n, `a plain value cannot hold ": "; quote it: ${t}`);
  if (t === '' || t === '~' || t === 'null') return null;
  if (t === 'true' || t === 'false') return t === 'true';
  return NUMBER.test(t) ? Number(t) : t;
}

// A value after "key:" or "- ", comment included.
function scalar(s, n) {
  s = s.trim();
  if (s[0] === '"' || s[0] === "'") {
    const [v, tail] = quoted(s, n);
    rest(tail, n);
    return v;
  }
  if (s[0] === '[') {
    const items = [];
    let body = s.slice(1);
    for (;;) {
      body = body.trimStart();
      if (body[0] === ']') break;
      let v;
      if (body[0] === '"' || body[0] === "'") [v, body] = quoted(body, n);
      else {
        const m = body.match(/^[^,\]]*/);
        if (!m[0].trim()) fail(n, body ? 'empty item in flow list' : 'unclosed flow list');
        if (/(^|\s)#/.test(m[0])) fail(n, 'a comment inside a flow list is not supported; quote the #');
        if (/[[{]/.test(m[0])) fail(n, 'nested flow collections are not supported');
        v = plain(m[0].trim(), n);
        body = body.slice(m[0].length);
      }
      items.push(v);
      body = body.trimStart();
      if (body[0] === ']') break;
      if (body[0] !== ',') fail(n, 'unclosed flow list');
      body = body.slice(1);
    }
    rest(body.slice(1), n);
    return items;
  }
  const hash = s.search(/(^|\s)#/); // linear; a regex with \s+ is quadratic on long blanks
  return plain(hash < 0 ? s : s.slice(0, hash).trimEnd(), n);
}

// Folded body: lines join with spaces, a blank line becomes a newline.
function fold(body, strip, first) {
  while (body.length && !body[body.length - 1].trim()) body.pop();
  const ind = body.find(l => l.trim()).match(/^ */)[0].length;
  const text = body.map((l, j) => {
    if (!l.trim()) return '';
    const w = l.match(/^ */)[0].length;
    if (w < ind) fail(first + j, 'folded line is less indented than the first');
    if (w > ind) fail(first + j, 'more-indented lines in a folded scalar are not supported');
    return l.slice(ind).trimEnd();
  }).join('\n').replace(/([^\n])\n(?=[^\n])/g, '$1 ').replace(/\n\n/g, '\n');
  return strip ? text : text + '\n';
}

function parse(text) {
  const lines = [];
  const source = String(text).split(/\r?\n/);
  source.forEach((raw, i) => {
    const n = i + 1;
    if (/^\s*(#.*)?$/.test(raw)) return;
    const indent = raw.match(/^ */)[0].length;
    if (raw[indent] === '\t') fail(n, 'tabs are not allowed; indent with spaces');
    if (/^(---|\.\.\.)(\s|$)/.test(raw)) fail(n, 'document markers are not supported');
    lines.push({ n, indent, text: raw.slice(indent).trimEnd() });
  });
  let i = 0;

  function block(indent) {
    const isList = l => l.text === '-' || l.text.startsWith('- ');
    return isList(lines[i]) ? list(indent, isList) : map(indent);
  }

  function map(indent) {
    const out = {};
    while (i < lines.length && lines[i].indent === indent) {
      const { n, text } = lines[i];
      const m = text.match(KEY);
      if (!m) fail(n, `expected "key: value": ${text}`);
      if (m[1] === '__proto__') fail(n, 'key __proto__ is not allowed');
      if (Object.prototype.hasOwnProperty.call(out, m[1])) fail(n, `duplicate key ${m[1]}`);
      i++;
      const value = m[2] && !/^#/.test(m[2]) ? m[2] : '';
      const folded = value.match(/^>(-?)\s*(#.*)?$/);
      if (folded) {
        // Raw lines, so blank and #-led lines inside the block are text, not skipped.
        let k = n;
        while (k < source.length && (!source[k].trim() || source[k].match(/^ */)[0].length > indent)) k++;
        const body = source.slice(n, k);
        if (!body.some(l => l.trim())) fail(n, 'empty folded scalar');
        out[m[1]] = fold(body, folded[1], n + 1);
        while (i < lines.length && lines[i].n <= k) i++;
      } else if (value) out[m[1]] = scalar(value, n);
      else if (i < lines.length && lines[i].indent > indent) out[m[1]] = block(lines[i].indent);
      else if (i < lines.length && lines[i].indent === indent && /^-( |$)/.test(lines[i].text)) out[m[1]] = block(indent);
      else out[m[1]] = null;
    }
    if (i < lines.length && lines[i].indent > indent) fail(lines[i].n, 'unexpected indent (multi-line values are not supported)');
    return out;
  }

  function list(indent, isList) {
    const out = [];
    while (i < lines.length && lines[i].indent === indent && isList(lines[i])) {
      const line = lines[i];
      const item = line.text.slice(1);
      const pad = item.match(/^ */)[0].length;
      if (!item.trim() || item.trim().startsWith('#')) fail(line.n, 'empty list item');
      if (KEY.test(item.trim())) {
        // "- key: v": the item is a map whose keys sit at the column after "- ".
        lines[i] = { n: line.n, indent: indent + 1 + pad, text: item.trim() };
        out.push(map(indent + 1 + pad));
      } else {
        out.push(scalar(item, line.n));
        i++;
      }
    }
    if (i < lines.length && lines[i].indent > indent) fail(lines[i].n, 'unexpected indent');
    return out;
  }

  if (!lines.length) return null;
  if (lines[0].indent) fail(lines[0].n, 'the top level must not be indented');
  const doc = block(0);
  if (i < lines.length) fail(lines[i].n, 'unexpected line');
  return doc;
}

// Scalars, flow lists of scalars, nested maps. Quotes any string that would not read back as itself.
function stringify(value, indent = '') {
  const one = (v, flow) => {
    if (v === null) return 'null';
    if (typeof v === 'boolean' || (typeof v === 'number' && Number.isFinite(v))) return String(v);
    if (typeof v !== 'string') throw new Error(`cannot write ${typeof v === 'object' ? 'a nested value here' : String(v)}`);
    if (/[\0-\x08\x0a-\x1f\x7f\x85\u2028\u2029]/.test(v)) throw new Error('value has a line break or control character');
    if (flow && /[,[\]{}]/.test(v)) return JSON.stringify(v);
    try { if (scalar(v, 0) === v && v.trim() === v) return v; } catch (e) {}
    return JSON.stringify(v);
  };
  return Object.entries(value).map(([k, v]) => {
    if (Array.isArray(v)) return `${indent}${k}: [${v.map(x => one(x, true)).join(', ')}]\n`;
    if (v && typeof v === 'object') return `${indent}${k}:\n` + stringify(v, indent + '  ');
    return `${indent}${k}: ${one(v)}\n`;
  }).join('');
}

module.exports = { parse, stringify };
