#!/usr/bin/env bash
# Two blinded LLM judges score every run's script (PLAN.md 8; proxy, not human
# evidence). Scripts go in verbatim under neutral ids (make-packet.js judge),
# one script per call, from an empty cwd, with all tools off.
# Writes <ts>/judge/, <ts>/judge-key.tsv, raw replies in <ts>/judge-raw/ and
# the parsed <ts>/judges.tsv.
# Usage: eval/judge.sh eval/results/<timestamp>
# Env: JUDGE_BUDGET in USD per call (default 0.5); JUDGE_CLAUDE (tests only).
set -euo pipefail

EVAL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JUDGES=(claude-opus-5-5 claude-sonnet-5-5)
JUDGE_BUDGET="${JUDGE_BUDGET:-0.5}"
JUDGE_CLAUDE="${JUDGE_CLAUDE:-claude}"

die() {
  echo "judge.sh: $*" >&2
  exit 1
}

(($# == 1)) && [[ -d "$1" ]] || die "usage: eval/judge.sh eval/results/<timestamp>"
ts="$(cd "$1" && pwd)"
raw="$ts/judge-raw"
[[ ! -e "$raw" ]] || die "$raw exists; judges already ran for $ts"

node "$EVAL_DIR/make-packet.js" "$ts" judge
mkdir "$raw"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
instructions="$(<"$EVAL_DIR/judge-prompt.md")"

# make-packet.js numbered the scripts in shuffled order, so J01, J02, ... is that order.
shopt -s nullglob
for script in "$ts"/judge/*.R "$ts"/judge/*.py; do
  id="$(basename "$script")"
  id="${id%.*}"
  lang=r
  [[ "$script" == *.py ]] && lang=python
  prompt="$instructions"$'\n\n```'"$lang"$'\n'"$(<"$script")"$'\n```\n'
  for model in "${JUDGES[@]}"; do
    echo "judge.sh: $id $model" >&2
    (
      cd "$work"
      # Same isolation as run.sh: the user's login, no CLAUDE.md, plugins or hooks.
      env -u ANTHROPIC_API_KEY \
        "$JUDGE_CLAUDE" --safe-mode --no-session-persistence -p "$prompt" --model "$model" --output-format json \
        --max-budget-usd "$JUDGE_BUDGET" --tools "" \
        </dev/null >"$raw/$id.$model.json" 2>"$raw/$id.$model.stderr"
    ) || echo "judge.sh: $id $model exited non-zero; parse will record NA" >&2
  done
done

node "$EVAL_DIR/report.js" --judges "$ts"
