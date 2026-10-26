#!/usr/bin/env bash
# Run the eval grid (PLAN.md 8): every task x arm x rep, one at a time.
# Usage: eval/run.sh --tasks 01,02,03,04 --arms A,B,C,D --reps 2 [--out NAME] [--dry-run]
# --out NAME writes into results/NAME instead of a new timestamp, so several
# run.sh processes (one per arm) can share one results folder and run in parallel.
# Env: MODEL (default claude-opus-5-5), BUDGET in USD per run (default 2),
#      EVAL_TASKS_DIR and EVAL_RESULTS_DIR (tests only).
# Runs use the user's Claude Code login. --safe-mode turns off their CLAUDE.md,
# plugins, hooks and MCP servers, so arms differ only by the system-prompt file.
# Arms B-D reach the model through --append-system-prompt-file, because plugin
# hooks do not run in an isolated session (step 0, 2026-10-26).
set -euo pipefail

EVAL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TASKS_DIR="${EVAL_TASKS_DIR:-$EVAL_DIR/tasks}"
RESULTS_DIR="${EVAL_RESULTS_DIR:-$EVAL_DIR/results}"
MODEL="${MODEL:-claude-opus-5-5}"
BUDGET="${BUDGET:-2}"
# Both paths come from `module show statistical/R/4.5.2/gcc.8.5.0`.
R_ROOT="/share/software/statistical/R/R-4.5.2/gcc-8.5.0/R-4.5.2"
ALLOWED_TOOLS=(Read Write Edit "Bash(Rscript:*)" "Bash(ls:*)" "Bash(mkdir:*)" "Bash(cat:*)" "Bash(head:*)")

die() {
  echo "run.sh: $*" >&2
  exit 1
}

usage() {
  die "usage: eval/run.sh --tasks 01,02 --arms A,B,C,D --reps 2 [--out NAME] [--dry-run]"
}

# "01" becomes the task dir named 01_*; an unknown id is kept as given.
resolve_task() {
  local matches=("$TASKS_DIR/$1"_*/)
  if [[ -d "${matches[0]}" ]]; then
    basename "${matches[0]}"
    return
  fi
  echo "$1"
}

write_meta() {
  node -e '
    const a = process.argv;
    const meta = { task: a[2], arm: a[3], rep: Number(a[4]), model: a[5], budget: Number(a[6]),
      start: a[7], end: a[8], exit_code: Number(a[9]), task_dir: a[10], arm_channel: a[11] };
    require("fs").writeFileSync(a[1], JSON.stringify(meta, null, 2) + "\n");
  ' "$@"
}

tasks=""
arms=""
reps=""
out_name=""
dry_run=0
while (($#)); do
  case "$1" in
    --tasks | --arms | --reps)
      (($# >= 2)) || usage
      printf -v "${1#--}" '%s' "$2"
      shift 2
      ;;
    --out)
      (($# >= 2)) || usage
      out_name="$2"
      shift 2
      ;;
    --dry-run)
      dry_run=1
      shift
      ;;
    *) usage ;;
  esac
done
[[ -n "$tasks" && -n "$arms" && "$reps" =~ ^[1-9][0-9]*$ ]] || usage
IFS=, read -r -a task_ids <<<"$tasks"
IFS=, read -r -a arm_ids <<<"$arms"
for arm in "${arm_ids[@]}"; do
  [[ "$arm" =~ ^[ABCD]$ ]] || die "unknown arm $arm; use A, B, C or D"
done

[[ -z "$out_name" || "$out_name" =~ ^[A-Za-z0-9_-]+$ ]] || die "--out takes a plain folder name, not a path"
out="$RESULTS_DIR/${out_name:-$(date -u +%Y%m%dT%H%M%SZ)}"
if ((!dry_run)); then
  command -v claude >/dev/null || die "claude is not on PATH"
  mkdir -p "$RESULTS_DIR"
  # A shared --out folder may already exist; each run dir inside is still refused if taken.
  if [[ -n "$out_name" ]]; then
    mkdir -p "$out"
  else
    mkdir "$out" 2>/dev/null || die "$out already exists; wait a second and rerun"
  fi
fi

for task_id in "${task_ids[@]}"; do
  task="$(resolve_task "$task_id")"
  task_dir="$TASKS_DIR/$task"
  prompt_file="$task_dir/prompt.md"
  for arm in "${arm_ids[@]}"; do
    for ((rep = 1; rep <= reps; rep++)); do
      run_dir="$out/${task}_${arm}_r${rep}"
      arm_file="$run_dir/arm-prompt.md"
      arm_args=()
      channel="none"
      if [[ "$arm" != "A" ]]; then
        arm_args=(--append-system-prompt-file "$arm_file")
        channel="append-system-prompt-file"
      fi
      claude_tail=(--model "$MODEL" --max-budget-usd "$BUDGET" --output-format json
        --permission-mode acceptEdits "${arm_args[@]}" --allowedTools "${ALLOWED_TOOLS[@]}")

      if ((dry_run)); then
        line=""
        if [[ "$arm" != "A" ]]; then
          printf -v line 'node %q %s > %q && ' "$EVAL_DIR/arm-prompt.js" "$arm" "$arm_file"
        fi
        printf '%scd %q && claude --safe-mode --no-session-persistence -p "$(cat %q)"' \
          "$line" "$run_dir/repo" "$prompt_file"
        printf ' %q' "${claude_tail[@]}"
        printf '\n'
        continue
      fi

      [[ -f "$prompt_file" && -d "$task_dir/inputs" ]] || die "$task_dir needs prompt.md and inputs/"
      mkdir "$run_dir" 2>/dev/null || die "$run_dir already exists"
      mkdir "$run_dir/repo"
      cp -R "$task_dir/inputs/." "$run_dir/repo/"
      git -C "$run_dir/repo" init -q
      if [[ "$arm" != "A" ]]; then
        node "$EVAL_DIR/arm-prompt.js" "$arm" >"$arm_file"
      fi

      echo "run.sh: $task arm $arm rep $rep" >&2
      start="$(date -u +%FT%TZ)"
      rc=0
      (
        cd "$run_dir/repo"
        # An exported ANTHROPIC_API_KEY would override the login, so drop it.
        env -u ANTHROPIC_API_KEY \
          PATH="$R_ROOT/bin:$PATH:$HOME/.local/bin" \
          LD_LIBRARY_PATH="$R_ROOT/lib64/R/lib:${LD_LIBRARY_PATH:-}" \
          claude --safe-mode --no-session-persistence -p "$(<"$prompt_file")" "${claude_tail[@]}" \
          </dev/null >"$run_dir/claude.json" 2>"$run_dir/claude.stderr"
      ) || rc=$?
      write_meta "$run_dir/meta.json" "$task" "$arm" "$rep" "$MODEL" "$BUDGET" \
        "$start" "$(date -u +%FT%TZ)" "$rc" "$task_dir" "$channel"
      # An account-level API error (credits, auth) fails every later run too; stop here.
      api_error="$(node -e 'try { const j = require(process.argv[1]); if (j.terminal_reason === "api_error") console.log(j.result) } catch (e) {}' "$run_dir/claude.json")"
      [[ -z "$api_error" ]] || die "API error in $task arm $arm rep $rep: $api_error; stopping"
    done
  done
done

((dry_run)) || echo "run.sh: done; score with: node eval/score.js $out" >&2
