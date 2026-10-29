# Analysis contracts: how a pipeline skill asks its grill questions

Every pipeline skill asks its "Parameters" rows through an analysis contract.
The contract is `contracts/<analysis>.yaml` in the project. It records the
agreed analysis before any code runs. Its fields are `skill`, `engine`,
`question`, `estimand`, `interpretation` (what a positive sign means), `unit`,
`observation_unit`, `replicate`, `columns` (the skill card's required inputs),
`design` (`formula`, `reference_levels`, `covariates`, `blocking`), `contrasts`,
`test`, `multiple_testing` (`method`, `family`), `thresholds`, `filters`,
`seed`, `outputs` and `params`. Any field may hold the literal `ask`. The
`approved` block holds `date`, `by` and `sha256`.

The CLI is `node "${CLAUDE_PLUGIN_ROOT}/scripts/contract.js" <command>`:

| Command | Does |
|---|---|
| `validate <file>` | checks the schema |
| `asks <file>` | prints the dotted paths still `ask` |
| `status <file>` | prints `approved` only when the approval matches the file |
| `approve <file> --by <name>` | records approval; refuses while any field is `ask` |
| `mycelium <file> --table <path> [--sample-col <col>] [--counts <h5ad>]` | prints a mycelium-extra `data_contract.v1` JSON |
| `card <skill>` | validates and summarises the skill's `card.yaml` |
| `gates <file>` | prints the card and variant report gates as a checklist |

## Procedure

1. **Find the contract.** Look for `contracts/*.yaml` whose `skill:` names
   this skill. Several match, or the user's request is unclear: ask which
   analysis. None match: offer to draft one from
   `${CLAUDE_PLUGIN_ROOT}/templates/contract.yaml`. Copy it on a yes, set
   `skill:`, fill fields only from the user's answers, then go to step 2.
   Nothing runs before step 4.
2. **Ask only what is open.** Run `validate`, then `asks`. Ask the user every
   path `asks` prints. Ask `estimand`, `interpretation`, `columns` and `engine`
   explicitly when they are `ask` or missing. Also ask every row of this
   skill's parameter table that has no value in the contract. Ask nothing the
   contract already answers. Write each answer into the contract as the user
   gave it.

   Optional and conditional rows: an omitted `params` key means the optional
   flag is not passed. The user confirms the omission once. Ask only the rows
   that apply to the chosen engine or mode.
3. **Approve.** Show the user the full contract. Ask for approval in words.
   On a yes, run `approve <file> --by <name>` with the name the user gives.
4. **Run only when approved.** Run `status` before the skill's script. Run the
   script only when it prints `approved`. A `--demo` run needs no contract and
   takes nothing from one. A stale status means the file
   changed after approval: show the user what changed (`git diff` on the
   file, or the fields you edited), then repeat step 3.
5. **Map fields to arguments.** Pass contract fields as the skill's
   arguments: `design.formula` to `--formula` or `--design`, `seed` to
   `--seed`, `multiple_testing.method` to `--adjust`, a cutoff in
   `thresholds` to its flag (`--padj-max`), `contrasts` to the skill's
   contrast flags (`--contrast`, `--contrast-term`, `--term`, as the skill
   names them), `filters` to the `--min-*` flags, `outputs` to `--outdir`,
   `engine` to `--engine`, `test` to the skill's method or engine flag. A
   skill argument with no contract field lives under `params`, keyed by the
   flag name without dashes (`params.cell_type_col` for `--cell-type-col`).
   Never pass a value that is not in the approved contract. A contract value
   with no flag in the skill means stop and tell the user. Do not approve or
   run.
6. **Walk the gates.** After `status` prints `approved` and before reporting
   results, run `gates <file>`. The report gates are a checklist the agent
   walks with the user. They do not block the run.

## Prior answers from memory

The project's Claude Code memory may record answers from earlier sessions.
Show one only as a question that names its source session or project: "In
the CD14 project you chose padj < 0.05; same here?" The user confirms or
restates it. Write the user's answer into the contract. Nothing is pre-filled,
and a remembered value never fills a field silently.

## Thresholds stay the user's answers

A contract only records answers. Every threshold, cutoff, filter and seed in
it is still a grill question the user answers. The agent never proposes a
value as a default, and the "Demo-only value" column stays demo-only.

## In a mycelium-extra repo

`.living/` marks core Mycelium. Mycelium-extra is active where
`.mycelium-extra/gate.json` exists. There, after approval, offer to run
`contract.js mycelium <file> --table <path>` on the skill's input table. Pass
the JSON by process substitution and never write it into the repo. Hand it to
the `data-contract-check` skill only on a yes.
