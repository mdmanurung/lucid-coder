# Form review: shaping the scientific-coder persona

## 1. Identity

Ponytail spends three sentences: role, lived pain, creed. "You are a lazy senior developer. Lazy means efficient, not careless. You have seen every over-engineered codebase and been paged at 3am for one. The best code is the code never written." Optimises for the shortest working diff; fears the 3am page. ub-mode has no identity: "Every time you are asked to make a change follow:" then three phases. Caveman: one sentence. Karpathy: none. i-have-adhd puts the identity on the reader: "The reader has ADHD. Output is not just brief. It is shaped so an ADHD brain can act on it."

## 2. Mechanism

Ponytail: a seven-rung ladder, "Stop at the first rung that holds", plus seven rules, each with its why inline: "re-implementing what's a few files over is the most common slop". ub-mode: three phases holding ~40 bare imperatives, no why ("Take complexity out before you add anything."). Karpathy: four principles, a bold thesis each, one test per section. i-have-adhd: ten rules, but the why is paid once up front, "Five facts drive every rule below", so each rule costs one line. That is the cheapest why-device of the set.

## 3. Voice

Ponytail: second person, 8-15 word sentences, concrete nouns (`<input type="date">`, PCA9685), aphorisms that compress a rule: "clever is what someone decodes at 3am", "later can scaffold for itself", "Read fully, then be lazy." ub-mode's triads stick despite no voice: "Pin what the code does. Change the structure. Prove the behavior did not change." Caveman and i-have-adhd teach by Bad/Good pairs, expensive but exact.

## 4. Persistence and scope

Ponytail and caveman share a block: "ACTIVE EVERY RESPONSE. No drift back to over-building. Still active if unsure." Scope by exclusion: "Ponytail governs what you build, not how you talk." Override: "User insists on the full version → build it, no re-arguing." i-have-adhd: "When a rule would delete the answer itself, the task wins; the shape stays." ub-mode has none.

## 5. Output contract

Ponytail: "Code first. Then at most three short lines", pattern `[code] → skipped: [X], add when [Y].` Caveman: `[thing] [action] [reason]. [next step].` ub-mode: "Keep a log ... Report that log upon completion", the opposite. core.md has none.

## 6. ub-mode against ponytail

Weaker: no identity, no why, no persistence, no output shape, tool names the agent may lack ("crapper", "dryer", "mutator"). Stronger: its gates are checkable (a CRAP threshold; "shortest diff" is not), and its boundary triad is the line this persona should steal. "Check the data where it enters. Inside, trust it." is core.md's "Only a check that stops silent data loss, mislabel or overwrite earns its place" in nine words.

## 7. core.md today

As a spec it is strong: an authority stack, a readability ladder, an explicit check policy. As a persona it fails every device above. Line 1 is a tagline, not a character. No why on any rule ("Bars start at 0"). No output contract. Its one memorable line ("Corpus code is evidence, not a template") is about provenance, not analysis. The sections are a style guide's, so the model reads it as lint, which drifts after three turns.

## Blueprint (3,350 bytes)

Order and budget:

1. Identity, 300 B. Three sentences: senior computational biologist; has rerun someone else's pipeline and found a figure whose labels no longer matched its data; writes for the colleague who reruns this in six months. Fear: an untraceable result.
2. Five facts about where analysis code fails, 450 B. Donors not cells are the replicates; positional alignment permutes silently; a filter without a count hides loss; a panel without a denominator misleads; a number without a script cannot rerun. Written as the five candidate lines below, ~90 B each; each habit points at one, so none needs its own why.
3. Six habits as a first-rung-that-holds ladder, 1,000 B. Inline once-or-twice logic; name what a later step reuses; parameters at the top; one map per repeated analysis; a function at 3+ uses; the package's getters. One clause of why each, borrowed from the facts.
4. Trust boundary, 250 B. The ub-mode triad reworded for data: check where it enters, trust inside, at most one extra check per script, none for what fails loudly.
5. Conventions that stay rules, 950 B. Layout, comments, data and figures as flat imperatives, no why, under "House rules (user-decided)". The persona does not argue with the user's picks.
6. Output contract, 200 B. `[code] → removed: n of N per group; check: <the one line that can fail>; decision: <what the code cannot show>`. At most three lines after the code.
7. Authority and persistence, 200 B. Keep "scientific correctness > user's request > learned rules > this style" as the override device, plus "active every response; the user's request wins, the shape stays". Learning: one sentence.

Sum: 3,350.

Paid for by cutting from core.md: "Checks" 600 to 250 B (the identical-IDs, unique-join-key and run-dir trio stays, the user asked for it; the "Not checks: input audits, type checks..." enumeration goes), the tagline (~120 B), "Control flow" folded into the habits (~200 B), "Comments" cut to two imperatives (~250 B), "Learning" to one sentence (~150 B), figure mechanics moved to r.md and python.md (~150 B).

Candidate lines, each encoding an existing core.md rule:

- "Donors are the replicates; cells are the measurements." (replication unit; show each donor as a point)
- "Position is a permutation waiting for a reorder; align by name." (select by name; IDs identical before a positional copy)
- "A silent filter is a result nobody chose." (report removals per filter step)
- "A panel is a claim; put the denominator on it." (state every denominator)
- "A number without a script is a rumour." (a script for every reported number)

Leave to the language files: formatter commands, idioms, check syntax, save calls, session info. Lookups, not character.

Skip: an intensity table and per-habit Bad/Good pairs; one pair at most, under the output contract.
