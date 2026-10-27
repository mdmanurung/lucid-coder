"""Shared pipeline-skill contract; see docs/pipeline-skills.md.

Every skills/<name>/scripts/<name>.py vendors this file next to plots.py and uses it
so all pipeline skills fail, write and record the same way (ClawBio F-05..F-09,
F-16..F-20, F-33). The analysis itself stays in the skill script.
"""

import argparse
import hashlib
import json
import os
import random
import shlex
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np


class SkillError(ValueError):
    """Bad input or a refused run: printed as `ERROR: <message>`, exit 2."""


def fail(errors, key, **fields):
    raise SkillError(errors[key].format(**fields))


def flag(key):
    return f"--{key.replace('_', '-')}"


def parse_required(parser, required, demo, argv, optional=(), demo_extra=("outdir",)):
    """Parse argv; every analysis arg is required unless --demo pins them all.

    `required` lists dest names with no default (each one is a grill question).
    `demo` maps every required arg (and any optional) to its demo-only value.
    Beside --demo only the `demo_extra` args may be given.
    """
    parser.add_argument("--outdir")
    parser.add_argument("--demo", action="store_true")
    args = vars(parser.parse_args(argv))
    if args.pop("demo"):
        given = [k for k in [*required, *optional] if args[k] is not None]
        conflicts = [k for k in given if k not in demo_extra]
        if conflicts:
            raise SkillError(
                "--demo fixes the analysis arguments; drop: "
                + " ".join(flag(k) for k in conflicts)
            )
        return {**demo, **{k: v for k, v in args.items() if v is not None}}
    missing = [k for k in [*required, "outdir"] if args[k] is None]
    if missing:
        raise SkillError(
            "missing required arguments: "
            + " ".join(flag(k) for k in missing)
            + "; ask the user, or use --demo"
        )
    return args


def seed_everything(seed):
    random.seed(seed)
    np.random.seed(seed)


def make_outdir(path):
    """Create <outdir>/{tables,figures}; refuse an existing outdir."""
    outdir = Path(path)
    if outdir.exists():
        raise SkillError(f"outdir already exists: {outdir}; pick a new one")
    for sub in ("tables", "figures"):
        (outdir / sub).mkdir(parents=True)
    return outdir


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def package_versions(packages):
    lines = [f"python {sys.version.split()[0]}"]
    for pkg in packages:
        try:
            lines.append(f"{pkg} {version(pkg)}")
        except PackageNotFoundError:
            lines.append(f"{pkg} not installed")
    return lines


def write_reproducibility(outdir, script, params, keys, argv, packages, extra=()):
    """Write commands.sh (resolved argv), environment.txt and checksums.sha256."""
    repro = Path(outdir) / "reproducibility"
    repro.mkdir()
    resolved = [sys.executable, str(Path(script).resolve())]
    for key in [*keys, "outdir"]:
        value = params.get(key)
        if value is None:
            continue
        # A repeatable flag (argparse action="append") holds a list: one flag per item.
        for item in value if isinstance(value, (list, tuple)) else [value]:
            resolved += [flag(key), str(item)]
    (repro / "commands.sh").write_text(
        "#!/usr/bin/env bash\n"
        f"# argv as run: {shlex.join(argv)}\n"
        f"{shlex.join(resolved)}\n"
    )
    env = package_versions(packages) + list(extra)
    (repro / "environment.txt").write_text("\n".join(env) + "\n")
    outdir = Path(outdir)
    hashed = sorted(outdir.glob("tables/*")) + sorted(outdir.glob("figures/*.pdf"))
    lines = [f"{sha256(p)}  {p.relative_to(outdir)}" for p in hashed if p.is_file()]
    (repro / "checksums.sha256").write_text("\n".join(lines) + "\n")


def write_report(outdir, title, params, keys, summary, captions):
    """report.md: parameter table, summary bullets, one caption per figure.

    Captions follow house-style fig.text: open with what is shown, <= 60 words.
    """
    rows = "\n".join(f"| `{flag(k)}` | `{params.get(k)}` |" for k in keys)
    figures = "\n\n".join(
        f"![{name}](figures/{name}.png)\n\n{caption}"
        for name, caption in captions.items()
    )
    bullets = "\n".join(f"- {line}" for line in summary)
    (Path(outdir) / "report.md").write_text(
        f"# {title}\n\n| parameter | value |\n|---|---|\n{rows}\n\n{bullets}\n\n"
        f"## Figures\n\n{figures}\n"
    )


def check_unique(names, what):
    """Stop on duplicated identifiers (gene names, sample ids) before any join."""
    seen = {}
    for name in names:
        seen[name] = seen.get(name, 0) + 1
    dupes = sorted(k for k, n in seen.items() if n > 1)
    if dupes:
        raise SkillError(
            f"duplicate {what}: {', '.join(map(str, dupes[:5]))}"
            + (f" and {len(dupes) - 5} more" if len(dupes) > 5 else "")
            + "; make them unique before running"
        )


def check_design_rank(meta, terms):
    """Stop when fixed-effect terms are aliased (batch confounded with condition).

    `meta` is a per-sample DataFrame; `terms` are its fixed-effect columns. Categorical
    columns become treatment-coded dummies, numeric ones stay as they are.
    """
    import pandas as pd

    design = pd.get_dummies(meta[list(terms)], drop_first=True, dtype=float)
    design.insert(0, "intercept", 1.0)
    rank = np.linalg.matrix_rank(design.to_numpy())
    if rank < design.shape[1]:
        raise SkillError(
            f"design is rank-deficient ({rank} < {design.shape[1]} columns): terms "
            f"{', '.join(terms)} are aliased, so no effect is identifiable; "
            "drop or merge the confounded term"
        )


def check_replication(obs, sample_col, donor_col, design_terms):
    """Stop on pseudoreplication: a donor with several samples but no donor term.

    Returns (n_donors, n_samples) for the report; `donor_col=None` skips the check.
    """
    n_samples = obs[sample_col].nunique()
    if donor_col is None:
        return None, n_samples
    per_donor = obs.groupby(donor_col, observed=True)[sample_col].nunique()
    repeated = per_donor[per_donor > 1]
    if len(repeated) and donor_col not in design_terms:
        raise SkillError(
            f"{len(repeated)} donor(s) in {donor_col!r} have more than one sample, but "
            f"the design has no {donor_col!r} term; samples are not independent, "
            f"so add {donor_col} to the design (fixed or (1|{donor_col}))"
        )
    return int(per_donor.size), int(n_samples)


def write_result(outdir, result, summary, artifacts, reference=None):
    """result.json with the ClawBio handoff fields (F-33).

    `reference` is the user-stated genome build and annotation release (e.g.
    "GRCh38, GENCODE 44"); a gene-level skill requires it and never guesses it.
    """
    out = {**result, "chat_summary_lines": summary, "preferred_artifacts": artifacts}
    if reference is not None:
        out["reference"] = reference
    (Path(outdir) / "result.json").write_text(
        json.dumps(out, indent=2, default=str) + "\n"
    )


def run(name, analyse, argv=None):
    """CLI wrapper: fixed PDF dates, `ERROR:` + exit 2 on SkillError, no traceback."""
    argv = sys.argv[1:] if argv is None else argv
    # Fixed PDF CreationDate so checksums match across reruns.
    os.environ.setdefault("SOURCE_DATE_EPOCH", "0")
    try:
        outdir = analyse(argv, [sys.argv[0], *argv])
    except SkillError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(f"{name} wrote {Path(outdir) / 'report.md'}")
    return 0


if __name__ == "__main__":
    import tempfile

    errors = {"bad": "bad value {value}"}
    try:
        fail(errors, "bad", value=3)
    except SkillError as exc:
        assert str(exc) == "bad value 3"
    parser = argparse.ArgumentParser()
    parser.add_argument("--alpha", type=float)
    try:
        parse_required(parser, ["alpha"], {"alpha": 0.1}, [])
    except SkillError as exc:
        assert "missing required arguments: --alpha --outdir" in str(exc)
    parser = argparse.ArgumentParser()
    parser.add_argument("--alpha", type=float)
    assert parse_required(parser, ["alpha"], {"alpha": 0.1}, ["--demo"])["alpha"] == 0.1
    with tempfile.TemporaryDirectory() as tmp:
        out = make_outdir(Path(tmp) / "run")
        (out / "tables" / "t.tsv").write_text("a\n")
        write_reproducibility(out, __file__, {"alpha": 1}, ["alpha"], ["x"], ["numpy"])
        assert "--alpha 1" in (out / "reproducibility" / "commands.sh").read_text()
        assert "tables/t.tsv" in (out / "reproducibility/checksums.sha256").read_text()
        try:
            make_outdir(out)
        except SkillError:
            pass
        else:
            raise SystemExit("make_outdir accepted an existing outdir")
    import pandas as pd

    try:
        check_unique(["A", "B", "A"], "gene names")
    except SkillError as exc:
        assert "duplicate gene names: A" in str(exc)
    else:
        raise SystemExit("check_unique accepted a duplicate")
    meta = pd.DataFrame({"cond": ["a", "a", "b", "b"], "batch": ["x", "x", "y", "y"]})
    try:
        check_design_rank(meta, ["cond", "batch"])
    except SkillError as exc:
        assert "rank-deficient" in str(exc)
    else:
        raise SystemExit("check_design_rank accepted aliased terms")
    check_design_rank(meta.assign(batch=["x", "y", "x", "y"]), ["cond", "batch"])
    obs = pd.DataFrame({"sample": ["s1", "s2", "s3"], "donor": ["d1", "d1", "d2"]})
    try:
        check_replication(obs, "sample", "donor", ["cond"])
    except SkillError as exc:
        assert "more than one sample" in str(exc)
    else:
        raise SystemExit("check_replication accepted pseudoreplication")
    assert check_replication(obs, "sample", "donor", ["cond", "donor"]) == (2, 3)
    with tempfile.TemporaryDirectory() as tmp:
        out = make_outdir(Path(tmp) / "run")
        write_reproducibility(out, __file__, {"de": ["a=1", "b=2"]}, ["de"], [], [])
        assert "--de a=1 --de b=2" in (out / "reproducibility/commands.sh").read_text()
