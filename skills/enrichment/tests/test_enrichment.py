"""CLI tests for scripts/enrichment.py: demo goldens, output contract, fail strings.

Run from skills/enrichment with pytest;
see references/maintainer.md for the exact command.
fgsea tests need Rscript with fgsea (RSCRIPT or PATH) and are skipped without it.
"""

import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "scripts" / "enrichment.py"
EXAMPLES = SKILL_DIR / "examples"
TEMPLATES = SKILL_DIR.parents[1] / "templates"
RSCRIPT = os.environ.get("RSCRIPT") or shutil.which("Rscript")

DE = [f"T={EXAMPLES / 'demo_T.tsv'}", f"Mono={EXAMPLES / 'demo_Mono.tsv'}"]
COMMON = {
    "--de": DE,
    "--gene-sets": str(EXAMPLES / "demo.gmt"),
    "--min-size": "10",
    "--max-size": "500",
    "--padj-max-sets": "0.05",
    "--limit": "2.5",
    "--seed": "0",
    "--reference": "synthetic demo, no genome build",
}
ORA = {
    **COMMON,
    "--method": "ora",
    "--universe": "tested",
    "--padj-max": "0.05",
    "--lfc-min": "0.5",
    "--direction": "both",
}
FGSEA = {
    **COMMON,
    "--method": "fgsea",
    "--rank-col": "stat",
    "--n-perm-simple": "10000",
}


def run(*argv, env=None):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *argv],
        capture_output=True,
        text=True,
        check=False,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1", **(env or {})},
    )


def run_with(tmp_path, base=ORA, drop=(), env=None, **changes):
    """`base` arguments with `changes` ({"limit": "0"}) and `drop` applied."""
    args = {**base, **{f"--{k.replace('_', '-')}": v for k, v in changes.items()}}
    argv = []
    for key, value in args.items():
        if key not in drop:
            for v in value if isinstance(value, list) else [value]:
                argv += [key, v]
    return run(*argv, "--outdir", str(tmp_path / "out"), env=env)


def enrichment(outdir):
    return pd.read_csv(outdir / "tables" / "enrichment.tsv", sep="\t")


def significant(outdir, contrast):
    table = enrichment(outdir)
    return set(table.loc[table["significant"] & (table["contrast"] == contrast), "set"])


def contract_paths():
    """File paths listed in the SKILL.md output tree, relative to outdir."""
    text = (SKILL_DIR / "SKILL.md").read_text()
    tree = text.split("## Output structure")[1].split("```")[1]
    paths, stack = [], []
    for line in tree.splitlines():
        match = re.search(r"[├└]── (\S+)", line)
        if not match:
            continue
        depth = match.start() // 4
        stack[depth:] = [match.group(1)]
        if not match.group(1).endswith("/"):
            paths.append("".join(stack))
    return paths


def is_error(completed, message):
    return (
        completed.returncode == 2
        and f"ERROR: {message}" in completed.stderr
        and "Traceback" not in completed.stderr
    )


@pytest.fixture(scope="module")
def demo_out(tmp_path_factory):
    if not RSCRIPT:
        pytest.skip("Rscript not found")
    out = tmp_path_factory.mktemp("demo") / "out"
    completed = run("--demo", "--outdir", str(out))
    assert completed.returncode == 0, completed.stderr
    return out


@pytest.fixture(scope="module")
def ora_out(tmp_path_factory):
    out = tmp_path_factory.mktemp("ora")
    completed = run_with(out)
    assert completed.returncode == 0, completed.stderr
    return out / "out"


def test_demo_golden_values(demo_out):
    result = json.loads((demo_out / "result.json").read_text())
    assert result["method"] == "fgsea" and result["seed"] == 0
    assert result["n_sets_in_file"] == 10
    t, mono = result["contrasts"]
    # SET_TINY (8 genes) is below --min-size 10 in both contrasts.
    assert t["n_sets_tested"] == mono["n_sets_tested"] == 9
    assert t["n_genes"] == mono["n_genes"] == 1000
    # All 317 distinct set genes are in both demo DE tables; no NA row is dropped.
    assert (t["n_matched"], t["n_in_resource"]) == (317, 317)
    assert result["dropped_per_group"] == {"T": 0, "Mono": 0}
    assert result["reference"] == "synthetic demo, no genome build"
    assert set(t["significant_sets"]) == {"SET_UP", "SET_DOWN"}
    assert (
        t["significant_sets"]["SET_UP"] > 2 and t["significant_sets"]["SET_DOWN"] < -2
    )
    assert mono["n_significant"] == 0
    table = enrichment(demo_out)
    assert table["set"].head(2).tolist() == ["SET_UP", "SET_DOWN"]
    assert table["contrast"].unique().tolist() == ["T", "Mono"]
    curves = pd.read_csv(demo_out / "figures" / "running_es.csv")
    assert curves.groupby("contrast")["set"].first().to_dict() == {
        "Mono": "SET_NULL5",
        "T": "SET_UP",
    }
    assert 2 <= len(result["chat_summary_lines"]) <= 3


def test_ora_golden_values(ora_out):
    result = json.loads((ora_out / "result.json").read_text())
    t, mono = result["contrasts"]
    assert (t["n_universe"], t["n_foreground"]) == (980, 31)
    assert mono["n_foreground"] == 0
    assert result["dropped_per_group"] == {"T": 20, "Mono": 20}
    assert significant(ora_out, "T") == {"SET_UP", "SET_DOWN"}
    assert significant(ora_out, "Mono") == set()
    up = enrichment(ora_out).query("contrast == 'T' and set == 'SET_UP'").iloc[0]
    assert (up["set_size"], up["size"], up["overlap"]) == (40, 40, 15)


def test_na_rows_dropped_per_contrast(tmp_path):
    de = pd.read_csv(EXAMPLES / "demo_T.tsv", sep="\t")
    de.loc[:2, "padj"] = None
    path = tmp_path / "na.tsv"
    de.to_csv(path, sep="\t", index=False)
    completed = run_with(tmp_path, de=[f"T={path}", DE[1]])
    assert completed.returncode == 0, completed.stderr
    result = json.loads((tmp_path / "out" / "result.json").read_text())
    # The demo tables already hold 20 NA-padj rows each (independent filtering).
    assert result["dropped_per_group"] == {"T": 23, "Mono": 20}
    assert result["contrasts"][0]["n_genes"] == 977


@pytest.mark.parametrize("fixture", ["ora_out", "demo_out"])
def test_commands_sh_reproduces_run(request, tmp_path, fixture):
    out = request.getfixturevalue(fixture)
    argv = shlex.split(
        (out / "reproducibility/commands.sh").read_text().splitlines()[-1]
    )
    assert argv.count("--de") == 2
    argv[argv.index("--outdir") + 1] = str(tmp_path / "again")
    completed = subprocess.run(argv, capture_output=True, text=True, check=False)
    assert completed.returncode == 0, completed.stderr
    again = tmp_path / "again" / "tables" / "enrichment.tsv"
    assert again.read_bytes() == (out / "tables" / "enrichment.tsv").read_bytes()


def test_tsv_gene_sets_match_gmt(tmp_path, ora_out):
    rows = [
        (line.split("\t")[0], gene)
        for line in (EXAMPLES / "demo.gmt").read_text().splitlines()
        for gene in line.split("\t")[2:]
    ]
    path = tmp_path / "sets.tsv"
    pd.DataFrame(rows, columns=["set", "gene"]).to_csv(path, sep="\t", index=False)
    completed = run_with(tmp_path, gene_sets=str(path))
    assert completed.returncode == 0, completed.stderr
    pd.testing.assert_frame_equal(enrichment(tmp_path / "out"), enrichment(ora_out))


def test_output_contract(demo_out):
    paths = contract_paths()
    assert len(paths) == 15
    for rel in paths:
        assert (demo_out / rel).stat().st_size > 0, rel
    check = subprocess.run(
        ["sha256sum", "-c", "reproducibility/checksums.sha256"],
        cwd=demo_out,
        capture_output=True,
        text=True,
        check=False,
    )
    assert check.returncode == 0, check.stdout + check.stderr
    commands = (demo_out / "reproducibility" / "commands.sh").read_text()
    assert "--de T=" in commands and "--de Mono=" in commands
    assert "--seed 0" in commands
    assert "--reference 'synthetic demo, no genome build'" in commands
    assert "SerialParam" in (demo_out / "reproducibility/environment.txt").read_text()


def test_rerun_is_identical(demo_out, tmp_path):
    completed = run("--demo", "--outdir", str(tmp_path / "again"))
    assert completed.returncode == 0, completed.stderr
    check = subprocess.run(
        ["sha256sum", "-c", str(demo_out / "reproducibility" / "checksums.sha256")],
        cwd=tmp_path / "again",
        capture_output=True,
        text=True,
        check=False,
    )
    assert check.returncode == 0, check.stdout


@pytest.mark.parametrize("arg", ["--seed", "--reference"])
def test_missing_argument_fails_without_demo(tmp_path, arg):
    assert is_error(
        run_with(tmp_path, drop=[arg]), f"missing required arguments: {arg}"
    )


def test_demo_conflict(tmp_path):
    completed = run("--demo", "--limit", "2", "--outdir", str(tmp_path / "out"))
    assert is_error(completed, "--demo fixes the analysis arguments; drop: --limit")


def test_outdir_exists(tmp_path):
    (tmp_path / "out").mkdir()
    assert is_error(run_with(tmp_path), "outdir already exists")


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


@pytest.mark.parametrize(
    ("base", "drop", "changes", "message"),
    [
        (ORA, (), {"limit": "0"}, "--limit must be > 0, got 0.0"),
        (ORA, (), {"max_size": "5"}, "--max-size must be >= --min-size, got 5"),
        (FGSEA, ("--rank-col",), {}, "--method fgsea needs --rank-col"),
        (ORA, (), {"rank_col": "stat"}, "--method ora does not use --rank-col"),
        (ORA, (), {"de": ["T cells=x.tsv"]}, "--de must be name=path"),
        (ORA, (), {"de": [DE[0], DE[0]]}, "contrast 'T' is given more than once"),
        (ORA, (), {"de": ["T=/no/such.tsv"]}, "input file not found: /no/such.tsv"),
        (FGSEA, (), {"rank_col": "wald"}, "column 'wald' not found in"),
        (ORA, (), {"gene_sets": "/no/sets.gmt"}, "input file not found: /no/sets.gmt"),
        (
            ORA,
            (),
            {"padj_max": "1e-30"},
            "no contrast has a both gene with padj < 1e-30",
        ),
        (ORA, (), {"min_size": "100", "max_size": "200"}, "no gene set has 100 to 200"),
    ],
)
def test_argument_errors(tmp_path, base, drop, changes, message):
    assert is_error(run_with(tmp_path, base=base, drop=drop, **changes), message)


@pytest.mark.parametrize(
    ("key", "name", "text", "message"),
    [
        (
            "de",
            "de.tsv",
            "gene\tpadj\tlog2_fc\nA\tlow\t1\n",
            "column 'padj' in {path} must be numeric and finite",
        ),
        (
            "de",
            "de.tsv",
            "gene\tpadj\tlog2_fc\nA\t0.1\t1\nA\t0.2\t1\n",
            "duplicate gene names in {path}: A; make them unique before running",
        ),
        (
            "gene_sets",
            "sets.txt",
            "S\tx\tA\n",
            "--gene-sets must be a .gmt file or a .tsv with set and gene columns",
        ),
        (
            "gene_sets",
            "sets.gmt",
            "S\tx\tA\nS\tx\tB\n",
            "duplicate gene set names in {path}: S",
        ),
        (
            "gene_sets",
            "sets.gmt",
            "S\tx\tENSG01\n",
            "no gene of --gene-sets is in contrast 'T'",
        ),
    ],
)
def test_input_errors(tmp_path, key, name, text, message):
    path = write(tmp_path, name, text)
    value = [f"T={path}"] if key == "de" else path
    completed = run_with(tmp_path, **{key: value})
    assert is_error(completed, message.format(path=path))


def test_no_rscript(tmp_path):
    env = {"RSCRIPT": "", "PATH": str(tmp_path)}
    completed = run_with(tmp_path, base=FGSEA, env=env)
    assert is_error(completed, "Rscript not found; set RSCRIPT or put Rscript on PATH")


def test_missing_r_package(tmp_path):
    fake = Path(write(tmp_path, "Rscript", "#!/bin/sh\nexit 3\n"))
    fake.chmod(0o755)
    completed = run_with(tmp_path, base=FGSEA, env={"RSCRIPT": str(fake)})
    assert is_error(completed, f"R package fgsea is not installed for {fake}")


def test_r_failure(tmp_path):
    env = {"RSCRIPT": shutil.which("false")}
    completed = run_with(tmp_path, base=FGSEA, env=env)
    assert is_error(completed, "fgsea.R failed with exit 1")


@pytest.mark.skipif(not TEMPLATES.is_dir(), reason="no templates/ beside the skill")
@pytest.mark.parametrize("name", ["skillkit.py", "plots.py", "theme.py", "palettes.py"])
def test_vendored_copies_match_templates(name):
    vendored = (SKILL_DIR / "scripts" / name).read_bytes()
    assert vendored == (TEMPLATES / name).read_bytes()
