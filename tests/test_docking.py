"""Vina and OpenBabel command construction (mocked), plus real end-to-end checks when the tools are installed."""
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest

import docking_pipeline as dp

GRID = {"center_x": 1.0, "center_y": 2.0, "center_z": 3.0, "size_x": 20, "size_y": 20, "size_z": 20}
HAVE_TOOLS = bool(shutil.which("vina") and shutil.which("obabel"))


@pytest.fixture
def fake_vina(monkeypatch):
    """Replace subprocess.run with a fake Vina that records the command and writes an output file."""
    state = {"cmd": None, "write_output": True}

    def fake_run(cmd, **kwargs):
        state["cmd"] = cmd
        if state["write_output"]:
            Path(cmd[cmd.index("--out") + 1]).write_text("MODEL 1\nENDMDL\n")
        return SimpleNamespace(stdout="   1       -7.512          0          0\n", stderr="")

    monkeypatch.setattr(dp, "find_tool", lambda path, name: name)
    monkeypatch.setattr(dp.subprocess, "run", fake_run)
    return state


def test_run_vina_builds_expected_command_and_parses_affinity(fake_vina, tmp_path):
    lig = tmp_path / "lig.pdbqt"
    lig.write_text("x")
    affinity, out = dp.run_vina(lig, "c1", tmp_path, "rec.pdbqt", GRID, exhaustiveness=16)
    cmd = fake_vina["cmd"]
    assert affinity == -7.512 and out == tmp_path / "c1_docked.pdbqt"
    assert cmd[cmd.index("--exhaustiveness") + 1] == "16"
    assert cmd[cmd.index("--center_y") + 1] == "2.0"
    assert cmd[cmd.index("--size_z") + 1] == "20"
    assert (tmp_path / "c1_vina.log").exists()


def test_run_vina_passes_seed_only_when_given(fake_vina, tmp_path):
    lig = tmp_path / "lig.pdbqt"
    lig.write_text("x")
    dp.run_vina(lig, "c1", tmp_path, "rec.pdbqt", GRID)
    assert "--seed" not in fake_vina["cmd"]
    dp.run_vina(lig, "c1", tmp_path, "rec.pdbqt", GRID, seed=42)
    cmd = fake_vina["cmd"]
    assert cmd[cmd.index("--seed") + 1] == "42"


def test_run_vina_failure_returns_none_pair(fake_vina, tmp_path):
    fake_vina["write_output"] = False
    lig = tmp_path / "lig.pdbqt"
    lig.write_text("x")
    assert dp.run_vina(lig, "c1", tmp_path, "rec.pdbqt", GRID) == (None, None)


def test_smiles_to_pdbqt_returns_none_when_obabel_makes_nothing(monkeypatch, tmp_path):
    monkeypatch.setattr(dp, "find_tool", lambda path, name: name)
    monkeypatch.setattr(dp.subprocess, "run", lambda *a, **k: SimpleNamespace(stdout="", stderr=""))
    assert dp.smiles_to_pdbqt("CCO", "x", tmp_path) is None


def test_find_tool_raises_when_missing(monkeypatch):
    monkeypatch.setattr(dp.shutil, "which", lambda name: None)
    with pytest.raises(FileNotFoundError, match="Cannot find"):
        dp.find_tool("/definitely/not/here", "vina")


# ── Real tools ───────────────────────────────────────────────────────────────
@pytest.fixture(scope="module")
def tiny_receptor(tmp_path_factory):
    """A small rigid 'receptor' made from a short peptide, so Vina has something to dock into."""
    import subprocess
    d = tmp_path_factory.mktemp("receptor")
    (d / "pep.smi").write_text("NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)O")
    subprocess.run(["obabel", str(d / "pep.smi"), "-O", str(d / "pep.pdb"), "--gen3d", "-h"],
                   capture_output=True)
    subprocess.run(["obabel", str(d / "pep.pdb"), "-xr", "-O", str(d / "receptor.pdbqt")],
                   capture_output=True)
    assert (d / "receptor.pdbqt").stat().st_size > 0
    return str(d / "receptor.pdbqt")


@pytest.fixture
def tiny_grid(tiny_receptor):
    pts = [(float(l[30:38]), float(l[38:46]), float(l[46:54]))
           for l in open(tiny_receptor) if l.startswith(("ATOM", "HETATM"))]
    c = [round((max(a) + min(a)) / 2, 2) for a in zip(*pts)]
    return {"center_x": c[0], "center_y": c[1], "center_z": c[2], "size_x": 22, "size_y": 22, "size_z": 22}


@pytest.mark.skipif(not HAVE_TOOLS, reason="needs vina and obabel on PATH")
def test_real_ligand_preparation_makes_a_pdbqt(tmp_path):
    p = dp.smiles_to_pdbqt("CCO", "eth", tmp_path)
    assert p is not None and "ROOT" in p.read_text()


@pytest.mark.skipif(not HAVE_TOOLS, reason="needs vina and obabel on PATH")
def test_real_docking_scores_match_the_log_the_parser_reads(tmp_path, tiny_receptor, tiny_grid):
    lig = dp.smiles_to_pdbqt("CCO", "eth", tmp_path)
    affinity, docked = dp.run_vina(lig, "eth", tmp_path, tiny_receptor, tiny_grid, exhaustiveness=8, seed=42)
    parsed = dp.parse_vina_log(str(tmp_path / "eth_vina.log"))
    assert docked.exists() and affinity is not None and affinity < 0
    assert parsed["best_affinity"] == affinity
    assert parsed["ligand_name"] == "eth"


@pytest.mark.skipif(not HAVE_TOOLS, reason="needs vina and obabel on PATH")
def test_real_docking_same_prepared_ligand_and_seed_gives_same_score(tmp_path, tiny_receptor, tiny_grid):
    """The seed makes Vina's search repeatable for ONE prepared ligand file.
    (Preparing the ligand again from SMILES is a separate, unseeded step: see NEXT_STEPS.md.)"""
    lig = dp.smiles_to_pdbqt("CCO", "eth", tmp_path)
    scores = []
    for k in range(3):
        work = tmp_path / f"run{k}"
        work.mkdir()
        scores.append(dp.run_vina(lig, "eth", work, tiny_receptor, tiny_grid, 8, seed=42)[0])
    assert len(set(scores)) == 1
