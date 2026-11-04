#!/usr/bin/env python3
"""Build a corpus of ligand PDBQT files for the viewer eval.

Each compound in compounds.tsv is turned into 3D coordinates with Open Babel,
then translated so it sits inside a real receptor pocket and wrapped in
MODEL/ENDMDL blocks with a VINA RESULT remark, so the files look like genuine
Vina output (several poses, best pose first) without needing to run docking.

Usage:
    python eval/make_corpus.py --receptor path/to/receptor.pdbqt [--out eval/corpus]
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def atom_lines(text: str) -> list[str]:
    return [l for l in text.splitlines() if l.startswith(("ATOM", "HETATM"))]


def coords(lines: list[str]) -> list[tuple[float, float, float]]:
    out = []
    for l in lines:
        try:
            out.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
        except (ValueError, IndexError):
            pass
    return out


def centroid(pts):
    n = len(pts) or 1
    return (sum(p[0] for p in pts) / n,
            sum(p[1] for p in pts) / n,
            sum(p[2] for p in pts) / n)


def shift(lines: list[str], dx: float, dy: float, dz: float) -> list[str]:
    """Translate the x/y/z columns of PDBQT atom lines, keeping fixed widths."""
    out = []
    for l in lines:
        try:
            x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
        except (ValueError, IndexError):
            out.append(l)
            continue
        out.append(f"{l[:30]}{x + dx:8.3f}{y + dy:8.3f}{z + dz:8.3f}{l[54:]}")
    return out


def smiles_to_atoms(name: str, smiles: str, work: Path) -> list[str] | None:
    obabel = shutil.which("obabel")
    if not obabel:
        sys.exit("obabel not found on PATH")
    out = work / f"{name}.pdbqt"
    proc = subprocess.run(
        [obabel, f"-:{smiles}", "-opdbqt", "-O", str(out),
         "--gen3d", "-h", "--partialcharge", "gasteiger"],
        capture_output=True, text=True, timeout=180,
    )
    if not out.exists() or out.stat().st_size == 0:
        print(f"  ! {name}: obabel failed: {proc.stderr.strip()[:120]}")
        return None
    return atom_lines(out.read_text())


def as_docked(lines: list[str], scores: list[float]) -> str:
    """Wrap atom lines as a multi-pose docked PDBQT, best pose first."""
    blocks = []
    for i, score in enumerate(scores, 1):
        # Later poses are nudged, the way alternative Vina modes differ.
        posed = shift(lines, (i - 1) * 1.4, (i - 1) * -0.9, (i - 1) * 0.7)
        blocks.append("\n".join([
            f"MODEL {i}",
            f"REMARK VINA RESULT:    {score:.3f}      0.000      0.000",
            *posed,
            "ENDMDL",
        ]))
    return "\n".join(blocks) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--out", default=str(HERE / "corpus"))
    ap.add_argument("--compounds", default=str(HERE / "compounds.tsv"))
    args = ap.parse_args()

    rec_lines = atom_lines(Path(args.receptor).read_text())
    if not rec_lines:
        sys.exit(f"no atoms in receptor {args.receptor}")
    site = centroid(coords(rec_lines))
    print(f"receptor: {len(rec_lines)} atoms, pocket target {site[0]:.1f} {site[1]:.1f} {site[2]:.1f}")

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    work = out_dir / ".work"
    work.mkdir(exist_ok=True)

    rows = [l.split("\t") for l in
            Path(args.compounds).read_text().strip().splitlines() if l.strip()]
    made = 0
    for name, smiles in rows:
        lines = smiles_to_atoms(name, smiles, work)
        if not lines:
            continue
        cx, cy, cz = centroid(coords(lines))
        posed = shift(lines, site[0] - cx, site[1] - cy, site[2] - cz)
        scores = [-9.4, -8.7, -8.1]
        (out_dir / f"{name}_docked.pdbqt").write_text(as_docked(posed, scores))
        print(f"  + {name}: {len(lines)} atoms")
        made += 1

    # One compound deliberately parked far outside the protein: the viewer
    # should still draw it, and report zero contact residues.
    lines = smiles_to_atoms("offsite_benzene", "c1ccccc1", work)
    if lines:
        cx, cy, cz = centroid(coords(lines))
        posed = shift(lines, site[0] - cx + 200, site[1] - cy, site[2] - cz)
        (out_dir / "offsite_benzene_docked.pdbqt").write_text(as_docked(posed, [-3.1]))
        print("  + offsite_benzene (far from receptor)")
        made += 1

    shutil.rmtree(work, ignore_errors=True)
    print(f"\n{made} ligand files in {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
