#!/usr/bin/env python3
"""Download receptor structures from the RCSB PDB for the viewer eval.

For each PDB id in targets.tsv this writes:

    targets/<id>.pdb       receptor atoms (protein/nucleic acid + metals and
                           cofactors, with waters and the native ligand removed)
    targets/sites.json     where the binding site is, per target

The site is the centroid of the structure's own co-crystallised ligand where
there is one, so eval ligands land in a real pocket rather than an arbitrary
point. Targets without a bound ligand fall back to the centre of the structure.

Downloaded files are treated as untrusted input: they are parsed by column, and
nothing from them is executed.

Usage:
    python eval/make_targets.py [--out eval/targets] [--force]
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
RCSB = "https://files.rcsb.org/download/{pdb_id}.pdb"

# Kept as part of the receptor rather than treated as the bound ligand.
SOLVENT = {"HOH", "WAT", "DOD"}
IONS_COFACTORS = {
    "ZN", "MG", "CA", "FE", "MN", "CU", "NA", "K", "CL", "SO4", "PO4",
    "HEM", "HEC", "NAD", "NAP", "FAD", "FMN", "SF4", "FES",
}
# Crystallisation additives: not the biological ligand, not part of the receptor.
BUFFERS = {
    "GOL", "EDO", "PEG", "PG4", "MPD", "DMS", "ACT", "ACY", "TRS", "MES",
    "EPE", "IMD", "FMT", "NO3", "AZI", "BME", "CIT", "TLA", "SCN", "IOD",
    "BR", "MRD", "P6G", "1PE", "PGE", "OLC", "OLA", "PLM", "CLR", "NH2",
}


def fetch(pdb_id: str, dest: Path, force: bool) -> str | None:
    raw = dest / f"{pdb_id}_raw.pdb"
    if raw.exists() and not force and raw.stat().st_size > 0:
        return raw.read_text(errors="replace")
    url = RCSB.format(pdb_id=pdb_id)
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            text = r.read().decode("utf-8", errors="replace")
    except (urllib.error.URLError, TimeoutError, OSError) as e:
        print(f"  ! {pdb_id}: download failed ({e})")
        return None
    raw.write_text(text)
    return text


def first_model(lines: list[str]) -> list[str]:
    """NMR entries hold many models; the viewer only ever shows one."""
    out, started = [], False
    for l in lines:
        if l.startswith("MODEL"):
            if started:
                break
            started = True
            continue
        if l.startswith("ENDMDL"):
            break
        out.append(l)
    return out


def xyz(l: str):
    try:
        return float(l[30:38]), float(l[38:46]), float(l[46:54])
    except (ValueError, IndexError):
        return None


def centroid(pts):
    n = len(pts) or 1
    return [round(sum(p[i] for p in pts) / n, 3) for i in range(3)]


def split_structure(text: str):
    """Return (receptor lines, native ligand lines, ligand resname)."""
    lines = first_model(text.splitlines())
    receptor, het_groups = [], {}

    for l in lines:
        if l.startswith("ATOM"):
            # Protein and nucleic acid backbone: always receptor.
            receptor.append(l)
        elif l.startswith("HETATM"):
            resn = l[17:20].strip().upper()
            if resn in SOLVENT:
                continue
            if resn in IONS_COFACTORS:
                receptor.append(l)
                continue
            if resn in BUFFERS:
                continue
            key = (resn, l[21:22], l[22:27])
            het_groups.setdefault(key, []).append(l)

    # The biggest remaining HETATM group is the bound ligand.
    ligand_lines, ligand_resn = [], ""
    if het_groups:
        key, ligand_lines = max(het_groups.items(), key=lambda kv: len(kv[1]))
        if len(ligand_lines) < 6:        # too small to be a real ligand
            ligand_lines, ligand_resn = [], ""
        else:
            ligand_resn = key[0]

    return receptor, ligand_lines, ligand_resn


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "targets"))
    ap.add_argument("--targets", default=str(HERE / "targets.tsv"))
    ap.add_argument("--force", action="store_true", help="re-download")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    for line in Path(args.targets).read_text().strip().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        pdb_id, _, desc = line.partition("\t")
        rows.append((pdb_id.strip().upper(), desc.strip()))

    sites = {}
    for pdb_id, desc in rows:
        text = fetch(pdb_id, out, args.force)
        if not text:
            continue
        receptor, ligand, ligand_resn = split_structure(text)
        if not receptor:
            print(f"  ! {pdb_id}: no receptor atoms found")
            continue

        rec_pts = [p for p in (xyz(l) for l in receptor) if p]
        if ligand:
            lig_pts = [p for p in (xyz(l) for l in ligand) if p]
            site, source = centroid(lig_pts), f"ligand {ligand_resn}"
        else:
            site, source = centroid(rec_pts), "structure centre"

        (out / f"{pdb_id}.pdb").write_text("\n".join(receptor) + "\n")
        chains = sorted({l[21:22] for l in receptor if l[21:22].strip()})
        sites[pdb_id] = {
            "description": desc,
            "site": site,
            "site_from": source,
            "atoms": len(receptor),
            "chains": chains,
        }
        print(f"  + {pdb_id}: {len(receptor):>6} atoms, "
              f"{len(chains)} chain(s), site from {source}")

    (out / "sites.json").write_text(json.dumps(sites, indent=2))
    # <id>_raw.pdb is kept as a download cache; --force re-fetches.
    print(f"\n{len(sites)} targets in {out}")
    return 0 if sites else 1


if __name__ == "__main__":
    raise SystemExit(main())
