#!/usr/bin/env python3
"""Eval for render_3dmol: does the viewer hold up across many ligands?

Two layers:

1. Python-side checks (this script). For every case, render_3dmol is called and
   the returned HTML is checked for things that would silently break the viewer:
   unreplaced placeholders, invalid embedded JSON, a ligand that lost all its
   atoms, nonsense contact counts, or ligand text that escapes the <script>.

2. Browser-side checks (run separately). Each case is also written to
   eval/out/<case>.html as a standalone page. Loading one sets
   window.__v3dReady = {ligAtoms, recAtoms, contacts} once 3Dmol has drawn it,
   so a browser can confirm the pose really rendered rather than just that the
   HTML looked plausible. eval/out/manifest.json lists every case for that pass.

Usage:
    python eval/eval_viewer.py --receptor path/to/receptor.pdbqt
"""
from __future__ import annotations

import argparse
import ast
import json
import re
import sys
from collections import Counter
from pathlib import Path

# Running under python -I drops the script's own directory from sys.path.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_corpus as mc  # noqa: E402  (same directory: geometry helpers)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def load_render_3dmol(app_path: Path):
    """Load render_3dmol alone, without running the Streamlit script."""
    tree = ast.parse(app_path.read_text(encoding="utf-8"))
    node = next(n for n in tree.body
                if isinstance(n, ast.FunctionDef) and n.name == "render_3dmol")
    ns: dict = {}
    exec(compile(ast.Module([node], []), str(app_path), "exec"), ns)
    return ns["render_3dmol"]


def atom_lines(text: str) -> list[str]:
    return [l for l in text.splitlines() if l.startswith(("ATOM", "HETATM"))]


# ── cases ────────────────────────────────────────────────────────────────────

AD_ELEM = {
    "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N", "O": "O", "OA": "O",
    "OS": "O", "S": "S", "SA": "S", "H": "H", "HD": "H", "HS": "H", "F": "F",
    "P": "P", "CL": "Cl", "BR": "Br", "I": "I", "MG": "Mg", "MN": "Mn",
    "ZN": "Zn", "CA": "Ca", "FE": "Fe", "CU": "Cu", "SE": "Se",
}


def reference_elements(text: str, fmt: str) -> Counter:
    """Work out each atom's element independently of the viewer's own code.

    Deliberately a separate, simpler implementation: if this and the viewer
    agree, the viewer is reading the right columns for the format. PDBQT keeps
    its atom type in 78-79, plain PDB keeps the element in 77-78.
    """
    counts: Counter = Counter()
    for l in atom_lines(text):
        if fmt == "pdbqt":
            counts[AD_ELEM.get(l[77:79].strip().upper(), "?")] += 1
        else:
            counts[(l[76:78].strip() or "?").capitalize()] += 1
    return counts


def place_in_site(ligand_pdbqt: str, site: list[float]) -> str:
    """Move a corpus ligand so its best pose sits in this target's pocket."""
    lines = atom_lines(ligand_pdbqt.split("ENDMDL")[0])
    pts = mc.coords(lines)
    if not pts:
        return ligand_pdbqt
    cx, cy, cz = mc.centroid(pts)
    moved = mc.shift(lines, site[0] - cx, site[1] - cy, site[2] - cz)
    return "MODEL 1\nREMARK VINA RESULT:    -8.500      0.000      0.000\n" \
           + "\n".join(moved) + "\nENDMDL\n"


def build_target_cases(targets_dir: Path, per_target: int = 3) -> list[dict]:
    """Every downloaded target, each with a few different ligands in its pocket.

    The ligands rotate through the corpus so that across the targets a wide
    range of chemistry is seen, without rendering every pair.
    """
    sites_file = targets_dir / "sites.json"
    if not sites_file.exists():
        return []
    sites = json.loads(sites_file.read_text())

    corpus = sorted(f for f in (HERE / "corpus").glob("*_docked.pdbqt")
                    if not f.name.startswith("offsite"))
    if not corpus:
        return []

    cases, k = [], 0
    for pdb_id in sorted(sites):
        meta = sites[pdb_id]
        rec_file = targets_dir / f"{pdb_id}.pdb"
        if not rec_file.exists():
            continue
        receptor = rec_file.read_text()
        for _ in range(per_target):
            lig_file = corpus[k % len(corpus)]
            k += 1
            lig_name = lig_file.name.replace("_docked.pdbqt", "")
            cases.append({
                "name": f"target:{pdb_id}/{lig_name}",
                "ligand": place_in_site(lig_file.read_text(), meta["site"]),
                "receptor": receptor,
                "receptor_fmt": "pdb",
                "expect_atoms": True,
                # Placed at the pocket centre, so there must be contacts.
                "expect_contacts": "some",
                "note": meta["description"],
            })
    return cases


def build_cases(receptor: str) -> list[dict]:
    """Each case: name, ligand text, receptor text, and what we expect."""
    cases: list[dict] = []

    # Generated ligands placed in the pocket (see make_corpus.py).
    for f in sorted((HERE / "corpus").glob("*_docked.pdbqt")):
        name = f.name.replace("_docked.pdbqt", "")
        offsite = name.startswith("offsite")
        cases.append({
            "name": f"ligand:{name}",
            "ligand": f.read_text(),
            "receptor": receptor,
            "expect_atoms": True,
            "expect_contacts": "zero" if offsite else "some",
        })

    # Real Vina output committed in the repo.
    for f in sorted((ROOT / "data" / "poses").rglob("*.pdbqt")):
        cases.append({
            "name": f"real:{f.parent.name}/{f.stem}",
            "ligand": f.read_text(),
            "receptor": receptor,
            "expect_atoms": True,
            "expect_contacts": "any",  # these were docked to another receptor
        })

    # Edge cases: the viewer must degrade gracefully, never crash the page.
    one = atom_lines((HERE / "corpus" / "aspirin_docked.pdbqt").read_text())
    aspirin_plain = "\n".join(one)

    cases += [
        {"name": "edge:empty-ligand", "ligand": "", "receptor": receptor,
         "expect_atoms": False, "expect_contacts": "zero"},
        {"name": "edge:no-model-block", "ligand": aspirin_plain, "receptor": receptor,
         "expect_atoms": True, "expect_contacts": "some"},
        {"name": "edge:no-receptor", "ligand": aspirin_plain, "receptor": "",
         "expect_atoms": True, "expect_contacts": "zero"},
        {"name": "edge:remarks-only", "ligand": "REMARK nothing here\nMODEL 1\nENDMDL\n",
         "receptor": receptor, "expect_atoms": False, "expect_contacts": "zero"},
        {"name": "edge:truncated-coords",
         "ligand": "MODEL 1\nATOM      1  C   LIG A   1      BAD\nENDMDL\n",
         "receptor": receptor, "expect_atoms": True, "expect_contacts": "zero"},
        {"name": "edge:blank-both", "ligand": "", "receptor": "",
         "expect_atoms": False, "expect_contacts": "zero"},
        # A </script> inside the data must not break out of the script tag.
        # It has to sit inside an ATOM line within the first 78 columns,
        # because render_3dmol keeps only ATOM/HETATM lines and truncates them.
        {"name": "edge:script-injection",
         "ligand": ("MODEL 1\n"
                    + "ATOM      1 </script><h1>x LIG A   1      "
                    + "  19.600  24.100  29.900  1.00  0.00    +0.000 C\n"
                    + aspirin_plain + "\nENDMDL\n"),
         "receptor": receptor, "expect_atoms": True, "expect_contacts": "some",
         # The injected line is deliberately malformed, so 3Dmol drops it.
         # Python still counts it, hence one fewer atom in the browser.
         "browser_atom_delta": -1},
        # Backslashes and quotes in the text must survive JSON embedding.
        {"name": "edge:quotes-backslash",
         "ligand": 'MODEL 1\n' + aspirin_plain + '\nREMARK "quoted" C:\\path\\x\nENDMDL\n',
         "receptor": receptor, "expect_atoms": True, "expect_contacts": "some"},
    ]
    return cases


# ── python-side checks ───────────────────────────────────────────────────────

def extract(html: str, var: str):
    """Pull a JSON literal back out of the generated page."""
    m = re.search(rf"var {var}\s*=\s*(.*?);\n", html, re.S)
    if not m:
        return None, f"{var} assignment not found"
    try:
        return json.loads(m.group(1)), None
    except json.JSONDecodeError as e:
        return None, f"{var} is not valid JSON: {e}"


def check(case: dict, html: str) -> list[str]:
    problems = []

    for ph in ("__H__", "__CUT__", "__REC__", "__LIG__", "__CON__"):
        if ph in html:
            problems.append(f"placeholder {ph} left in output")

    lig, err = extract(html, "LIGAND")
    if err:
        problems.append(err)
    rec, err = extract(html, "RECEPTOR")
    if err:
        problems.append(err)
    con, err = extract(html, "CONTACTS")
    if err:
        problems.append(err)

    # The data must not terminate the script element early. The page has
    # exactly two </script> of its own (the CDN tag and the inline block), so
    # any extra one came out of the embedded ligand/receptor text.
    if html.count("</script>") != 2:
        problems.append(f"{html.count('</script>')} </script> in page, expected 2: "
                        "embedded data breaks out of the script tag")
    # Checked on the raw page, not the decoded JSON: after decoding, a '<'
    # that was correctly escaped as \\u003c reappears as a plain '<', which is
    # fine. What matters is that the literal text sent to the browser is inert.
    for label in ("LIGAND", "RECEPTOR"):
        m = re.search(rf"var {label}\s*=\s*(.*?);\n", html, re.S)
        if m and "</script" in m.group(1).lower():
            problems.append(f"{label} assignment holds a raw </script")

    # Elements must survive the rewrite. A receptor read with the wrong
    # format's columns still passes every other check: the atom counts match,
    # it just quietly turns calcium into carbon.
    if rec and case["receptor"]:
        want = reference_elements(case["receptor"], case.get("receptor_fmt", "pdbqt"))
        got = reference_elements(rec, "pdb")
        for elem in sorted(set(want) | set(got)):
            if elem == "?":
                continue
            if want.get(elem, 0) != got.get(elem, 0):
                problems.append(
                    f"receptor element {elem}: {got.get(elem, 0)} in output, "
                    f"{want.get(elem, 0)} in source")
        if got.get("?", 0):
            problems.append(f"{got['?']} receptor atoms have no element")

    n_lig = len(atom_lines(lig or ""))
    if case["expect_atoms"] and n_lig == 0:
        problems.append("expected ligand atoms, got none")
    if not case["expect_atoms"] and n_lig > 0:
        problems.append(f"expected no ligand atoms, got {n_lig}")

    # Only the best pose (first MODEL) should be shown.
    src_models = case["ligand"].count("\nMODEL") + case["ligand"].startswith("MODEL")
    if src_models > 1 and n_lig:
        per_pose = len(atom_lines(case["ligand"])) / src_models
        if n_lig > per_pose * 1.5:
            problems.append(f"shows {n_lig} atoms, more than one pose ({per_pose:.0f})")

    if con is not None:
        want = case["expect_contacts"]
        if want == "some" and len(con) == 0:
            problems.append("expected contact residues, got none")
        if want == "zero" and len(con) > 0:
            problems.append(f"expected no contact residues, got {len(con)}")
        for c in con:
            if not str(c.get("resi", "")).strip():
                problems.append("contact with empty resi")
                break

    if rec is not None and case["receptor"] and len(atom_lines(rec)) == 0:
        problems.append("receptor given but no atoms survived")

    if "3Dmol-min.js" not in html:
        problems.append("3Dmol script tag missing")

    return problems


PAGE = """<!doctype html>
<html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>viewer eval: __NAME__</title>
<style>body{margin:0;background:#111;font-family:system-ui,sans-serif}
h1{color:#94a3b8;font:600 13px system-ui;margin:8px 12px}</style>
</head><body><h1>__NAME__</h1>__FRAGMENT__</body></html>
"""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True,
                    help="PDBQT receptor for the ligand sweep and edge cases")
    ap.add_argument("--targets", default=str(HERE / "targets"),
                    help="directory from make_targets.py (skipped if missing)")
    ap.add_argument("--per-target", type=int, default=3,
                    help="how many different ligands to try per target")
    ap.add_argument("--out", default=str(HERE / "out"))
    args = ap.parse_args()

    render = load_render_3dmol(ROOT / "docking_app_v2.py")
    receptor = Path(args.receptor).read_text()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    for old in out_dir.glob("*.html"):
        old.unlink()

    cases = build_cases(receptor)
    target_cases = build_target_cases(Path(args.targets), args.per_target)
    if target_cases:
        cases += target_cases
    else:
        print("no targets found; run make_targets.py for the multi-target pass\n")
    manifest, failed = [], 0

    print(f"{len(cases)} cases\n" + "-" * 62)
    for case in cases:
        try:
            html = render(case["ligand"], case["receptor"])
        except Exception as e:                                  # noqa: BLE001
            print(f"FAIL  {case['name']}: raised {type(e).__name__}: {e}")
            failed += 1
            manifest.append({"name": case["name"], "python": "raised"})
            continue

        problems = check(case, html)
        lig, _ = extract(html, "LIGAND")
        con, _ = extract(html, "CONTACTS")
        n_lig = len(atom_lines(lig or ""))
        n_con = len(con or [])

        slug = re.sub(r"[^a-z0-9]+", "_", case["name"].lower()).strip("_")
        (out_dir / f"{slug}.html").write_text(
            PAGE.replace("__NAME__", case["name"]).replace("__FRAGMENT__", html))

        manifest.append({
            "name": case["name"], "file": f"{slug}.html",
            "lig_atoms": n_lig, "contacts": n_con,
            "browser_atoms": n_lig + case.get("browser_atom_delta", 0),
            "expect_atoms": case["expect_atoms"],
            "python": "ok" if not problems else "problems",
            "kb": round(len(html) / 1024, 1),
        })

        if problems:
            failed += 1
            print(f"FAIL  {case['name']}")
            for p in problems:
                print(f"        - {p}")
        else:
            print(f"ok    {case['name']:<34} {n_lig:>3} atoms  "
                  f"{n_con:>3} contacts  {len(html)/1024:>6.1f} kB")

    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print("-" * 62)
    print(f"{len(cases) - failed}/{len(cases)} passed python checks")
    print(f"pages for the browser pass: {out_dir}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
