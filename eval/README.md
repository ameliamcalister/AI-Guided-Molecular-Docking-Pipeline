# Viewer eval

Checks that `render_3dmol` in `docking_app_v2.py` holds up across many different
ligands, not just the one compound that happened to be on screen.

There are two passes, because each catches things the other cannot. The Python
pass sees the HTML that gets generated; only a real browser can tell you whether
3Dmol actually drew the pose.

## 1. Build the ligands

Turns each SMILES in `compounds.tsv` into 3D coordinates with Open Babel, then
places it inside a real receptor pocket and wraps it in `MODEL`/`ENDMDL` blocks
with a `VINA RESULT` remark, so the files look like genuine Vina output without
needing to run docking.

```bash
python eval/make_corpus.py --receptor /path/to/receptor.pdbqt
```

Writes `eval/corpus/*.pdbqt`. Add compounds by adding `name<TAB>SMILES` rows to
`compounds.tsv`. One compound, `offsite_benzene`, is parked far outside the
protein on purpose: the viewer should still draw it and report zero contacts.

## 2. Download the targets

The ligand corpus only exercises one receptor. This downloads ten different
structures from the RCSB PDB so the viewer is tested against a range of
proteins, not just MRP1:

```bash
python eval/make_targets.py
```

Edit `targets.tsv` to change the list. The current set is chosen to vary the
things the viewer has to cope with:

| PDB  | Why it is in the set                                  |
|------|-------------------------------------------------------|
| 1HSG | HIV-1 protease, a homodimer and the classic benchmark |
| 3PTB | Trypsin, with a bound calcium ion                     |
| 4DFR | Dihydrofolate reductase, two chains                   |
| 6LU7 | SARS-CoV-2 main protease                              |
| 2RH1 | A membrane GPCR                                       |
| 1AKE | Adenylate kinase, no drug-like ligand                 |
| 4HHB | Haemoglobin: four chains, four haem irons             |
| 1CRN | Crambin, only 46 residues                             |
| 1BNA | B-DNA, nucleic acid rather than protein               |
| 5UJA | MRP1 itself, a large ABC transporter                  |

Each structure is split into receptor and native ligand. Where there is a
co-crystallised ligand its centroid becomes the binding site, so eval ligands
land in a real pocket rather than an arbitrary point.

Note these are plain PDB files, while the pipeline's own receptors are PDBQT.
That difference matters: the two formats put the element in different columns.

## 3. Python pass

```bash
python eval/eval_viewer.py --receptor /path/to/receptor.pdbqt
```

For every case it calls `render_3dmol` and checks the generated page for things
that would break the viewer silently: unreplaced placeholders, invalid embedded
JSON, a ligand that lost all its atoms, more atoms than one pose, nonsense
contact counts, and ligand text escaping the `<script>` block.

It covers the generated ligands against the PDBQT receptor, the real Vina
output in `data/poses/`, every downloaded target with a rotating selection of
ligands in its pocket, and edge cases: empty ligand, no `MODEL` block, no
receptor, remarks only, truncated coordinates, and deliberately hostile
content (`</script>`, quotes, backslashes).

For each receptor it also recounts the elements with a separate, simpler
parser and compares. A receptor read with the wrong format's columns passes
every other check — the atom counts match, it just quietly turns calcium into
carbon — so this is the only check that catches it.

It also writes each case to `eval/out/<case>.html` for the browser pass, plus
`eval/out/manifest.json` with the counts Python computed.

## 4. Browser pass

The generated pages set `window.__v3dReady = {ligAtoms, recAtoms, contacts}`
once 3Dmol has finished drawing, so a browser can confirm the scene really
rendered and matches what Python expected.

```bash
cd eval && python -m http.server 8765
```

Then open <http://127.0.0.1:8765/browser_check.html>. It loads every page in
turn and prints a pass/fail table. It has to be served over HTTP; opening it as
a `file://` URL will not work.

## What this has caught

- The whole function returned `None`, so the viewer never rendered at all.
- `json.dumps` leaves `<` alone, so a `</script>` in a PDBQT line closed the
  script element early and broke the page.
- Truncating PDBQT lines to 78 columns keeps the AutoDock type (`A`, `OA`,
  `HD`) in the column PDB uses for the element, so 3Dmol read `A` as an
  element and dropped atoms.
- 3Dmol's PDB parser drops hydrogens unless you pass `keepH`, silently losing
  the polar H that make H-bonds worth looking at.
- Reading the AutoDock type column on a plain PDB file turned trypsin's
  calcium ion and haemoglobin's haem irons into carbon, because PDB keeps the
  element in columns 77-78 and PDBQT keeps its own type in 78-79. Found only
  once the eval covered targets other than MRP1.
