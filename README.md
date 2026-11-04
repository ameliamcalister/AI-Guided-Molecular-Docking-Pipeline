# AI-Guided Docking Pipeline

A generalised, AI-guided structure-based drug design pipeline. Works with any protein
target, any PDB structure, and any binding site (initially built for but not limited to MRP1). 
Parses AutoDock Vina logs, retrieves live ChEMBL reference data, auto-detects residue contacts,
builds a target-specific system prompt, calls Claude (or a local Ollama model) for scaffold
modification suggestions, validates returned SMILES with RDKit, and can auto-dock
candidates across multiple autonomous optimisation cycles providing an binding affinity output 
(kcal/mol).

<img width="1422" height="773" alt="Screenshot 2026-11-04 at 10 54 16" src="https://github.com/user-attachments/assets/3edabc83-6f75-4b8f-9bdd-f4c3c4d35848" />
<img width="1422" height="773" alt="Screenshot 2026-11-04 at 10 55 14" src="https://github.com/user-attachments/assets/01a061ba-313c-46de-a2c7-c00a15871880" />
<img width="1422" height="773" alt="Screenshot 2026-11-04 at 11 06 21" src="https://github.com/user-attachments/assets/78c6fba1-ee07-4c67-974e-75e538ecbdce" />
<img width="1422" height="773" alt="Screenshot 2026-11-04 at 11 06 27" src="https://github.com/user-attachments/assets/76323b6b-2a49-4f93-8203-78602059eef0" />
<img width="1422" height="773" alt="Screenshot 2026-11-04 at 11 06 34" src="https://github.com/user-attachments/assets/7d5bcd60-92c8-40ba-a956-7c88f47fc052" />
<img width="1422" height="773" alt="Screenshot 2026-11-04 at 11 08 21" src="https://github.com/user-attachments/assets/6dcc1792-fce8-4a20-9829-4eef7829a5ad" />



Two entry points are included:
- `docking_pipeline.py` — command-line pipeline
- `docking_app_v2.py` — Streamlit web interface (autonomous iterative loop, in-browser
  3D pose viewer via 3Dmol.js, plateau-detection auto-stop)

## Setup

```bash
pip install -r requirements.txt
```

You also need [AutoDock Vina](http://vina.scripps.edu/) and [OpenBabel](https://openbabel.org/)
installed locally, and a receptor PDBQT file for your target.

## Configuration

Set these environment variables before running (defaults assume `vina`/`obabel` are on
your `PATH`):

```bash
export ANTHROPIC_API_KEY="your-key-here"      # required — never commit this
export VINA_PATH="/path/to/vina"              # optional override
export OBABEL_PATH="/path/to/obabel"          # optional override
```

## Command-line usage

```bash
python docking_pipeline.py \
    --log compound_vina.log \
    --smiles "CCO..." \
    --receptor path/to/receptor.pdbqt \
    --center_x 45.0 --center_y 32.0 --center_z 18.0 \
    --size_x 25.0 --size_y 25.0 --size_z 25.0 \
    --target "EGFR kinase" \
    --pdb_id 1IEP \
    --cavity "ATP binding pocket, hinge region" \
    --key_residues "Thr790,Met793,Lys745,Glu762" \
    --chembl_id CHEMBL203 \
    --mode affinity \
    --chembl \
    --n 3 \
    --autodock \
    --yes
```

Example from the MRP1 project specifically:

```bash
python docking_pipeline.py \
    --log smiles4_log.txt --smiles "..." \
    --receptor 5UJA_receptor.pdbqt \
    --center_x 94.806 --center_y 59.056 --center_z 56.653 \
    --size_x 21.920 --size_y 21.297 --size_z 23.806 \
    --target "MRP1 (ABCC1)" --pdb_id 5UJA \
    --key_residues "Thr1241,Lys332,Gln450,Phe594,Met1092,Trp1245,Tyr1032" \
    --mode affinity --chembl --n 3 --autodock --yes
```

Output: `<name>_report.txt` with AI suggestions, RDKit validation, and docking results.

## Streamlit interface

```bash
streamlit run docking_app_v2.py
```

Upload a receptor PDBQT, enter target/grid details (pre-filled with MRP1/5UJA defaults),
optionally upload a Vina log and docked ligand PDBQT, and run one or more autonomous
optimisation cycles interactively, with 3D pose visualisation and per-cycle downloads.
Requires `docking_pipeline.py` in the same folder.

## Limitations
- Validated primarily on one target family (ABC transporters); testing across more diverse
  binding-site geometries would strengthen the generality claim.
- LLM-proposed modifications are checked for structural validity (RDKit) but not yet screened
  against synthesisability.
  
## Notes

- Receptor/ligand structure files (`*.pdb`, `*.pdbqt`) and generated reports are excluded
  from version control by default — see `.gitignore`.
- Built as part of an BSc final year project on computational MRP1 (ABCC1) inhibitor
  discovery at Imperial College London, generalised to support arbitrary targets.

