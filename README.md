# AI-Guided Molecular Docking Pipeline

A structure-based drug design pipeline that puts an LLM in the optimisation loop. It docks a ligand, reads the pose, asks a model for scaffold modifications, validates what comes back, and docks the result; repeating until scores plateau.

It works with any protein target, PDB structure and binding site. It was built for MRP1 (ABCC1) as part of a BSc final year project at Imperial College London, but is now generalised to any protein target.

Live app: https://ai-guided-molecular-docking-pipeline-ntbanardtcxdmqtdk6yina.streamlit.app/

Two entry points:

docking_pipeline.py — command-line pipeline
docking_app_v2.py — Streamlit web interface, with an autonomous iterative loop, an in-browser 3D pose viewer and plateau-detection auto-stop

<img width="1422" height="773" alt="Screenshot 2026-11-04 at 10 54 16" src="https://github.com/user-attachments/assets/3edabc83-6f75-4b8f-9bdd-f4c3c4d35848" />
<img width="1422" height="773" alt="Screenshot 2026-11-04 at 11 06 21" src="https://github.com/user-attachments/assets/78c6fba1-ee07-4c67-974e-75e538ecbdce" />
<img width="1422" height="773" alt="Screenshot 2026-11-04 at 11 08 21" src="https://github.com/user-attachments/assets/6dcc1792-fce8-4a20-9829-4eef7829a5ad" />

## Quick start

### With Docker

```bash
git clone https://github.com/USERNAME/REPO-NAME.git
cd REPO-NAME
docker build -t docking-app .
export ANTHROPIC_API_KEY="your-key-here"       # don't commit this
docker run -p 8501:8501 -e ANTHROPIC_API_KEY -v "$(pwd)/data:/app/data" docking-app
# open http://localhost:8501
```

The image includes AutoDock Vina and OpenBabel. `./data` keeps your run history and saved poses.

### Without Docker

```bash
pip install -r requirements.txt
```

Install [AutoDock Vina](http://vina.scripps.edu/) and [OpenBabel](https://openbabel.org/)
locally (on Ubuntu: `sudo apt install openbabel autodock-vina`), and have a receptor PDBQT file
for your target.

```bash
export ANTHROPIC_API_KEY="your-key-here"       # required, never commit this
export VINA_PATH="/path/to/vina"               # optional override
export OBABEL_PATH="/path/to/obabel"           # optional override
streamlit run docking_app_v2.py
```

Upload a receptor PDBQT, enter target and grid details (pre-filled with MRP1/5UJA defaults), optionally upload a Vina log and docked ligand PDBQT, then run one or more optimisation cycles.

# Repository layout
```bash
docking_pipeline.py      # CLI pipeline and shared logic
docking_app_v2.py        # Streamlit app (imports docking_pipeline)
provenance.py            # SQLite run store
requirements.txt
requirements-dev.txt
packages.txt             # apt packages for Streamlit Community Cloud
Dockerfile
data/
  runs.db                # run history
  poses/<run_id>/        # saved poses
docs/
  DEPLOYMENT.md          # Streamlit Community Cloud notes
NEXT_STEPS.md            # known issues
.github/workflows/ci.yml
```
docking_app_v2.py imports docking_pipeline, so the two must stay in the same directory.

# Command-line usage
```bash
python docking_pipeline.py \
    --log compound_vina.log \
    --smiles "CC(=O)Oc1ccccc1C(=O)O" \
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
Output: <name>_report.txt with the model's suggestions, RDKit validation and docking results.

## Run history and reproducibility
 
Every web-app run is saved to `data/runs.db` (SQLite) with:
 
- the target, grid, Vina settings, LLM backend and model, and a SHA-256 of the receptor file
- Python and package versions, the Vina version, and the git commit
- start time, end time and status (`completed`, `failed` or `cancelled`)
- one row per compound per cycle: SMILES, Vina score, and the pose in `data/poses/<run_id>/`
```python
from provenance import RunStore
 
store = RunStore("data/runs.db")
run_id = store.list_runs()[0]["run_id"]
print(store.get_run(run_id)["params"])
print(store.best_results(run_id, n=5))
```
 
**What "Fix random seed" does:** it passes `--seed 42` to Vina, making the search repeatable for a
given prepared ligand file. OpenBabel generates each ligand's 3D structure from SMILES and is not
seeded, so scores for the same SMILES can still differ slightly between runs. See `NEXT_STEPS.md`.
 
## Deployment
 
The app shells out to `vina` and `obabel`, so it cannot be hosted anywhere that serves only static
files. It runs on Streamlit Community Cloud, which installs the Debian packages listed in
`packages.txt` alongside the pip requirements.
 
One trap worth knowing: `packages.txt` must contain bare package names and nothing else. Community
Cloud passes the file's words straight to `apt-get install` without stripping comments, so a `#`
line mentioning a flag such as `--ligand` makes apt fail with `Command line option --ligand is not
understood` and the whole install aborts.
 
Full setup, secrets configuration and free-tier constraints (2 CPU cores, non-persistent disk,
sleeps after 12 hours idle) are in [`docs/DEPLOYMENT.md`](docs/DEPLOYMENT.md). `Dockerfile` is the
route for self-hosting with a persistent volume.

 
## Tests
 
```bash
pip install -r requirements-dev.txt
pytest
```
 
CI runs on every push (`.github/workflows/ci.yml`). Tests cover log and PDBQT parsing, contact
detection, SMILES validation, ADMET thresholds, prompt building, iteration-mode logic, ChEMBL and
LLM calls (faked, so no key or network needed), Vina command construction, the auto-grid helper,
the run store, and real Vina/OpenBabel docking when those tools are installed. Three tests are
marked `xfail` deliberately: they document known issues in `NEXT_STEPS.md`.
 
## Limitations
 
- **Docking scores are not binding affinities.** Vina's scoring function has a reported standard
  error of roughly 2–3 kcal/mol, which exceeds most cycle-to-cycle improvements this pipeline
  reports. Treat the output as a ranking heuristic for prioritising compounds, not a prediction.
- **Validated primarily on one target family** (ABC transporters). Testing across more diverse
  binding-site geometries would strengthen the generality claim.
- **No synthesisability screening.** Proposed modifications are checked for structural validity
  (RDKit) and ADMET limits, but nothing checks whether a molecule can actually be made.
- **Rigid receptor.** Side-chain and backbone flexibility are not modelled, so induced-fit effects
  are missed.
- **The LLM has no structural model.** It reasons over a text description of contacts and geometry,
  not the coordinates themselves.
