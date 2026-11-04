import ast
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def pdbqt_line(serial, resname, resnum, x, y, z, name="CA", record="ATOM  "):
    """One fixed-width PDBQT atom line, the way the pipeline's parser expects it."""
    return (f"{record}{serial:5d} {name:<4s} {resname:>3s} A{resnum:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00    +0.000 C ")


def write_pdbqt(path, atoms, models=False):
    """atoms: list of (resname, resnum, x, y, z). If models=True, atoms is a list of lists (one per MODEL)."""
    lines = []
    groups = atoms if models else [atoms]
    for m, group in enumerate(groups, 1):
        if models:
            lines.append(f"MODEL {m}")
        for i, (res, num, x, y, z) in enumerate(group, 1):
            lines.append(pdbqt_line(i, res, num, x, y, z))
        if models:
            lines.append("ENDMDL")
    Path(path).write_text("\n".join(lines) + "\n")
    return str(path)


def load_function(path, name):
    """Load one top-level function from a file without running the file.

    docking_app_v2.py is a Streamlit script that runs UI code on import, so its helper
    functions cannot be imported normally. Moving the pure helpers into a module would
    remove the need for this (see NEXT_STEPS.md).
    """
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
    namespace = {}
    exec(compile(ast.Module([node], []), str(path), "exec"), namespace)
    return namespace[name]


@pytest.fixture
def root():
    return ROOT
