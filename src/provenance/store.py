"""SQLite store that records each docking run and its results.

Design goals:
- Reproducible: every run stores its parameters and software environment.
- Simple: standard library only (sqlite3), no server needed.
- Safe: foreign keys are enforced, and a crashed run is marked 'failed'.

Vina scores are binding energies in kcal/mol, so LOWER (more negative) is better.
"""
from __future__ import annotations

import json
import platform
import shutil
import sqlite3
import subprocess
import sys
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Any, Iterator

# Packages whose versions are recorded with every run. Edit to match requirements.txt.
TRACKED_PACKAGES = ["numpy", "pandas", "rdkit", "vina", "streamlit", "openbabel-wheel"]

SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
    run_id      TEXT PRIMARY KEY,
    started_at  TEXT NOT NULL,
    finished_at TEXT,
    status      TEXT NOT NULL DEFAULT 'running',
    target      TEXT NOT NULL,
    params_json TEXT NOT NULL,
    env_json    TEXT NOT NULL,
    notes       TEXT
);
CREATE TABLE IF NOT EXISTS results (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id     TEXT NOT NULL REFERENCES runs(run_id) ON DELETE CASCADE,
    ligand_id  TEXT NOT NULL,
    smiles     TEXT NOT NULL,
    iteration  INTEGER NOT NULL DEFAULT 0,
    score      REAL,
    pose_path  TEXT,
    extra_json TEXT
);
CREATE INDEX IF NOT EXISTS idx_results_run ON results(run_id);
"""


def _utcnow() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _git_commit() -> str | None:
    if shutil.which("git") is None:
        return None
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, timeout=5
        )
        return (out.stdout.strip() or None) if out.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        return None


def _vina_cli_version() -> str | None:
    exe = shutil.which("vina")
    if exe is None:
        return None
    try:
        out = subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=5)
        text = (out.stdout or out.stderr).strip()
        return text.splitlines()[0] if text else None
    except (OSError, subprocess.SubprocessError):
        return None


def capture_environment() -> dict[str, Any]:
    """Snapshot of the software environment, stored with each run."""
    packages: dict[str, str] = {}
    for name in TRACKED_PACKAGES:
        try:
            packages[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            continue
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": packages,
        "vina_cli": _vina_cli_version(),
        "git_commit": _git_commit(),
    }


class RunStore:
    """Create with a file path (e.g. 'runs.db') or ':memory:' for tests."""

    def __init__(self, path: str | Path = "runs.db") -> None:
        self.path = str(path)
        if self.path != ":memory:":
            Path(self.path).parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(self.path)
        self._conn.row_factory = sqlite3.Row
        self._conn.execute("PRAGMA foreign_keys = ON")
        self._conn.executescript(SCHEMA)
        self._conn.commit()

    def close(self) -> None:
        self._conn.close()

    # ---- runs -----------------------------------------------------------
    def start_run(
        self, target: str, params: dict[str, Any], notes: str | None = None
    ) -> str:
        run_id = uuid.uuid4().hex[:12]
        self._conn.execute(
            "INSERT INTO runs (run_id, started_at, target, params_json, env_json, notes) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (run_id, _utcnow(), target, json.dumps(params, sort_keys=True),
             json.dumps(capture_environment(), sort_keys=True), notes),
        )
        self._conn.commit()
        return run_id

    def finish_run(self, run_id: str, status: str = "completed") -> None:
        if status not in {"completed", "failed", "cancelled"}:
            raise ValueError(f"unknown status: {status!r}")
        cur = self._conn.execute(
            "UPDATE runs SET status = ?, finished_at = ? WHERE run_id = ?",
            (status, _utcnow(), run_id),
        )
        if cur.rowcount == 0:
            raise KeyError(f"no such run: {run_id}")
        self._conn.commit()

    @contextmanager
    def run(
        self, target: str, params: dict[str, Any], notes: str | None = None
    ) -> Iterator[str]:
        """Context manager: marks the run completed, or failed if an error escapes."""
        run_id = self.start_run(target, params, notes)
        try:
            yield run_id
        except BaseException:
            self.finish_run(run_id, "failed")
            raise
        else:
            self.finish_run(run_id, "completed")

    def get_run(self, run_id: str) -> dict[str, Any]:
        row = self._conn.execute("SELECT * FROM runs WHERE run_id = ?", (run_id,)).fetchone()
        if row is None:
            raise KeyError(f"no such run: {run_id}")
        d = dict(row)
        d["params"] = json.loads(d.pop("params_json"))
        d["env"] = json.loads(d.pop("env_json"))
        return d

    def list_runs(self) -> list[dict[str, Any]]:
        rows = self._conn.execute(
            "SELECT run_id, started_at, status, target FROM runs ORDER BY started_at DESC, rowid DESC"
        ).fetchall()
        return [dict(r) for r in rows]

    # ---- results --------------------------------------------------------
    def add_result(
        self,
        run_id: str,
        ligand_id: str,
        smiles: str,
        score: float | None,
        iteration: int = 0,
        pose_path: str | None = None,
        extra: dict[str, Any] | None = None,
    ) -> None:
        try:
            self._conn.execute(
                "INSERT INTO results (run_id, ligand_id, smiles, iteration, score, pose_path, extra_json) "
                "VALUES (?, ?, ?, ?, ?, ?, ?)",
                (run_id, ligand_id, smiles, iteration, score, pose_path,
                 json.dumps(extra, sort_keys=True) if extra else None),
            )
        except sqlite3.IntegrityError as exc:
            raise KeyError(f"no such run: {run_id}") from exc
        self._conn.commit()

    def get_results(self, run_id: str) -> list[dict[str, Any]]:
        rows = self._conn.execute(
            "SELECT * FROM results WHERE run_id = ? ORDER BY iteration, id", (run_id,)
        ).fetchall()
        out = []
        for r in rows:
            d = dict(r)
            extra = d.pop("extra_json")
            d["extra"] = json.loads(extra) if extra else {}
            out.append(d)
        return out

    def best_results(self, run_id: str, n: int = 10) -> list[dict[str, Any]]:
        """Top-n by score. Lower (more negative) kcal/mol is better; NULL scores are excluded."""
        rows = self._conn.execute(
            "SELECT * FROM results WHERE run_id = ? AND score IS NOT NULL "
            "ORDER BY score ASC, id ASC LIMIT ?",
            (run_id, n),
        ).fetchall()
        return [dict(r) for r in rows]
