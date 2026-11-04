import json
import pytest

from provenance import RunStore, capture_environment


@pytest.fixture
def store():
    s = RunStore(":memory:")
    yield s
    s.close()


PARAMS = {"exhaustiveness": 8, "num_modes": 9, "center": [1.0, 2.0, 3.0]}


def test_start_run_records_params_and_environment(store):
    run_id = store.start_run("MRP1", PARAMS, notes="smoke")
    run = store.get_run(run_id)
    assert run["target"] == "MRP1"
    assert run["params"] == PARAMS
    assert run["status"] == "running"
    assert run["finished_at"] is None
    assert run["notes"] == "smoke"
    assert "python" in run["env"] and "packages" in run["env"]


def test_environment_snapshot_has_expected_keys():
    env = capture_environment()
    assert set(env) == {"python", "platform", "packages", "vina_cli", "git_commit"}
    assert isinstance(env["packages"], dict)


def test_finish_run_sets_status_and_timestamp(store):
    run_id = store.start_run("MRP1", PARAMS)
    store.finish_run(run_id, "completed")
    run = store.get_run(run_id)
    assert run["status"] == "completed"
    assert run["finished_at"] is not None


def test_finish_run_rejects_unknown_status_and_unknown_run(store):
    run_id = store.start_run("MRP1", PARAMS)
    with pytest.raises(ValueError):
        store.finish_run(run_id, "banana")
    with pytest.raises(KeyError):
        store.finish_run("doesnotexist", "completed")


def test_results_round_trip_in_iteration_order(store):
    run_id = store.start_run("MRP1", PARAMS)
    store.add_result(run_id, "L2", "CCO", -5.1, iteration=1, extra={"note": "b"})
    store.add_result(run_id, "L1", "CCN", -6.3, iteration=0, pose_path="poses/L1.pdbqt")
    rows = store.get_results(run_id)
    assert [r["ligand_id"] for r in rows] == ["L1", "L2"]
    assert rows[0]["pose_path"] == "poses/L1.pdbqt"
    assert rows[1]["extra"] == {"note": "b"}
    assert rows[0]["extra"] == {}


def test_best_results_lower_score_is_better_and_nulls_excluded(store):
    run_id = store.start_run("MRP1", PARAMS)
    store.add_result(run_id, "A", "C", -4.0)
    store.add_result(run_id, "B", "CC", -9.2)
    store.add_result(run_id, "C", "CCC", None)
    store.add_result(run_id, "D", "CCCC", -7.5)
    best = store.best_results(run_id, n=2)
    assert [r["ligand_id"] for r in best] == ["B", "D"]


def test_add_result_requires_existing_run(store):
    with pytest.raises(KeyError):
        store.add_result("nope", "L1", "C", -1.0)


def test_context_manager_marks_completed(store):
    with store.run("MRP1", PARAMS) as run_id:
        store.add_result(run_id, "L1", "C", -3.0)
    assert store.get_run(run_id)["status"] == "completed"


def test_context_manager_marks_failed_and_reraises(store):
    captured = {}
    with pytest.raises(RuntimeError, match="vina crashed"):
        with store.run("MRP1", PARAMS) as run_id:
            captured["id"] = run_id
            raise RuntimeError("vina crashed")
    assert store.get_run(captured["id"])["status"] == "failed"


def test_list_runs_newest_first(store):
    first = store.start_run("T1", PARAMS)
    second = store.start_run("T2", PARAMS)
    ids = [r["run_id"] for r in store.list_runs()]
    assert ids == [second, first]


def test_params_are_stored_as_sorted_json(store):
    run_id = store.start_run("MRP1", {"b": 1, "a": 2})
    raw = store._conn.execute(
        "SELECT params_json FROM runs WHERE run_id = ?", (run_id,)
    ).fetchone()[0]
    assert raw == json.dumps({"a": 2, "b": 1})


def test_file_database_creates_folder_and_persists(tmp_path):
    db = tmp_path / "nested" / "data" / "runs.db"
    first = RunStore(db)
    run_id = first.start_run("MRP1", PARAMS)
    first.add_result(run_id, "L1", "CCO", -5.5)
    first.close()
    second = RunStore(db)
    assert second.get_run(run_id)["target"] == "MRP1"
    assert second.get_results(run_id)[0]["score"] == -5.5
    second.close()
