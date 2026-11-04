import json

import docking_pipeline as dp

GRID = dict(grid_centre=(94.8, 59.0, 56.6), grid_size=(21.9, 21.3, 23.8))


def test_system_prompt_contains_target_grid_and_residues():
    p = dp.build_system_prompt("MRP1 (ABCC1)", "5UJA", "substrate cavity",
                               ["Thr1241", "Lys332"], reference_smiles="CCO", chembl_id="CHEMBL4523", **GRID)
    for text in ["MRP1 (ABCC1)", "5UJA", "substrate cavity", "Thr1241", "Lys332",
                 "centre_x = 94.8", "size_z   = 23.8", "CCO", "CHEMBL4523"]:
        assert text in p


def test_system_prompt_without_residues_says_so():
    p = dp.build_system_prompt("T", "1ABC", "site", [], **GRID)
    assert "not specified" in p


def test_system_prompt_states_scoring_convention():
    p = dp.build_system_prompt("T", "1ABC", "site", [], **GRID)
    assert "more negative = stronger" in p


def _prompt(**kw):
    base = dict(smiles="CCO", docking_score=-9.0, residues=[], contact_detail={},
                ligand_centre=None, **GRID)
    base.update(kw)
    return dp.build_prompt(**base)


def test_build_prompt_grades_contacts_by_distance():
    p = _prompt(contact_detail={"Ala1": 2.0, "Gly2": 2.7, "Leu3": 3.3})
    assert "STRONG (<2.5 Å)" in p and "Ala1: 2.0" in p
    assert "MEDIUM (2.5–3.0 Å)" in p and "Gly2: 2.7" in p
    assert "WEAK (3.0–3.5 Å)" in p and "Leu3: 3.3" in p


def test_build_prompt_without_contacts_targets_full_cavity():
    assert "Not detected — target full cavity" in _prompt()


def test_build_prompt_reports_ligand_displacement_from_grid_centre():
    p = _prompt(ligand_centre=(104.8, 59.0, 56.6))
    assert "ΔX=10.0" in p and "Expand toward: negative X" in p


def test_build_prompt_modes_change_the_task():
    assert "improve binding affinity" in _prompt(admet_mode=False)
    admet = _prompt(admet_mode=True)
    assert "ADMET optimisation mode" in admet and "-9.00 kcal/mol" in admet


def test_build_prompt_includes_name_notes_and_count():
    p = _prompt(compound_name="cmpd1", notes="watch the hinge", n_suggestions=5)
    assert "Compound name: cmpd1" in p and "watch the hinge" in p and "Suggest 5 structural" in p


# ── iteration mode resolution ────────────────────────────────────────────────
def test_resolve_mode_manual_modes():
    assert dp.resolve_mode("affinity", [9.0, 10.0], 13.0)[0] == "affinity"
    assert dp.resolve_mode("admet", [9.0, 10.0], 13.0)[0] == "admet"


def test_resolve_mode_auto_starts_with_affinity():
    assert dp.resolve_mode("auto", [], 13.0)[0] == "affinity"


def test_resolve_mode_auto_switches_when_target_reached():
    mode, msg = dp.resolve_mode("auto", [9.0, 13.5], 13.0)
    assert mode == "admet" and "reached" in msg


def test_resolve_mode_auto_switches_on_plateau():
    mode, msg = dp.resolve_mode("auto", [9.0, 9.1, 9.2], 13.0)
    assert mode == "admet" and "Plateau" in msg


def test_resolve_mode_auto_keeps_going_while_improving():
    assert dp.resolve_mode("auto", [9.0, 10.0, 11.5], 13.0)[0] == "affinity"


def test_resolve_mode_auto_needs_enough_history_to_call_a_plateau():
    assert dp.resolve_mode("auto", [9.0, 9.0], 13.0)[0] == "affinity"


def test_resolve_mode_uses_absolute_values():
    assert dp.resolve_mode("auto", [-9.0, -13.5], 13.0)[0] == "admet"


def test_history_round_trip_and_missing_file(tmp_path):
    f = tmp_path / "history.json"
    assert dp.load_history(str(f)) == []
    dp.save_history(str(f), [9.1, 10.2])
    assert dp.load_history(str(f)) == [9.1, 10.2]
    assert json.loads(f.read_text()) == [9.1, 10.2]
