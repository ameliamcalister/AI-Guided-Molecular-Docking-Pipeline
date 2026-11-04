import pytest

import docking_pipeline as dp

pytestmark = pytest.mark.skipif(not dp.RDKIT_AVAILABLE, reason="RDKit not installed")

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


def test_validate_smiles_accepts_valid_and_rejects_invalid():
    ok, mol = dp.validate_smiles(ASPIRIN)
    assert ok and mol is not None
    ok, mol = dp.validate_smiles("C1CC(=O)N(C")       # unclosed ring and branch
    assert not ok and mol is None


def test_validate_smiles_strips_whitespace():
    assert dp.validate_smiles("  CCO \n")[0]


def test_calculate_properties_for_ethanol():
    _, mol = dp.validate_smiles("CCO")
    p = dp.calculate_properties(mol)
    assert p["MW"] == pytest.approx(46.04, abs=0.01)
    assert (p["HBD"], p["HBA"], p["RotBonds"]) == (1, 1, 0)
    assert p["TPSA"] == pytest.approx(20.2, abs=0.1)


def test_calculate_properties_none_returns_empty():
    assert dp.calculate_properties(None) == {}


GOOD = {"MW": 400, "cLogP": 3.0, "HBD": 2, "HBA": 6, "TPSA": 90, "RotBonds": 5}


def test_check_admet_clean_molecule_has_no_issues():
    assert dp.check_admet(GOOD) == []


def test_check_admet_empty_props_has_no_issues():
    assert dp.check_admet({}) == []


@pytest.mark.parametrize("change, level, label", [
    ({"MW": 610}, "WARN", "MW"), ({"MW": 650}, "FLAG", "MW"),
    ({"cLogP": 5.8}, "WARN", "cLogP"), ({"cLogP": 6.5}, "FLAG", "cLogP"),
    ({"HBD": 4}, "FLAG", "HBD"), ({"HBA": 11}, "FLAG", "HBA"),
    ({"TPSA": 150}, "FLAG", "TPSA"),
    ({"RotBonds": 11}, "WARN", "RotBonds"), ({"RotBonds": 13}, "FLAG", "RotBonds"),
])
def test_check_admet_thresholds(change, level, label):
    issues = dp.check_admet({**GOOD, **change})
    assert len(issues) == 1
    assert issues[0].startswith(level) and label in issues[0]


def test_check_admet_limits_are_inclusive():
    edge = {"MW": 600, "cLogP": 5.5, "HBD": 3, "HBA": 10, "TPSA": 140, "RotBonds": 10}
    assert dp.check_admet(edge) == []


def test_extract_smiles_finds_backticked_and_labelled_and_dedupes():
    text = (f"1. Add acetyl: `{ASPIRIN}`\n"
            f"**SMILES**: {ASPIRIN}\n"
            "2. Another: **SMILES:** `CCN(CC)C(=O)c1ccccc1`\n")
    found = [f["smiles"] for f in dp.extract_smiles_from_response(text)]
    assert found.count(ASPIRIN) == 1
    assert "CCN(CC)C(=O)c1ccccc1" in found


def test_extract_smiles_ignores_short_or_plain_words():
    assert dp.extract_smiles_from_response("Use `CCO` or `benzene` here.") == []


def test_validate_response_smiles_reports_valid_and_invalid():
    text = f"Valid: `{ASPIRIN}`\nBroken: `C1CC(=O)N(CC(C`\n"
    annotated, results = dp.validate_response_smiles(text)
    assert [r["valid"] for r in results] == [True, False]
    assert "1/2 valid" in annotated and "1 invalid" in annotated


def test_validate_response_smiles_without_smiles_returns_empty_list():
    annotated, results = dp.validate_response_smiles("No molecules here.")
    assert results == [] and "No SMILES detected" in annotated
