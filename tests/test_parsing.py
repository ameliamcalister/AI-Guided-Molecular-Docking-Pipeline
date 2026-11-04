import pytest

import docking_pipeline as dp
from conftest import write_pdbqt

VINA_LOG = """\
AutoDock Vina v1.2.5
Rigid receptor: rec/5UJA_receptor.pdbqt
Ligand: /tmp/work/compound_7.pdbqt
Grid center: X 94.806 Y 59.056 Z 56.653
Grid size  : X 21.92 Y 21.297 Z 23.806
Exhaustiveness: 8
Performing docking (random seed: 42) ...

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -9.123          0          0
   2       -8.900      1.234      2.345
   3       -8.410      2.001      3.110
"""


def test_parse_vina_log_reads_all_fields(tmp_path):
    log = tmp_path / "run.log"
    log.write_text(VINA_LOG)
    d = dp.parse_vina_log(str(log))
    assert d["best_affinity"] == -9.123
    assert d["n_modes"] == 3
    assert d["ligand_name"] == "compound_7"
    assert d["receptor_name"] == "5UJA_receptor"
    assert d["grid_center"] == (94.806, 59.056, 56.653)
    assert d["exhaustiveness"] == 8
    assert d["random_seed"] == 42
    assert d["all_modes"][1] == {"mode": 2, "affinity": -8.9, "rmsd_lb": 1.234, "rmsd_ub": 2.345}


def test_parse_vina_log_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        dp.parse_vina_log(str(tmp_path / "nope.log"))


def test_parse_vina_log_without_modes_raises(tmp_path):
    log = tmp_path / "empty.log"
    log.write_text("AutoDock Vina v1.2.5\nLigand: x.pdbqt\n")
    with pytest.raises(ValueError, match="No docking modes"):
        dp.parse_vina_log(str(log))


@pytest.mark.xfail(strict=True, reason="grid-centre regex only matches positive numbers, so a negative "
                                       "centre is silently dropped (grid_center is not used downstream)")
def test_parse_vina_log_negative_grid_centre(tmp_path):
    log = tmp_path / "neg.log"
    log.write_text(VINA_LOG.replace("X 94.806 Y 59.056 Z 56.653", "X -12.5 Y 3.0 Z -7.25"))
    assert dp.parse_vina_log(str(log))["grid_center"] == (-12.5, 3.0, -7.25)


def test_parse_pdbqt_atoms_reads_fields_and_skips_other_lines(tmp_path):
    f = write_pdbqt(tmp_path / "r.pdbqt", [("ALA", 10, 1.5, -2.25, 3.0)])
    with open(f, "a") as fh:
        fh.write("REMARK not an atom\nTORSDOF 0\n")
    atoms = dp.parse_pdbqt_atoms(f)
    assert atoms == [{"record": "ATOM", "resname": "ALA", "resnum": 10, "x": 1.5, "y": -2.25, "z": 3.0}]


RECEPTOR = [("ALA", 10, 0.0, 0.0, 0.0), ("GLY", 11, 3.0, 0.0, 0.0), ("LEU", 12, 10.0, 0.0, 0.0)]


def test_detect_contacts_within_cutoff_sorted_by_distance(tmp_path):
    rec = write_pdbqt(tmp_path / "rec.pdbqt", RECEPTOR)
    lig = write_pdbqt(tmp_path / "lig.pdbqt", [("LIG", 1, 0.0, 0.0, 1.0)])
    residues, detail = dp.detect_contacts(lig, rec, cutoff=3.5)
    assert residues == ["Ala10", "Gly11"]
    assert detail == {"Ala10": 1.0, "Gly11": 3.16}


def test_detect_contacts_cutoff_changes_result(tmp_path):
    rec = write_pdbqt(tmp_path / "rec.pdbqt", RECEPTOR)
    lig = write_pdbqt(tmp_path / "lig.pdbqt", [("LIG", 1, 0.0, 0.0, 1.0)])
    assert dp.detect_contacts(lig, rec, cutoff=2.0)[0] == ["Ala10"]


def test_detect_contacts_applies_numbering_offset(tmp_path):
    rec = write_pdbqt(tmp_path / "rec.pdbqt", RECEPTOR)
    lig = write_pdbqt(tmp_path / "lig.pdbqt", [("LIG", 1, 0.0, 0.0, 1.0)])
    residues, _ = dp.detect_contacts(lig, rec, numbering_offset=100)
    assert residues == ["Ala110", "Gly111"]


def test_detect_contacts_keeps_nearest_distance_per_residue(tmp_path):
    rec = write_pdbqt(tmp_path / "rec.pdbqt", [("ALA", 10, 0.0, 0.0, 0.0), ("ALA", 10, 0.0, 0.0, 2.0)])
    lig = write_pdbqt(tmp_path / "lig.pdbqt", [("LIG", 1, 0.0, 0.0, 3.0)])
    assert dp.detect_contacts(lig, rec)[1] == {"Ala10": 1.0}


def test_detect_contacts_unknown_residue_name_is_capitalised(tmp_path):
    rec = write_pdbqt(tmp_path / "rec.pdbqt", [("XYZ", 5, 0.0, 0.0, 0.0)])
    lig = write_pdbqt(tmp_path / "lig.pdbqt", [("LIG", 1, 0.0, 0.0, 1.0)])
    assert dp.detect_contacts(lig, rec)[0] == ["Xyz5"]


def test_get_ligand_centre_is_mean_of_atoms(tmp_path):
    lig = write_pdbqt(tmp_path / "lig.pdbqt", [("LIG", 1, 0.0, 0.0, 0.0), ("LIG", 1, 2.0, 4.0, 6.0)])
    assert dp.get_ligand_centre(lig) == (1.0, 2.0, 3.0)


def test_get_ligand_centre_empty_file_is_none(tmp_path):
    f = tmp_path / "empty.pdbqt"
    f.write_text("REMARK nothing\n")
    assert dp.get_ligand_centre(str(f)) is None


# ── Known issue: a Vina output file holds up to 9 poses (MODEL blocks) and the pipeline reads them all.
POSES = [[("LIG", 1, 0.0, 0.0, 1.0)], [("LIG", 1, 10.0, 0.0, 1.0)]]


@pytest.mark.xfail(strict=True, reason="detect_contacts reads every pose in a docked file, so residues touching "
                                       "poses 2-9 are reported as contacts of the best pose")
def test_detect_contacts_uses_best_pose_only(tmp_path):
    rec = write_pdbqt(tmp_path / "rec.pdbqt", RECEPTOR)
    lig = write_pdbqt(tmp_path / "lig.pdbqt", POSES, models=True)
    assert dp.detect_contacts(lig, rec)[0] == ["Ala10", "Gly11"]


@pytest.mark.xfail(strict=True, reason="get_ligand_centre averages atoms over all poses, not just the best one")
def test_ligand_centre_uses_best_pose_only(tmp_path):
    lig = write_pdbqt(tmp_path / "lig.pdbqt", POSES, models=True)
    assert dp.get_ligand_centre(lig) == (0.0, 0.0, 1.0)
