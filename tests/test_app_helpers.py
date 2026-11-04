import pytest

from conftest import ROOT, load_function, write_pdbqt

auto_grid = load_function(ROOT / "docking_app_v2.py", "auto_grid_from_receptor")


def test_auto_grid_covers_receptor_with_padding(tmp_path):
    rec = write_pdbqt(tmp_path / "r.pdbqt", [("ALA", 1, 0.0, 0.0, 0.0), ("GLY", 2, 10.0, 20.0, 30.0)])
    g = auto_grid(rec)
    assert (g["center_x"], g["center_y"], g["center_z"]) == (5.0, 10.0, 15.0)
    assert (g["size_x"], g["size_y"], g["size_z"]) == (18.0, 28.0, 38.0)


def test_auto_grid_padding_is_adjustable(tmp_path):
    rec = write_pdbqt(tmp_path / "r.pdbqt", [("ALA", 1, 0.0, 0.0, 0.0), ("GLY", 2, 10.0, 10.0, 10.0)])
    assert auto_grid(rec, padding=0.0)["size_x"] == 10.0


def test_auto_grid_handles_negative_coordinates(tmp_path):
    rec = write_pdbqt(tmp_path / "r.pdbqt", [("ALA", 1, -10.0, -4.0, -2.0), ("GLY", 2, 10.0, 4.0, 2.0)])
    g = auto_grid(rec)
    assert (g["center_x"], g["center_y"], g["center_z"]) == (0.0, 0.0, 0.0)


def test_auto_grid_counts_hetatm_and_skips_malformed_lines(tmp_path):
    f = tmp_path / "r.pdbqt"
    f.write_text("REMARK x\nATOM      1  CA  ALA A   1        BAD     0.000   0.000\n")
    with open(f, "a") as fh:
        fh.write("HETATM    2  O   HOH A   2       2.000   2.000   2.000  1.00  0.00    +0.000 O \n")
    g = auto_grid(str(f))
    assert g["center_x"] == 2.0


def test_auto_grid_empty_receptor_raises(tmp_path):
    f = tmp_path / "empty.pdbqt"
    f.write_text("REMARK nothing\n")
    with pytest.raises(ValueError, match="No atoms"):
        auto_grid(str(f))
