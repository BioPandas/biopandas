"""Utility function for reading Tripos MOL2 files from files"""

# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import sys

if sys.version_info >= (3, 9):
    import importlib.resources as pkg_resources
else:
    import importlib_resources as pkg_resources

import tests.mol2.data
from biopandas.mol2 import PandasMol2
from biopandas.mol2.mol2_io import split_multimol2
from tests.testutils import assert_raises

TEST_DATA = pkg_resources.files(tests.mol2.data)


def test_read_mol2():
    data_path_1 = str(TEST_DATA.joinpath("40_mol2_files.mol2"))
    data_path_2 = str(TEST_DATA.joinpath("40_mol2_files.mol2.gz"))
    data_path_3 = str(TEST_DATA.joinpath("empty_line.mol2"))

    for data_path in (data_path_1, data_path_2):
        pdmol = PandasMol2().read_mol2(data_path)
        assert pdmol.df.shape == (65, 9)
        assert pdmol.code == "ZINC38611810"

        expect = [
            "atom_id",
            "atom_name",
            "x",
            "y",
            "z",
            "atom_type",
            "subst_id",
            "subst_name",
            "charge",
        ]
        assert expect == list(pdmol.df.columns)
        assert len(pdmol.mol2_text) == 6469
        assert pdmol.mol2_path == data_path

    pdmol = PandasMol2().read_mol2(data_path_3)
    assert pdmol.df.shape == (59, 9)


def test_read_mol2_from_list():
    data_path = str(TEST_DATA.joinpath("40_mol2_files.mol2"))
    mol2 = next(split_multimol2(data_path))

    pdmol = PandasMol2().read_mol2_from_list(
        mol2_lines=mol2[1], mol2_code=mol2[0]
    )
    assert pdmol.df.shape == (65, 9)
    assert pdmol.code == "ZINC38611810"


def test_rmsd():
    data_path_1 = str(TEST_DATA.joinpath("1b5e_1.mol2"))
    data_path_2 = str(TEST_DATA.joinpath("1b5e_2.mol2"))

    pdmol_1 = PandasMol2().read_mol2(data_path_1)
    pdmol_2 = PandasMol2().read_mol2(data_path_2)

    assert pdmol_1.rmsd(pdmol_1.df, pdmol_2.df, heavy_only=False) == 1.5523
    assert pdmol_1.rmsd(pdmol_1.df, pdmol_2.df) == 1.1609


def test_distance():
    data_path = str(TEST_DATA.joinpath("1b5e_1.mol2"))

    pdmol = PandasMol2().read_mol2(data_path)
    assert round(pdmol.distance().values[0], 3) == 31.185


def test_distance_external_df():
    data_path = str(TEST_DATA.joinpath("1b5e_1.mol2"))

    pdmol = PandasMol2().read_mol2(data_path)
    new_df = pdmol.df.iloc[1:, :].copy()
    assert round(PandasMol2.distance_df(df=new_df).values[0], 3) == 31.165


def test_overwrite_df():
    data_path = str(TEST_DATA.joinpath("1b5e_1.mol2"))
    pdmol = PandasMol2().read_mol2(data_path)

    def overwrite():
        pdmol.df = pdmol.df[(pdmol.df["atom_type"] != "H")]

    expect = (
        "Please use `PandasMol2._df = ... `"
        " instead\nof `PandasMol2.df = ... `"
        " if you are sure that\nyou want"
        " to overwrite the `df` attribute."
    )

    assert_raises(AttributeError, expect, overwrite)


def test__get_atomsection_raises():
    """Test if ValueError is raised if input list is not in the mol2 format."""

    expect = (
        "Structural data could not be loaded. "
        "Is the input file/text in the mol2 format?"
    )

    def run_code():
        PandasMol2()._get_atomsection(["", ""])

    assert_raises(ValueError, expect, run_code)
