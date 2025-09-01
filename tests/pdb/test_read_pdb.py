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

from urllib.error import HTTPError

import numpy as np
import pandas as pd
import pytest

import tests.pdb.data
from biopandas.pdb import PandasPdb
from tests.testutils import assert_raises

TEST_DATA = pkg_resources.files(tests.pdb.data)

TESTDATA_FILENAME = str(TEST_DATA.joinpath("3eiy.pdb"))
TESTDATA_FILENAME2 = str(TEST_DATA.joinpath("4eiy_anisouchunk.pdb"))
TESTDATA_FILENAME_GZ = str(TEST_DATA.joinpath("3eiy.pdb.gz"))
TESTDATA_FILENAME_AF2_V4 = str(TEST_DATA.joinpath("AF-Q5VSL9-F1-model_v4.pdb"))

TESTDATA_FILENAME_AF2_V3 = str(TEST_DATA.joinpath("AF-Q5VSL9-F1-model_v3.pdb"))

ATOM_DF_COLUMNS = [
    "record_name",
    "atom_number",
    "blank_1",
    "atom_name",
    "alt_loc",
    "residue_name",
    "blank_2",
    "chain_id",
    "residue_number",
    "insertion",
    "blank_3",
    "x_coord",
    "y_coord",
    "z_coord",
    "occupancy",
    "b_factor",
    "blank_4",
    "segment_id",
    "element_symbol",
    "charge",
    "line_idx",
]

ANISOU_DF_COLUMNS = [
    "record_name",
    "atom_number",
    "blank_1",
    "atom_name",
    "alt_loc",
    "residue_name",
    "blank_2",
    "chain_id",
    "residue_number",
    "insertion",
    "blank_3",
    "U(1,1)",
    "U(2,2)",
    "U(3,3)",
    "U(1,2)",
    "U(1,3)",
    "U(2,3)",
    "blank_4",
    "element_symbol",
    "charge",
    "line_idx",
]

with open(TESTDATA_FILENAME, "r") as f:
    three_eiy = f.read()

with open(TESTDATA_FILENAME2, "r") as f:
    four_eiy = f.read()

with open(TESTDATA_FILENAME_AF2_V4, "r") as f:
    af_test_struct_v4 = f.read()

with open(TESTDATA_FILENAME_AF2_V3, "r") as f:
    af_test_struct_v3 = f.read()


def test__read_pdb():
    """Test private _read_pdb"""
    ppdb = PandasPdb()
    _, txt = ppdb._read_pdb(TESTDATA_FILENAME)
    assert txt == three_eiy


def test__read_pdb_raises():
    """Test private _read_pdb:
    Test if ValueError is raised for wrong file formats."""

    expect = (
        "Wrong file format; allowed file formats are " ".pdb, .pdb.gz, .ent, .ent.gz"
    )

    def run_code_1():
        PandasPdb()._read_pdb("protein.mol2")

    assert_raises(ValueError, expect, run_code_1)

    def run_code_2():
        PandasPdb()._read_pdb("protein.mol2.gz")

    assert_raises(ValueError, expect, run_code_2)


def test_fetch_pdb():
    """Test fetch_pdb"""

    try:
        ppdb = PandasPdb()
        _, txt = ppdb._fetch_pdb("3eiy")
    except HTTPError:
        _, txt = None, None
    except ConnectionResetError:
        _, txt = None, None

    if txt:  # skip if PDB down
        txt[:100] == three_eiy[:100]
        ppdb.fetch_pdb("3eiy")
        assert ppdb.pdb_text == txt
        assert ppdb.pdb_path == "https://files.rcsb.org/download/3eiy.pdb"


def test_fetch_af2():
    """Test fetch_pdb"""
    # Check latest release
    try:
        ppdb = PandasPdb()
        _, txt = ppdb._fetch_af2("Q5VSL9", af2_version=4)
    except HTTPError:
        _, txt = None, None
    except ConnectionResetError:
        _, txt = None, None

    if txt:  # skip if AF2 DB down
        txt[:100] == af_test_struct_v4[:100]
        ppdb.fetch_pdb(uniprot_id="Q5VSL9", source="alphafold2-v4")
        assert ppdb.pdb_text == txt
        assert (
            ppdb.pdb_path
            == "https://alphafold.ebi.ac.uk/files/AF-Q5VSL9-F1-model_v4.pdb"
        )

    # Check legacy release
    try:
        ppdb = PandasPdb()
        _, txt = ppdb._fetch_af2("Q5VSL9", af2_version=3)
    except HTTPError:
        _, txt = None, None
    except ConnectionResetError:
        _, txt = None, None

    if txt:  # skip if AF2 DB down
        txt[:100] == af_test_struct_v3[:100]
        ppdb.fetch_pdb(uniprot_id="Q5VSL9", source="alphafold2-v3")
        assert ppdb.pdb_text == txt
        assert (
            ppdb.pdb_path
            == "https://alphafold.ebi.ac.uk/files/AF-Q5VSL9-F1-model_v3.pdb"
        )


def test__read_pdb_gz():
    """Test public _read_pdb with gzip files"""
    ppdb = PandasPdb()
    _, txt = ppdb._read_pdb(TESTDATA_FILENAME_GZ)
    assert txt == three_eiy


def test__construct_df():
    """Test pandas dataframe construction"""
    ppdb = PandasPdb()
    dfs = ppdb._construct_df(three_eiy.splitlines())
    assert set(dfs.keys()) == {"OTHERS", "ATOM", "ANISOU", "HETATM"}
    assert set(dfs["ATOM"].columns) == set(ATOM_DF_COLUMNS)
    assert set(dfs["HETATM"].columns) == set(ATOM_DF_COLUMNS)
    assert set(dfs["ANISOU"].columns) == set(ANISOU_DF_COLUMNS)
    exp = pd.Series(
        np.array(
            [
                "ATOM",
                1,
                "",
                "N",
                "",
                "SER",
                "",
                "A",
                2,
                "",
                "",
                2.527,
                54.656,
                -1.667,
                1.0,
                52.73,
                "",
                "",
                "N",
                None,
                609,
            ]
        ),
        index=[
            "record_name",
            "atom_number",
            "blank_1",
            "atom_name",
            "alt_loc",
            "residue_name",
            "blank_2",
            "chain_id",
            "residue_number",
            "insertion",
            "blank_3",
            "x_coord",
            "y_coord",
            "z_coord",
            "occupancy",
            "b_factor",
            "blank_4",
            "segment_id",
            "element_symbol",
            "charge",
            "line_idx",
        ],
    )
    assert exp.equals(dfs["ATOM"].loc[0, :])


def test_read_pdb():
    """Test public read_pdb"""
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    assert ppdb.pdb_text == three_eiy
    assert ppdb.code == "3eiy", ppdb.code
    assert ppdb.pdb_path == TESTDATA_FILENAME


def test_read_pdb_from_list():
    """Test public read_pdb_from_list"""

    for pdb_text, code in zip([three_eiy, four_eiy], ["3eiy", "4eiy"]):
        ppdb = PandasPdb()
        ppdb.read_pdb_from_list(pdb_text.splitlines(True))
        assert ppdb.pdb_text == pdb_text
        assert ppdb.code == code
        assert ppdb.pdb_path == ""


def test_anisou_input_handling():
    """Test public read_pdb"""
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME2)
    assert ppdb.pdb_text == four_eiy
    assert ppdb.code == "4eiy", ppdb.code


@pytest.mark.xfail(raises=AttributeError)
def test_get_exceptions():
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    ppdb.get("main-chai")


def test_get_all():
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    for i in ["c-alpha", "hydrogen", "main chain"]:
        ppdb.get(i)


def test_get_df():
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)

    shape = ppdb.get("c-alpha").shape
    assert shape == (174, 21), shape

    shape = ppdb.get("hydrogen", invert=True, records=("ATOM",)).shape
    assert shape == (1330, 21), shape

    # deprecated use of string
    shape = ppdb.get("hydrogen", invert=True, records="ATOM").shape
    assert shape == (1330, 21), shape

    shape = ppdb.get("hydrogen").shape
    assert shape == (0, 21), shape

    shape = ppdb.get("main chain", records=("ATOM",)).shape
    assert shape == (696, 21), shape

    shape = ppdb.get("heavy", records=("ATOM",)).shape
    assert shape == (1330, 21), shape

    shape = ppdb.get("carbon", records=("ATOM",)).shape
    assert shape == (857, 21), shape
