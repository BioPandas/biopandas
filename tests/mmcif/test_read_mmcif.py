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

from pathlib import Path
from urllib.error import HTTPError

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import tests.mmcif.data
from biopandas.mmcif import PandasMmcif
from biopandas.pdb import PandasPdb
from tests.testutils import assert_raises

TEST_DATA = pkg_resources.files(tests.mmcif.data)

TESTDATA_FILENAME = str(TEST_DATA.joinpath("3eiy.cif"))

# Not clear on how ANISOU records are handled in mmCIF files so skipping
# TESTDATA_FILENAME2 = str(TEST_DATA.joinpath("4eiy_anisouchunk.cif"))
TESTDATA_FILENAME2 = str(TEST_DATA.joinpath("4eiy.cif"))
TESTDATA_FILENAME_GZ = str(TEST_DATA.joinpath("3eiy.cif.gz"))

TESTDATA_FILENAME_AF2_V4 = str(TEST_DATA.joinpath("AF-Q5VSL9-F1-model_v4.cif"))
TESTDATA_FILENAME_AF2_V3 = str(TEST_DATA.joinpath("AF-Q5VSL9-F1-model_v3.cif"))

ATOM_DF_COLUMNS = [
    "B_iso_or_equiv",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "auth_asym_id",
    "auth_atom_id",
    "auth_comp_id",
    "auth_seq_id",
    "group_PDB",
    "id",
    "label_alt_id",
    "label_asym_id",
    "label_atom_id",
    "label_comp_id",
    "label_entity_id",
    "label_seq_id",
    "occupancy",
    "pdbx_PDB_ins_code",
    "pdbx_PDB_model_num",
    "pdbx_formal_charge",
    "type_symbol",
]

ANISOU_DF_COLUMNS = [
    "id",
    "type_symbol",
    "pdbx_label_atom_id",
    "pdbx_label_alt_id",
    "pdbx_label_comp_id",
    "pdbx_label_asym_id",
    "pdbx_label_seq_id",
    "pdbx_PDB_ins_code",
    "U[1][1]",
    "U[2][2]",
    "U[3][3]",
    "U[1][2]",
    "U[1][3]",
    "U[2][3]",
    "pdbx_auth_seq_id",
    "pdbx_auth_comp_id",
    "pdbx_auth_asym_id",
    "pdbx_auth_atom_id",
]

with open(TESTDATA_FILENAME, "r") as f:
    three_eiy = f.read()

with open(TESTDATA_FILENAME2, "r") as f:
    four_eiy = f.read()

with open(TESTDATA_FILENAME_AF2_V4, "r") as f:
    af2_test_struct_v4 = f.read()

with open(TESTDATA_FILENAME_AF2_V3, "r") as f:
    af2_test_struct_v3 = f.read()


def test__read_pdb():
    """Test private _read_pdb"""
    ppdb = PandasMmcif()
    _, txt = ppdb._read_mmcif(TESTDATA_FILENAME)
    assert txt == three_eiy


def test__read_pdb_raises():
    """Test private _read_pdb:
    Test if ValueError is raised for wrong file formats."""

    expect = (
        "Wrong file format; allowed file formats are "
        ".cif, .cif.gz, .mmcif, .mmcif.gz"
    )

    def run_code_1():
        PandasMmcif()._read_mmcif("protein.mol2")

    assert_raises(ValueError, expect, run_code_1)

    def run_code_2():
        PandasMmcif()._read_mmcif("protein.mol2.gz")

    assert_raises(ValueError, expect, run_code_2)


def test_fetch_pdb():
    """Test fetch_pdb"""

    try:
        ppdb = PandasMmcif()
        _, txt = ppdb._fetch_mmcif("3eiy")
    except (HTTPError, ConnectionResetError):
        _, txt = None, None
    if txt:  # skip if PDB down
        txt[:100] == three_eiy[:100]
        ppdb.fetch_mmcif("3eiy")
        assert ppdb.mmcif_text == txt
        assert ppdb.mmcif_path == "https://files.rcsb.org/download/3eiy.cif"


def test_fetch_af2():
    """Test fetch_af2"""
    # Test latest release
    try:
        ppdb = PandasMmcif()
        _, txt = ppdb._fetch_af2("Q5VSL9", af2_version=4)
    except (HTTPError, ConnectionResetError):
        _, txt = None, None
    if txt:  # skip if AF DB down
        txt[:100] == af2_test_struct_v4[:100]
        ppdb.fetch_mmcif(uniprot_id="Q5VSL9", source="alphafold2-v4")
        assert ppdb.mmcif_text == txt
        assert (
            ppdb.mmcif_path
            == "https://alphafold.ebi.ac.uk/files/AF-Q5VSL9-F1-model_v4.cif"
        )

    # Test legacy release
    try:
        ppdb = PandasMmcif()
        _, txt = ppdb._fetch_af2("Q5VSL9", af2_version=3)
    except (HTTPError, ConnectionResetError):
        _, txt = None, None
    if txt:  # skip if AF DB down
        txt[:100] == af2_test_struct_v3[:100]
        ppdb.fetch_mmcif(uniprot_id="Q5VSL9", source="alphafold2-v3")
        assert ppdb.mmcif_text == txt
        assert (
            ppdb.mmcif_path
            == "https://alphafold.ebi.ac.uk/files/AF-Q5VSL9-F1-model_v3.cif"
        )


def test__read_pdb_gz():
    """Test public _read_pdb with gzip files"""
    ppdb = PandasMmcif()
    _, txt = ppdb._read_mmcif(TESTDATA_FILENAME_GZ)
    assert txt == three_eiy


def test__construct_df():
    """Test pandas dataframe construction"""
    ppdb = PandasMmcif()
    dfs = ppdb._construct_df(three_eiy)
    # assert set(dfs.keys()) == {"OTHERS", "ATOM", "ANISOU", "HETATM"}
    # Currently don't parse OTHERS records as I'm not sure where they're located in mmCIF files
    assert set(dfs.keys()) == {"ATOM", "ANISOU", "HETATM"}
    assert set(dfs["ATOM"].columns) == set(ATOM_DF_COLUMNS)
    assert set(dfs["HETATM"].columns) == set(ATOM_DF_COLUMNS)
    assert set(dfs["ANISOU"].columns) == set(ANISOU_DF_COLUMNS)
    exp = pd.Series(
        [
            52.73,
            2.527,
            54.656,
            -1.667,
            "A",
            "N",
            "SER",
            2,
            "ATOM",
            1,
            None,
            "A",
            "N",
            "SER",
            1,
            "23",
            1.0,
            "None",
            1,
            None,
            "N",
        ],
        index=[
            "B_iso_or_equiv",
            "Cartn_x",
            "Cartn_y",
            "Cartn_z",
            "auth_asym_id",
            "auth_atom_id",
            "auth_comp_id",
            "auth_seq_id",
            "group_PDB",
            "id",
            "label_alt_id",
            "label_asym_id",
            "label_atom_id",
            "label_comp_id",
            "label_entity_id",
            "label_seq_id",
            "occupancy",
            "pdbx_PDB_ins_code",
            "pdbx_PDB_model_num",
            "pdbx_formal_charge",
            "type_symbol",
        ],
    )
    # There's some weird behaviour here I can't quite figure out.
    # Comparing the two series fails on exp.label_alt_id and exp.pdbx_formal_charge
    # However, if I compare these values directly, they are equal.
    # assert exp.equals(dfs["ATOM"].loc[0, :])
    for k, v in exp.items():
        assert dfs["ATOM"].loc[0, :][k] == v, k


def test_read_pdb():
    """Test public read_pdb"""
    ppdb = PandasMmcif()
    ppdb.read_mmcif(TESTDATA_FILENAME)
    assert ppdb.pdb_text == three_eiy
    assert ppdb.code == "3eiy", ppdb.code
    assert ppdb.mmcif_path == TESTDATA_FILENAME


def test_read_pdb_from_list():
    """Test public read_pdb_from_list"""

    for pdb_text, code in zip([three_eiy, four_eiy], ["3eiy", "4eiy"]):
        ppdb = PandasMmcif()
        ppdb.read_mmcif_from_list(pdb_text)
        assert ppdb.pdb_text == pdb_text
        assert ppdb.code == code
        assert ppdb.mmcif_path == ""


def test_read_pdb_with_pathlib():
    """Test public read_pdb with pathlib.Path object as input"""
    ppdb = PandasMmcif()
    ppdb.read_mmcif(Path(TESTDATA_FILENAME))
    assert ppdb.pdb_text == three_eiy
    assert ppdb.code == "3eiy", ppdb.code
    assert ppdb.mmcif_path == TESTDATA_FILENAME


# Again, not sure how ANISOU are handled, so skipping
# It seems like the whole PDB ATOM record is duplicated - whereas
# previously ANISOU records were only present for a subset of atoms
# def test_anisou_input_handling():
#    """Test public read_pdb"""
#    ppdb = PandasMmcif()
#    ppdb.read_pdb(TESTDATA_FILENAME2)
#    assert ppdb.pdb_text == four_eiy
#    assert ppdb.code == "4eiy", ppdb.code


@pytest.mark.xfail(raises=AttributeError)
def test_get_exceptions():
    ppdb = PandasMmcif()
    ppdb.read_mmcif(TESTDATA_FILENAME)
    ppdb.get("main-chai")


def test_get_all():
    ppdb = PandasMmcif()
    ppdb.read_mmcif(TESTDATA_FILENAME)
    for i in ["c-alpha", "hydrogen", "main chain"]:
        ppdb.get(i)


def test_get_df():
    ppdb = PandasMmcif()
    ppdb.read_mmcif(TESTDATA_FILENAME)

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


def test_mmcif_pdb_conversion():
    """Tests conversion from mmCIF df to PDB df"""
    # Multichain test
    pdb = PandasPdb().fetch_pdb("4hhb")
    mmcif = PandasMmcif().fetch_mmcif("4hhb")
    mmcif_pdb = mmcif.convert_to_pandas_pdb()

    assert_frame_equal(
        pdb.df["ATOM"].drop(columns=["line_idx"]),
        mmcif_pdb.df["ATOM"].drop(columns=["line_idx"]),
    )
    assert_frame_equal(
        pdb.df["HETATM"].drop(columns=["line_idx"]),
        mmcif_pdb.df["HETATM"].drop(columns=["line_idx"]).reset_index(drop=True),
    )

    # single chain test
    pdb = PandasPdb().fetch_pdb("3eiy")
    mmcif = PandasMmcif().fetch_mmcif("3eiy")
    mmcif_pdb = mmcif.convert_to_pandas_pdb()

    assert_frame_equal(
        pdb.df["ATOM"].drop(columns=["line_idx"]),
        mmcif_pdb.df["ATOM"].drop(columns=["line_idx"]),
    )
    assert_frame_equal(
        pdb.df["HETATM"].drop(columns=["line_idx"]),
        mmcif_pdb.df["HETATM"].drop(columns=["line_idx"]).reset_index(drop=True),
    )
