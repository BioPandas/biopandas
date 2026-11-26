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


@pytest.mark.skipif(
    sys.platform == "win32",
    reason="Windows CI cert store triggers ssl.SSLError: [ASN1: NOT_ENOUGH_DATA]",
)
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


@pytest.mark.skipif(
    sys.platform == "win32",
    reason="Windows CI cert store triggers ssl.SSLError: [ASN1: NOT_ENOUGH_DATA]",
)
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
            float("nan"),
            "A",
            "N",
            "SER",
            1,
            "23",
            1.0,
            float("nan"),
            1,
            float("nan"),
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
        actual = dfs["ATOM"].loc[0, :][k]
        if pd.isna(v):
            assert pd.isna(actual), k
        else:
            assert actual == v, k


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


@pytest.mark.skipif(
    sys.platform == "win32",
    reason="Windows CI cert store triggers ssl.SSLError: [ASN1: NOT_ENOUGH_DATA]",
)
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


def test_convert_to_pandas_pdb_preserves_alt_loc(tmp_path):
    """Test convert_to_pandas_pdb preserves alternate locations."""
    mmcif = PandasMmcif().read_mmcif(TESTDATA_FILENAME2)
    expected_alt = mmcif.df["ATOM"]["label_alt_id"].fillna("").astype(str)
    assert (expected_alt != "").sum() == 186

    conv = mmcif.convert_to_pandas_pdb()
    assert (conv.df["ATOM"]["alt_loc"] == expected_alt).all()
    assert set(conv.df["ATOM"]["alt_loc"].unique()) == {"", "A", "B"}

    out_pdb = tmp_path / "4eiy_converted.pdb"
    conv.to_pdb(str(out_pdb))
    rt = PandasPdb().read_pdb(str(out_pdb))
    assert (rt.df["ATOM"]["alt_loc"].values == expected_alt.values).all()


def test_convert_to_pandas_pdb_preserves_insertion_codes(tmp_path):
    """Test convert_to_pandas_pdb preserves insertion codes and amino3to1."""
    cif_text = (
        "data_TEST\n"
        "_entry.id TEST\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_alt_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_entity_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n"
        "_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n"
        "_atom_site.occupancy\n"
        "_atom_site.B_iso_or_equiv\n"
        "_atom_site.pdbx_formal_charge\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.auth_comp_id\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.auth_atom_id\n"
        "_atom_site.pdbx_PDB_model_num\n"
        "ATOM 1 C CA . ASN H 1 1 ? 1.0 2.0 3.0 1.0 10.0 ? 52 ASN H CA 1\n"
        "ATOM 2 C CA . PRO H 1 2 A 4.0 5.0 6.0 1.0 11.0 ? 52 PRO H CA 1\n"
        "ATOM 3 C CA . LYS H 1 3 B 7.0 8.0 9.0 1.0 12.0 ? 52 LYS H CA 1\n"
        "ATOM 4 C CA . SER H 1 4 ? 10.0 11.0 12.0 1.0 13.0 ? 53 SER H CA 1\n"
    )
    mmcif = PandasMmcif().read_mmcif_from_list(cif_text)
    conv = mmcif.convert_to_pandas_pdb()

    assert list(conv.df["ATOM"]["insertion"]) == ["", "A", "B", ""]
    assert "".join(conv.amino3to1()["residue_name"].values) == "NPKS"
    assert "".join(conv.amino3to1()["residue_name"].values) == "".join(
        mmcif.amino3to1()["auth_comp_id"].values
    )

    out_pdb = tmp_path / "ins_converted.pdb"
    conv.to_pdb(str(out_pdb))
    rt = PandasPdb().read_pdb(str(out_pdb))
    assert list(rt.df["ATOM"]["insertion"]) == ["", "A", "B", ""]
    assert "".join(rt.amino3to1()["residue_name"].values) == "NPKS"


def test_convert_to_pandas_pdb_preserves_formal_charges(tmp_path):
    """Test convert_to_pandas_pdb preserves signed formal charges."""
    cif_text = (
        "data_CHG\n"
        "_entry.id CHG\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_alt_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_entity_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n"
        "_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n"
        "_atom_site.occupancy\n"
        "_atom_site.B_iso_or_equiv\n"
        "_atom_site.pdbx_formal_charge\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.auth_comp_id\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.auth_atom_id\n"
        "_atom_site.pdbx_PDB_model_num\n"
        "ATOM 1 N N 1 GLU A 1 1 ? 1.0 2.0 3.0 0.5 10.0 ? 1 GLU A N 1\n"
        "ATOM 2 O OE2 . GLU A 1 1 ? 3.0 4.0 5.0 1.0 12.0 -1 1 GLU A OE2 1\n"
        "ATOM 3 N NZ . LYS A 1 2 ? 4.0 5.0 6.0 1.0 13.0 +1 2 LYS A NZ 1\n"
        "HETATM 4 FE FE . HEM A 2 . ? 5.0 6.0 7.0 1.0 14.0 2 3 HEM A FE 1\n"
    )
    mmcif = PandasMmcif().read_mmcif_from_list(cif_text)
    conv = mmcif.convert_to_pandas_pdb()

    # Numeric label_alt_id ('1') must remain a string
    assert conv.df["ATOM"]["alt_loc"].tolist() == ["1", "", ""]
    assert pd.isna(conv.df["ATOM"]["charge"].iloc[0])
    assert conv.df["ATOM"]["charge"].iloc[1] == -1.0
    assert conv.df["ATOM"]["charge"].iloc[2] == 1.0
    assert conv.df["HETATM"]["charge"].iloc[0] == 2.0

    out_pdb = tmp_path / "charged_converted.pdb"
    conv.to_pdb(str(out_pdb))
    rt = PandasPdb().read_pdb(str(out_pdb))
    assert pd.isna(rt.df["ATOM"]["charge"].iloc[0])
    assert rt.df["ATOM"]["charge"].iloc[1] == -1.0
    assert rt.df["ATOM"]["charge"].iloc[2] == 1.0
    assert rt.df["HETATM"]["charge"].iloc[0] == 2.0


def test_convert_to_pandas_pdb_encounter_order_chain_offsets():
    """Test encounter-order TER offsets for non-alphabetical chains and MSE."""
    cif_text = (
        "data_ORD\n"
        "_entry.id ORD\n"
        "loop_\n"
        "_atom_site.group_PDB\n"
        "_atom_site.id\n"
        "_atom_site.type_symbol\n"
        "_atom_site.label_atom_id\n"
        "_atom_site.label_alt_id\n"
        "_atom_site.label_comp_id\n"
        "_atom_site.label_asym_id\n"
        "_atom_site.label_entity_id\n"
        "_atom_site.label_seq_id\n"
        "_atom_site.pdbx_PDB_ins_code\n"
        "_atom_site.Cartn_x\n"
        "_atom_site.Cartn_y\n"
        "_atom_site.Cartn_z\n"
        "_atom_site.occupancy\n"
        "_atom_site.B_iso_or_equiv\n"
        "_atom_site.pdbx_formal_charge\n"
        "_atom_site.auth_seq_id\n"
        "_atom_site.auth_comp_id\n"
        "_atom_site.auth_asym_id\n"
        "_atom_site.auth_atom_id\n"
        "_atom_site.pdbx_PDB_model_num\n"
        "ATOM 1 N N . ALA L 1 1 ? 1.0 0.0 0.0 1.0 10.0 ? 1 ALA L N 1\n"
        "HETATM 2 SE SE . MSE L 1 2 ? 2.0 0.0 0.0 1.0 10.0 ? 2 MSE L SE 1\n"
        "ATOM 3 C OXT . CYS L 1 3 ? 3.0 0.0 0.0 1.0 10.0 ? 3 CYS L OXT 1\n"
        "ATOM 4 N N . GLN H 2 1 ? 4.0 0.0 0.0 1.0 10.0 ? 1 GLN H N 1\n"
        "ATOM 5 N N . ASP A 3 1 ? 5.0 0.0 0.0 1.0 10.0 ? 1 ASP A N 1\n"
        "HETATM 6 O O . HOH W 4 . ? 6.0 0.0 0.0 1.0 10.0 ? 101 HOH W O 1\n"
    )
    mmcif = PandasMmcif().read_mmcif_from_list(cif_text)
    conv = mmcif.convert_to_pandas_pdb(offset_chains=True)

    # Chain L (0 TERs before): ATOM ids 1, 3 -> 1, 3; in-chain MSE id 2 -> 2
    # Chain H (1 TER before): ATOM id 4 -> 5
    # Chain A (2 TERs before): ATOM id 5 -> 7
    # Trailing HOH (3 TERs before): HETATM id 6 -> 9
    assert conv.df["ATOM"]["atom_number"].tolist() == [1, 3, 5, 7]
    assert conv.df["HETATM"]["atom_number"].tolist() == [2, 9]

    # Single-record selection and missing optional columns must not fail
    mmcif_opt = PandasMmcif().read_mmcif_from_list(cif_text)
    mmcif_opt.df["ATOM"] = mmcif_opt.df["ATOM"].drop(
        columns=["label_alt_id", "pdbx_PDB_ins_code", "pdbx_formal_charge"]
    )
    conv_atom_only = mmcif_opt.convert_to_pandas_pdb(records=["ATOM"])
    assert set(conv_atom_only.df.keys()) == {"ATOM"}
    assert conv_atom_only.df["ATOM"]["atom_number"].tolist() == [1, 3, 5, 7]
    assert (conv_atom_only.df["ATOM"]["alt_loc"] == "").all()
    assert (conv_atom_only.df["ATOM"]["insertion"] == "").all()
    assert conv_atom_only.df["ATOM"]["charge"].isna().all()
