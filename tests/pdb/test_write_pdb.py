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

import os
import warnings

import pandas as pd

import tests.pdb.data
from biopandas.pdb import PandasPdb

TEST_DATA = pkg_resources.files(tests.pdb.data)

TESTDATA_FILENAME = str(TEST_DATA.joinpath("3eiy.pdb"))
TESTDATA_FILENAME2 = str(TEST_DATA.joinpath("4eiy_anisouchunk.pdb"))
TESTDATA_FILENAME3 = str(TEST_DATA.joinpath("5mtn_multichain.pdb"))
OUTFILE = str(TEST_DATA.joinpath("tmp.pdb"))
OUTFILE_GZ = str(TEST_DATA.joinpath("tmp.pdb.gz"))

hetatm = ""
with open(TESTDATA_FILENAME, "r") as f:
    for line in f:
        if line.startswith("HETATM"):
            hetatm += line

with open(TESTDATA_FILENAME2, "r") as f:
    four_eiy = f.read()


def test_defaults():
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    ppdb.to_pdb(path=OUTFILE, records=None)
    with open(TESTDATA_FILENAME, "r") as f:
        f1 = f.read()
    with open(OUTFILE, "r") as f:
        f2 = f.read()
    assert f1 == f2
    os.remove(OUTFILE)


def test_nonexpected_column():
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    ppdb.df["HETATM"]["test"] = pd.Series("test", index=ppdb.df["HETATM"].index)
    with warnings.catch_warnings(record=True) as w:
        ppdb.to_pdb(path=OUTFILE, records=["HETATM"])
    with open(OUTFILE, "r") as f:
        f1 = f.read()
    os.remove(OUTFILE)
    assert f1 == hetatm


def test_records():
    """Test private _read_pdb."""
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    ppdb.to_pdb(path=OUTFILE, records=["HETATM"])
    with open(OUTFILE, "r") as f:
        f1 = f.read()
    os.remove(OUTFILE)
    assert f1 == hetatm


def test_anisou():
    """Test writing ANISOU entries."""
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME2)
    ppdb.to_pdb(path=OUTFILE, records=None)
    with open(OUTFILE, "r") as f:
        f1 = f.read()
    os.remove(OUTFILE)
    assert f1 == four_eiy


def test_add_remark():
    """Test adding a REMARK entry."""
    # Add remark
    code = 3
    remark1 = "THIS IS A HIGHLY IMPORTANT FREE-TEXT REMARK WHICH IS EXACTLY 80 CHARACTERS LONG."
    remark2 = ""
    remark3 = "THIS IS A NEXT MULTI-LINE INDENTED REMARK\n FOLLOWING THE BLANK REMARK."
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    n_atoms = len(ppdb.df["ATOM"])
    ppdb.add_remark(code, remark1)
    ppdb.add_remark(code, remark2)
    ppdb.add_remark(code, remark3, 5)
    ppdb.to_pdb(path=OUTFILE)

    # Test modified file contains remarks
    with open(OUTFILE, "r") as f:
        f1 = f.read()
    expected_substr = (
        "REMARK   3  OTHER REFINEMENT REMARKS: HYDROGENS HAVE BEEN ADDED IN THE RIDING   \n"
        "REMARK   3  POSITIONS                                                           \n"
        "REMARK   3 THIS IS A HIGHLY IMPORTANT FREE-TEXT REMARK WHICH IS EXACTLY 80      \n"
        "REMARK   3 CHARACTERS LONG.                                                     \n"
        "REMARK   3                                                                      \n"
        "REMARK   3      THIS IS A NEXT MULTI-LINE INDENTED REMARK                       \n"
        "REMARK   3      FOLLOWING THE BLANK REMARK.                                     \n"
        "REMARK   4                                                                      \n"
    )
    assert expected_substr in f1

    # Test number of atoms remained the same
    ppdb = PandasPdb()
    ppdb.read_pdb(OUTFILE)
    os.remove(OUTFILE)
    assert len(ppdb.df["ATOM"]) == n_atoms


def test_introduce_remark():
    """Test introducing a REMARK entry to the file with no remarks."""
    # Add remark
    code = 3
    remark = "THIS IS A HIGHLY IMPORTANT FREE-TEXT REMARK WHICH IS EXACTLY 80 CHARACTERS LONG."
    indent = 1
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME3)
    n_atoms = len(ppdb.df["ATOM"])
    ppdb.add_remark(code, remark, indent)
    ppdb.to_pdb(path=OUTFILE)

    # Test modified file starts with new remark
    with open(OUTFILE, "r") as f:
        f1 = f.read()
    expected_prefix = (
        "REMARK   3  THIS IS A HIGHLY IMPORTANT FREE-TEXT REMARK WHICH IS EXACTLY 80     \n"
        "REMARK   3  CHARACTERS LONG.                                                    \n"
    )
    assert f1.startswith(expected_prefix)

    # Test number of atoms remained the same
    ppdb = PandasPdb()
    ppdb.read_pdb(OUTFILE)
    os.remove(OUTFILE)
    assert len(ppdb.df["ATOM"]) == n_atoms


def test_b_factor_shift():
    """Test b_factor shifting one white space when saving the fetched pdb."""
    ppdb = PandasPdb()
    ppdb.fetch_pdb("2e28")
    ppdb.to_pdb(path=OUTFILE, records=None)
    tmp_df = ppdb.read_pdb(path=OUTFILE).df["ATOM"]
    os.remove(OUTFILE)
    assert tmp_df[
        tmp_df["element_symbol"].isnull() | (tmp_df["element_symbol"] == "")
    ].empty
    assert not tmp_df[tmp_df["blank_4"].isnull() | (tmp_df["blank_4"] == "")].empty


def test_to_pdb_stream():
    """Test public write_pdb_stream"""
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    stream = ppdb.to_pdb_stream()

    lines_to_check = open(TESTDATA_FILENAME).read().split("\n")
    lines_to_check = [
        line for line in lines_to_check if line.startswith(("ATOM", "HETATM"))
    ]
    lines_to_check.append("")

    source_pdb = "\n".join(lines_to_check)
    assert stream.read() == source_pdb
