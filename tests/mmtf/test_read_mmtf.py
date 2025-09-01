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
import unittest

import pandas as pd

import tests.mmtf.data
import tests.pdb.data
from biopandas.mmtf import PandasMmtf
from biopandas.pdb import PandasPdb

TEST_DATA = pkg_resources.files(tests.mmtf.data)
PDB_TEST_DATA = pkg_resources.files(tests.pdb.data)

MMTF_TESTDATA_FILENAME = str(TEST_DATA.joinpath("3eiy.mmtf"))
MMTF_TESTDATA_FILENAME_GZ = str(TEST_DATA.joinpath("3eiy.mmtf.gz"))

PDB_TESTDATA_FILENAME = str(PDB_TEST_DATA.joinpath("3eiy.pdb"))
PDB_TESTDATA_FILENAME_GZ = str(PDB_TEST_DATA.joinpath("3eiy.pdb.gz"))


ATOM_DF_COLUMNS = [
    "record_name",
    "atom_number",
    "atom_name",
    # "alt_loc",
    "residue_name",
    "chain_id",
    "residue_number",
    # "insertion",
    "x_coord",
    "y_coord",
    "z_coord",
    "occupancy",
    "b_factor",
    "element_symbol",
    # "charge",
]


@unittest.skip(reason="PDB No longer serves MMTF files.")
def test_fetch_pdb():
    """Test fetch_pdb"""
    ppdb = PandasMmtf()
    ppdb.fetch_mmtf("3eiy")
    assert max(ppdb.df["ATOM"].residue_number) == 175


def test__read_mmtf():
    """Test public _read_pdb with gzip files"""
    pmmtf = PandasMmtf()
    ppdb = PandasPdb()
    pmmtf.read_mmtf(MMTF_TESTDATA_FILENAME)

    ppdb = ppdb.read_pdb(PDB_TESTDATA_FILENAME)

    pd.testing.assert_frame_equal(
        pmmtf.df["ATOM"][ATOM_DF_COLUMNS].reset_index(drop=True),
        ppdb.df["ATOM"][ATOM_DF_COLUMNS].reset_index(drop=True),
    )

    ATOM_DF_COLUMNS.remove("atom_number")
    ATOM_DF_COLUMNS.remove("element_symbol")
    pd.testing.assert_frame_equal(
        pmmtf.df["HETATM"][ATOM_DF_COLUMNS].reset_index(drop=True),
        ppdb.df["HETATM"][ATOM_DF_COLUMNS].reset_index(drop=True),
    )


def test__read_mmtf_gz():
    """Test public _read_pdb with gzip files"""
    pmmtf = PandasMmtf()
    ppdb = PandasPdb()
    pmmtf.read_mmtf(MMTF_TESTDATA_FILENAME_GZ)
    ppdb = ppdb.read_pdb(PDB_TESTDATA_FILENAME_GZ)

    pmmtf.df["ATOM"].alt_loc.replace("\x00", "", inplace=True)
    pmmtf.df["HETATM"].alt_loc.replace("\x00", "", inplace=True)

    pd.testing.assert_frame_equal(
        pmmtf.df["ATOM"][ATOM_DF_COLUMNS].reset_index(drop=True),
        ppdb.df["ATOM"][ATOM_DF_COLUMNS].reset_index(drop=True),
    )
    # pd.testing.assert_frame_equal(
    #    pmmtf.df["HETATM"][ATOM_DF_COLUMNS].reset_index(drop=True),
    #    ppdb.df["HETATM"][ATOM_DF_COLUMNS].reset_index(drop=True),
    #    )


def test_read_mmtf():
    """Test public read_pdb"""
    ppdb = PandasMmtf()
    ppdb.read_mmtf(MMTF_TESTDATA_FILENAME)
    assert ppdb.mmtf_path == MMTF_TESTDATA_FILENAME
