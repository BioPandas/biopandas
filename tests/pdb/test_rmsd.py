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

import pytest

import tests.pdb.data
from biopandas.pdb import PandasPdb

TEST_DATA = pkg_resources.files(tests.pdb.data)

TESTDATA_1t48 = str(TEST_DATA.joinpath("1t48_995.pdb"))
TESTDATA_1t49 = str(TEST_DATA.joinpath("1t49_995.pdb"))
TESTDATA_lig1 = str(TEST_DATA.joinpath("lig_conf_1.pdb"))
TESTDATA_lig2 = str(TEST_DATA.joinpath("lig_conf_2.pdb"))

TESTDATA_rna = str(TEST_DATA.joinpath("1ehz-rna_short.pdb"))

p1t48 = PandasPdb()
p1t48.read_pdb(TESTDATA_1t48)
p1t49 = PandasPdb()
p1t49.read_pdb(TESTDATA_1t49)

pl1 = PandasPdb()
pl1.read_pdb(TESTDATA_lig1)
pl2 = PandasPdb()
pl2.read_pdb(TESTDATA_lig2)


def test_equal():
    r = PandasPdb.rmsd(p1t48.df["ATOM"], p1t48.df["ATOM"], s=None)
    assert r == 0.000, r


@pytest.mark.xfail(raises=AttributeError)
def test_wrong_arg():
    PandasPdb.rmsd(p1t48.df["ATOM"].loc[1:, :], p1t48.df["ATOM"], s="bla")


@pytest.mark.xfail(raises=AttributeError)
def test_incompatible():
    PandasPdb.rmsd(p1t48.df["ATOM"].loc[1:, :], p1t48.df["ATOM"], s=None)


@pytest.mark.xfail(raises=AttributeError)
def test_invalid_query():
    PandasPdb.rmsd(p1t48.df["ATOM"].loc[1:, :], p1t48.df["ATOM"], s="bla")


def test_protein():
    r = PandasPdb.rmsd(
        p1t48.df["ATOM"], p1t49.df["ATOM"], s="c-alpha", invert=False
    )
    assert r == 0.4785, r


def test_rna_and_nonmatching_indices():
    ehz = PandasPdb().read_pdb(TESTDATA_rna)
    at = ehz.df["ATOM"]
    a64 = at[at["residue_number"] == 64]
    a66 = at[at["residue_number"] == 66]
    r = PandasPdb.rmsd(a64, a66)
    assert r == 10.2007, r


def test_ligand():
    r = PandasPdb.rmsd(
        pl1.df["HETATM"], pl2.df["HETATM"], s="hydrogen", invert=True
    )
    assert r == 1.9959, r


def test_ligand_default():
    r = PandasPdb.rmsd(pl1.df["HETATM"], pl2.df["HETATM"], s=None)
    assert r == 2.6444, r


def test_decimals():
    r = PandasPdb.rmsd(p1t48.df["ATOM"], p1t49.df["ATOM"], decimals=3)
    assert r == 0.738, r
