# BioPandas
# Author: Arian Jamasb <arian@jamasb.io> Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import sys

if sys.version_info >= (3, 9):
    import importlib.resources as pkg_resources
else:
    import importlib_resources as pkg_resources

import pytest

import tests.mmtf.data
from biopandas.mmtf import PandasMmtf

TEST_DATA = pkg_resources.files(tests.mmtf.data)

TESTDATA_1t48 = str(TEST_DATA.joinpath("1t48.mmtf"))
TESTDATA_1t49 = str(TEST_DATA.joinpath("1t49.mmtf"))
# TESTDATA_lig1 = str(TEST_DATA.joinpath("lig_conf_1.pdb"))
# TESTDATA_lig2 = str(TEST_DATA.joinpath("lig_conf_2.pdb"))

TESTDATA_rna = str(TEST_DATA.joinpath("1ehz.mmtf"))

p1t48 = PandasMmtf()
p1t48.read_mmtf(TESTDATA_1t48)
p1t49 = PandasMmtf()
p1t49.read_mmtf(TESTDATA_1t49)
# Subset to the first 995 atoms for consistency with PDB tests
p1t48.df["ATOM"] = p1t48.df["ATOM"].loc[p1t48.df["ATOM"].atom_number <= 995]
p1t49.df["ATOM"] = p1t49.df["ATOM"].loc[p1t49.df["ATOM"].atom_number <= 995]

# pl1 = PandasPdb()
# pl1.read_pdb(TESTDATA_lig1)
# pl2 = PandasPdb()
# pl2.read_pdb(TESTDATA_lig2)


def test_equal():
    r = PandasMmtf.rmsd(p1t48.df["ATOM"], p1t48.df["ATOM"], s=None)
    assert r == 0.000, r


@pytest.mark.xfail(raises=AttributeError)
def test_wrong_arg():
    PandasMmtf.rmsd(p1t48.df["ATOM"].loc[1:, :], p1t48.df["ATOM"], s="bla")


@pytest.mark.xfail(raises=AttributeError)
def test_incompatible():
    PandasMmtf.rmsd(p1t48.df["ATOM"].loc[1:, :], p1t48.df["ATOM"], s=None)


@pytest.mark.xfail(raises=AttributeError)
def test_invalid_query():
    PandasMmtf.rmsd(p1t48.df["ATOM"].loc[1:, :], p1t48.df["ATOM"], s="bla")


def test_protein():
    r = PandasMmtf.rmsd(
        p1t48.df["ATOM"], p1t49.df["ATOM"], s="c-alpha", invert=False
    )
    assert r == 0.4785, r


def test_rna_and_nonmatching_indices():
    ehz = PandasMmtf().read_mmtf(TESTDATA_rna)
    at = ehz.df["ATOM"]
    a64 = at[at["residue_number"] == 64]
    a66 = at[at["residue_number"] == 66]
    r = PandasMmtf.rmsd(a64, a66)
    assert r == 10.2007, r


# Skipping ligand tests as the ligand coordinates used in PDB test are unknown
# def test_ligand():
#    r = PandasMmtf.rmsd(pl1.df["HETATM"], pl2.df["HETATM"], s="hydrogen", invert=True)
#    assert r == 1.9959, r


# def test_ligand_default():
#    r = PandasMmtf.rmsd(pl1.df["HETATM"], pl2.df["HETATM"], s=None)
#    assert r == 2.6444, r
