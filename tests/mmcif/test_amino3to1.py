# BioPandas
# Author: Arian Jamasb <arian@jamasb.io>, Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import sys

if sys.version_info >= (3, 9):
    import importlib.resources as pkg_resources
else:
    import importlib_resources as pkg_resources

import numpy as np

import tests.mmcif.data
from biopandas.mmcif import PandasMmcif

TEST_DATA = pkg_resources.files(tests.mmcif.data)


def test_defaults():
    TESTDATA_1t48 = str(TEST_DATA.joinpath("1t48.cif"))
    p1t48 = PandasMmcif()
    p1t48.read_mmcif(TESTDATA_1t48)
    expect_res = [
        "M",
        "E",
        "M",
        "E",
        "K",
        "E",
        "F",
        "E",
        "Q",
        "I",
        "D",
        "K",
        "S",
        "G",
        "S",
        "W",
        "A",
        "A",
        "I",
        "Y",
        "Q",
        "D",
        "I",
        "R",
        "H",
        "E",
        "A",
        "S",
        "D",
        "F",
        "P",
        "C",
        "R",
        "V",
        "A",
        "K",
        "L",
        "P",
        "K",
        "N",
        "K",
        "N",
        "R",
        "N",
        "R",
        "Y",
        "R",
        "D",
        "V",
        "S",
        "P",
        "F",
        "D",
        "H",
        "S",
        "R",
        "I",
        "K",
        "L",
        "H",
        "Q",
        "E",
        "D",
        "N",
        "D",
        "Y",
        "I",
        "N",
        "A",
        "S",
        "L",
        "I",
        "K",
        "M",
        "E",
        "E",
        "A",
        "Q",
        "R",
        "S",
        "Y",
        "I",
        "L",
        "T",
        "Q",
        "G",
        "P",
        "L",
        "P",
        "N",
        "T",
        "C",
        "G",
        "H",
        "F",
        "W",
        "E",
        "M",
        "V",
        "W",
        "E",
        "Q",
        "K",
        "S",
        "R",
        "G",
        "V",
        "V",
        "M",
        "L",
        "N",
        "R",
        "V",
        "M",
        "E",
        "K",
        "G",
        "S",
        "L",
        "K",
        "C",
        "A",
        "Q",
        "Y",
        "W",
        "P",
        "Q",
        "K",
        "E",
        "E",
        "K",
        "E",
        "M",
        "I",
        "F",
        "E",
        "D",
        "T",
        "N",
        "L",
        "K",
        "L",
        "T",
        "L",
        "I",
        "S",
        "E",
        "D",
        "I",
        "K",
        "S",
        "Y",
        "Y",
        "T",
        "V",
        "R",
        "Q",
        "L",
        "E",
        "L",
        "E",
        "N",
        "L",
        "T",
        "T",
        "Q",
        "E",
        "T",
        "R",
        "E",
        "I",
        "L",
        "H",
        "F",
        "H",
        "Y",
        "T",
        "T",
        "W",
        "P",
        "D",
        "F",
        "G",
        "V",
        "P",
        "E",
        "S",
        "P",
        "A",
        "S",
        "F",
        "L",
        "N",
        "F",
        "L",
        "F",
        "K",
        "V",
        "R",
        "E",
        "S",
        "G",
        "S",
        "L",
        "S",
        "P",
        "E",
        "H",
        "G",
        "P",
        "V",
        "V",
        "V",
        "H",
        "C",
        "S",
        "A",
        "G",
        "I",
        "G",
        "R",
        "S",
        "G",
        "T",
        "F",
        "C",
        "L",
        "A",
        "D",
        "T",
        "C",
        "L",
        "L",
        "L",
        "M",
        "D",
        "K",
        "R",
        "K",
        "D",
        "P",
        "S",
        "S",
        "V",
        "D",
        "I",
        "K",
        "K",
        "V",
        "L",
        "L",
        "E",
        "M",
        "R",
        "K",
        "F",
        "R",
        "M",
        "G",
        "L",
        "I",
        "Q",
        "T",
        "A",
        "D",
        "Q",
        "L",
        "R",
        "F",
        "S",
        "Y",
        "L",
        "A",
        "V",
        "I",
        "E",
        "G",
        "A",
        "K",
        "F",
        "I",
        "M",
    ]

    transl = p1t48.amino3to1()
    expect_chain = ["A" for _ in range(transl.shape[0])]
    got_chain = list(transl["auth_asym_id"].values)
    got_res = list(transl["auth_comp_id"].values)

    assert expect_chain == got_chain
    assert expect_res == got_res


def test_sameindex():
    TESTDATA_1t48 = str(TEST_DATA.joinpath("1t48.cif"))
    p1t48 = PandasMmcif()
    p1t48.read_mmcif(TESTDATA_1t48)
    p1t48.df["ATOM"].index = np.zeros(p1t48.df["ATOM"].shape[0], dtype=int)

    expect_res = [
        "M",
        "E",
        "M",
        "E",
        "K",
        "E",
        "F",
        "E",
        "Q",
        "I",
        "D",
        "K",
        "S",
        "G",
        "S",
        "W",
        "A",
        "A",
        "I",
        "Y",
        "Q",
        "D",
        "I",
        "R",
        "H",
        "E",
        "A",
        "S",
        "D",
        "F",
        "P",
        "C",
        "R",
        "V",
        "A",
        "K",
        "L",
        "P",
        "K",
        "N",
        "K",
        "N",
        "R",
        "N",
        "R",
        "Y",
        "R",
        "D",
        "V",
        "S",
        "P",
        "F",
        "D",
        "H",
        "S",
        "R",
        "I",
        "K",
        "L",
        "H",
        "Q",
        "E",
        "D",
        "N",
        "D",
        "Y",
        "I",
        "N",
        "A",
        "S",
        "L",
        "I",
        "K",
        "M",
        "E",
        "E",
        "A",
        "Q",
        "R",
        "S",
        "Y",
        "I",
        "L",
        "T",
        "Q",
        "G",
        "P",
        "L",
        "P",
        "N",
        "T",
        "C",
        "G",
        "H",
        "F",
        "W",
        "E",
        "M",
        "V",
        "W",
        "E",
        "Q",
        "K",
        "S",
        "R",
        "G",
        "V",
        "V",
        "M",
        "L",
        "N",
        "R",
        "V",
        "M",
        "E",
        "K",
        "G",
        "S",
        "L",
        "K",
        "C",
        "A",
        "Q",
        "Y",
        "W",
        "P",
        "Q",
        "K",
        "E",
        "E",
        "K",
        "E",
        "M",
        "I",
        "F",
        "E",
        "D",
        "T",
        "N",
        "L",
        "K",
        "L",
        "T",
        "L",
        "I",
        "S",
        "E",
        "D",
        "I",
        "K",
        "S",
        "Y",
        "Y",
        "T",
        "V",
        "R",
        "Q",
        "L",
        "E",
        "L",
        "E",
        "N",
        "L",
        "T",
        "T",
        "Q",
        "E",
        "T",
        "R",
        "E",
        "I",
        "L",
        "H",
        "F",
        "H",
        "Y",
        "T",
        "T",
        "W",
        "P",
        "D",
        "F",
        "G",
        "V",
        "P",
        "E",
        "S",
        "P",
        "A",
        "S",
        "F",
        "L",
        "N",
        "F",
        "L",
        "F",
        "K",
        "V",
        "R",
        "E",
        "S",
        "G",
        "S",
        "L",
        "S",
        "P",
        "E",
        "H",
        "G",
        "P",
        "V",
        "V",
        "V",
        "H",
        "C",
        "S",
        "A",
        "G",
        "I",
        "G",
        "R",
        "S",
        "G",
        "T",
        "F",
        "C",
        "L",
        "A",
        "D",
        "T",
        "C",
        "L",
        "L",
        "L",
        "M",
        "D",
        "K",
        "R",
        "K",
        "D",
        "P",
        "S",
        "S",
        "V",
        "D",
        "I",
        "K",
        "K",
        "V",
        "L",
        "L",
        "E",
        "M",
        "R",
        "K",
        "F",
        "R",
        "M",
        "G",
        "L",
        "I",
        "Q",
        "T",
        "A",
        "D",
        "Q",
        "L",
        "R",
        "F",
        "S",
        "Y",
        "L",
        "A",
        "V",
        "I",
        "E",
        "G",
        "A",
        "K",
        "F",
        "I",
        "M",
    ]
    transl = p1t48.amino3to1()
    expect_chain = ["A" for _ in range(transl.shape[0])]
    got_chain = list(transl["auth_asym_id"].values)
    got_res = list(transl["auth_comp_id"].values)

    assert expect_chain == got_chain
    assert expect_res == got_res


def test_multichain():
    TESTDATA_5mtn = str(TEST_DATA.joinpath("5mtn_multichain.cif"))
    mtn = PandasMmcif()
    mtn.read_mmcif(TESTDATA_5mtn)
    expect_res_a = [
        "S",
        "L",
        "E",
        "P",
        "E",
        "P",
        "W",
        "F",
        "F",
        "K",
        "N",
        "L",
        "S",
        "R",
        "K",
        "D",
        "A",
        "E",
        "R",
        "Q",
        "L",
        "L",
        "A",
        "P",
        "G",
        "N",
        "T",
        "H",
        "G",
        "S",
        "F",
        "L",
        "I",
        "R",
        "E",
        "S",
        "E",
        "S",
        "T",
        "A",
        "G",
        "S",
        "F",
        "S",
        "L",
        "S",
        "V",
        "R",
        "D",
        "F",
        "D",
        "Q",
        "G",
        "E",
        "V",
        "V",
        "K",
        "H",
        "Y",
        "K",
        "I",
        "R",
        "N",
        "L",
        "D",
        "N",
        "G",
        "G",
        "F",
        "Y",
        "I",
        "S",
        "P",
        "R",
        "I",
        "T",
        "F",
        "P",
        "G",
        "L",
        "H",
        "E",
        "L",
        "V",
        "R",
        "H",
        "Y",
        "T",
    ]
    expect_res_b = [
        "S",
        "V",
        "S",
        "S",
        "V",
        "P",
        "T",
        "K",
        "L",
        "E",
        "V",
        "V",
        "A",
        "A",
        "T",
        "P",
        "T",
        "S",
        "L",
        "L",
        "I",
        "S",
        "W",
        "D",
        "A",
        "P",
        "A",
        "V",
        "T",
        "V",
        "V",
        "Y",
        "Y",
        "L",
        "I",
        "T",
        "Y",
        "G",
        "E",
        "T",
        "G",
        "S",
        "P",
        "W",
        "P",
        "G",
        "G",
        "Q",
        "A",
        "F",
        "E",
        "V",
        "P",
        "G",
        "S",
        "K",
        "S",
        "T",
        "A",
        "T",
        "I",
        "S",
        "G",
        "L",
        "K",
        "P",
        "G",
        "V",
        "D",
        "Y",
        "T",
        "I",
        "T",
        "V",
        "Y",
        "A",
        "H",
        "R",
        "S",
        "S",
        "Y",
        "G",
        "Y",
        "S",
        "E",
        "N",
        "P",
        "I",
        "S",
        "I",
        "N",
        "Y",
        "R",
        "T",
    ]

    transl = mtn.amino3to1()

    expect_chain = ["A" for _ in range(88)] + ["B" for _ in range(94)]
    got_chain = list(transl["auth_asym_id"].values)

    got_res_a = list(
        transl.loc[transl["auth_asym_id"] == "A", "auth_comp_id"].values
    )
    got_res_b = list(
        transl.loc[transl["auth_asym_id"] == "B", "auth_comp_id"].values
    )

    assert expect_chain == got_chain
    assert expect_res_a == got_res_a
    assert expect_res_b == got_res_b


def test_pdb_with_insertion_codes():
    PDB_2D7T_PATH = str(TEST_DATA.joinpath("2d7t.cif"))

    ppdb = PandasMmcif().read_mmcif(PDB_2D7T_PATH)
    sequence = ppdb.amino3to1()
    assert "".join(sequence[50:60]["auth_comp_id"].values) == "INPKSGDTNY"
