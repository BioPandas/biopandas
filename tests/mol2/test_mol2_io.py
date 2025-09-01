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
from biopandas.mol2.mol2_io import split_multimol2
from tests.testutils import assert_raises

TEST_DATA = pkg_resources.files(tests.mol2.data)


def test_split_multimol2():
    all_mol2 = []
    for i in split_multimol2(str(TEST_DATA.joinpath("40_mol2_files.mol2"))):
        all_mol2.append(i[0])
    assert all_mol2[1] == "ZINC04084113"
    assert len(all_mol2) == 40


def test_split_multimol2_wrong_format():
    expect = (
        "Wrong file format;" "allowed file formats are .mol2 and .mol2.gz."
    )

    def run_code():
        next(split_multimol2("40_mol2_files.pdb"))

    assert_raises(ValueError, expect, run_code)


def test_split_multimol2_gz():
    all_mol2 = []
    for i in split_multimol2(str(TEST_DATA.joinpath("40_mol2_files.mol2.gz"))):
        all_mol2.append(i[0])
    assert all_mol2[1].decode() == "ZINC04084113"
    assert len(all_mol2) == 40
