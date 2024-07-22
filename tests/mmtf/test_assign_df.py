# BioPandas
# Author: Arian Jamasb <arian"jamasb.io>, Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import sys

if sys.version_info >= (3, 9):
    import importlib.resources as pkg_resources
else:
    import importlib_resources as pkg_resources

import tests.mmtf.data
from biopandas.mmtf import PandasMmtf
from tests.testutils import assert_raises

TEST_DATA = pkg_resources.files(tests.mmtf.data)

TESTDATA_FILENAME = str(TEST_DATA.joinpath("3eiy.mmtf"))


def test_overwrite_df():
    pdb = PandasMmtf().read_mmtf(TESTDATA_FILENAME)

    def overwrite():
        pdb.df = "bla"

    expect = (
        "Please use `PandasMmtf._df = ... ` instead\n"
        "of `PandasMmtf.df = ... ` if you are sure that\n"
        "you want to overwrite the `df` attribute."
    )

    assert_raises(AttributeError, expect, overwrite)
