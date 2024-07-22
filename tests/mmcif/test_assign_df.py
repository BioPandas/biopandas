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

import tests.mmcif.data
from biopandas.mmcif import PandasMmcif
from tests.testutils import assert_raises

TEST_DATA = pkg_resources.files(tests.mmcif.data)


def test_overwrite_df():
    data_path = str(TEST_DATA.joinpath("3eiy.cif"))
    pdb = PandasMmcif().read_mmcif(data_path)

    def overwrite():
        pdb.df = "bla"

    expect = (
        "Please use `PandasMmcif._df = ... ` instead\n"
        "of `PandasMmcif.df = ... ` if you are sure that\n"
        "you want to overwrite the `df` attribute."
    )

    assert_raises(AttributeError, expect, overwrite)
