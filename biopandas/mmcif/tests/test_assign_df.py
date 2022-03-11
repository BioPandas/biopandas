# BioPandas
# Author: Arian Jamasb <arian"jamasb.io>, Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import os

from biopandas.mmcif import PandasMMCIF
from biopandas.testutils import assert_raises

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), "data", "3eiy.cif")


def test_overwrite_df():
    data_path = os.path.join(os.path.dirname(__file__), "data", "3eiy.cif")
    pdb = PandasMMCIF().read_mmcif(data_path)

    def overwrite():
        pdb.df = "bla"

    expect = (
        "Please use `PandasMMCIF._df = ... ` instead\n"
        "of `PandasMMCIF.df = ... ` if you are sure that\n"
        "you want to overwrite the `df` attribute."
    )

    assert_raises(AttributeError, expect, overwrite)
