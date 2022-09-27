# BioPandas
# Author: Arian Jamasb <arian"jamasb.io>, Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import os

from biopandas.mmtf import PandasMmtf
from biopandas.testutils import assert_raises

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), "data", "3eiy.mmtf")


def test_overwrite_df():
    data_path = os.path.join(os.path.dirname(__file__), "data", "3eiy.mmtf")
    pdb = PandasMmtf().read_mmtf(data_path)

    def overwrite():
        pdb.df = "bla"

    expect = (
        "Please use `PandasMmtf._df = ... ` instead\n"
        "of `PandasMmtf.df = ... ` if you are sure that\n"
        "you want to overwrite the `df` attribute."
    )

    assert_raises(AttributeError, expect, overwrite)
