# BioPandas
# Author: Arian Jamasb <arian@jamasb.io>, Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import os

import pandas as pd
from biopandas.mmtf import PandasMmtf


def test_equal():
    TESTDATA_1t48 = os.path.join(os.path.dirname(__file__), "data", "1t48.mmtf")

    p1t48 = PandasMmtf()
    p1t48.read_mmtf(TESTDATA_1t48)
    dist = p1t48.distance(xyz=(70.785, 15.477, 23.359), records=("ATOM",))

    expect = pd.Series(
        [2.533259, 1.520502, 0.000000, 1.257597, 1.252510], index=[12, 13, 14, 15, 16]
    )
    assert dist[dist < 3].all() == expect.all()


def test_deprecated_str_arg():
    TESTDATA_1t48 = os.path.join(os.path.dirname(__file__), "data", "1t48.mmtf")

    p1t48 = PandasMmtf()
    p1t48.read_mmtf(TESTDATA_1t48)
    dist = p1t48.distance(xyz=(70.785, 15.477, 23.359), records="ATOM")

    expect = pd.Series(
        [2.533259, 1.520502, 0.000000, 1.257597, 1.252510], index=[12, 13, 14, 15, 16]
    )
    assert dist[dist < 3].all() == expect.all()


def test_use_external_df():
    TESTDATA_1t48 = os.path.join(os.path.dirname(__file__), "data", "1t48.mmtf")

    p1t48 = PandasMmtf()
    p1t48.read_mmtf(TESTDATA_1t48)
    new_df = p1t48.df["ATOM"].iloc[:-1, :].copy()
    dist = PandasMmtf.distance_df(df=new_df, xyz=(70.785, 15.477, 23.359))

    expect = pd.Series([2.533259, 1.520502, 0.000000, 1.257597], index=[12, 13, 14, 15])
    assert dist[dist < 3].all() == expect.all()
