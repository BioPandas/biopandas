# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# Author: Arian Jamasb <arian@jamasb.io>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import sys

if sys.version_info >= (3, 9):
    import importlib.resources as pkg_resources
else:
    import importlib_resources as pkg_resources
import os

from pandas.testing import assert_frame_equal

import tests.mmtf.data
from biopandas.mmtf import PandasMmtf

TEST_DATA = pkg_resources.files(tests.mmtf.data)

TESTDATA_FILENAME = str(TEST_DATA.joinpath("2jyf.mmtf"))


def test_label_models():
    df = PandasMmtf().read_mmtf(TESTDATA_FILENAME)
    assert "model_id" in df.df["ATOM"].columns


def test_get_model():
    df = PandasMmtf().read_mmtf(TESTDATA_FILENAME)
    MODEL_INDEX = 1
    new_df = df.get_model(MODEL_INDEX)
    assert new_df.df["ATOM"]["model_id"].all() == MODEL_INDEX


def test_get_models():
    df = PandasMmtf().read_mmtf(TESTDATA_FILENAME)
    MODEL_INDICES = [1, 3, 5]

    df = df.get_models(MODEL_INDICES)
    assert df.df["ATOM"]["model_id"].all() in MODEL_INDICES

    # Test a round trip through biopandas
    df.to_mmtf("test.mmtf")
    written = PandasMmtf().read_mmtf("test.mmtf")

    # Note: No way to preserve model ID as far as I can tell
    assert_frame_equal(
        df.df["ATOM"].drop("model_id", axis=1).reset_index(drop=True),
        written.df["ATOM"].drop("model_id", axis=1).reset_index(drop=True),
    )
    assert_frame_equal(
        df.df["HETATM"].drop("model_id", axis=1).reset_index(drop=True),
        written.df["HETATM"].drop("model_id", axis=1).reset_index(drop=True),
    )

    # Clean
    os.remove("test.mmtf")
