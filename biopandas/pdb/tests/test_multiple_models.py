# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# Author: Arian Jamasb <arian@jamasb.io>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas
import os

from biopandas.pdb import PandasPdb

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), "data", "2jyf.pdb")


def test_get_index_df():
    df = PandasPdb().read_pdb(TESTDATA_FILENAME)
    idxs = df.get_model_start_end()["model_idx"]
    assert len(idxs) == max(idxs.astype(int))


def test_label_models():
    df = PandasPdb().read_pdb(TESTDATA_FILENAME)
    df.label_models()
    assert "model_id" in df.df["ATOM"].columns


def test_get_model():
    df = PandasPdb().read_pdb(TESTDATA_FILENAME)
    MODEL_INDEX = 1
    new_df = df.get_model(MODEL_INDEX)
    assert new_df.df["ATOM"]["model_id"].all() == MODEL_INDEX


def test_get_models():
    df = PandasPdb().read_pdb(TESTDATA_FILENAME)
    MODEL_INDICES = [1, 3, 5]

    new_df = df.get_models(MODEL_INDICES)
    assert new_df.df["ATOM"]["model_id"].all() in MODEL_INDICES
