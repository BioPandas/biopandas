import os
import unittest

import pandas as pd
from pandas.testing import assert_frame_equal

from biopandas.mmtf.pandas_mmtf import PandasMmtf, write_mmtf


@unittest.skip(reason="PDB No longer serves MMTF files.")
def test_write_mmtf_bp():
    PDB_CODES = [
        "4hhb",
        "3eiy",
        "1t48",
        "1ehz",
        "4ggb",
        "1bxa",
        "1cbn",
        "1rcf",
    ]
    for pdb in PDB_CODES:
        pm1 = PandasMmtf().fetch_mmtf(pdb)
        pm1.to_mmtf("test.mmtf")
        assert os.path.exists("test.mmtf")

        pm2 = PandasMmtf().read_mmtf("test.mmtf")
        assert_frame_equal(
            pm1.df["ATOM"].reset_index(drop=True),
            pm2.df["ATOM"].reset_index(drop=True),
        )
        assert_frame_equal(
            pm1.df["HETATM"].reset_index(drop=True),
            pm2.df["HETATM"].reset_index(drop=True),
        )

    os.remove("test.mmtf")


@unittest.skip(reason="PDB No longer serves MMTF files.")
def test_write_mmtf():
    PDB_CODES = [
        "4hhb",
        "3eiy",
        "1t48",
        "1ehz",
        "4ggb",
        "1bxa",
        "1cbn",
        "1rcf",
    ]
    for pdb in PDB_CODES:
        pm1 = PandasMmtf().fetch_mmtf(pdb)
        write_mmtf(pd.concat([pm1.df["ATOM"], pm1.df["HETATM"]]), "test.mmtf")
        assert os.path.exists("test.mmtf")

        pm2 = PandasMmtf().read_mmtf("test.mmtf")
        assert_frame_equal(
            pm1.df["ATOM"].reset_index(drop=True),
            pm2.df["ATOM"].reset_index(drop=True),
        )
        assert_frame_equal(
            pm1.df["HETATM"].reset_index(drop=True),
            pm2.df["HETATM"].reset_index(drop=True),
        )

    os.remove("test.mmtf")
