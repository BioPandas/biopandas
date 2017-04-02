# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import pandas as pd
from biopandas.pdb import PandasPdb
import os
from nose.tools import raises


def test_equal():
    TESTDATA_1t48 = os.path.join(os.path.dirname(__file__), 'data',
                                                            '1t48_995.pdb')

    p1t48 = PandasPdb()
    p1t48.read_pdb(TESTDATA_1t48)
    dist = p1t48.distance(xyz=(70.785, 15.477, 23.359), record='ATOM')

    expect = pd.Series([2.533259, 1.520502, 0.000000, 1.257597, 1.252510],
                       index=[12, 13, 14, 15, 16])
    assert dist[dist < 3].all() == expect.all()
