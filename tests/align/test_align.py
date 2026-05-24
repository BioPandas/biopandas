# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import pytest
from biopandas.pdb import PandasPdb
from biopandas.align import Align
import numpy as np
import os

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), "data", "3eiy.pdb")
TESTDATA_FILENAME2 = os.path.join(
    os.path.dirname(__file__), "data", "4eiy_anisouchunk.pdb"
)

OUTFILE = os.path.join(os.path.dirname(__file__), "data", "tmp.pdb")

def test_write_pdb_to_temp_file():
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)

    align = Align()

    with align.write_pdb_to_temp_file(ppdb) as temp_file:
        assert os.path.exists(temp_file.name)
        assert temp_file.name.endswith('.pdb')

        temp_content = open(temp_file.name, 'r').read()
        assert temp_content == open(TESTDATA_FILENAME, 'r').read()

def test_transform():
    align = Align()
    coords = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    translation = np.array([3, 2, 1])

    transformed_coords = align.transform(coords, matrix, translation)
    target_coords = np.array([[4, 4, 4], [7, 7, 7], [10, 10, 10]])

    assert np.array_equal(transformed_coords, target_coords)

def test_filter_and_validate_chain():
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    print(TESTDATA_FILENAME)
    align = Align()

    filtered_pdb = align.filter_and_validate_chain(ppdb, 'A')
    assert filtered_pdb.df['ATOM']['chain_id'].unique() == ['A']

    with pytest.raises(ValueError):
        align.filter_and_validate_chain(ppdb, 'B')


