# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas


from biopandas.pdb import PandasPDB
import os
import numpy as np
import pandas as pd
from nose.tools import raises


TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'data', '3eiy.pdb')
TESTDATA_FILENAME_GZ = os.path.join(os.path.dirname(__file__), 'data', '3eiy.pdb.gz')

ATOM_DF_COLUMNS = ['record_name', 'atom_number', 'blank_1',
                 'atom_name', 'alt_loc', 'resi_name',
                 'blank_2', 'chain_id', 'resi_number',
                 'resi_insert_code', 'blank_3',
                 'x_coord', 'y_coord', 'z_coord',
                 'occupancy', 'b_factor', 'blank_4',
                 'segment_id', 'element_symbol',
                 'charge', 'line_idx']

with open(TESTDATA_FILENAME, 'r') as f:
    three_eiy = f.read()

def test__read_pdb():
    """Test private _read_pdb"""
    ppdb = PandasPDB()
    txt = ppdb._read_pdb(TESTDATA_FILENAME)
    print(txt)
    assert txt == three_eiy

def test_fetch_pdb():
    """Test fetch_pdb"""
    ppdb = PandasPDB()
    txt = ppdb._fetch_pdb('3eiy')
    with open('./text.pdb', 'w') as f:
        f.write(txt)
    txt[:100] == three_eiy[:100]
    ppdb.fetch_pdb('3eiy')
    assert ppdb.pdb_text == txt
    txt = ppdb._fetch_pdb('3ey')
    err = "We're sorry, but the requested file is not available"
    assert err in txt

def test__read_pdb_gz():
    """Test public _read_pdb with gzip files"""
    ppdb = PandasPDB()
    txt = ppdb._read_pdb(TESTDATA_FILENAME_GZ)
    assert txt == three_eiy

def test__construct_df():
    """Test pandas dataframe construction"""
    ppdb = PandasPDB()
    dfs = ppdb._construct_df(three_eiy.splitlines())
    assert set(dfs.keys()) == {'OTHERS', 'ATOM', 'ANISOU', 'HETATM'}
    assert set(dfs['ATOM'].columns) == set(ATOM_DF_COLUMNS)
    assert set(dfs['HETATM'].columns) == set(ATOM_DF_COLUMNS)
    assert set(dfs['ANISOU'].columns) == set(ATOM_DF_COLUMNS)
    exp = pd.Series(np.array(['ATOM', 1, '', 'N', '', 'SER', '', 'A', 2, '', '',
              2.527, 54.656, -1.667, 1.0, 52.73, '', '', 'N', None, 609]),
          index=['record_name', 'atom_number', 'blank_1',
                 'atom_name', 'alt_loc', 'resi_name',
                 'blank_2', 'chain_id', 'resi_number',
                 'resi_insert_code', 'blank_3',
                 'x_coord', 'y_coord', 'z_coord',
                 'occupancy', 'b_factor', 'blank_4',
                 'segment_id', 'element_symbol',
                 'charge', 'line_idx'])
    assert exp.equals(dfs['ATOM'].loc[0, :])

def test_read_pdb():
    """Test public read_pdb"""
    ppdb = PandasPDB()
    ppdb.read_pdb(TESTDATA_FILENAME)
    assert ppdb.pdb_text == three_eiy
    assert ppdb.code == '3eiy', ppdb.code

@raises(AttributeError)
def test_get_exceptions():
    ppdb = PandasPDB()
    ppdb.read_pdb(TESTDATA_FILENAME)
    ppdb.get('main-chai')

def test_get_all():
    ppdb = PandasPDB()
    ppdb.read_pdb(TESTDATA_FILENAME)
    for i in ['c-alpha', 'hydrogen', 'main chain']:
        ppdb.get(i)

def test_get_df():
    ppdb = PandasPDB()
    ppdb.read_pdb(TESTDATA_FILENAME)

    shape = ppdb.get('c-alpha').shape
    assert shape == (174, 21), shape

    shape = ppdb.get('hydrogen', invert=True).shape
    assert shape == (1330, 21), shape

    shape = ppdb.get('hydrogen').shape
    assert shape == (0, 21), shape

    shape = ppdb.get('main chain').shape
    assert shape == (696, 21), shape
