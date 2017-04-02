# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

from biopandas.pdb import PandasPdb
import warnings
import pandas as pd
import os


TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'data', '3eiy.pdb')
TESTDATA_FILENAME2 = os.path.join(os.path.dirname(__file__), 'data',
                                  '4eiy_anisouchunk.pdb')
OUTFILE = os.path.join(os.path.dirname(__file__), 'data', 'tmp.pdb')
OUTFILE_GZ = os.path.join(os.path.dirname(__file__), 'data', 'tmp.pdb.gz')

hetatm = ''
with open(TESTDATA_FILENAME, 'r') as f:
    for line in f:
        if line.startswith('HETATM'):
            hetatm += line

with open(TESTDATA_FILENAME2, 'r') as f:
    four_eiy = f.read()


def test_defaults():
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    ppdb.to_pdb(path=OUTFILE, records=None)
    with open(TESTDATA_FILENAME, 'r') as f:
        f1 = f.read()
    with open(OUTFILE, 'r') as f:
        f2 = f.read()
    assert f1 == f2
    os.remove(OUTFILE)


def test_nonexpected_column():
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    ppdb.df['HETATM']['test'] = pd.Series('test',
                                          index=ppdb.df['HETATM'].index)
    with warnings.catch_warnings(record=True) as w:
        ppdb.to_pdb(path=OUTFILE, records=['HETATM'])
    with open(OUTFILE, 'r') as f:
        f1 = f.read()
    os.remove(OUTFILE)
    assert f1 == hetatm


def test_records():
    """Test private _read_pdb."""
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    ppdb.to_pdb(path=OUTFILE, records=['HETATM'])
    with open(OUTFILE, 'r') as f:
        f1 = f.read()
    os.remove(OUTFILE)
    assert f1 == hetatm


def test_anisou():
    """Test writing ANISOU entries."""
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME2)
    ppdb.to_pdb(path=OUTFILE, records=None)
    with open(OUTFILE, 'r') as f:
        f1 = f.read()
    os.remove(OUTFILE)
    assert f1 == four_eiy
