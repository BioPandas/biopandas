# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

from biopandas.pdb import PandasPdb
from biopandas.align.tmalign import TMAlign
import warnings
import pandas as pd
import os

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), "data", "3eiy.pdb")
TESTDATA_FILENAME2 = os.path.join(
    os.path.dirname(__file__), "data", "4eiy_anisouchunk.pdb"
)

OUTFILE = os.path.join(os.path.dirname(__file__), "data", "tmp.pdb")

def test_write_pdb_to_temp_file():
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)

    tmalign = TMAlign(TESTDATA_FILENAME, TESTDATA_FILENAME2, 'C:/Users/jvarg/Downloads/Usalign')

    with tmalign.write_pdb_to_temp_file(ppdb) as temp_file:
        assert os.path.exists(temp_file.name)
        assert temp_file.name.endswith('.pdb')

        temp_content = open(temp_file.name, 'r').read()
        assert temp_content == open(TESTDATA_FILENAME, 'r').read()