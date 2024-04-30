# BioPandas extension for
# Author: Sebastian Raschka, Julia Varga
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

from biopandas.stack.stack import PandasPdbStack
import os
import difflib


TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), "data", "3eiy.pdb")
TESTDATA_FILENAME2 = os.path.join(
    os.path.dirname(__file__), "data", "4eiy_anisouchunk.pdb"
)

TESTDATA_FILENAME_GZ = os.path.join(os.path.dirname(__file__), "data", "3eiy.pdb.gz")
TESTDATA_FILENAME_AF2_V4 = os.path.join(
    os.path.dirname(__file__), "data", "AF-Q5VSL9-F1-model_v4.pdb"
)

def test_write_pdb():
    stack = PandasPdbStack()
    stack.add_pdbs([TESTDATA_FILENAME, TESTDATA_FILENAME2,
                   TESTDATA_FILENAME_AF2_V4])
    stack.write_entries(".")

    for key, pdb in stack.pdbs.items():
        outfile = f'{key}.pdb'
        infile = f'data/{key}.pdb'
        assert os.path.exists(f'{key}.pdb')
        print(infile, outfile)
        with open(infile, "r") as f:
            f1 = f.read().replace('\r\n', '\n').rstrip('\n')
        with open(outfile, "r") as f:
            f2 = f.read().replace('\r\n', '\n').rstrip('\n')
        assert f1 == f2

        os.remove(outfile)