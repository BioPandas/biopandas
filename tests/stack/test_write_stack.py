# BioPandas extension for
# Author: Sebastian Raschka, Julia Varga
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

from biopandas.stack import PandasPdbStack
import os
import difflib


TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), "data", "3eiy.pdb")
TESTDATA_FILENAME2 = os.path.join(
    os.path.dirname(__file__), "data", "4eiy_anisouchunk.pdb"
)

TESTDATA_FILENAME_GZ = os.path.join(os.path.dirname(__file__), "data", "3eiy.pdb.gz")
TESTDATA_FILENAME_AF2_V6 = os.path.join(
    os.path.dirname(__file__), "data", "AF-Q5VSL9-F1-model_v6.pdb"
)

def test_write_pdb():
    stack = PandasPdbStack()
    stack.add_pdbs([TESTDATA_FILENAME, TESTDATA_FILENAME2,
                   TESTDATA_FILENAME_AF2_V6])
    base_dir = os.path.dirname(__file__)
    stack.write_entries(base_dir)

    for key, pdb in stack.pdbs.items():
        outfile = os.path.join(base_dir, f'{key}.pdb')
        infile = os.path.join(base_dir, 'data', f'{key}.pdb')
        assert os.path.exists(outfile)

        with open(infile, "r") as f:
            f1 = [line.rstrip() for line in f.read().replace('\r', '\n').split('\n')]
        with open(outfile, "r") as f:
            f2 = [line.rstrip() for line in f.read().replace('\r', '\n').split('\n')]
        assert f1 == f2

        os.remove(outfile)

def test_write_pdb_no_dir_exists():
    stack = PandasPdbStack()
    stack.add_pdbs([TESTDATA_FILENAME, TESTDATA_FILENAME2,
                   TESTDATA_FILENAME_AF2_V6])
    base_dir = os.path.dirname(__file__)
    test_dir = os.path.join(base_dir, "test_dir")
    stack.write_entries(test_dir)

    for key, pdb in stack.pdbs.items():
        outfile = os.path.join(test_dir, f'{key}.pdb')
        infile = os.path.join(base_dir, 'data', f'{key}.pdb')
        assert os.path.exists(outfile)

        with open(infile, "r") as f:
            f1 = [line.rstrip() for line in f.read().replace('\r', '\n').split('\n')]
        with open(outfile, "r") as f:
            f2 = [line.rstrip() for line in f.read().replace('\r', '\n').split('\n')]
        assert f1 == f2

        os.remove(outfile)

    os.rmdir(test_dir)