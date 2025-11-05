# BioPandas extension for
# Author: Sebastian Raschka
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

from biopandas.stack.stack import PandasPdbStack
from nose.tools import assert_raises
import os

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), "data", "3eiy.pdb")
TESTDATA_FILENAME2 = os.path.join(
    os.path.dirname(__file__), "data", "4eiy_anisouchunk.pdb"
)
TESTDATA_FILENAME_GZ = os.path.join(os.path.dirname(__file__), "data", "3eiy.pdb.gz")
TESTDATA_FILENAME_AF2_V6 = os.path.join(
    os.path.dirname(__file__), "data", "AF-Q5VSL9-F1-model_v6.pdb"
)


TESTDATA_FILENAME_MMCIF = os.path.join(os.path.dirname(__file__), "data", "3eiy.cif")
TESTDATA_FILENAME_MMCIF_GZ = os.path.join(os.path.dirname(__file__), "data", "3eiy.cif.gz")

with open(TESTDATA_FILENAME, "r") as f:
    three_eiy = f.read()

with open(TESTDATA_FILENAME2, "r") as f:
    four_eiy = f.read()

def test_add_pdb():
    stack = PandasPdbStack()
    stack.add_pdb(TESTDATA_FILENAME)
    assert len(stack.pdbs.keys()) == 1
    assert list(stack.pdbs.keys())[0] == '3eiy'

def test_add_pdb_gz():
    stack = PandasPdbStack()
    stack.add_pdb(TESTDATA_FILENAME_GZ)
    assert len(stack.pdbs.keys()) == 1
    assert list(stack.pdbs.keys())[0] == '3eiy'

def test_add_multiple_pdbs():
    stack = PandasPdbStack()
    stack.add_pdbs([TESTDATA_FILENAME, TESTDATA_FILENAME_AF2_V6])
    assert len(stack.pdbs.keys()) == 2
    assert '3eiy' in stack.pdbs
    assert 'AF-Q5VSL9-F1-model_v6' in stack.pdbs

def test_add_mmcif():
    stack = PandasPdbStack()
    stack.add_pdb(TESTDATA_FILENAME_MMCIF)
    assert len(stack.pdbs.keys()) == 1
    assert list(stack.pdbs.keys())[0] == '3eiy'

def test_add_mmcif_gz():
    stack = PandasPdbStack()
    stack.add_pdb(TESTDATA_FILENAME_MMCIF_GZ)
    assert len(stack.pdbs.keys()) == 1
    assert list(stack.pdbs.keys())[0] == '3eiy'

def test_fetch_multiple_pdbs_from_pdb():
    stack = PandasPdbStack()
    stack.add_pdbs(['3eiy', '4eiy'])
    assert len(stack.pdbs.keys()) == 2
    assert '3eiy' in stack.pdbs
    assert '4eiy' in stack.pdbs

def test_fetch_multiple_pdbs_from_af2():
    stack = PandasPdbStack()
    stack.add_pdbs(['Q5VSL9', 'P99999'])
    assert len(stack.pdbs.keys()) == 2
    assert 'Q5VSL9' in stack.pdbs
    assert 'P99999' in stack.pdbs

def test_fetch_multiple_pdbs_mixed():
    stack = PandasPdbStack()
    stack.add_pdbs(['3eiy', 'Q5VSL9'])
    assert len(stack.pdbs.keys()) == 2
    assert '3eiy' in stack.pdbs
    assert 'Q5VSL9' in stack.pdbs

def test_add_multiple_structures_mixed():
    stack = PandasPdbStack()
    stack.add_pdbs([TESTDATA_FILENAME, TESTDATA_FILENAME_GZ,
                    TESTDATA_FILENAME_AF2_V6,
                    TESTDATA_FILENAME_MMCIF, TESTDATA_FILENAME_MMCIF_GZ,
                    'Q5VSL9', 'P99999', '1YCR', '6V4E'])
    # There are multiple 3eiy, they overwrite each other, hence the length of dict is not 9
    assert len(stack.pdbs.keys()) == 6
    assert '3eiy' in stack.pdbs
    assert 'AF-Q5VSL9-F1-model_v6' in stack.pdbs
    assert 'Q5VSL9' in stack.pdbs
    assert 'P99999' in stack.pdbs

def test_read_pdb_from_list():
    """Test public read_pdb_from_list"""

    line_dict = dict(zip(["3eiy", "4eiy"], [three_eiy, four_eiy]))
    stack = PandasPdbStack()
    stack.add_pdbs(line_dict)
    assert len(stack.pdbs.keys()) == 2
    assert '3eiy' in stack.pdbs
    assert '4eiy' in stack.pdbs

def test_add_to_key():
    stack = PandasPdbStack()
    stack.add_pdb(TESTDATA_FILENAME, key='3eiy')
    assert len(stack.pdbs.keys()) == 1
    assert '3eiy' in stack.pdbs
    stack.add_pdb(TESTDATA_FILENAME, key='1ycr')
    assert len(stack.pdbs.keys()) == 2
    assert '3eiy' in stack.pdbs
    assert '1ycr' in stack.pdbs

def test_add_pdb_nonexistent():
    stack = PandasPdbStack()
    assert_raises(FileNotFoundError, stack.add_pdb, 'nonexistent.pdb')

def test_add_mmcifnonexistent():
    stack = PandasPdbStack()
    assert_raises(FileNotFoundError, stack.add_pdbs, ['nonexistent.cif'])

def test_add_random_id ():
    stack = PandasPdbStack()
    assert_raises(ValueError, stack.add_pdb, 'AGDHK')

def test_wrong_input():
    stack = PandasPdbStack()
    assert_raises(TypeError, stack.add_pdb, 1234)

def test_fetch_pdb_uniprot_id():
    stack = PandasPdbStack()
    assert_raises(ValueError, stack.fetch_pdb, uniprot_id='Q5VSL9', pdb_id='1YCR')

def test_fetch_no_id():
    stack = PandasPdbStack()
    assert_raises(ValueError, stack.fetch_pdb)