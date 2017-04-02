# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import gzip


def split_multimol2(mol2_path):
    r"""
    Splits a multi-mol2 file into individual Mol2 file contents.

    Parameters
    -----------
    mol2_path : str
      Path to the multi-mol2 file. Parses gzip files if the filepath
      ends on .gz.

    Returns
    -----------
    A generator object for lists for every extracted mol2-file. Lists contain
        the molecule ID and the mol2 file contents.
        e.g., ['ID1234', ['@<TRIPOS>MOLECULE\n', '...']]. Note that bytestrings
        are returned (for reasons of efficieny) if the Mol2 content is read
        from a gzip (.gz) file.

    """
    if mol2_path.endswith('.gz'):
        open_file = gzip.open
        read_mode = 'rb'
    else:
        open_file = open
        read_mode = 'r'
    check = {'rb': b'@<TRIPOS>MOLECULE', 'r': '@<TRIPOS>MOLECULE'}

    with open_file(mol2_path, read_mode) as f:
        mol2 = ['', []]
        while True:
            try:
                line = next(f)
                if line.startswith(check[read_mode]):
                    if mol2[0]:
                        yield(mol2)
                    mol2 = ['', []]
                    mol2_id = next(f)
                    mol2[0] = mol2_id.rstrip()
                    mol2[1].append(line)
                    mol2[1].append(mol2_id)
                else:
                    mol2[1].append(line)
            except StopIteration:
                yield(mol2)
                return
