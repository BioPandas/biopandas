## split_multimol2

*split_multimol2(mol2_path)*

Splits a multi-mol2 file into individual Mol2 file contents.

**Parameters**

- `mol2_path` : str

    Path to the multi-mol2 file. Parses gzip files if the filepath
    ends on .gz.

**Returns**

A generator object for lists for every extracted mol2-file. Lists contain
    the molecule ID and the mol2 file contents.
    e.g., ['ID1234', ['@<TRIPOS>MOLECULE\n', '...']]. Note that bytestrings
    are returned (for reasons of efficieny) if the Mol2 content is read
    from a gzip (.gz) file.

