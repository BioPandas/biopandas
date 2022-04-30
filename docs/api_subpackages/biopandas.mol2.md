biopandas version: 0.3.0
## PandasMol2

*PandasMol2()*

Object for working with Tripos Mol2 structure files.

**Attributes**

- `df` : pandas.DataFrame

    DataFrame of a Mol2's ATOM section


- `mol2_text` : str

    Mol2 file contents in string format


- `code` : str

    ID, code, or name of the molecule stored


- `pdb_path` : str

    Location of the MOL2 file that was read in via `read_mol2`

### Methods

<hr>

*distance(xyz=(0.0, 0.0, 0.0))*

Computes Euclidean distance between atoms in
    self.df and a 3D point.

**Parameters**

- `xyz` : tuple (0.00, 0.00, 0.00)

    X, Y, and Z coordinate of the reference center for the distance
    computation

**Returns**

- `pandas.Series` : Pandas Series object containing the Euclidean

    distance between the atoms in the atom section and `xyz`.

<hr>

*distance_df(df, xyz=(0.0, 0.0, 0.0))*

Computes Euclidean distance between atoms and a 3D point.

**Parameters**

- `df` : DataFrame

    DataFrame containing entries similar to the PandasMol2.df
    format for the
    the distance computation to the `xyz` reference coordinates.

- `xyz` : tuple (0.00, 0.00, 0.00)

    X, Y, and Z coordinate of the reference center for the distance
    computation

**Returns**

- `pandas.Series` : Pandas Series object containing the Euclidean

    distance between the atoms in the atom section and `xyz`.

<hr>

*read_mol2(path, columns=None)*

Reads Mol2 files (unzipped or gzipped) from local drive

    Note that if your mol2 file contains more than one molecule,
    only the first molecule is loaded into the DataFrame

**Attributes**

- `path` : str

    Path to the Mol2 file in .mol2 format or gzipped format (.mol2.gz)


- `columns` : dict or None (default: None)

    If None, this methods expects a 9-column ATOM section that contains
    the following columns:

    {0:('atom_id', int), 1:('atom_name', str),
    2:('x', float), 3:('y', float), 4:('z', float),
    5:('atom_type', str), 6:('subst_id', int),
    7:('subst_name', str), 8:('charge', float)}

    If your Mol2 files are formatted differently, you can provide your
    own column_mapping dictionary in a format similar to the one above.
    However, note that not all assert_raise_message methods
    may be supported then.

**Returns**

self

<hr>

*read_mol2_from_list(mol2_lines, mol2_code, columns=None)*

Reads Mol2 file from a list into DataFrames

**Attributes**

- `mol2_lines` : list

    A list of lines containing the mol2 file contents. For example,
    ['@<TRIPOS>MOLECULE\n',
    'ZINC38611810\n',
    '   65    68     0     0     0\n',
    'SMALL\n',
    'NO_CHARGES\n',
    '\n',
    '@<TRIPOS>ATOM\n',
    '      1 C1  -1.1786  2.7011  -4.0323 C.3  1 <0>   -0.1537\n',
    '      2 C2  -1.2950  1.2442  -3.5798 C.3  1 <0>   -0.1156\n',
    ...]


- `mol2_code` : str or None

    Name or ID of the molecule.


- `columns` : dict or None (default: None)

    If None, this methods expects a 9-column ATOM section that contains
    the following columns:
    {0:('atom_id', int), 1:('atom_name', str),
    2:('x', float), 3:('y', float), 4:('z', float),
    5:('atom_type', str), 6:('subst_id', int),
    7:('subst_name', str), 8:('charge', float)}
    If your Mol2 files are formatted differently, you can provide your
    own column_mapping dictionary in a format similar to the one above.
    However, note that not all assert_raise_message methods may be
    supported then.

**Returns**

self

<hr>

*rmsd(df1, df2, heavy_only=True)*

Compute the Root Mean Square Deviation between molecules

**Parameters**

- `df1` : pandas.DataFrame

    DataFrame with HETATM, ATOM, and/or ANISOU entries

- `df2` : pandas.DataFrame

    Second DataFrame for RMSD computation against df1. Must have the
    same number of entries as df1

- `heavy_only` : bool (default: True)

    Which atoms to compare to compute the RMSD. If `True` (default),
    computes the RMSD between non-hydrogen atoms only.

**Returns**

- `rmsd` : float

    Root Mean Square Deviation between df1 and df2

### Properties

<hr>

*df*

Acccesses the pandas DataFrame

## split_multimol2

*split_multimol2(mol2_path)*

Generator function that
    splits a multi-mol2 file into individual Mol2 file contents.

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

