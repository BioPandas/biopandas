biopandas version: 0.2.4
## PandasPdb

*PandasPdb()*

Object for working with Protein Databank structure files.

**Attributes**

- `df` : dict

    Dictionary storing pandas DataFrames for PDB record sections.
    The dictionary keys are {'ATOM', 'HETATM', 'ANISOU', 'OTHERS'}
    where 'OTHERS' contains all entries that are not parsed as
    'ATOM', 'HETATM', or 'ANISOU'.


- `pdb_text` : str

    PDB file contents in raw text format.


- `pdb_path` : str

    Location of the PDB file that was read in via `read_pdb`
    or URL of the page where the PDB content was fetched from
    if `fetch_pdb` was called.


- `header` : str

    PDB file description.


- `code` : str

    PDB code

### Methods

<hr>

*amino3to1(record='ATOM', residue_col='residue_name', fillna='?')*

Creates 1-letter amino acid codes from DataFrame

    Non-canonical amino-acids are converted as follows:
    ASH (protonated ASP) => D
    CYX (disulfide-bonded CYS) => C
    GLH (protonated GLU) => E
    HID/HIE/HIP (different protonation states of HIS) = H
    HYP (hydroxyproline) => P
    MSE (selenomethionine) => M

**Parameters**

- `record` : str, default: 'ATOM'

    Specfies the record DataFrame.

- `residue_col` : str,  default: 'residue_name'

    Column in `record` DataFrame to look for 3-letter amino acid
    codes for the conversion.

- `fillna` : str, default: '?'

    Placeholder string to use for unknown amino acids.

**Returns**

- `pandas.DataFrame` : Pandas DataFrame object consisting of two columns,

    `'chain_id'` and `'residue_name'`, where the former contains
    the chain ID of the amino acid and the latter
    contains the 1-letter amino acid code, respectively.

<hr>

*distance(xyz=(0.0, 0.0, 0.0), records=('ATOM', 'HETATM'))*

Computes Euclidean distance between atoms and a 3D point.

**Parameters**

- `xyz` : tuple, default: (0.00, 0.00, 0.00)

    X, Y, and Z coordinate of the reference center for the distance
    computation.

- `records` : iterable, default: ('ATOM', 'HETATM')

    Specify which record sections to consider. For example, to consider
    both protein and ligand atoms, set `records=('ATOM', 'HETATM')`.
    This setting is ignored if `df` is not set to None.
    For downward compatibility, a string argument is still supported
    but deprecated and will be removed in future versions.

**Returns**

- `pandas.Series` : Pandas Series object containing the Euclidean

    distance between the atoms in the record section and `xyz`.

<hr>

*distance_df(df, xyz=(0.0, 0.0, 0.0))*

Computes Euclidean distance between atoms and a 3D point.

**Parameters**

- `df` : DataFrame

    DataFrame containing entries in the `PandasPdb.df['ATOM']`
    or `PandasPdb.df['HETATM']` format for the
    the distance computation to the `xyz` reference coordinates.

- `xyz` : tuple, default: (0.00, 0.00, 0.00)

    X, Y, and Z coordinate of the reference center for the distance
    computation.

**Returns**

- `pandas.Series` : Pandas Series object containing the Euclidean

    distance between the atoms in the record section and `xyz`.

<hr>

*fetch_pdb(pdb_code)*

Fetches PDB file contents from the Protein Databank at rcsb.org.

**Parameters**

- `pdb_code` : str

    A 4-letter PDB code, e.g., "3eiy".

**Returns**

self

<hr>

*get(s, df=None, invert=False, records=('ATOM', 'HETATM'))*

Filter PDB DataFrames by properties

**Parameters**

- `s` : str  in {'main chain', 'hydrogen', 'c-alpha', 'heavy'}

    String to specify which entries to return.


- `df` : pandas.DataFrame, default: None

    Optional DataFrame to perform the filter operation on.
    If df=None, filters on self.df['ATOM'].


- `invert` : bool, default: True

    Inverts the search query. For example if s='hydrogen' and
    invert=True, all but hydrogen entries are returned.


- `records` : iterable, default: ('ATOM', 'HETATM')

    Specify which record sections to consider. For example, to consider
    both protein and ligand atoms, set `records=('ATOM', 'HETATM')`.
    This setting is ignored if `df` is not set to None.
    For downward compatibility, a string argument is still supported
    but deprecated and will be removed in future versions.

**Returns**

- `df` : pandas.DataFrame

    Returns a DataFrame view on the filtered entries.

<hr>

*impute_element(records=('ATOM', 'HETATM'), inplace=False)*

Impute element_symbol from atom_name section.

**Parameters**

- `records` : iterable, default: ('ATOM', 'HETATM')

    Coordinate sections for which the element symbols should be
    imputed.


- `inplace` : bool, (default: False

    Performs the operation in-place if True and returns a copy of the
    PDB DataFrame otherwise.

**Returns**

DataFrame

<hr>

*parse_sse()*

Parse secondary structure elements

<hr>

*read_pdb(path)*

Read PDB files (unzipped or gzipped) from local drive

**Attributes**

- `path` : str

    Path to the PDB file in .pdb format or gzipped format (.pdb.gz).

**Returns**

self

<hr>

*rmsd(df1, df2, s=None, invert=False)*

Compute the Root Mean Square Deviation between molecules.

**Parameters**

- `df1` : pandas.DataFrame

    DataFrame with HETATM, ATOM, and/or ANISOU entries.


- `df2` : pandas.DataFrame

    Second DataFrame for RMSD computation against df1. Must have the
    same number of entries as df1.


- `s` : {'main chain', 'hydrogen', 'c-alpha', 'heavy', 'carbon'} or None,

    default: None
    String to specify which entries to consider. If None, considers
    all atoms for comparison.


- `invert` : bool, default: False

    Inverts the string query if true. For example, the setting
    `s='hydrogen', invert=True` computes the RMSD based on all
    but hydrogen atoms.

**Returns**

- `rmsd` : float

    Root Mean Square Deviation between df1 and df2

<hr>

*to_pdb(path, records=None, gz=False, append_newline=True)*

Write record DataFrames to a PDB file or gzipped PDB file.

**Parameters**

- `path` : str

    A valid output path for the pdb file


- `records` : iterable, default: None

    A list of PDB record sections in
    {'ATOM', 'HETATM', 'ANISOU', 'OTHERS'} that are to be written.
    Writes all lines to PDB if `records=None`.


- `gz` : bool, default: False

    Writes a gzipped PDB file if True.


- `append_newline` : bool, default: True

    Appends a new line at the end of the PDB file if True

### Properties

<hr>

*df*

Acccess dictionary of pandas DataFrames for PDB record sections.

