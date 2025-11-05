biopandas version: 0.6.0dev
## fetch_mmtf

### fetch_mmtf

*fetch_mmtf(pdb_code: 'str') -> 'pd.DataFrame'*

Returns a dataframe from a PDB code.

    :param pdb_code: 4-letter PDB accession code
    :type pdb_code: str
    :return: Dataframe of protein structure
    :rtype: pd.DataFrame

## PandasMmtf

### PandasMmtf

*PandasMmtf()*

None

### Methods

<hr>

### amino3to1

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

### distance

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

### distance_df

*distance_df(df, xyz=(0.0, 0.0, 0.0))*

Computes Euclidean distance between atoms and a 3D point.

**Parameters**

- `df` : DataFrame

    DataFrame containing entries in the `PandasMmtf.df['ATOM']`
    or `PandasMmtf.df['HETATM']` format for the
    the distance computation to the `xyz` reference coordinates.

- `xyz` : tuple, default: (0.00, 0.00, 0.00)

    X, Y, and Z coordinate of the reference center for the distance
    computation.

**Returns**

- `pandas.Series` : Pandas Series object containing the Euclidean

    distance between the atoms in the record section and `xyz`.

<hr>

### fetch_mmtf

*fetch_mmtf(pdb_code: 'str')*

None

<hr>

### get_model

*get_model(model_index: 'int') -> 'PandasMmtf'*

Returns a new PandasPDB object with the dataframes subset to the
    given model index.

**Parameters**

- `model_index` : int

    An integer representing the model index to subset to.

**Returns**

- `pandas_pdb.PandasPdb` : A new PandasPdb object containing the

    structure subsetted to the given model.

<hr>

### get_models

*get_models(model_indices: 'List[int]') -> 'PandasMmtf'*

Returns a new PandasMmtf object with the dataframes subset to the
    given model index.

**Parameters**

- `model_indices` : List[int]

    A list representing the model indexes to subset to.

**Returns**

- `pandas_pdb.PandasMmtf` : A new PandasPdb object

    containing the structure subsetted to the given model.

<hr>

### impute_element

*impute_element(records=('ATOM', 'HETATM'), inplace=False)*

Impute element_symbol from atom_name section.

**Parameters**

- `records` : iterable, default: ('ATOM', 'HETATM')

    Coordinate sections for which the element symbols should be
    imputed.


- `inplace` : bool, (default: False

    Performs the operation in-place if True and returns a copy of the
    MMTF DataFrame otherwise.

**Returns**

DataFrame

<hr>

### parse_sse

*parse_sse()*

Parse secondary structure elements

<hr>

### read_mmtf

*read_mmtf(filename: 'Union[str, os.PathLike]')*

None

<hr>

### rmsd

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

### to_mmtf

*to_mmtf(path, records=('ATOM', 'HETATM'))*

Write record DataFrames to an MMTF file.

**Parameters**

- `path` : str

    A valid output path for the mmtf file

    records: tuple(str):
    A tuple of records to write. Defaults to ("ATOM". "HETATM")

<hr>

### to_pdb

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

### df

Access dictionary of pandas DataFrames for MMTF record sections.

## parse_mmtf

### parse_mmtf

*parse_mmtf(file_path: 'str') -> 'pd.DataFrame'*

Parses an MMTF file from disk and returns a pandas DataFrame.

    :param file_path: Path to mmtf file.
    :type file_path: str
    :return: Dataframe of protein structure.
    :rtype: pd.DataFrame

