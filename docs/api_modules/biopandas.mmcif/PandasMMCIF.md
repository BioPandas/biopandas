## PandasMmcif

### PandasMmcif

*PandasMmcif(use_auth: 'bool' = True)*

None

### Methods

<hr>

### amino3to1

*amino3to1(record: 'str' = 'ATOM', residue_col: 'str' = 'auth_comp_id', residue_number_col: 'str' = 'auth_seq_id', chain_col: 'str' = 'auth_asym_id', fillna: 'str' = '?')*

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

### convert_to_pandas_pdb

*convert_to_pandas_pdb(offset_chains: 'bool' = True, records: 'List[str]' = ['ATOM', 'HETATM']) -> 'PandasPdb'*

Returns a PandasPdb object with the same data as the PandasMmcif
    object.

**Attributes**

offset_chains: bool
    Whether or not to offset atom numbering based on number of chains.
    This can arise due to the presence of TER records in PDBs which are
    not found in mmCIFs.
    records: List[str]
    List of record types to save. Any of ["ATOM", "HETATM", "OTHERS"].
    Defaults to ["ATOM", "HETATM"].

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

### fetch_mmcif

*fetch_mmcif(pdb_code: 'Optional[str]' = None, uniprot_id: 'Optional[str]' = None, source: 'str' = 'pdb')*

Fetches mmCIF file contents from the Protein Databank at rcsb.org or AlphaFold database at https://alphafold.ebi.ac.uk/.
    .

**Parameters**

- `pdb_code` : str, optional

    A 4-letter PDB code, e.g., `"3eiy"` to retrieve structures from the PDB. Defaults to `None`.


- `uniprot_id` : str, optional

    A UniProt Identifier, e.g., `"Q5VSL9"` to retrieve structures from the AF2 database. Defaults to `None`.


- `source` : str

    The source to retrieve the structure from
    (`"pdb"`, `"alphafold2-v3"` or `"alphafold2-v4"`). Defaults to `"pdb"`.

**Returns**

self

<hr>

### get

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

### get_model

*get_model(model_index: 'int') -> 'PandasMmcif'*

Returns a new PandasMmcif object with the dataframes subset to the
    given model index.

**Parameters**

- `model_index` : int

    An integer representing the model index to subset to.

**Returns**

- `pandas_pdb.PandasPdb` : A new PandasMMcif object containing the

    structure subsetted to the given model.

<hr>

### get_models

*get_models(model_indices: 'List[int]') -> 'PandasMmcif'*

Returns a new PandasMmcif object with the dataframes subset to the
    given model index.

**Parameters**

- `model_indices` : List[int]

    A list representing the model indexes to subset to.

**Returns**

- `pandas_pdb.PandasMmtf` : A new PandasMmcif object

    containing the structure subsetted to the given model.

<hr>

### label_models

*label_models()*

Adds a column ("model_id") to the underlying
    DataFrames containing the model number.

<hr>

### read_mmcif

*read_mmcif(path)*

Read MMCIF files (unzipped or gzipped) from local drive

**Attributes**

- `path` : Union[str, os.PathLike]

    Path to the MMCIF file in .cif format or gzipped format (.cif.gz).

**Returns**

self

<hr>

### read_mmcif_from_list

*read_mmcif_from_list(mmcif_lines)*

Reads mmCIF file from a list into DataFrames

**Attributes**

- `pdb_lines` : list

    A list of lines containing the mmCIF file contents.

**Returns**

self

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

### to_mmcif

*to_mmcif(path, records=None, gz=False)*

Write record DataFrames to an mmCIF file or gzipped mmCIF file.

**Parameters**

- `path` : str

    A valid output path for the mmcif file

- `records` : iterable, default: None

    A list of record sections in {'ATOM', 'HETATM', 'ANISOU'}
    that are to be written. Writes all sections if `records=None`.

- `gz` : bool, default: False

    Writes a gzipped mmCIF file if True.

### Properties

<hr>

### df

Acccess dictionary of pandas DataFrames for PDB record sections.

