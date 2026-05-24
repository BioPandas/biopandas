biopandas version: 0.6.0dev
## PandasPdb

### PandasPdb

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


- `pdb_path` : Union[str, os.PathLike]

    Location of the PDB file that was read in via `read_pdb`
    or URL of the page where the PDB content was fetched from
    if `fetch_pdb` was called.


- `header` : str

    PDB file description.


- `code` : str

    PDB code

### Methods

<hr>

### add_remark

*add_remark(code, text='', indent=0)*

Add custom REMARK entry.

    The remark will be inserted to preserve the ordering of REMARK codes, i.e. if the code is
    `n` it will be added after all remarks with codes less or equal to `n`. If the object does
    not store any remarks the remark will be inserted right before the first of ATOM, HETATM or
    ANISOU records.

**Parameters**

- `code` : int

    REMARK code according to PDB standards.


- `text` : str

    The text of the remark. If the text does not fit into a single line it will be wrapped
    into multiple lines of REMARK entries. Likewise, if the text contains new line
    characters it will be split accordingly.


- `indent` : int, default: 0

    Number of white spaces inserted before the text of the remark.

**Returns**

Nothing

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

    Specifies the record DataFrame.

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

### fetch_pdb

*fetch_pdb(pdb_code: 'Optional[str]' = None, uniprot_id: 'Optional[str]' = None, source: 'str' = 'pdb')*

Fetches PDB file contents from the Protein Databank at rcsb.org or AlphaFold database
    at https://alphafold.ebi.ac.uk/.
    .

**Parameters**

- `pdb_code` : str, optional

    A 4-letter PDB code, e.g., `"3eiy"` to retrieve structures from the PDB.
    Defaults to `None`.


- `uniprot_id` : str, optional

    A UniProt Identifier, e.g., `"Q5VSL9"` to retrieve structures from the AF2 database.
    Defaults to `None`.


- `source` : str

    The source to retrieve the structure from
    # (`"pdb"`, `"alphafold2-v3"`, `"alphafold2-v4"`(latest)). #deprecated
    (`"pdb"`, `"alphafold2-v6"`(latest)).
    Defaults to `"pdb"`.

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

*get_model(model_index: 'int') -> 'PandasPdb'*

Returns a new PandasPDB object with the dataframes subset to the given model index.

**Parameters**

- `model_index` : int

    An integer representing the model index to subset to.

**Returns**

- `pandas_pdb.PandasPdb` : A new PandasPdb object containing the

    structure subsetted to the given model.

<hr>

### get_model_start_end

*get_model_start_end() -> 'pd.DataFrame'*

Get the start and end of the models contained in the PDB file.

    Extracts model start and end line indexes based
    on lines labelled 'OTHERS' during parsing.

**Returns**

- `pandas.DataFrame` : Pandas DataFrame object containing

    the start and end line indexes of the models.

<hr>

### get_models

*get_models(model_indices: 'List[int]') -> 'PandasPdb'*

Returns a new PandasPDB object with the dataframes subset to the given model index.

**Parameters**

- `model_indices` : List[int]

    A list representing the model indexes to subset to.

**Returns**

- `pandas_pdb.PandasPdb` : A new PandasPdb object

    containing the structure subsetted to the given model.

<hr>

### gyradius

*gyradius(records: 'tuple[str]' = ('ATOM',), decimals: 'int' = 4) -> 'float'*

Compute the Radius of Gyration of a molecule

**Parameters**

- `records` : iterable, default: ("ATOM",)

    Records from PandasPdb object for which to calculate the radius of gyration.
    Any of `("ATOM", "HETATM")`.


- `decimals` : int, default: 4

    Specifies the number of decimal places to round the final value to.

**Returns**

- `rg` : float

    Radius of Gyration of df in Angstrom

<hr>

### impute_element

*impute_element(records=('ATOM', 'HETATM'), inplace=False)*

Impute element_symbol from atom_name section.

**Parameters**

- `records` : iterable, default: ('ATOM', 'HETATM')

    Coordinate sections for which the element symbols should be
    imputed.


- `inplace` : bool, default: False

    Performs the operation in-place if True and returns a copy of the
    PDB DataFrame otherwise.

**Returns**

DataFrame

<hr>

### label_models

*label_models()*

Adds a column (`"model_id"`) to the underlying
    DataFrames containing the model number.

<hr>

### parse_sse

*parse_sse()*

Parse secondary structure elements

<hr>

### read_pdb

*read_pdb(path)*

Read PDB files (unzipped or gzipped) from local drive

**Attributes**

- `path` : str

    Path to the PDB file in .pdb format or gzipped format (.pdb.gz).

**Returns**

self

<hr>

### read_pdb_from_list

*read_pdb_from_list(pdb_lines)*

Reads PDB file from a list into DataFrames

**Attributes**

- `pdb_lines` : list

    A list of lines containing the pdb file contents.

**Returns**

self

<hr>

### rmsd

*rmsd(df1, df2, s=None, invert=False, decimals=4)*

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


- `decimals` : int, default: 4

    Specifies the number of decimal places to round the final value to.

**Returns**

- `rmsd` : float

    Root Mean Square Deviation between df1 and df2

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

<hr>

### to_pdb_stream

*to_pdb_stream(records: 'tuple[str]' = ('ATOM', 'HETATM')) -> 'StringIO'*

Writes a PDB dataframe to a stream.

**Parameters**

- `records` : iterable, default: ('ATOM', 'HETATM')

    Iterable of record names to save to stream. Any of `["ATOM", "HETATM", "OTHERS"]`.

**Returns**

- `io.StringIO` : Filestream of PDB file.


### Properties

<hr>

### df

Access dictionary of pandas DataFrames for PDB record sections.

