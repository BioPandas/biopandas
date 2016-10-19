biopandas version: 0.1.5.dev0
## PandasPDB

*PandasPDB()*

Object for working with Protein Databank structure files.

**Attributes**

- `df` : dict

    Dictionary storing pandas DataFrames for PDB record sections.
    The dictionary keys are {'ATOM', 'HETATM', 'ANISOU', 'OTHERS'}
    where 'OTHERS' contains all entries that are not parsed as
    'ATOM', 'HETATM', or 'ANISOU'


- `pdb_text` : str

    PDB file contents in raw text format


- `header` : str

    PDB file description


- `code` : str

    PDB code

### Methods

<hr>

*fetch_pdb(pdb_code)*

Fetches PDB file contents from the Protein Databank at rcsb.org.

**Parameters**

- `pdb_code` : str

    A 4-letter PDB code, e.g., "3eiy"

**Returns**

self

<hr>

*get(s, df=None, invert=False)*

Filter PDB DataFrames by properties

**Parameters**

- `s` : str  in {'main chain', 'hydrogen', 'c-alpha', 'heavy'}

    String to specify which entries to return


- `df` : pandas.DataFrame, default: None

    Optional DataFrame to perform the filter operation on.
    If df=None, filters on self.df['ATOM']


- `invert` : bool, default: True

    Inverts the search query. For example if s='hydrogen' and
    invert=True, all but hydrogen entries are returned

**Returns**

- `df` : pandas.DataFrame

    Returns a DataFrame view on the filtered entries.

<hr>

*impute_element(sections=['ATOM', 'HETATM'], inplace=False)*

Impute element_symbol from atom_name section.

**Parameters**

- `sections` : iterable (default: ['ATOM', 'HETATM'])

    Coordinate sections for which the element symbols should be
    imputed.


- `inplace` : bool (default: False)

    Performs the operation in-place if True and returns a copy of the
    PDB DataFrame otherwise.

**Returns**

DataFrame

<hr>

*read_pdb(path)*

Read PDB files (unzipped or gzipped) from local drive

**Attributes**

- `path` : str

    Path to the PDB file in .pdb format or gzipped format (.pdb.gz)

**Returns**

self

<hr>

*rmsd(df1, df2, s=None, invert=False)*

Compute the Root Mean Square Deviation between molecules.

**Parameters**

- `df1` : pandas.DataFrame

    DataFrame with HETATM, ATOM, and/or ANISOU entries


- `df2` : pandas.DataFrame

    Second DataFrame for RMSD computation against df1. Must have the
    same number of entries as df1


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
    Writes all lines to PDB if records=None


- `gz` : bool, default: False

    Writes a gzipped PDB file if True


- `append_newline` : bool, default: True

    Appends a new line at the end of the PDB file if True

### Properties

<hr>

*df*

Acccess dictionary of pandas DataFrames for PDB record sections.

