biopandas version: 0.6.0dev
## PandasPdbStack

### PandasPdbStack

*PandasPdbStack()*

None

### Methods

<hr>

### add_pdb

*add_pdb(source: 'Union[str, Dict[str, List[str]]]', key=None, convert_mmcif=True)*

Adds a single PDB to the stack with automatic or explicit keying.
    :param source: a string which defines a filename, a PDB or UniProt ID or a dictionary with lists of strings.
    :param method: a method to process the source. Auto will automatically determine the source type and processes accordingly.

    :return: None

<hr>

### add_pdbs

*add_pdbs(sources: 'List[Union[str, List[str], Dict[str, List[str]]]]')*

Adds multiple PDBs to the BioPandas collection from a list of sources.
    :param sources: a list of multiple PDB files or PDB/UniProt ID-s, that needs to be processed the input types can be mixed.
    :param method: a method to process the sources.

    :return: None

<hr>

### apply_calculation

*apply_calculation(calculation_func: 'Callable') -> 'dict'*

Applies a calculation across all PDBs in the stack and returns a dictionary with the calculated values.
    :param calculation_func: a function that processes a pandaspdb objects and returns objects (values, arrays or dictionaries, etc.)

    :return: a dictionary with the calculated values.

<hr>

### apply_filter

*apply_filter(filter_func: 'Callable', keep_null=True, **kwargs) -> 'PandasPdbStack'*

Applies a filter across all PDBs in the stack and returns a new stack with the filtered PDBs.
    Mandatory inputs for the filter function are the key and the PandasPdb object.
    :param filter_func: a function that processes a pandaspdb objects and returns a modified pandaspdb object.
    :param keep_null: a boolean to modify whether to keep empty structures.

    :return: a new stack with the filtered PDBs.

<hr>

### delete_entry

*delete_entry(key: 'str') -> 'None'*

Deletes a PDB entry from the stack.
    :param key: the key of the PDB entry to delete.

    :return: None

<hr>

### fetch_pdb

*fetch_pdb(key: 'str' = None, pdb_id: 'str' = None, uniprot_id: 'str' = None)*

Fetches a PDB file from the RCSB PDB  repository or AF2 database and assigns it to the keyof the same ID.
    :param pdb_id: the PDB ID to fetch.
    :param uniprot_id: the UniProt ID to fetch from the AF2 database.

    :return: None

<hr>

### read_mmcif

*read_mmcif(file_path: 'str', key: 'str' = None, convert=True)*

Reads mmCIF file from disk or URL.
    :param file_path: the path to the mmCIF file. Reads the file and assigns it to the key of the filename.
    :param key: the key to associate the PDB with.

    :return: None

<hr>

### read_pdb

*read_pdb(file_path: 'str', key: 'str' = None)*

Reads PDB file from disk or URL.
    :param file_path: the path to the PDB file. Reads the file and assigns it to the key of the filename.
    :param key: the key to associate the PDB with.

    :return: None

<hr>

### read_pdb_from_list

*read_pdb_from_list(pdb_lines: 'List[str]', key: 'str')*

Reads PDB data from a list of lines and assigns it to the specified key.
    :param pdb_lines: a list of PDB lines.
    :param key: a key to associate the PDB with.

    :return: None

<hr>

### update_entry

*update_entry(key: 'str', new_pdb: 'Union[str, Dict[str, List[str]]]') -> 'None'*

Updates a PDB entry in the stack.
    :param key: the key of the PDB entry to update.
    :param new_pdb: the new PandasPdb object to associate with the key.

    :return: None

<hr>

### write_entries

*write_entries(outdir: 'str', outfmt='pdb', records=['ATOM']) -> 'None'*

Writes all PDB entries in the stack to a directory.
    :param outdir: the directory to write the PDB files to.

    :return: None

### Properties

