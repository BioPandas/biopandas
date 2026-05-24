biopandas version: 0.6.0dev
## Align

### Align

*Align()*

None

### Methods

<hr>

### filter_and_validate_chain

*filter_and_validate_chain(pdb, chain_id)*

Filter the PandasPdb by chain_id and validate the presence of the chain.
    :param pdb: the PandasPdb object to filter.
    :param chain_id: the chain ID to filter by.

    :return: filtered_pdb

<hr>

### transform

*transform(coords, matrix, translation)*

Apply the rotation matrix and translation vector to the structure.
    :param coords: the coordinates to transform.
    :param matrix: the rotation matrix.
    :param translation: the translation vector.

    :return: transformed coordinates as a numpy array.

<hr>

### write_pdb_to_temp_file

*write_pdb_to_temp_file(pdb)*

Write a PandasPdb/PandasMmcif object's data to a temporary structure file and return the file handle.
    :param pdb: the Pandas object to write to the file.

    :return: file handle

### Properties

## TMAlign

### TMAlign

*TMAlign(tmalign_path: str = None)*

Class to align structures using TMalign and transform the mobile structure(s) while extracting TM-scores.
    TODO: extend to handle multiple chains in multiple structures.

### Methods

<hr>

### filter_and_validate_chain

*filter_and_validate_chain(pdb, chain_id)*

Filter the PandasPdb by chain_id and validate the presence of the chain.
    :param pdb: the PandasPdb object to filter.
    :param chain_id: the chain ID to filter by.

    :return: filtered_pdb

<hr>

### parse_tmalign_rotation_matrix

*parse_tmalign_rotation_matrix(file_path: str) -> (<built-in function array>, <built-in function array>)*

Parse the rotation matrix of TMalign and translation vector from the TMalign output file.
    :param file_path: the path to the TMalign output file.

    :return: matrix, translation

<hr>

### process_structure_for_tmalign

*process_structure_for_tmalign(target_file, mobile_pdb: biopandas.pdb.pandas_pdb.PandasPdb, mobile_chain: str) -> (<class 'biopandas.pdb.pandas_pdb.PandasPdb'>, <class 'float'>)*

Handle the TMalign execution and transformation for a given mobile structure and return the transformed mobile structure and TM-score
    :param target_file: the target structure's filepath
    :param mobile_pdb: the mobile structure

    :return: transformed_mobile, tm_score

<hr>

### run_tmalign

*run_tmalign(target: str, mobile: str, matrix_file_path: str = None) -> (<class 'str'>, <class 'float'>)*

Function to execute TMalign with a rotation matrix output for one target-mobile pair.
    :param target: the structure to align to, a filepath
    :param mobile: the structure to align, a filepath

    :return: matrix_file_path, tm_score

<hr>

### tmalign_in_stack

*tmalign_in_stack(stack: biopandas.stack.stack.PandasPdbStack, mobile_chains: dict, target: str = None) -> (<class 'biopandas.stack.stack.PandasPdbStack'>, <class 'dict'>)*

For doing TMalign inside a stack, with one of its entries
    :param stack: PandasPdbStack with the structures to align. All of them must have only one chain!
    :param target: the target structure to align to. If not provided, the first structure in the stack will be used.

    :return: matrix_file_path, tm_score

<hr>

### tmalign_to

*tmalign_to(target: biopandas.pdb.pandas_pdb.PandasPdb, mobiles: [<class 'biopandas.pdb.pandas_pdb.PandasPdb'>, <class 'biopandas.stack.stack.PandasPdbStack'>], target_chain: str, mobile_chains: [<class 'str'>, <class 'dict'>]) -> ([<class 'biopandas.pdb.pandas_pdb.PandasPdb'>, <class 'biopandas.stack.stack.PandasPdbStack'>], [<class 'float'>, <class 'dict'>])*

Run TMalign and transform the mobile structure(s) while extracting TM-scores, specifying chains to align.
    :param target: the target structure to align to, a PandasPdb object.
    :param mobiles: the structure(s) to align, either a PandasPdb object or a PandasPdbStack.
    :param target_chain: the chain of the target structure to align to.
    :param mobile_chains: the chain(s) to align. A dictionary for each structure in the stack or a single chain ID.

    :return: return the transformed structures and the corresponding TM-scores

<hr>

### transform

*transform(coords, matrix, translation)*

Apply the rotation matrix and translation vector to the structure.
    :param coords: the coordinates to transform.
    :param matrix: the rotation matrix.
    :param translation: the translation vector.

    :return: transformed coordinates as a numpy array.

<hr>

### transform_coords

*transform_coords(pdb, matrix, translation, type='ATOM')*

Apply the rotation matrix and translation vector to the structure.
    :param pdb: the PandasPdb object to transform.
    :param type: the record type to transform.
    :param matrix: the rotation matrix.
    :param translation: the translation vector.

    :return: transformed_pdb

<hr>

### write_pdb_to_temp_file

*write_pdb_to_temp_file(pdb)*

Write a PandasPdb/PandasMmcif object's data to a temporary structure file and return the file handle.
    :param pdb: the Pandas object to write to the file.

    :return: file handle

### Properties

