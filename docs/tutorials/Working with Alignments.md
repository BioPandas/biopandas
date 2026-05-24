BioPandas

Author: Sebastian Raschka <mail@sebastianraschka.com>  
License: BSD 3 clause  
Project Website: http://rasbt.github.io/biopandas/  
Code Repository: https://github.com/rasbt/biopandas 


```python
%load_ext watermark
%watermark -d -u -p pandas,biopandas
```

    Last updated: 2024-05-06
    
    pandas   : 2.2.1
    biopandas: 0.5.1.dev0
    


```python
# auto-reload biopandas
%load_ext autoreload
%autoreload 2
```


```python
import pandas as pd
pd.set_option('display.width', 600)
pd.set_option('display.max_columns', 8)
```

# Align class

The Align class is an extandable class with a collection of functions frequently used for structural alignment.


### Apply rotation and translation to a set of coordinates
A set of coordinates can be transformed by a rotation matrix and a translation vector. The `transform` function takes a set of coordinates, a rotation matrix, and a translation vector as input and returns the transformed coordinates.


```python
from biopandas.align import Align
import numpy as np

align = Align()
coords = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
translation = np.array([3, 2, 1])

transformed_coords = align.transform(coords, matrix, translation)
target_coords = np.array([[4, 4, 4], [7, 7, 7], [10, 10, 10]])

transformed_coords == target_coords
```




    array([[ True,  True,  True],
           [ True,  True,  True],
           [ True,  True,  True]])



## TMAlign subclass

TMalign (or its newest version, USalign, which is installed with this package) can be used to superimpose two structures or a set of structures onto a target structures.

### Align two or more structures loaded into biopandas onto a specified target structure

The `tmalign_to` structure aligns a structure or structures to a target structure.

In the simplest case, two structure objects are aligned. The chains of the target and mobile structure(s) needs to be specified. The function returns the transformed structure and the TM-score.


```python
from biopandas.align import TMAlign
from biopandas.pdb import PandasPdb
tmalign = TMAlign()

# Align two structures
tmalign = TMAlign()
ppdb = PandasPdb()
ppdb.read_pdb('data/1ycr.pdb')

ppdb_mobile = PandasPdb()
ppdb_mobile.read_pdb('data/2d7t.pdb')

transformed_structure, tm_score = tmalign.tmalign_to(ppdb, ppdb_mobile, 'A', 'H')

print('TM-scores:', tm_score)
```


    ---------------------------------------------------------------------------

    AttributeError                            Traceback (most recent call last)

    Cell In[5], line 13
         10 ppdb_mobile = PandasPdb()
         11 ppdb_mobile.read_pdb('data/2d7t.pdb')
    ---> 13 transformed_structure, tm_score = tmalign.tmalign_to(ppdb, ppdb_mobile, 'A', 'H')
         15 print('TM-scores:', tm_score)
    

    File D:\work\huji\biopandas\biopandas\align\tmalign.py:120, in TMAlign.tmalign_to(self, target, mobiles, target_chain, mobile_chains)
        118 if isinstance(mobiles, PandasPdb):
        119     mobile_atoms = self.filter_and_validate_chain(mobiles, mobile_chains)
    --> 120     transformed_mobile, tm_score = self.process_structure_for_tmalign(target_file, mobile_atoms, mobile_chains)
        121     return transformed_mobile, tm_score
        122 elif  isinstance(mobiles, PandasPdbStack):
    

    File D:\work\huji\biopandas\biopandas\align\tmalign.py:98, in TMAlign.process_structure_for_tmalign(self, target_file, mobile_pdb, mobile_chain)
         95     transformed_mobile = self.transform_coords(transformed_mobile, type='HETATM', matrix=matrix,
         96                                                translation=translation)
         97 # clean up
    ---> 98 os.remove(matrix_file_path.name) if os.path.exists(matrix_file_path.name) else None
         99 os.remove(mobile_file.name) if os.path.exists(mobile_file.name) else None
        101 return transformed_mobile, tm_score
    

    AttributeError: 'str' object has no attribute 'name'


 In case of several mobile structures, either a single chain or a dictionary with chain identifiers can be provided. The returned values are a stack of transformed structures and a dictionary with TM-scores. 


```python
from biopandas.align import TMAlign
from biopandas.pdb import PandasPdb
from biopandas.stack import PandasPdbStack
tmalign = TMAlign()

# Align a stack of structures to a target structure
tmalign = TMAlign()
ppdb = PandasPdb()
ppdb.read_pdb('data/1ycr.pdb')

ppdb_stack = PandasPdbStack()
ppdb_stack.add_pdbs(['data/3eiy.cif.gz', 'data/2d7t.pdb', '4LWV'])

transformed_structures, tm_scores = tmalign.tmalign_to(ppdb, ppdb_stack, 'A', {'2d7t': 'H', '3eiy': 'A', '4LWV': 'A'})

print('TM-scores:', tm_scores)
```

### Align two structures using the `run_tmalign` function


The `run_tmalign` function takes two PDB file paths as input and returns the path to the matrix file and the TM-score. If no matrix file is written, the path is an automatically generated temporary file. If the directory of the file path does not exist, it is created.
The matrix file can be parsed and the rotation-translation matrix can be extracted with the `parse_matrix_file` function.

For convenience, this cell exposes the functions of the TMAlign class for custom usage. Running `tmalign_to` is the recommended way to align structures.


```python
from biopandas.align import TMAlign
tmalign = TMAlign()
matrix_file_path, tm_score = tmalign.run_tmalign('data/1ycr.pdb', 'data/2d7t.pdb', 'data/tmalign_output/tmalign_matrix.txt')
matrix, translation = tmalign.parse_tmalign_rotation_matrix(matrix_file_path)

print('Matrix:\n', matrix)
print('Translation:', translation)
print('TM-score:', tm_score)

# To apply this transformation to the coordinates of the second structure, the `transform` function can be used
coords = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
transformed_coords = align.transform(coords, matrix, translation)

# Check new coordinates
print('Transformed coordinates:\n', transformed_coords)
```

## Alignment of all PDB structures in a stack

A wrapper around `tmalign_to` is provided to align all PDB structures in a PandasPdbStack object to a member of the stack. For performing this alignment, all structures must have only one chain. The chains can be filtered with the stack's `apply_filter` function if necessary.


```python
from biopandas.stack import PandasPdbStack
stack = PandasPdbStack()
stack.add_pdbs(['1ycr', '2d7t', '3eiy', '4LWV'])

chains = {'1ycr': 'A', '2d7t': 'H', '3eiy': 'A', '4LWV': 'A'}

# Align all structures to the first one and show TM-scores
tmalign = TMAlign()
target_pdb_id, transformed_structures, tm_scores = tmalign.tmalign_in_stack(stack, chains)
print(f'TM_scores to {target_pdb_id}: {tm_scores}')

# Write out the transformed structures
transformed_structures.write_entries('data/aligned_structures/')
```


```python
# Align all structures to a specified one and show TM-scores
target_pdb_id, transformed_structures, tm_scores = tmalign.tmalign_in_stack(stack, chains, target='2d7t')
print(f'TM_scores to {target_pdb_id}: {tm_scores}')
```
