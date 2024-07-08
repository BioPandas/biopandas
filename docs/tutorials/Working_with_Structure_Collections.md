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
import pandas as pd
pd.set_option('display.width', 600)
pd.set_option('display.max_columns', 8)
```

# Working with Structure Collections

## Loading files
The `PandasPdbStack` class is a wrapper around the `PandasPdb` class that allows for loading and handling of multiple structures at once.
The stack can currently handle PDB and mmCIF file formats. Structures from files, PDB and UniProt can be loaded with the same command, at the same time. 

### 1 -- Loading files from file or database


```python
from biopandas.stack import PandasPdbStack
ppdb_stack = PandasPdbStack()
ppdb_stack.add_pdbs(['1l2y', 'data/1ycr.pdb', 'P99999', 'data/3eiy.cif.gz'])

# Check the content
ppdb_stack.pdbs.keys()
```




    dict_keys(['1l2y', '1ycr', 'P99999', '3eiy'])



If needed, further structures can be added to the stack with the `add_pdbs` or `add_pdb` method.


```python
# By specifying a different key, a structure can also be added twice to the stack.
ppdb_stack.add_pdb('data/1ycr.pdb', '1YCR')

# Check the content
ppdb_stack.pdbs.keys()
```




    dict_keys(['1l2y', '1ycr', 'P99999', '3eiy', '1YCR'])



### 2 -- Loading files from list of lines

Loading stuctures from list of lines can be from a dictionary where the keys will be used to refer to the structure.
Alternatively, it can also be done using the `add_pdb` method done one-by with specifying a key.




```python
from biopandas.stack import PandasPdbStack
ppdb_stack = PandasPdbStack()

with open('./data/3eiy.pdb', 'r') as f:
    three_eiy = f.read()

with open('./data/1ycr.pdb', 'r') as f:
    one_ycr = f.read()
    
with open('./data/2d7t.pdb', 'r') as f:
    two_d7t = f.read()

lines_input = {'3EIY': three_eiy, '1YCR': one_ycr}    
ppdb_stack.add_pdbs(lines_input)

ppdb_stack.add_pdb(two_d7t, '2D7T')

# Check the content
ppdb_stack.pdbs.keys()
```

    D:\work\huji\biopandas\biopandas\pdb\pandas_pdb.py:545: UserWarning: No ATOM/HETATM entries have been loaded. Is the input file/text in the pdb format?
      warnings.warn(
    D:\work\huji\biopandas\biopandas\pdb\pandas_pdb.py:545: UserWarning: No ATOM/HETATM entries have been loaded. Is the input file/text in the pdb format?
      warnings.warn(
    D:\work\huji\biopandas\biopandas\pdb\pandas_pdb.py:545: UserWarning: No ATOM/HETATM entries have been loaded. Is the input file/text in the pdb format?
      warnings.warn(
    




    dict_keys(['3EIY', '1YCR', '2D7T'])



## Applying functions to all members of the collection

Two apply methods are available to apply a function to all members of the collection. The function can do manipulation on the members (`apply_filter`) or calculation on the members (`apply_calculation`). apply_filter always returns a new PandasPdbStack object, while apply_calculation returns a dictionary with the results for each structure. 

### 1 -- apply_filter

Mandatory arguments for functions to be applied with this are the key (name of the structure, as string) and the PandasPdb object. Optional arguments can be passed as a dictionary.



```python
from biopandas.stack import PandasPdbStack
ppdb_stack = PandasPdbStack()
ppdb_stack.add_pdbs(['1l2y', 'data/1ycr.pdb', 'P99999', 'data/3eiy.cif.gz'])

def filter_by_chain(key, pdb, chain):
    # Example function for applying filtering
    pdb.df['ATOM'] = pdb.df['ATOM'].query('chain_id==@chain')
    return pdb


ppdb_stack_filtered = ppdb_stack.apply_filter(filter_by_chain, chain='H')
ppdb_stack_filtered.pdbs.keys()
```




    dict_keys(['1l2y', '1ycr', 'P99999', '3eiy'])



If a filtering would return an empty structure, the keep_nulls argument can be set to True to keep the structure in the collection or False to remove it. By default, keep_nulls is set to True.


```python
ppdb_stack_filtered_no_null = ppdb_stack.apply_filter(filter_by_chain, chain='B', keep_null=False)
ppdb_stack_filtered_no_null.pdbs.keys()
```




    dict_keys(['1ycr'])



In case, different structures are filtered for different attributes, a dictionary defining the attributes can be passed to the apply_filter method.


```python
from biopandas.stack import PandasPdbStack
ppdb_stack = PandasPdbStack()
ppdb_stack.add_pdbs(['1l2y', 'data/1ycr.pdb', 'P99999', 'data/3eiy.cif.gz'])

def filter_by_chains(key, pdb, chains):
    # Example function for applying filtering
    chain = chains[key]
    pdb.df['ATOM'] = pdb.df['ATOM'].query('chain_id==@chain')
    return pdb

args = {'chains': {'1ycr': 'A', '1l2y': 'H', '3eiy': 'A', 'P99999': 'B'}}
filtered_stack = ppdb_stack.apply_filter(filter_by_chains, keep_null=False, **args)

# Check content
filtered_stack.pdbs.keys()
```




    dict_keys(['1ycr', '3eiy'])



### 2 -- apply_calculation

The `apply_calculation` method is used to apply a function that calculates an output(s) for each structure in the collection.


```python
from biopandas.stack import PandasPdbStack
ppdb_stack = PandasPdbStack()
ppdb_stack.add_pdbs(['1l2y', 'data/1ycr.pdb', 'P99999', 'data/3eiy.cif.gz'])

def calculate_chain_lengths(key, pdb):
    # Assuming `pdb` is a PandasPdb object
    lengths = {}
    for ch in pdb.df['ATOM']['chain_id'].unique():
      ch_len = len(pdb.df['ATOM'].query('chain_id==@ch and atom_name=="CA"'))
      lengths[ch] = ch_len
    return lengths

chain_lengths = ppdb_stack.apply_calculation(calculate_chain_lengths)
chain_lengths
```




    {'1l2y': {'A': 760},
     '1ycr': {'A': 85, 'B': 13},
     'P99999': {'A': 105},
     '3eiy': {'A': 174}}



## Saving files

Using stack, it is also possible to write out the structure files in PDB format into an output directory. The directory will be created if does not exist.



```python
from biopandas.stack import PandasPdbStack
import glob

ppdb_stack = PandasPdbStack()
ppdb_stack.add_pdbs(['1l2y', 'data/1ycr.pdb', 'P99999', 'data/3eiy.cif.gz'])

ppdb_stack.write_entries('data/stack_output')

# Check the presence of the output files
glob.glob('data/stack_output/*')
```




    ['data/stack_output\\1l2y.pdb',
     'data/stack_output\\1ycr.pdb',
     'data/stack_output\\3eiy.pdb',
     'data/stack_output\\P99999.pdb']


