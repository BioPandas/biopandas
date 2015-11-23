
![Logo](./img/logos/logo_size_1.png)

**Working with molecular structures in pandas DataFrames**


[![Continuous Integration](https://travis-ci.org/rasbt/biopandas.svg?branch=master)](https://travis-ci.org/rasbt/biopandas)
[![PyPI Version](https://img.shields.io/pypi/v/biopandas.svg)](https://pypi.python.org/pypi/biopandas/)
[![License](https://img.shields.io/badge/license-new%20BSD-blue.svg)](https://github.com/rasbt/biopandas/blob/master/LICENSE)
![Python 2.7](https://img.shields.io/badge/python-2.7-blue.svg)
![Python 3.5](https://img.shields.io/badge/python-3.5-blue.svg)

<br>

<hr>

#### Links

- Documentation: [http://rasbt.github.io/biopandas/](http://rasbt.github.io/biopandas/)
- Source code repository: [https://github.com/rasbt/biopandas](https://github.com/rasbt/biopandas)
- PyPI: [https://pypi.python.org/pypi/biopandas](https://pypi.python.org/pypi/biopandas)

<hr>

#### About  


If you are a computational biologist, chances are that you cursed one too many times about protein structure files. Yes, I am talking about ye Goode Olde Protein Data Bank format, aka "PDB files." Nothing against PDB, it's a neatly structured format (if deployed correctly); yet, it is a bit cumbersome to work with PDB files in "modern" programming languages -- I am pretty sure we all agree on this.

As machine learning and "data science" person, I fell in love with [pandas](http://pandas.pydata.org) DataFrames for handling just about everything that can be loaded into memory.  
So, why don't we take pandas to the structural biology world? Working with molecular structures of biological macromolecules in pandas DataFrames is what BioPandas is all about!



<hr>

#### Examples

![3eiy](./img/index/3eiy.png)

```python
# Initialize a new PandasPDB object
# and fetch the PDB file from rcsb.org
>>> from biopandas.pdb import PandasPDB
>>> ppdb = PandasPDB().fetch_pdb('3eiy')
>>> ppdb.df['ATOM'].head()
```

![3eiy head](./img/index/3eiy_head.png)

<br><br>


![3eiy head](./img/index/ligand_rmsd.png)

```python
# Load structures from your drive and compute the
# Root Mean Square Deviation
>>> from biopandas.pdb import PandasPDB
>>> pl1 = PandasPDB().read_pdb('./docking_pose_1.pdb')
>>> pl2 = PandasPDB().read_pdb('./docking_pose_2.pdb')
>>> r = PandasPDB.rmsd(pl1.df['HETATM'], pl2.df['HETATM'],
                       s='hydrogen', invert=True)
>>> print('RMSD: %.4f Angstrom' % r)

RMSD: 2.6444 Angstrom
```
