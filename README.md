![Logo](./docs/sources/img/logos/logo.png)

**Working with molecular structures in pandas DataFrames**


[![Continuous Integration](https://travis-ci.org/rasbt/biopandas.svg?branch=master)](https://travis-ci.org/rasbt/biopandas)
[![Build status](https://ci.appveyor.com/api/projects/status/jcp91fvbgmqws30p/branch/master?svg=true)](https://ci.appveyor.com/project/rasbt/biopandas/branch/master)
[![Code Coverage](https://coveralls.io/repos/rasbt/biopandas/badge.svg?branch=master&service=github)](https://coveralls.io/github/rasbt/biopandas?branch=master)
[![Code Health](https://landscape.io/github/rasbt/biopandas/master/landscape.svg?style=flat)](https://landscape.io/github/rasbt/biopandas/master)
[![PyPI Version](https://img.shields.io/pypi/v/biopandas.svg)](https://pypi.python.org/pypi/biopandas/)
[![License](https://img.shields.io/badge/license-new%20BSD-blue.svg)](https://github.com/rasbt/biopandas/blob/master/LICENSE)
![Python 2.7](https://img.shields.io/badge/python-2.7-blue.svg)
![Python 3.5](https://img.shields.io/badge/python-3.5-blue.svg)

<hr>

## Links
- Documentation: [http://rasbt.github.io/biopandas/](http://rasbt.github.io/biopandas/)
- Source code repository: [https://github.com/rasbt/biopandas](https://github.com/rasbt/biopandas)
- PyPI: [https://pypi.python.org/pypi/biopandas](https://pypi.python.org/pypi/biopandas)
- How to contribute: [http://rasbt.github.io/biopandas/contributing/](http://rasbt.github.io/biopandas/contributing/)
- Changelog: [./docs/sources/CHANGELOG.md](./docs/sources/CHANGELOG.md)

<br>

If you are a computational biologist, chances are that you cursed one too many times about protein structure files. Yes, I am talking about ye Goode Olde Protein Data Bank format, aka "PDB files." Nothing against PDB, it's a neatly structured format (if deployed correctly); yet, it is a bit cumbersome to work with PDB files in "modern" programming languages -- I am pretty sure we all agree on this.

As machine learning and "data science" person, I fell in love with [pandas](http://pandas.pydata.org) DataFrames for handling just about everything that can be loaded into memory.  
So, why don't we take pandas to the structural biology world? Working with molecular structures of biological macromolecules in pandas DataFrames is what BioPandas is all about!

<br>

## Examples

![3eiy](./docs/sources/img/index/3eiy.png)

```python
# Initialize a new PandasPDB object
# and fetch the PDB file from rcsb.org
>>> from biopandas.pdb import PandasPDB
>>> ppdb = PandasPDB().fetch_pdb('3eiy')
>>> ppdb.df['ATOM'].head()
```

![3eiy head](./docs/sources/img/index/3eiy_head.png)

<br><br>
<br><br>


![3eiy head](./docs/sources/img/index/ligand_rmsd.png)

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

<br><br>
<br><br>


## Quick Install

- latest version (from GitHub): `pip install git+git://github.com/rasbt/biopandas.git#egg=biopandas`
- latest PyPI version: `pip install biopandas`

For more information, please see [http://rasbt.github.io/biopandas/installation/](http://rasbt.github.io/biopandas/installation/).
