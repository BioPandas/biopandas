![Logo](./docs/img/logos/logo.png#gh-light-mode-only)![Logo](./docs/img/logos/logo_dark.png#gh-dark-mode-only)

**Working with molecular structures in pandas DataFrames**


[![Continuous Integration](https://travis-ci.org/rasbt/biopandas.svg?branch=master)](https://travis-ci.org/rasbt/biopandas)
[![Build status](https://ci.appveyor.com/api/projects/status/jcp91fvbgmqws30p/branch/master?svg=true)](https://ci.appveyor.com/project/rasbt/biopandas/branch/master)
[![Code Coverage](https://coveralls.io/repos/rasbt/biopandas/badge.svg?branch=master&service=github)](https://coveralls.io/github/rasbt/biopandas?branch=master)
[![PyPI Version](https://img.shields.io/pypi/v/biopandas.svg)](https://pypi.python.org/pypi/biopandas/)
[![License](https://img.shields.io/badge/license-new%20BSD-blue.svg)](https://github.com/rasbt/biopandas/blob/master/LICENSE)
![Python 3](https://img.shields.io/badge/python-3-blue.svg)
[![JOSS](http://joss.theoj.org/papers/10.21105/joss.00279/status.svg)](http://joss.theoj.org/papers/10.21105/joss.00279)
[![Discuss](https://img.shields.io/badge/discuss-github-blue.svg)](https://github.com/rasbt/biopandas/discussions)

<hr>

## Links
- Documentation: [https://BioPandas.github.io/biopandas/](https://BioPandas.github.io/biopandas/)
- Source code repository: [https://github.com/rasbt/biopandas](https://github.com/rasbt/biopandas)
- PyPI: [https://pypi.python.org/pypi/biopandas](https://pypi.python.org/pypi/biopandas)
- How to contribute: [https://biopandas.github.io/biopandas/CONTRIBUTING/](https://biopandas.github.io/biopandas/CONTRIBUTING/)
- Changelog: [./docs/sources/CHANGELOG.md](./docs/sources/CHANGELOG.md)

<br>

If you are a computational biologist, chances are that you cursed one too many times about protein structure files. Yes, I am talking about ye Goode Olde Protein Data Bank format, aka "PDB files." Nothing against PDB, it's a neatly structured format (if deployed correctly); yet, it is a bit cumbersome to work with PDB files in "modern" programming languages -- I am pretty sure we all agree on this.

As machine learning and "data science" person, I fell in love with [pandas](http://pandas.pydata.org) DataFrames for handling just about everything that can be loaded into memory.  
So, why don't we take pandas to the structural biology world? Working with molecular structures of biological macromolecules (from PDB and MOL2 files) in pandas DataFrames is what BioPandas is all about!

<br>

## Examples

![3eiy](./docs/img/index/3eiy.png#gh-light-mode-only)![3eiy](./docs/img/index/3eiy_dark.png#gh-dark-mode-only)

```python
# Initialize a new PandasPdb object
# and fetch the PDB file from rcsb.org
>>> from biopandas.pdb import PandasPdb
>>> ppdb = PandasPdb().fetch_pdb('3eiy')
>>> ppdb.df['ATOM'].head()
```

![3eiy head](./docs/img/index/3eiy_head.png#gh-light-mode-only)![3eiy head](./docs/img/index/3eiy_head_dark.png#gh-dark-mode-only)

<br><br>
<br><br>


![3eiy head](./docs/img/index/ligand_rmsd.png#gh-light-mode-only)![3eiy head](./docs/img/index/ligand_rmsd_dark.png#gh-dark-mode-only)

```python
# Load structures from your drive and compute the
# Root Mean Square Deviation
>>> from biopandas.pdb import PandasPdb
>>> pl1 = PandasPdb().read_pdb('./docking_pose_1.pdb')
>>> pl2 = PandasPdb().read_pdb('./docking_pose_2.pdb')
>>> r = PandasPdb.rmsd(pl1.df['HETATM'], pl2.df['HETATM'],
                       s='hydrogen', invert=True)
>>> print('RMSD: %.4f Angstrom' % r)

RMSD: 2.6444 Angstrom
```

<br><br>
<br><br>


## Quick Install

- install the latest version (from GitHub): `pip install git+git://github.com/rasbt/biopandas.git#egg=biopandas`
- install the latest PyPI version: `pip install biopandas`
- install biopandas via conda-forge: `conda install biopandas -c conda-forge`

#### Requirements

- [Python](https://www.python.org) >=3.7
- [NumPy](http://www.numpy.org) >= 1.11.2
- [SciPy](https://www.scipy.org/scipylib/index.html) >= 0.18.1
- [Pandas](http://pandas.pydata.org) >= 0.19.1


For more information, please see [https://BioPandas.github.io/biopandas/installation/](https://BioPandas.github.io/biopandas/installation/).

<br><br>
<br><br>


### Cite as

If you use BioPandas as part of your workflow in a scientific publication, please consider citing the BioPandas repository with the following DOI:

- Sebastian Raschka. Biopandas: Working with molecular structures in pandas dataframes. *The Journal of Open Source Software*, 2(14), jun 2017. doi: 10.21105/joss.00279. URL http://dx.doi.org/10.21105/joss.00279.

```
@article{raschkas2017biopandas,
  doi = {10.21105/joss.00279},
  url = {http://dx.doi.org/10.21105/joss.00279},
  year  = {2017},
  month = {jun},
  publisher = {The Open Journal},
  volume = {2},
  number = {14},
  author = {Sebastian Raschka},
  title = {BioPandas: Working with molecular structures in pandas DataFrames},
  journal = {The Journal of Open Source Software}
}
```
