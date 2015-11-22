![Logo](./docs/sources/img/logos/logo.png)

**Working with molecular structures in pandas DataFrames**


[![Continuous Integration](https://travis-ci.org/rasbt/biopandas.svg?branch=master)](https://travis-ci.org/rasbt/biopandas)
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

<br><br>

## Examples

![3eiy](./docs/sources/img/index/3eiy.png)

```python
>>> ppdb = PandasPDB()
>>> ppdb.fetch_pdb('3eiy')
>>> ppdb.df['ATOM'].head()
```

![3eiy head](./docs/sources/img/index/3eiy_head.png)

<br><br>
<br><br>


![3eiy head](./docs/sources/img/index/ligand_rmsd.png)

```python
>>> pl1 = PandasPDB()
>>> pl1.read_pdb('./docking_pose_1.pdb')
>>> pl2 = PandasPDB()
>>> pl2.read_pdb('./docking_pose_2.pdb')

>>> r = PandasPDB.rmsd(pl1.df['HETATM'], pl2.df['HETATM'], s='no hydrogen')
>>> print('RMSD: %f' % r)
RMSD: 2.6444
```

<br><br>
<br><br>


## Quick Install

- latest version (from GitHub): `pip install git+git://github.com/rasbt/biopandas.git#egg=biopandas`
- latest PyPI version: `pip install biopandas`

For more information, please see [http://rasbt.github.io/biopandas/installation/](http://rasbt.github.io/biopandas/installation/).
