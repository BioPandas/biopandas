# Installing BioPandas ![](img/logos/1j1v_120.png)


## Requirements

BioPandas requires the following software and packages:

- [Python](https://www.python.org) 2.7, 3.5, or 3.6
- [NumPy](http://www.numpy.org) >= 1.11.2
- [SciPy](https://www.scipy.org/scipylib/index.html) >= 0.18.1
- [Pandas](http://pandas.pydata.org) >= 0.19.1


## PyPI

You can install the latest stable release of `biopandas` directly from Python's package index via `pip` by executing the following code from your command line:  

```bash
pip install biopandas  
```


## Conda-forge

Versions of `biopandas` are now also available via [conda-forge](https://github.com/conda-forge/biopandas-feedstock); you can install it via


```bash
conda install biopandas -c conda-forge
```

 or simply

```bash
conda install biopandas
```

if you have `conda-forge` already [added to your channels](https://github.com/conda-forge/biopandas-feedstock).



## Latest GitHub Source Code

<br>

You want to try out the latest features before they go live on PyPI? Install the `biopandas` dev-version latest development version from the GitHub repository by executing

```bash
pip install git+git://github.com/rasbt/biopandas.git
```

<br>


Alternatively, you download the package manually from [PYPI](https://pypi.python.org/pypi/biopandas) or [GitHub](https://github.com/rasbt/biopandas), unzip it, navigate into the package, and execute the command:

```bash
python setup.py install
```


