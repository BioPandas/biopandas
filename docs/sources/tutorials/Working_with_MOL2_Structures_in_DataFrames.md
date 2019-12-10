# Working with MOL2 Structures in DataFrames

The Tripos MOL2 format is a common format for working with small molecules. In this tutorial, we will go over some examples that illustrate how we can use Biopandas' MOL2 DataFrames to analyze molecules conveniently.

## Loading MOL2 Files

Using the `read_mol2` method, we can read MOL2 files from standard .mol2 text files:


```python
from biopandas.mol2 import PandasMol2

pmol = PandasMol2().read_mol2('./data/1b5e_1.mol2')
```

[File link: [1b5e_1.mol2](https://raw.githubusercontent.com/rasbt/biopandas/master/docs/sources/tutorials/data/1b5e_1.mol2)]

The `read_mol2` method can also load structures from `.mol2.gz` files, but if you have a multi-mol2 file, keep in mind that it will only fetch the first molecule in this file. In the section "[Parsing Multi-MOL2 files](#parsing-multi-mol2-files)," we will see how we can parse files that contain multiple structures.


```python
pmol = PandasMol2().read_mol2('./data/40_mol2_files.mol2.gz')
```

[File link: [40_mol2_files.mol2.gz](https://github.com/rasbt/biopandas/blob/master/docs/sources/tutorials/data/40_mol2_files.mol2.gz?raw=true)]

After the file was succesfully loaded, we have access to the following basic `PandasMol2` attributes:


```python
print('Molecule ID: %s' % pmol.code)
print('\nRaw MOL2 file contents:\n\n%s\n...' % pmol.mol2_text[:500])
```

    Molecule ID: ZINC38611810
    
    Raw MOL2 file contents:
    
    @<TRIPOS>MOLECULE
    ZINC38611810
       65    68     0     0     0
    SMALL
    NO_CHARGES
    
    @<TRIPOS>ATOM
          1 C1         -1.1786    2.7011   -4.0323 C.3       1 <0>        -0.1537
          2 C2         -1.2950    1.2442   -3.5798 C.3       1 <0>        -0.1156
          3 C3         -0.1742    0.4209   -4.2178 C.3       1 <0>        -0.1141
          4 C4         -0.2887   -1.0141   -3.7721 C.2       1 <0>         0.4504
          5 O1         -1.1758   -1.3445   -3.0212 O.2       1 <0>        -0.4896
          6 O2       
    ...


The most interesting and useful attribute, however, is the [`PandasMol2.df`](../api_subpackages/biopandas.mol2#pandasmol2df) DataFrame, which contains the ATOM section of the MOL2 structure. Let's print the first 3 lines from the `ATOM` coordinate section to see how it looks like:


```python
pmol.df.head(3)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>atom_id</th>
      <th>atom_name</th>
      <th>x</th>
      <th>y</th>
      <th>...</th>
      <th>atom_type</th>
      <th>subst_id</th>
      <th>subst_name</th>
      <th>charge</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>C1</td>
      <td>-1.1786</td>
      <td>2.7011</td>
      <td>...</td>
      <td>C.3</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-0.1537</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>C2</td>
      <td>-1.2950</td>
      <td>1.2442</td>
      <td>...</td>
      <td>C.3</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-0.1156</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>C3</td>
      <td>-0.1742</td>
      <td>0.4209</td>
      <td>...</td>
      <td>C.3</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-0.1141</td>
    </tr>
  </tbody>
</table>
<p>3 rows × 9 columns</p>
</div>



## The MOL2 Data Format

`PandasMol2` expects the MOL2 file to be in the standard Tripos MOL2 format, and most importantly, that the "@<TRIPOS>ATOM" section is consistent with the following format convention:


> Format:
     **atom_id atom_name x y z atom_type [subst_id
        [subst_name [charge [status_bit]]]]**
       
> - atom_id (integer) = the ID number of the atom at the time the file was created. This is provided for reference only and is not used when the .mol2 file is read into SYBYL.
- atom_name (string) = the name of the atom.
- x (real) = the x coordinate of the atom.
- y (real) = the y coordinate of the atom.
- z (real) = the z coordinate of the atom.
- atom_type (string) = the SYBYL atom type for the atom.
- subst_id (integer) = the ID number of the substructure containing the atom.
- subst_name (string) = the name of the substructure containing the atom.
- charge (real) = the charge associated with the atom.
- status_bit (string) = the internal SYBYL status bits associated with the atom. These should never be set by the user. Valid status bits are DSPMOD, TYPECOL, CAP, BACKBONE, DICT, ESSENTIAL, WATER and DIRECT.

For example, the contents of a typical Tripos MOL2 file may look like this:

```
@<TRIPOS>MOLECULE
DCM Pose 1
   32    33     0     0     0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
      1 C1         18.8934    5.5819   24.1747 C.2       1 <0>       -0.1356 
      2 C2         18.1301    4.7642   24.8969 C.2       1 <0>       -0.0410 
      3 C3         18.2645    6.8544   23.7342 C.2       1 <0>        0.4856 
...
     31 H11        18.5977    8.5756   22.6932 H         1 <0>        0.4000 
     32 H12        14.2530    1.0535   27.4278 H         1 <0>        0.4000 
@<TRIPOS>BOND
    1     1     2 2
    2     1     3 1
    3     2    11 1
    4     3    10 2
    5     3    12 1
...
   28     8    27 1
   29     9    28 1
   30     9    29 1
   31    12    30 1
   32    12    31 1
   33    18    32 1
```

## Working with MOL2 DataFrames

In the previous sections, we've seen how to load MOL2 structures into DataFrames and how to access them. Once, we have the ATOM section of a MOL2 file in a DataFrame format, we can readily slice and dice the molecular structure and analyze it.
To demonstrate some typical use cases, let us load the structure of deoxycytidylate hydroxymethylase (DCM), which is shown in the figure below:

![](./img/1b5e_1.png)


```python
from biopandas.mol2 import PandasMol2

pmol = PandasMol2()
pmol.read_mol2('./data/1b5e_1.mol2')
pmol.df.tail(10)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>atom_id</th>
      <th>atom_name</th>
      <th>x</th>
      <th>y</th>
      <th>...</th>
      <th>atom_type</th>
      <th>subst_id</th>
      <th>subst_name</th>
      <th>charge</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>22</th>
      <td>23</td>
      <td>H3</td>
      <td>15.8520</td>
      <td>2.8983</td>
      <td>...</td>
      <td>H</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>23</th>
      <td>24</td>
      <td>H4</td>
      <td>14.3405</td>
      <td>3.3601</td>
      <td>...</td>
      <td>H</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>24</th>
      <td>25</td>
      <td>H5</td>
      <td>15.3663</td>
      <td>0.9351</td>
      <td>...</td>
      <td>H</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>25</th>
      <td>26</td>
      <td>H6</td>
      <td>16.6681</td>
      <td>1.6130</td>
      <td>...</td>
      <td>H</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>26</th>
      <td>27</td>
      <td>H7</td>
      <td>15.3483</td>
      <td>4.6961</td>
      <td>...</td>
      <td>H</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>27</th>
      <td>28</td>
      <td>H8</td>
      <td>18.8490</td>
      <td>1.8078</td>
      <td>...</td>
      <td>H</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>28</th>
      <td>29</td>
      <td>H9</td>
      <td>17.8303</td>
      <td>1.5497</td>
      <td>...</td>
      <td>H</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>29</th>
      <td>30</td>
      <td>H10</td>
      <td>19.9527</td>
      <td>7.4708</td>
      <td>...</td>
      <td>H</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>0.4</td>
    </tr>
    <tr>
      <th>30</th>
      <td>31</td>
      <td>H11</td>
      <td>18.5977</td>
      <td>8.5756</td>
      <td>...</td>
      <td>H</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>0.4</td>
    </tr>
    <tr>
      <th>31</th>
      <td>32</td>
      <td>H12</td>
      <td>14.2530</td>
      <td>1.0535</td>
      <td>...</td>
      <td>H</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>0.4</td>
    </tr>
  </tbody>
</table>
<p>10 rows × 9 columns</p>
</div>



[File link: [1b5e_1.mol2](https://raw.githubusercontent.com/rasbt/biopandas/master/docs/sources/tutorials/data/1b5e_1.mol2)]

For example, we can select all hydrogen atoms by filtering on the atom type column:


```python
pmol.df[pmol.df['atom_type'] != 'H'].tail(10)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>atom_id</th>
      <th>atom_name</th>
      <th>x</th>
      <th>y</th>
      <th>...</th>
      <th>atom_type</th>
      <th>subst_id</th>
      <th>subst_name</th>
      <th>charge</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>10</th>
      <td>11</td>
      <td>N2</td>
      <td>16.8196</td>
      <td>5.0644</td>
      <td>...</td>
      <td>N.am</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-0.4691</td>
    </tr>
    <tr>
      <th>11</th>
      <td>12</td>
      <td>N3</td>
      <td>19.0194</td>
      <td>7.7275</td>
      <td>...</td>
      <td>N.pl3</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-0.8500</td>
    </tr>
    <tr>
      <th>12</th>
      <td>13</td>
      <td>O1</td>
      <td>18.7676</td>
      <td>-2.3524</td>
      <td>...</td>
      <td>O.3</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-1.0333</td>
    </tr>
    <tr>
      <th>13</th>
      <td>14</td>
      <td>O2</td>
      <td>20.3972</td>
      <td>-0.3812</td>
      <td>...</td>
      <td>O.3</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-1.0333</td>
    </tr>
    <tr>
      <th>14</th>
      <td>15</td>
      <td>O3</td>
      <td>15.0888</td>
      <td>6.5824</td>
      <td>...</td>
      <td>O.2</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-0.5700</td>
    </tr>
    <tr>
      <th>15</th>
      <td>16</td>
      <td>O4</td>
      <td>18.9314</td>
      <td>-0.7527</td>
      <td>...</td>
      <td>O.2</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-1.0333</td>
    </tr>
    <tr>
      <th>16</th>
      <td>17</td>
      <td>O5</td>
      <td>16.9690</td>
      <td>3.4315</td>
      <td>...</td>
      <td>O.3</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-0.5600</td>
    </tr>
    <tr>
      <th>17</th>
      <td>18</td>
      <td>O6</td>
      <td>14.3223</td>
      <td>1.8946</td>
      <td>...</td>
      <td>O.3</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-0.6800</td>
    </tr>
    <tr>
      <th>18</th>
      <td>19</td>
      <td>O7</td>
      <td>17.9091</td>
      <td>-0.0135</td>
      <td>...</td>
      <td>O.3</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-0.5512</td>
    </tr>
    <tr>
      <th>19</th>
      <td>20</td>
      <td>P1</td>
      <td>19.0969</td>
      <td>-0.9440</td>
      <td>...</td>
      <td>P.3</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>1.3712</td>
    </tr>
  </tbody>
</table>
<p>10 rows × 9 columns</p>
</div>



Or, if we like  to count the number of keto-groups in this molecule, we can do the following:


```python
keto = pmol.df[pmol.df['atom_type'] == 'O.2']
print('number of keto groups: %d' % keto.shape[0])
keto
```

    number of keto groups: 2





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>atom_id</th>
      <th>atom_name</th>
      <th>x</th>
      <th>y</th>
      <th>...</th>
      <th>atom_type</th>
      <th>subst_id</th>
      <th>subst_name</th>
      <th>charge</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>14</th>
      <td>15</td>
      <td>O3</td>
      <td>15.0888</td>
      <td>6.5824</td>
      <td>...</td>
      <td>O.2</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-0.5700</td>
    </tr>
    <tr>
      <th>15</th>
      <td>16</td>
      <td>O4</td>
      <td>18.9314</td>
      <td>-0.7527</td>
      <td>...</td>
      <td>O.2</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-1.0333</td>
    </tr>
  </tbody>
</table>
<p>2 rows × 9 columns</p>
</div>



A list of all the allowed atom types that can be found in Tripos MOL2 files is provided below:

    Code       Definition
    C.3        carbon sp3
    C.2        carbon sp2
    C.1        carbon sp
    C.ar       carbon aromatic
    C.cat      cabocation (C+) used only in a guadinium group
    N.3        nitrogen sp3
    N.2        nitrogen sp2
    N.1        nitrogen sp
    N.ar       nitrogen aromatic
    N.am       nitrogen amide
    N.pl3      nitrogen trigonal planar
    N.4        nitrogen sp3 positively charged
    O.3        oxygen sp3
    O.2        oxygen sp2
    O.co2      oxygen in carboxylate and phosphate groups
    O.spc      oxygen in Single Point Charge (SPC) water model
    O.t3p      oxygen in Transferable Intermolecular Potential (TIP3P) water model
    S.3        sulfur sp3
    S.2        sulfur sp2
    S.O        sulfoxide sulfur
    S.O2/S.o2  sulfone sulfur
    P.3        phosphorous sp3
    F          fluorine
    H          hydrogen
    H.spc      hydrogen in Single Point Charge (SPC) water model
    H.t3p      hydrogen in Transferable Intermolecular Potential (TIP3P) water model
    LP         lone pair
    Du         dummy atom
    Du.C       dummy carbon
    Any        any atom
    Hal        halogen
    Het        heteroatom = N, O, S, P
    Hev        heavy atom (non hydrogen)
    Li         lithium
    Na         sodium
    Mg         magnesium
    Al         aluminum
    Si         silicon
    K          potassium
    Ca         calcium
    Cr.thm     chromium (tetrahedral)
    Cr.oh      chromium (octahedral)
    Mn         manganese
    Fe         iron
    Co.oh      cobalt (octahedral)
    Cu         copper


## Plotting

Since we are using pandas under the hood, which in turns uses matplotlib under the hood, we can produce quick summary plots of our MOL2 structures conveniently. Below are a few examples of how to visualize molecular properties.


```python
from biopandas.mol2 import PandasMol2

pmol = PandasMol2().read_mol2('./data/1b5e_1.mol2')
```

[File link: [1b5e_1.mol2](https://raw.githubusercontent.com/rasbt/biopandas/master/docs/sources/tutorials/data/1b5e_1.mol2)]


```python
%matplotlib inline
import matplotlib.pyplot as plt
from matplotlib import style
style.use('ggplot')
```

For instance, let's say we are interested in the counts of the different atom types that can be found in the MOL2 file; we could do the following:


```python
pmol.df['atom_type'].value_counts().plot(kind='bar')
plt.xlabel('atom type')
plt.ylabel('count')
plt.show()
```


![png](Working_with_MOL2_Structures_in_DataFrames_files/Working_with_MOL2_Structures_in_DataFrames_34_0.png)


If this is too fine-grained for our needs, we could summarize the different atom types by atomic elements:


```python
pmol.df['element_type'] = pmol.df['atom_type'].apply(lambda x: x.split('.')[0])

pmol.df['element_type'].value_counts().plot(kind='bar')
plt.xlabel('element type')
plt.ylabel('count')
plt.show()
```


![png](Working_with_MOL2_Structures_in_DataFrames_files/Working_with_MOL2_Structures_in_DataFrames_36_0.png)


One of the coolest features in pandas is the groupby method. Below is an example plotting the average charge of the different atom types with the standard deviation as error bars:


```python
groupby_charge = pmol.df.groupby(['atom_type'])['charge']
groupby_charge.mean().plot(kind='bar', yerr=groupby_charge.std())
plt.ylabel('charge')
plt.show()
```


![png](Working_with_MOL2_Structures_in_DataFrames_files/Working_with_MOL2_Structures_in_DataFrames_38_0.png)


## Computing the Root Mean Square Deviation

The Root-mean-square deviation (RMSD) is simply a measure of the average distance between atoms of 2 structures. This calculation of the Cartesian error follows the equation:

$$RMSD(a, b) = \sqrt{\frac{1}{n} \sum^{n}_{i=1} \big((a_{ix})^2 + (a_{iy})^2 + (a_{iz})^2 \big)} \\
= \sqrt{\frac{1}{n} \sum^{n}_{i=1} || a_i + b_i||_2^2}$$

So, assuming that the we have the following 2 conformations of a ligand molecule

![](./img/1b5e_poses.png)

we can compute the RMSD as follows:


```python
from biopandas.mol2 import PandasMol2

l_1 = PandasMol2().read_mol2('./data/1b5e_1.mol2')
l_2 = PandasMol2().read_mol2('./data/1b5e_2.mol2')

r_heavy = PandasMol2.rmsd(l_1.df, l_2.df)
r_all  = PandasMol2.rmsd(l_1.df, l_2.df, heavy_only=False)

print('Heavy-atom RMSD: %.4f Angstrom' % r_heavy)
print('All-atom RMSD: %.4f Angstrom' % r_all)
```

    Heavy-atom RMSD: 1.1609 Angstrom
    All-atom RMSD: 1.5523 Angstrom


[File links: [1b5e_1.mol2](https://raw.githubusercontent.com/rasbt/biopandas/master/docs/sources/tutorials/data/1b5e_1.mol2), [1b5e_2.mol2](https://raw.githubusercontent.com/rasbt/biopandas/master/docs/sources/tutorials/data/1b5e_2.mol2)]

<br>

## Filtering Atoms by Distance

We can use the `distance` method to compute the distance between each atom (or a subset of atoms) in our data frame and a three-dimensional reference point. For example, let's assume were are interested in computing the distance between a keto group in the DMC molecule, which we've seen earlier, and other atoms in the same molecule.

First, let's get the coordinates of all keto-groups in this molecule:


```python
from biopandas.mol2 import PandasMol2

pmol = PandasMol2().read_mol2('./data/1b5e_1.mol2')

keto_coord = pmol.df[pmol.df['atom_type'] == 'O.2'][['x', 'y', 'z']]
keto_coord
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>x</th>
      <th>y</th>
      <th>z</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>14</th>
      <td>15.0888</td>
      <td>6.5824</td>
      <td>25.0727</td>
    </tr>
    <tr>
      <th>15</th>
      <td>18.9314</td>
      <td>-0.7527</td>
      <td>24.1606</td>
    </tr>
  </tbody>
</table>
</div>



In the following example, we use `PandasMol2`'s `distance` method. The `distance` method returns a pandas `Series` object containing the Euclidean distance between an atom and all other atoms in the structure. In the following example, `keto_coord.values[0]` refers to the x, y, z coordinates of the first row (i.e., first keto group) in the array above:


```python
print('x, y, z coords:', keto_coord.values[0])
distances = pmol.distance(keto_coord.values[0])
```

    x, y, z coords: [15.0888  6.5824 25.0727]


For our convenience, we can add these `distances` to our MOL2 DataFrame:


```python
pmol.df['distances'] = distances
pmol.df.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>atom_id</th>
      <th>atom_name</th>
      <th>x</th>
      <th>y</th>
      <th>...</th>
      <th>subst_id</th>
      <th>subst_name</th>
      <th>charge</th>
      <th>distances</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>C1</td>
      <td>18.8934</td>
      <td>5.5819</td>
      <td>...</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-0.1356</td>
      <td>4.035144</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>C2</td>
      <td>18.1301</td>
      <td>4.7642</td>
      <td>...</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-0.0410</td>
      <td>3.547712</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>C3</td>
      <td>18.2645</td>
      <td>6.8544</td>
      <td>...</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>0.4856</td>
      <td>3.456969</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>C4</td>
      <td>16.2520</td>
      <td>6.2866</td>
      <td>...</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>0.8410</td>
      <td>1.232313</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>C5</td>
      <td>15.3820</td>
      <td>3.0682</td>
      <td>...</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>0.0000</td>
      <td>3.527546</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 10 columns</p>
</div>



Now, say we are interested in the Euclidean distance between the two keto groups in the molecule:


```python
pmol.df[pmol.df['atom_type'] == 'O.2']
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>atom_id</th>
      <th>atom_name</th>
      <th>x</th>
      <th>y</th>
      <th>...</th>
      <th>subst_id</th>
      <th>subst_name</th>
      <th>charge</th>
      <th>distances</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>14</th>
      <td>15</td>
      <td>O3</td>
      <td>15.0888</td>
      <td>6.5824</td>
      <td>...</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-0.5700</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>15</th>
      <td>16</td>
      <td>O4</td>
      <td>18.9314</td>
      <td>-0.7527</td>
      <td>...</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-1.0333</td>
      <td>8.330738</td>
    </tr>
  </tbody>
</table>
<p>2 rows × 10 columns</p>
</div>



In the example above, the distance between the two keto groups is 8 angstrom.

![](./img/1b5e_ketodist.png)

Another common task that we can perform using these atomic distances is to select only the neighboring atoms of the keto group (here: atoms within 3 angstrom). The code is as follows:


```python
all_within_3A = pmol.df[pmol.df['distances'] <= 3.0]
all_within_3A.tail()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>atom_id</th>
      <th>atom_name</th>
      <th>x</th>
      <th>y</th>
      <th>...</th>
      <th>subst_id</th>
      <th>subst_name</th>
      <th>charge</th>
      <th>distances</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>7</th>
      <td>8</td>
      <td>C8</td>
      <td>16.0764</td>
      <td>4.1199</td>
      <td>...</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>0.5801</td>
      <td>2.814490</td>
    </tr>
    <tr>
      <th>9</th>
      <td>10</td>
      <td>N1</td>
      <td>17.0289</td>
      <td>7.1510</td>
      <td>...</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-0.6610</td>
      <td>2.269690</td>
    </tr>
    <tr>
      <th>10</th>
      <td>11</td>
      <td>N2</td>
      <td>16.8196</td>
      <td>5.0644</td>
      <td>...</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-0.4691</td>
      <td>2.307553</td>
    </tr>
    <tr>
      <th>14</th>
      <td>15</td>
      <td>O3</td>
      <td>15.0888</td>
      <td>6.5824</td>
      <td>...</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>-0.5700</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>26</th>
      <td>27</td>
      <td>H7</td>
      <td>15.3483</td>
      <td>4.6961</td>
      <td>...</td>
      <td>1</td>
      <td>&lt;0&gt;</td>
      <td>0.0000</td>
      <td>2.446817</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 10 columns</p>
</div>



## Parsing Multi-MOL2 files

### Basic Multi-MOL2 File Parsing

As mentioned earlier, `PandasMol2.read_mol2` method only reads in the first molecule if it is given a multi-MOL2 file. However, if we want to create DataFrames from multiple structures in a MOL2 file, we can use the handy `split_multimol2` generator.

The `split_multimol2` generator yields tuples containing the molecule IDs and the MOL2 content as strings in a list -- each line in the MOL2 file is stored as a string in the list.


```python
from biopandas.mol2 import split_multimol2

mol2_id, mol2_cont = next(split_multimol2('./data/40_mol2_files.mol2'))

print('Molecule ID:\n', mol2_id)
print('First 10 lines:\n', mol2_cont[:10])
```

    Molecule ID:
     ZINC38611810
    First 10 lines:
     ['@<TRIPOS>MOLECULE\n', 'ZINC38611810\n', '   65    68     0     0     0\n', 'SMALL\n', 'NO_CHARGES\n', '\n', '@<TRIPOS>ATOM\n', '      1 C1         -1.1786    2.7011   -4.0323 C.3       1 <0>        -0.1537\n', '      2 C2         -1.2950    1.2442   -3.5798 C.3       1 <0>        -0.1156\n', '      3 C3         -0.1742    0.4209   -4.2178 C.3       1 <0>        -0.1141\n']


[File link: [40_mol2_files.mol2](https://raw.githubusercontent.com/rasbt/biopandas/master/docs/sources/tutorials/data/40_mol2_files.mol2)]

We can now use this generator to loop over all files in a multi-MOL2 file and create PandasMol2 DataFrames. A typical use case would be the filtering of mol2 files by certain properties:


```python
pdmol = PandasMol2()

with open('./data/filtered.mol2', 'w') as f:
    for mol2 in split_multimol2('./data/40_mol2_files.mol2'):
        pdmol.read_mol2_from_list(mol2_lines=mol2[1], mol2_code=mol2[0])
        
        # do some analysis
        keep_molecule = False
        
        # save molecule if it passes our filter criterion
        if keep_molecule: 
            # note that the mol2_text contains the original mol2 content
            f.write(pdmol.mol2_text) 
```

### Using Multiprocessing for Multi-MOL2 File Analysis

To improve the computational efficiency and throughput for multi-mol2 analyses, it is recommended to use the [`mputil`](https://github.com/rasbt/mputil) package, which evaluates Python generators lazily. The `lazy_imap` function from `mputil` is based on Python's standardlib multiprocessing `imap` function, but it doesn't consume the generator upfront. This lazy evaluation is important, for example, if we are parsing large (possibly Gigabyte- or Terabyte-large) multi-mol2 files for multiprocessing.

The following example provides a template for atom-type based molecule queries, but the `data_processor` function can be extended to do any kind of functional group queries (for example, involving the `'charge'` column and/or `PandasMol2.distance` method). 


```python
import pandas as pd
from mputil import lazy_imap
from biopandas.mol2 import PandasMol2
from biopandas.mol2 import split_multimol2
```


    ---------------------------------------------------------------------------

    ModuleNotFoundError                       Traceback (most recent call last)

    <ipython-input-23-1a655249f2ec> in <module>()
          1 import pandas as pd
    ----> 2 from mputil import lazy_imap
          3 from biopandas.mol2 import PandasMol2
          4 from biopandas.mol2 import split_multimol2


    ModuleNotFoundError: No module named 'mputil'



```python
# Selection strings to capture
# all molecules that contain at least one sp2 hybridized
# oxygen atom and at least one Fluorine atom
SELECTIONS = ["(pdmol.df.atom_type == 'O.2')",
              "(pdmol.df.atom_type == 'F')"]
 
# Path to the multi-mol2 input file
MOL2_FILE = "./data/40_mol2_files.mol2"

# Data processing function to be run in parallel
def data_processor(mol2):
    """Return molecule ID if there's a match and '' otherwise"""
    pdmol = PandasMol2().read_mol2_from_list(mol2_lines=mol2[1],
                                             mol2_code=mol2[0])
    match = mol2[0]
    for sub_sele in SELECTIONS:
        if not pd.eval(sub_sele).any():
            match = ''
            break
    return match

# Process molecules and save IDs of hits to disk
with open('./data/selected_ids.txt', 'w') as f:       
    searched, found = 0, 0
    for chunk in lazy_imap(data_processor=data_processor,
                           data_generator=split_multimol2(MOL2_FILE),
                           n_cpus=0): # means all available cpus

        for mol2_id in chunk:
            if mol2_id:
                # write IDs of matching molecules to disk
                f.write('%s\n' % mol2_id)
                found += 1
        searched += len(chunk)
            
        
print('Searched %d molecules. Got %d hits.' % (searched, found))
```

[Input File link: [40_mol2_files.mol2](https://raw.githubusercontent.com/rasbt/biopandas/master/docs/sources/tutorials/data/40_mol2_files.mol2)]

[Output File link: [selected_ids.txt](https://raw.githubusercontent.com/rasbt/biopandas/master/docs/sources/tutorials/data/selected_ids.txt)]
