# Working with PDB Structures in DataFrames

## Loading PDB Files

There are 2 1/2 ways to load a PDB structure into a `PandasPDB` object.


#### 1
PDB files can be directly fetched from The Protein Data Bank at [http://www.rcsb.org](http://www.rcsb.org) via its unique 4-letter after initializing a new [`PandasPDB`](../api/biopandas.pdb#pandaspdb) object and calling the [`fetch_pdb`](../api/biopandas.pdb#pandaspdbfetch_pdb) method:


```python
from biopandas.pdb import PandasPDB

# Initialize a new PandasPDB object
# and fetch the PDB file from rcsb.org
ppdb = PandasPDB().fetch_pdb('3eiy')
```

#### 2 a)

Alternatively, we can load PDB files from local directories as regular PDB files using [`read_pdb`](../api/biopandas.pdb#pandaspdbread_pdb):


```python
ppdb.read_pdb('./data/3eiy.pdb')
```




#### 2 b)

Or, we can load them from gzip archives like so (note that the file must end with a '.gz' suffix in order to be recognized as a gzip file):


```python
ppdb.read_pdb('./data/3eiy.pdb.gz')
```




After the file was succesfully loaded, we have access to the following attributes:


```python
print('PDB Code: %s' % ppdb.code)
print('PDB Header Line: %s' % ppdb.header)
print('\nRaw PDB file contents:\n\n%s\n...' % ppdb.pdb_text[:1000])
```

    PDB Code: 3eiy
    PDB Header Line:     HYDROLASE                               17-SEP-08   3EIY
    
    Raw PDB file contents:
    
    HEADER    HYDROLASE                               17-SEP-08   3EIY              
    TITLE     CRYSTAL STRUCTURE OF INORGANIC PYROPHOSPHATASE FROM BURKHOLDERIA      
    TITLE    2 PSEUDOMALLEI WITH BOUND PYROPHOSPHATE                                
    COMPND    MOL_ID: 1;                                                            
    COMPND   2 MOLECULE: INORGANIC PYROPHOSPHATASE;                                 
    COMPND   3 CHAIN: A;                                                            
    COMPND   4 EC: 3.6.1.1;                                                         
    COMPND   5 ENGINEERED: YES                                                      
    SOURCE    MOL_ID: 1;                                                            
    SOURCE   2 ORGANISM_SCIENTIFIC: BURKHOLDERIA PSEUDOMALLEI 1710B;                
    SOURCE   3 ORGANISM_TAXID: 320372;                                              
    SOURCE   4 GENE: PPA, BURPS1710B_1237;                                          
    SOURCE   5 EXPRESSION_SYSTEM
    ...


The most interesting / useful attribute is the [`PandasPDB.df`](../api/biopandas.pdb#pandaspdbdf) DataFrame dictionary though, which gives us access to the PDB files as pandas DataFrames. Let's print the first 3 lines from the `ATOM` coordinate section to see how it looks like:


```python
ppdb.df['ATOM'].head(3)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>record_name</th>
      <th>atom_number</th>
      <th>blank_1</th>
      <th>atom_name</th>
      <th>...</th>
      <th>segment_id</th>
      <th>element_symbol</th>
      <th>charge</th>
      <th>line_idx</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ATOM</td>
      <td>1</td>
      <td></td>
      <td>N</td>
      <td>...</td>
      <td></td>
      <td>N</td>
      <td>NaN</td>
      <td>609</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ATOM</td>
      <td>2</td>
      <td></td>
      <td>CA</td>
      <td>...</td>
      <td></td>
      <td>C</td>
      <td>NaN</td>
      <td>610</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ATOM</td>
      <td>3</td>
      <td></td>
      <td>C</td>
      <td>...</td>
      <td></td>
      <td>C</td>
      <td>NaN</td>
      <td>611</td>
    </tr>
  </tbody>
</table>
<p>3 rows × 21 columns</p>
</div>



But more on that in the next section.

## Looking at PDBs in DataFrames

PDB files are parsed according to the [PDB file format description](http://www.rcsb.org/pdb/static.do?p=file_formats/pdb/index.html). More specifically, BioPandas reads the columns of the ATOM, HETATM, and ANISOU section as shown in the following excerpt from [http://deposit.rcsb.org/adit/docs/pdb_atom_format.html#ATOM](http://deposit.rcsb.org/adit/docs/pdb_atom_format.html#ATOM); note that the ANISOU section looks a bit different: [http://deposit.rcsb.org/adit/docs/pdb_atom_format.html#ANISOU](http://deposit.rcsb.org/adit/docs/pdb_atom_format.html#ANISOU)

| COLUMNS | DATA TYPE    | CONTENTS                                   | biopandas column name |
|---------|--------------|--------------------------------------------|-----------------------|
| 1 - 6   | Record name  | "ATOM"                                     | record_name           |
| 7 - 11  | Integer      | Atom serial number.                        | atom_number           |
| 12      |              |                                            | blank_1               |
| 13 - 16 | Atom         | Atom name.                                 | atom_name             |
| 17      | Character    | Alternate location indicator.              | alt_loc               |
| 18 - 20 | Residue name | Residue name.                              | residue_name          |
| 21      |              |                                            | blank_2               |
| 22      | Character    | Chain identifier.                          | chain_id              |
| 23 - 26 | Integer      | Residue sequence number.                   | residue_number        |
| 27      | AChar        | Code for insertion of residues.            | insertion             |
| 28 - 30 |              |                                            | blank_3               |
| 31 - 38 | Real(8.3)    | Orthogonal coordinates for X in Angstroms. | x_coord               |
| 39 - 46 | Real(8.3)    | Orthogonal coordinates for Y in Angstroms. | y_coord               |
| 47 - 54 | Real(8.3)    | Orthogonal coordinates for Z in Angstroms. | z_coord               |
| 55 - 60 | Real(6.2)    | Occupancy.                                 | occupancy             |
| 61 - 66 | Real(6.2)    | Temperature factor (Default = 0.0).        | bfactor               |
| 67-72   |              |                                            | blank_4               |
| 73 - 76 | LString(4)   | Segment identifier, left-justified.        | segment_id            |
| 77 - 78 | LString(2)   | Element symbol, right-justified.           | element_symbol        |
| 79 - 80 | LString(2)   | Charge on the atom.                        | charge                |

    Example: 
             1         2         3         4         5         6         7         8
    12345678901234567890123456789012345678901234567890123456789012345678901234567890
    ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92      A1   N
    ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85      A1   C
    ATOM    147  C   VAL A  25      30.447  15.105  58.363  1.00 12.34      A1   C
    ATOM    148  O   VAL A  25      29.520  15.059  59.174  1.00 15.65      A1   O
    ATOM    149  CB AVAL A  25      30.385  17.437  57.230  0.28 13.88      A1   C
    ATOM    150  CB BVAL A  25      30.166  17.399  57.373  0.72 15.41      A1   C
    ATOM    151  CG1AVAL A  25      28.870  17.401  57.336  0.28 12.64      A1   C
    ATOM    152  CG1BVAL A  25      30.805  18.788  57.449  0.72 15.11      A1   C
    ATOM    153  CG2AVAL A  25      30.835  18.826  57.661  0.28 13.58      A1   C
    ATOM    154  CG2BVAL A  25      29.909  16.996  55.922  0.72 13.25      A1   C

After loading a PDB file from rcsb.org or our local drive, the [`PandasPDB.df`](../api/biopandas.pdb/#pandaspdbdf) attribute should contain the following 4 DataFrame objects:


```python
from biopandas.pdb import PandasPDB
ppdb = PandasPDB()
ppdb.read_pdb('./data/3eiy.pdb')
ppdb.df.keys()
```




    dict_keys(['ATOM', 'OTHERS', 'HETATM', 'ANISOU'])



- 'ATOM': contains the entries from the ATOM coordinate section
- 'ATOM':  ... entries from the "HETATM" coordinate section    
- 'ANISOU': ... entries from the "ANISOU" coordinate section 
- 'OTHERS': Everything else that is *not* a 'ATOM', 'HETATM', or 'ANISOU' entry

![](./img/df_dict.jpg)

The columns of the 'HETATM' DataFrame are indentical to the 'ATOM' DataFrame that we've seen earlier:


```python
ppdb.df['HETATM'].head(2)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>record_name</th>
      <th>atom_number</th>
      <th>blank_1</th>
      <th>atom_name</th>
      <th>...</th>
      <th>segment_id</th>
      <th>element_symbol</th>
      <th>charge</th>
      <th>line_idx</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>HETATM</td>
      <td>1332</td>
      <td></td>
      <td>K</td>
      <td>...</td>
      <td></td>
      <td>K</td>
      <td>NaN</td>
      <td>1940</td>
    </tr>
    <tr>
      <th>1</th>
      <td>HETATM</td>
      <td>1333</td>
      <td></td>
      <td>NA</td>
      <td>...</td>
      <td></td>
      <td>NA</td>
      <td>NaN</td>
      <td>1941</td>
    </tr>
  </tbody>
</table>
<p>2 rows × 21 columns</p>
</div>



And 'ANISOU' entries are handled a bit differently, ...


```python
ppdb.df['ANISOU'].head(2)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>record_name</th>
      <th>atom_number</th>
      <th>blank_1</th>
      <th>atom_name</th>
      <th>...</th>
      <th>blank_4</th>
      <th>element_symbol</th>
      <th>charge</th>
      <th>line_idx</th>
    </tr>
  </thead>
  <tbody>
  </tbody>
</table>
<p>0 rows × 21 columns</p>
</div>



Since the DataFrame is pretty broad, it'd probably better so show the column names separately:


```python
ppdb.df['ANISOU'].columns
```




    Index(['record_name', 'atom_number', 'blank_1', 'atom_name', 'alt_loc', 'residue_name', 'blank_2', 'chain_id', 'residue_number', 'insertion', 'blank_3', 'U(1,1)', 'U(2,2)', 'U(3,3)', 'U(1,2)', 'U(1,3)', 'U(2,3)', 'blank_4', 'element_symbol', 'charge', 'line_idx'], dtype='object')



Ah, another interesting thing to mention is that the columns already come with the types you'd expect (where `object` essentially "means" `str` here):


```python
ppdb.df['ATOM'].dtypes
```




    record_name        object
    atom_number         int64
    blank_1            object
    atom_name          object
    alt_loc            object
    residue_name       object
    blank_2            object
    chain_id           object
    residue_number      int64
    insertion          object
    blank_3            object
    x_coord           float64
    y_coord           float64
    z_coord           float64
    occupancy         float64
    b_factor          float64
    blank_4            object
    segment_id         object
    element_symbol     object
    charge            float64
    line_idx            int64
    dtype: object



As it appears, our '3eiy' structure does not contain any 'ANISOU' entries, hence, the DataFrame above remains empty. 

Typically, all good things come in threes, however, there is a 4th DataFrame, an'OTHER' DataFrame, which contains everything that wasn't parsed as 'ATOM', 'HETATM', or 'ANISOU' coordinate section:


```python
ppdb.df['OTHERS'].head(5)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>record_name</th>
      <th>entry</th>
      <th>line_idx</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>HEADER</td>
      <td>HYDROLASE                               17...</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>TITLE</td>
      <td>CRYSTAL STRUCTURE OF INORGANIC PYROPHOSPHA...</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>TITLE</td>
      <td>2 PSEUDOMALLEI WITH BOUND PYROPHOSPHATE</td>
      <td>2</td>
    </tr>
    <tr>
      <th>3</th>
      <td>COMPND</td>
      <td>MOL_ID: 1;</td>
      <td>3</td>
    </tr>
    <tr>
      <th>4</th>
      <td>COMPND</td>
      <td>2 MOLECULE: INORGANIC PYROPHOSPHATASE;</td>
      <td>4</td>
    </tr>
  </tbody>
</table>
</div>



Although these 'OTHER' entries are typically less useful for structure-related computations, you may still want to take a look at them to get a short summary of the PDB structure and learn about it's potential quirks and gotchas (typically listed in the REMARKs section). Lastly, the 'OTHERS' DataFrame comes in handy if we want to reconstruct the structure as PDB file as we will see later (note the `line_idx` columns in all of the DataFrames).

## Working with PDB DataFrames

In the previous sections, we've seen how to load PDB structures into DataFrames, and how to access them. Now, let's talk about manipulating PDB files in DataFrames.


```python
from biopandas.pdb import PandasPDB
ppdb = PandasPDB()
ppdb.read_pdb('./data/3eiy.pdb.gz')
ppdb.df['ATOM'].head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>record_name</th>
      <th>atom_number</th>
      <th>blank_1</th>
      <th>atom_name</th>
      <th>...</th>
      <th>segment_id</th>
      <th>element_symbol</th>
      <th>charge</th>
      <th>line_idx</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ATOM</td>
      <td>1</td>
      <td></td>
      <td>N</td>
      <td>...</td>
      <td></td>
      <td>N</td>
      <td>NaN</td>
      <td>609</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ATOM</td>
      <td>2</td>
      <td></td>
      <td>CA</td>
      <td>...</td>
      <td></td>
      <td>C</td>
      <td>NaN</td>
      <td>610</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ATOM</td>
      <td>3</td>
      <td></td>
      <td>C</td>
      <td>...</td>
      <td></td>
      <td>C</td>
      <td>NaN</td>
      <td>611</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ATOM</td>
      <td>4</td>
      <td></td>
      <td>O</td>
      <td>...</td>
      <td></td>
      <td>O</td>
      <td>NaN</td>
      <td>612</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ATOM</td>
      <td>5</td>
      <td></td>
      <td>CB</td>
      <td>...</td>
      <td></td>
      <td>C</td>
      <td>NaN</td>
      <td>613</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 21 columns</p>
</div>



Okay, there's actually not *that* much to say ...   
Once we have our PDB file in the DataFrame format, we have the whole convenience of [pandas](http://pandas.pydata.org) right there at our fingertips.

For example, let's get all Proline residues:


```python
ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name'] == 'PRO'].head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>record_name</th>
      <th>atom_number</th>
      <th>blank_1</th>
      <th>atom_name</th>
      <th>...</th>
      <th>segment_id</th>
      <th>element_symbol</th>
      <th>charge</th>
      <th>line_idx</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>38</th>
      <td>ATOM</td>
      <td>39</td>
      <td></td>
      <td>N</td>
      <td>...</td>
      <td></td>
      <td>N</td>
      <td>NaN</td>
      <td>647</td>
    </tr>
    <tr>
      <th>39</th>
      <td>ATOM</td>
      <td>40</td>
      <td></td>
      <td>CA</td>
      <td>...</td>
      <td></td>
      <td>C</td>
      <td>NaN</td>
      <td>648</td>
    </tr>
    <tr>
      <th>40</th>
      <td>ATOM</td>
      <td>41</td>
      <td></td>
      <td>C</td>
      <td>...</td>
      <td></td>
      <td>C</td>
      <td>NaN</td>
      <td>649</td>
    </tr>
    <tr>
      <th>41</th>
      <td>ATOM</td>
      <td>42</td>
      <td></td>
      <td>O</td>
      <td>...</td>
      <td></td>
      <td>O</td>
      <td>NaN</td>
      <td>650</td>
    </tr>
    <tr>
      <th>42</th>
      <td>ATOM</td>
      <td>43</td>
      <td></td>
      <td>CB</td>
      <td>...</td>
      <td></td>
      <td>C</td>
      <td>NaN</td>
      <td>651</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 21 columns</p>
</div>



Or main chain atoms:


```python
ppdb.df['ATOM'][ppdb.df['ATOM']['atom_name'] == 'C'].head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>record_name</th>
      <th>atom_number</th>
      <th>blank_1</th>
      <th>atom_name</th>
      <th>...</th>
      <th>segment_id</th>
      <th>element_symbol</th>
      <th>charge</th>
      <th>line_idx</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>2</th>
      <td>ATOM</td>
      <td>3</td>
      <td></td>
      <td>C</td>
      <td>...</td>
      <td></td>
      <td>C</td>
      <td>NaN</td>
      <td>611</td>
    </tr>
    <tr>
      <th>8</th>
      <td>ATOM</td>
      <td>9</td>
      <td></td>
      <td>C</td>
      <td>...</td>
      <td></td>
      <td>C</td>
      <td>NaN</td>
      <td>617</td>
    </tr>
    <tr>
      <th>19</th>
      <td>ATOM</td>
      <td>20</td>
      <td></td>
      <td>C</td>
      <td>...</td>
      <td></td>
      <td>C</td>
      <td>NaN</td>
      <td>628</td>
    </tr>
    <tr>
      <th>25</th>
      <td>ATOM</td>
      <td>26</td>
      <td></td>
      <td>C</td>
      <td>...</td>
      <td></td>
      <td>C</td>
      <td>NaN</td>
      <td>634</td>
    </tr>
    <tr>
      <th>33</th>
      <td>ATOM</td>
      <td>34</td>
      <td></td>
      <td>C</td>
      <td>...</td>
      <td></td>
      <td>C</td>
      <td>NaN</td>
      <td>642</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 21 columns</p>
</div>



It's also easy to strip our coordinate section from hydrogen atoms if there are any ...


```python
ppdb.df['ATOM'][ppdb.df['ATOM']['element_symbol'] != 'H'].head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>record_name</th>
      <th>atom_number</th>
      <th>blank_1</th>
      <th>atom_name</th>
      <th>...</th>
      <th>segment_id</th>
      <th>element_symbol</th>
      <th>charge</th>
      <th>line_idx</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ATOM</td>
      <td>1</td>
      <td></td>
      <td>N</td>
      <td>...</td>
      <td></td>
      <td>N</td>
      <td>NaN</td>
      <td>609</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ATOM</td>
      <td>2</td>
      <td></td>
      <td>CA</td>
      <td>...</td>
      <td></td>
      <td>C</td>
      <td>NaN</td>
      <td>610</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ATOM</td>
      <td>3</td>
      <td></td>
      <td>C</td>
      <td>...</td>
      <td></td>
      <td>C</td>
      <td>NaN</td>
      <td>611</td>
    </tr>
    <tr>
      <th>3</th>
      <td>ATOM</td>
      <td>4</td>
      <td></td>
      <td>O</td>
      <td>...</td>
      <td></td>
      <td>O</td>
      <td>NaN</td>
      <td>612</td>
    </tr>
    <tr>
      <th>4</th>
      <td>ATOM</td>
      <td>5</td>
      <td></td>
      <td>CB</td>
      <td>...</td>
      <td></td>
      <td>C</td>
      <td>NaN</td>
      <td>613</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 21 columns</p>
</div>



Or, let's compute the average temperature factor of our protein main chain:


```python
mainchain = ppdb.df['ATOM'][(ppdb.df['ATOM']['atom_name'] == 'C') | 
                            (ppdb.df['ATOM']['atom_name'] == 'O') | 
                            (ppdb.df['ATOM']['atom_name'] == 'N') | 
                            (ppdb.df['ATOM']['atom_name'] == 'CA')]

bfact_mc_avg = mainchain['b_factor'].mean()
print('Average B-Factor [Main Chain]: %.2f' % bfact_mc_avg)
```

    Average B-Factor [Main Chain]: 28.83


## Plotting

Since we are using pandas under the hood, which in turns uses matplotlib under the hood, we can produce quick summary plots of our PDB structures relatively conveniently:


```python
from biopandas.pdb import PandasPDB
ppdb = PandasPDB().read_pdb('./data/3eiy.pdb.gz')
```


```python
%matplotlib inline
import matplotlib.pyplot as plt
from matplotlib import style
style.use('ggplot')
```


```python
ppdb.df['ATOM']['b_factor'].plot(kind='hist')
plt.title('Distribution of B-Factors')
plt.xlabel('B-factor')
plt.ylabel('count')
plt.show()
```


![png](Working_with_PDB_Structures_in_DataFrames_files/Working_with_PDB_Structures_in_DataFrames_53_0.png)



```python
ppdb.df['ATOM']['b_factor'].plot(kind='line')
plt.title('B-Factors Along the Amino Acid Chain')
plt.xlabel('Residue Number')
plt.ylabel('B-factor in $A^2$')
plt.show()
```


![png](Working_with_PDB_Structures_in_DataFrames_files/Working_with_PDB_Structures_in_DataFrames_54_0.png)



```python
ppdb.df['ATOM']['element_symbol'].value_counts().plot(kind='bar')
plt.title('Distribution of Atom Types')
plt.xlabel('elements')
plt.ylabel('count')
plt.show()
```


![png](Working_with_PDB_Structures_in_DataFrames_files/Working_with_PDB_Structures_in_DataFrames_55_0.png)


## Computing the Root Mean Square Deviation

BioPandas also comes with certain convenience functions, for example, ...

The [Root-mean-square deviation] (RMSD) is simply a measure of the average distance between atoms of 2 protein or ligand structures. This calculation of the Cartesian error follows the equation:

$$RMSD(a, b) = \sqrt{\frac{1}{n} \sum^{n}_{i=1} \big((a_{ix})^2 + (a_{iy})^2 + (a_{iz})^2 \big)} \\
= \sqrt{\frac{1}{n} \sum^{n}_{i=1} || a_i + b_i||_2^2}$$

So, assuming that the we have the following 2 conformation of a ligand molecule

![](./img/ligand_rmsd.png)

we can compute the RMSD as follows:


```python
from biopandas.pdb import PandasPDB

l_1 = PandasPDB().read_pdb('./data/lig_conf_1.pdb')
l_2 = PandasPDB().read_pdb('./data/lig_conf_2.pdb')
r = PandasPDB.rmsd(l_1.df['HETATM'], l_2.df['HETATM'], s='hydrogen', invert=True)
print('RMSD: %.4f Angstrom' % r)
```

    RMSD: 2.6444 Angstrom


Note that the `s` parameter in [`PandasPDB.rmsd`](../api/biopandas.pdb#pandaspdbrmsd) specifies the search string "consider all hydrogen atoms" and via the `invert=True` option, we say "Compute the RMSD for everything BUT the hydrogen atoms" or "Compute the RSMD based on the heavy atoms only".

Similarly, we can compute the RMSD between 2 related protein structures:

![](./img/1t48_rmsd.png)

The hydrogen-free RMSD:


```python
p_1 = PandasPDB().read_pdb('./data/1t48_995.pdb')
p_2 = PandasPDB().read_pdb('./data/1t49_995.pdb')
r = PandasPDB.rmsd(p_1.df['ATOM'], p_2.df['ATOM'], s='hydrogen', invert=True)
print('RMSD: %.4f Angstrom' % r)
```

    RMSD: 0.7377 Angstrom


Or the RMSD between the main chains only:


```python
p_1 = PandasPDB().read_pdb('./data/1t48_995.pdb')
p_2 = PandasPDB().read_pdb('./data/1t49_995.pdb')
r = PandasPDB.rmsd(p_1.df['ATOM'], p_2.df['ATOM'], s='main chain', invert=False)
print('RMSD: %.4f Angstrom' % r)
```

    RMSD: 0.4781 Angstrom


<br>

[more to come]

## Wrapping it up - Saving PDB structures

Finally, let's talk about how to get the PDB structures out of the DataFrame format back into the beloved .pdb format.

Let's say we loaded a PDB structure, removed it from it's hydrogens:


```python
from biopandas.pdb import PandasPDB
ppdb = PandasPDB().read_pdb('./data/3eiy.pdb.gz')
ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['element_symbol'] != 'H']
```

We can save the file using the [`PandasPDB.to_pdb`](../api/biopandas.pdb#pandaspdbto_pdb) method:


```python
ppdb.to_pdb(path='./data/3eiy_stripped.pdb', records=None, gz=False, append_newline=True)
```

By default, all records (that is, 'ATOM', 'HETATM', 'OTHERS', 'ANISOU') are written if we set `records=None`. Alternatively, let's say we want to get rid of the 'ANISOU' entries and produce a compressed gzip archive of our PDB structure:


```python
ppdb.to_pdb(path='./data/3eiy_stripped.pdb.gz', records=['ATOM', 'HETATM', 'OTHERS'], gz=True, append_newline=True)
```
