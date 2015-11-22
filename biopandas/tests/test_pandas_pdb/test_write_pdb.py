"""
BioPandas
Author: Sebastian Raschka <mail@sebastianraschka.com>
License: BSD 3 clause
Project Website: http://rasbt.github.io/biopandas/
Code Repository: https://github.com/rasbt/biopandas
"""

from biopandas import PandasPDB
import os
import numpy as np
import pandas as pd
from biopandas.testutils import assertMultiLineEqual
from nose.tools import raises


TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'data', '3eiy.pdb')
OUTFILE = os.path.join(os.path.dirname(__file__), 'data', 'tmp.pdb')
OUTFILE_GZ = os.path.join(os.path.dirname(__file__), 'data', 'tmp.pdb.gz')

def test_defaults():
    pass
    ppdb = PandasPDB()
    ppdb.read_pdb(TESTDATA_FILENAME)
    ppdb.to_pdb(path=OUTFILE)
    with open(TESTDATA_FILENAME, 'r') as f:
        f1 = f.read()
    with open(OUTFILE, 'r') as f:
        f2 = f.read()
    print(ppdb.df['OTHERS'])
    assert f1 == f2
    os.remove(OUTFILE)


def test_defaults():
    """Test private _read_pdb"""
    ppdb = PandasPDB()
    ppdb.read_pdb(TESTDATA_FILENAME)
    ppdb.to_pdb(path=OUTFILE, records=['HETATM'])
    with open(OUTFILE, 'r') as f:
        f1 = f.read()
    os.remove(OUTFILE)
    assert f1 == hetatm


hetatm = """HETATM 1332  K     K A 176      24.990  43.276   0.005  0.50 24.45           K
HETATM 1333  NA   NA A 177       1.633  34.181  11.897  1.00 26.73          NA
HETATM 1334  NA   NA A 178       6.489  35.143   8.444  1.00 30.89          NA
HETATM 1335  P1  POP A 179       1.233  37.542  11.212  1.00 32.68           P
HETATM 1336  O1  POP A 179       1.910  38.831  11.612  1.00 32.62           O
HETATM 1337  O2  POP A 179       1.288  37.475   9.712  1.00 33.46           O
HETATM 1338  O3  POP A 179       1.948  36.362  11.841  1.00 30.47           O
HETATM 1339  O   POP A 179      -0.339  37.641  11.573  1.00 33.48           O
HETATM 1340  P2  POP A 179      -1.193  36.552  12.370  1.00 36.05           P
HETATM 1341  O4  POP A 179      -2.611  36.792  11.954  1.00 33.39           O
HETATM 1342  O5  POP A 179      -1.079  36.873  13.870  1.00 33.69           O
HETATM 1343  O6  POP A 179      -0.735  35.128  12.124  1.00 33.11           O
HETATM 1344  O1  PG4 A 180      25.225  36.156   6.596  1.00 42.63           O
HETATM 1345  C1  PG4 A 180      25.070  35.591   5.302  1.00 44.20           C
HETATM 1346  C2  PG4 A 180      25.728  36.597   4.386  1.00 46.03           C
HETATM 1347  O2  PG4 A 180      24.774  36.857   3.380  1.00 48.37           O
HETATM 1348  C3  PG4 A 180      25.201  37.981   2.614  1.00 47.94           C
HETATM 1349  C4  PG4 A 180      25.268  39.174   3.548  1.00 46.48           C
HETATM 1350  O3  PG4 A 180      26.617  39.507   3.811  1.00 45.59           O
HETATM 1351  C5  PG4 A 180      26.819  40.888   3.597  1.00 45.33           C
HETATM 1352  C6  PG4 A 180      28.266  41.161   3.954  1.00 47.26           C
HETATM 1353  O4  PG4 A 180      28.982  41.378   2.738  1.00 48.16           O
HETATM 1354  C7  PG4 A 180      30.163  40.616   2.604  1.00 47.41           C
HETATM 1355  C8  PG4 A 180      31.331  41.437   3.117  1.00 51.14           C
HETATM 1356  O5  PG4 A 180      32.081  40.607   4.011  1.00 53.41           O
HETATM 1357  C1  PEG A 181      13.856  28.176  23.722  1.00 59.47           C
HETATM 1358  O1  PEG A 181      14.416  29.487  23.735  1.00 59.79           O
HETATM 1359  C2  PEG A 181      13.605  27.741  22.291  1.00 59.20           C
HETATM 1360  O2  PEG A 181      12.908  26.499  22.379  1.00 59.79           O
HETATM 1361  C3  PEG A 181      11.848  26.393  21.440  1.00 59.65           C
HETATM 1362  C4  PEG A 181      10.883  25.309  21.906  1.00 59.90           C
HETATM 1363  O4  PEG A 181      10.045  25.835  22.942  1.00 58.85           O
HETATM 1364  C1  PEG A 182      19.721  56.811   7.363  1.00 56.17           C
HETATM 1365  O1  PEG A 182      18.409  56.724   7.966  1.00 55.97           O
HETATM 1366  C2  PEG A 182      20.835  56.288   8.271  1.00 54.49           C
HETATM 1367  O2  PEG A 182      21.824  57.291   8.516  1.00 56.31           O
HETATM 1368  C3  PEG A 182      23.176  56.838   8.332  1.00 57.82           C
HETATM 1369  C4  PEG A 182      23.815  57.479   7.094  1.00 58.58           C
HETATM 1370  O4  PEG A 182      24.621  56.536   6.357  1.00 59.04           O
HETATM 1371  O   HOH A 183      -6.268  35.297  21.123  1.00 22.39           O
HETATM 1372  O   HOH A 184      18.291  50.223   3.964  1.00 24.68           O
HETATM 1373  O   HOH A 185      22.029  42.728   8.623  1.00 20.47           O
HETATM 1374  O   HOH A 186      15.767  55.114   9.164  1.00 20.81           O
HETATM 1375  O   HOH A 187       0.754  45.724  20.847  1.00 21.83           O
HETATM 1376  O   HOH A 188      11.338  47.709  21.125  1.00 19.21           O
HETATM 1377  O   HOH A 189      19.695  33.152  21.711  1.00 38.29           O
HETATM 1378  O   HOH A 190      10.491  50.293  22.101  1.00 18.65           O
HETATM 1379  O   HOH A 191       4.344  42.176   9.606  1.00 23.68           O
HETATM 1380  O   HOH A 192      22.154  49.578  16.234  1.00 24.18           O
HETATM 1381  O   HOH A 193      19.540  33.045   0.367  0.50 15.02           O
HETATM 1382  O   HOH A 194      -0.281  31.958  19.977  1.00 27.04           O
HETATM 1383  O   HOH A 195       0.714  43.186  11.724  1.00 31.86           O
HETATM 1384  O   HOH A 196       3.567  51.947  22.586  1.00 23.51           O
HETATM 1385  O   HOH A 197      -1.180  48.644  17.742  1.00 15.53           O
HETATM 1386  O   HOH A 198       1.948  54.836  16.377  1.00 21.93           O
HETATM 1387  O   HOH A 199       2.509  44.733  28.229  1.00 21.51           O
HETATM 1388  O   HOH A 200       2.391  38.378  14.597  1.00 23.11           O
HETATM 1389  O   HOH A 201       1.535  33.854  14.371  1.00 27.83           O
HETATM 1390  O   HOH A 202      -6.330  29.284  19.591  1.00 32.67           O
HETATM 1391  O   HOH A 203      -2.980  54.585   7.752  1.00 25.04           O
HETATM 1392  O   HOH A 204      -1.529  52.313  21.859  1.00 28.95           O
HETATM 1393  O   HOH A 205      -4.219  37.228  18.172  1.00 29.05           O
HETATM 1394  O   HOH A 206      -9.725  38.390  17.936  1.00 25.82           O
HETATM 1395  O   HOH A 207      -0.251  50.231  15.527  1.00 41.88           O
HETATM 1396  O   HOH A 208      -7.548  43.260  22.219  1.00 19.57           O
HETATM 1397  O   HOH A 209      24.427  35.212  14.096  1.00 31.41           O
HETATM 1398  O   HOH A 210      17.913  38.657  31.156  1.00 36.76           O
HETATM 1399  O   HOH A 211      23.606  42.865  12.451  1.00 25.35           O
HETATM 1400  O   HOH A 212     -10.993  36.660  20.173  1.00 34.07           O
HETATM 1401  O   HOH A 213      12.726  35.081  31.604  1.00 41.53           O
HETATM 1402  O   HOH A 214      11.555  34.813   2.535  1.00 32.40           O
HETATM 1403  O   HOH A 215      26.985  45.074  13.907  1.00 33.93           O
HETATM 1404  O   HOH A 216      22.305  31.774  11.176  1.00 27.99           O
HETATM 1405  O   HOH A 217      -5.311  51.156  18.105  1.00 23.06           O
HETATM 1406  O   HOH A 218       4.813  38.197  26.560  1.00 30.39           O
HETATM 1407  O   HOH A 219      16.429  31.381   6.556  1.00 31.04           O
HETATM 1408  O   HOH A 220       2.094  31.244   8.621  1.00 39.67           O
HETATM 1409  O   HOH A 221       3.282  37.464   7.811  1.00 36.79           O
HETATM 1410  O   HOH A 222      12.081  42.928  -4.465  1.00 35.63           O
HETATM 1411  O   HOH A 223      10.457  52.853  26.328  1.00 29.83           O
HETATM 1412  O   HOH A 224     -10.834  37.817  15.537  1.00 32.13           O
HETATM 1413  O   HOH A 225      11.844  46.500  23.826  1.00 35.48           O
HETATM 1414  O   HOH A 226       4.136  40.856   6.472  1.00 40.41           O
HETATM 1415  O   HOH A 227      27.046  47.404  15.617  1.00 32.79           O
HETATM 1416  O   HOH A 228      27.357  54.793   8.601  1.00 45.22           O
HETATM 1417  O   HOH A 229      23.653  47.745  22.742  1.00 50.33           O
HETATM 1418  O   HOH A 230      13.464  46.940  -3.502  1.00 32.04           O
HETATM 1419  O   HOH A 231       3.976  32.286  14.525  1.00 31.56           O
HETATM 1420  O   HOH A 232      20.596  34.827  26.319  1.00 49.91           O
HETATM 1421  O   HOH A 233       0.268  53.053  17.314  1.00 42.31           O
HETATM 1422  O   HOH A 234       8.635  54.732  26.919  1.00 33.52           O
HETATM 1423  O   HOH A 235       8.723  45.643  25.830  1.00 31.33           O
HETATM 1424  O   HOH A 236       9.525  50.219  24.646  1.00 32.19           O
HETATM 1425  O   HOH A 237      29.165  40.813  21.487  1.00 56.55           O
HETATM 1426  O   HOH A 238       2.494  40.055   8.917  1.00 36.15           O
HETATM 1427  O   HOH A 239      -0.777  34.809   9.090  1.00 37.11           O
HETATM 1428  O   HOH A 240      23.455  48.536   9.084  1.00 28.44           O
HETATM 1429  O   HOH A 241      24.600  49.610  11.229  1.00 23.81           O
HETATM 1430  O   HOH A 242      21.015  49.121   7.957  1.00 22.66           O
HETATM 1431  O   HOH A 243      24.242  43.201   9.915  1.00 27.87           O
HETATM 1432  O   HOH A 244      26.271  42.655  13.711  1.00 28.13           O
HETATM 1433  O   HOH A 245       5.960  36.600   6.289  1.00 41.50           O
HETATM 1434  O   HOH A 246      18.246  32.111  23.785  1.00 56.31           O
HETATM 1435  O   HOH A 247       7.627  32.876   7.134  1.00 33.79           O
HETATM 1436  O   HOH A 248       1.165  51.316  -0.972  1.00 46.31           O
HETATM 1437  O   HOH A 249       1.154  41.748  26.990  1.00 28.03           O
HETATM 1438  O   HOH A 250       5.939  26.535  21.149  1.00 44.12           O
HETATM 1439  O   HOH A 251      -5.540  43.055  28.186  1.00 26.47           O
HETATM 1440  O   HOH A 252      11.949  30.493  27.231  1.00 52.78           O
HETATM 1441  O   HOH A 253      13.331  24.495  13.440  1.00 52.00           O
HETATM 1442  O   HOH A 254      -1.435  29.723  20.052  1.00 52.22           O
HETATM 1443  O   HOH A 255       2.115  43.165   6.250  1.00 34.91           O
HETATM 1444  O   HOH A 256      15.992  53.085   0.933  1.00 34.04           O
HETATM 1445  O   HOH A 257       2.292  33.719   9.368  1.00 44.24           O
HETATM 1446  O   HOH A 258       3.626  29.350  26.403  1.00 42.97           O
HETATM 1447  O   HOH A 259      -3.063  50.429  19.160  1.00 33.64           O
HETATM 1448  O   HOH A 260      -4.016  34.810  12.023  1.00 43.84           O
HETATM 1449  O   HOH A 261      -7.025  43.005   7.425  1.00 50.38           O
HETATM 1450  O   HOH A 262      20.664  37.966  29.449  1.00 39.40           O
HETATM 1451  O   HOH A 263       8.219  36.489   7.242  1.00 45.22           O
HETATM 1452  O   HOH A 264       7.603  30.314  26.808  1.00 49.95           O
HETATM 1453  O   HOH A 265      27.192  38.462  19.035  1.00 38.71           O
HETATM 1454  O   HOH A 266      -3.725  41.084  29.707  1.00 44.44           O
HETATM 1455  O   HOH A 267      -3.184  38.130   9.900  1.00 31.98           O
HETATM 1456  O   HOH A 268      -6.968  37.923  11.495  1.00 40.84           O
HETATM 1457  O   HOH A 269      15.958  47.852  -2.923  1.00 44.94           O
HETATM 1458  O   HOH A 270      13.106  37.129  -9.048  1.00 35.12           O
HETATM 1459  O   HOH A 271      23.944  34.080  16.627  1.00 56.20           O
HETATM 1460  O   HOH A 272      17.584  52.420  -1.888  1.00 38.11           O
HETATM 1461  O   HOH A 273       5.521  45.450  -0.564  1.00 47.04           O
HETATM 1462  O   HOH A 274      25.459  49.006   4.506  1.00 54.02           O
HETATM 1463  O   HOH A 275      16.595  27.943  19.928  1.00 49.01           O
HETATM 1464  O   HOH A 276      22.019  36.726  26.870  1.00 56.59           O
HETATM 1465  O   HOH A 277      10.474  36.654   8.412  1.00 45.47           O
HETATM 1466  O   HOH A 278      24.611  33.020  11.703  1.00 47.62           O
HETATM 1467  O   HOH A 279      26.899  46.740   4.551  1.00 52.92           O
HETATM 1468  O   HOH A 280       6.994  29.354  29.500  1.00 59.06           O
HETATM 1469  O   HOH A 281       1.016  41.471   9.152  1.00 46.63           O
HETATM 1470  O   HOH A 282      -3.738  28.411  20.723  1.00 58.68           O
HETATM 1471  O   HOH A 283      14.325  31.639   5.778  1.00 49.26           O
HETATM 1472  O   HOH A 284      27.777  37.727   7.983  1.00 49.01           O
HETATM 1473  O   HOH A 285      11.301  44.913  29.125  1.00 47.44           O
HETATM 1474  O   HOH A 286      -3.722  51.289  27.893  1.00 28.88           O
HETATM 1475  O   HOH A 287       5.173  46.827  -3.289  1.00 59.99           O
HETATM 1476  O   HOH A 288      -9.268  38.914  11.481  1.00 39.62           O
HETATM 1477  O   HOH A 289       9.346  31.934  27.060  1.00 45.27           O
HETATM 1478  O   HOH A 290      24.327  46.350   8.439  1.00 39.23           O
HETATM 1479  O   HOH A 291      17.766  26.234  14.240  1.00 58.80           O
HETATM 1480  O   HOH A 292      14.667  53.897  27.017  1.00 44.45           O
HETATM 1481  O   HOH A 293      13.540  30.919   2.872  1.00 48.68           O
HETATM 1482  O   HOH A 294      13.812  54.578  29.178  1.00 37.62           O
"""
