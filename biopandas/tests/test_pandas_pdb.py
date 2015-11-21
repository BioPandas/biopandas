
def read_pdb():
    ppdb = PandasPDB()
    ppdb.read_pdb('/Users/Sebastian/Desktop/3EIY.pdb')
    ppdb.to_pdb('/Users/Sebastian/Desktop/3EIY.txt')
