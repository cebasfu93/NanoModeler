import numpy as np

class Block:
    #A block is a sulphur atom, it's two neighbouring Au atoms, and the C atom
    def __init__(self, ndx_S=0, ndx_Au=0, ndx_C=0):
        self.S = ndx_S
        self.Au = ndx_Au
        self.C = ndx_C

class Gancho:
    def __init__(self, ndx_S=0, ndx_Au=0, ndx_C=0):
        self.S = np.array(ndx_S)
        self.Au = np.array(ndx_Au)
        self.C = np.array(ndx_C)

    def add(self, more_S, more_Au, more_C):
        self.S = np.append(self.S, more_S)
        self.Au = np.append(self.Au, more_Au)
        self.C = np.append(self.C, more_C)

class Staple:
    #A staple is as defined in paper 333
    def __init__(self, ndx_S=0, ndx_Au=0, ndx_C=0):
        #self.tipo = tipo
        self.S = np.array(ndx_S, dtype='int')
        #self.Au_l = np.array(ndx_Au_l, dtype='int')
        self.Au = np.array(ndx_Au, dtype='int')
        #self.Au_s = np.array(list(set(list(self.Au))-set(list(self.Au_l))), dtype='int')
        self.C = np.array(ndx_C, dtype='int')
        only = np.unique(self.Au, return_counts=True)
        self.Au_l = only[0][np.where(only[1]==2)[0]]
        self.Au_s = only[0][np.where(only[1]==1)[0]]
        self.Au = only[0]
        self.tipo = "UNK"

    def change_tipo(self, new_tipo):
        self.tipo = new_tipo

class Residue:
    #A residue is only defined by the residue type and atomnumbers that conform the residue
    def __init__(self, restype='UNK', ndx=[]):
        self.restype = restype
        self.ndx = ndx
