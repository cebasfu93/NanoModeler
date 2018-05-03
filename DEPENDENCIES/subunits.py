import numpy as np

class Block:
    #A staple is as defined in paper 333
    def __init__(self, ndx_S=0, ndx_Au=0, ndx_C=0, ndx_H=0, types_Au=[]):
        self.S = np.array(ndx_S, dtype='int')
        self.Au = np.array(ndx_Au, dtype='int')
        self.C = np.array(ndx_C, dtype='int')
        self.H = np.array(ndx_H, dtype='int')
        self.typesAu = np.array(types_Au, dtype='str')
