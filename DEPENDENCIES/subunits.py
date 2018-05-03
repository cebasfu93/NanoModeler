import numpy as np

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
