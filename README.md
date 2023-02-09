NanoModeler v1.0.4

Known bugs:
-None. NanoModeler is perfect.

-But ligands must consist of only one residue due to the saveoff .lib in tleap.
-When the mol2 is rewriten, the third term in the SUBSTRUCTURE section (the root atom) is always set to 1, i.e. there must be only one residue also because of this.
-There can't be an atom named ST in the input mol2 files
-In line approx. 91 in rewrite_mol2.py the charge of the cappint groups should be distributed in N atoms, not N-1:
---> charge_per_atom = np.sum(charge_cap)/N_at
-Ligands must consist of only one residue due to the saveoff .lib in tleap.
-When the mol2 is rewriten, the third term in the SUBSTRUCTURE section (the root atom) is always set to 1. It doesn't seem to bring any problem so far.

Perspectives:
