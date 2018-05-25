NanoModeler v0.3

Known bugs:
-None. NanoModeler is perfect.
-But ligands must consist of only one residue due to the saveoff .lib in tleap.
-When the mol2 is rewriten, the third term in the SUBSTRUCTURE section (the root atom) is always set to 1, i.e. there must be only one residue also because of this.
-There can't be an atom named ST in the input mol2 files
-log file is currently not part of the output, is part of the return of the function

Perspectives:
-Make a broader check of possible errors.
-Csearch must include the S atom.
-Clean log file.
-Comment code.
