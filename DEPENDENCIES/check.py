import numpy as np
import sys

def check_VAR(VAR):
    print("Checking presence of second ligand...")
    if VAR["LIG1_FRAC"] < 0 or VAR["LIG1_FRAC"] > 1.0:
        sys.exit("LIG1_FRAC must be between 0 and 1")
    if VAR["MORPHOLOGY"] != "random" and VAR["MORPHOLOGY"] != "janus" and VAR["MORPHOLOGY"] != "stripe" and VAR["LIG1_FRAC"] >= 0 and VAR["LIG1_FRAC"] <= 1.0:
        sys.exit("Unsupported morphology. So far we support 'random', 'janus', and 'stripe' coatings")
    if VAR["STRIPES"] < 1:
        sys.exit("The number of stripes must be at least one")
    if VAR["FIRST"] != 1 and VAR["FIRST"] != 2:
        sys.exit("The first ligand to place in the striped morphology bust be either 1 (for ligand 1) or 2 (for ligand 2)")

def check_mol2(fname):
    mol2 = np.genfromtxt(fname, delimiter='\n', dtype='str')

    found_ATOM = False
    found_BOND = False
    found_CONNECT = False
    for i in mol2:
        if "@<TRIPOS>ATOM" in i:
            found_ATOM = True
        elif "@<TRIPOS>BOND" in i:
            found_BOND = True
        elif "@<TRIPOS>RESIDUECONNECT" in i:
            found_CONNECT = True

    if not found_ATOM:
        sys.exit("Keyword '@<TRIPOS>ATOM' not found in mol2 file.")
    if not found_BOND:
        sys.exit("Keyword '@<TRIPOS>BOND' not found in mol2 file.")
    if not found_CONNECT:
        sys.exit("Keyword '@<TRIPOS>RESIDUECONNECT' not found in mol2 file.")

    N_lig_file=len(mol2)
    found_ATOM=False
    atoms = []
    names = []
    types = []
    for i in range(N_lig_file):
        if found_ATOM:
            if "@<TRIPOS>" in mol2[i]:
                break
            atoms.append(mol2[i].split())
            names.append(mol2[i].split()[1])
        elif "@<TRIPOS>ATOM" in mol2[i]:
            found_ATOM = True
    print("{} atoms were found in the mol2 file...".format(len(atoms)))

    print("Checking if columns 3, 4, and 5 correspond to floating numbers...")
    for i in range(len(atoms)):
        float(atoms[i][2]), float(atoms[i][3]), float(atoms[i][4])

    for i in range(N_lig_file):
        if "@<TRIPOS>RESIDUECONNECT" in mol2[i]:
            connect = mol2[i+1].split()[1]
            print("The name found for the connecting atom in the mol2 file is '{}'...".format(connect))
    names = np.array(names)
    ndx_con = np.where(names==connect)[0][0]
    print("The connecting atom was identified to be atom {}...".format(ndx_con+1))
