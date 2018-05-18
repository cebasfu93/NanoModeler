import numpy as np

def check_VAR(VAR, log):
    log += "Checking input options..."
    if float(VAR["LIG1_FRAC"]) < 0 or float(VAR["LIG1_FRAC"]) > 1.0:
        sys.exit("LIG1_FRAC must be between 0 and 1.")
    if VAR["MORPHOLOGY"] != "random" and VAR["MORPHOLOGY"] != "janus" and VAR["MORPHOLOGY"] != "stripe" and float(VAR["LIG1_FRAC"]) >= 0 and float(VAR["LIG1_FRAC"]) <= 1.0:
        sys.exit("Unsupported morphology. So far we support 'random', 'janus', and 'stripe' coatings.")
    if int(VAR["STRIPES"]) < 1:
        sys.exit("The number of stripes must be at least one.")

def check_mol2(fname, log):
    mol2 = np.genfromtxt(fname, delimiter='\n', dtype='str')

    found_MOLECULE = False
    found_ATOM = False
    found_BOND = False
    found_CONNECT = False
    for i in mol2:
        if "@<TRIPOS>MOLECULE" in i:
            found_MOLECULE = True
        elif "@<TRIPOS>ATOM" in i:
            found_ATOM = True
        elif "@<TRIPOS>BOND" in i:
            found_BOND = True
        elif "@<TRIPOS>RESIDUECONNECT" in i:
            found_CONNECT = True

    if not found_MOLECULE:
        sys.exit("Keyword '@<TRIPOS>MOLECULE' not found in mol2 file.")
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
    for i in range(N_lig_file):
        if found_ATOM:
            if "@<TRIPOS>" in mol2[i]:
                break
            atoms.append(mol2[i].split())
            names.append(mol2[i].split()[1])
        elif "@<TRIPOS>ATOM" in mol2[i]:
            found_ATOM = True
    log += "{} atoms were found in the mol2 file...".format(len(atoms))

    log += "Checking if columns 3, 4, and 5 correspond to floating numbers..."
    for i in range(len(atoms)):
        float(atoms[i][2]), float(atoms[i][3]), float(atoms[i][4])

    for i in range(N_lig_file):
        if "@<TRIPOS>RESIDUECONNECT" in mol2[i]:
            connect = mol2[i+1].split()[1]
            log += "The name found for the connecting atom in the mol2 file is '{}'...".format(connect)
    names = np.array(names)
    ndx_con = np.where(names==connect)[0][0]
    log += "The connecting atom was identified to be atom {}...".format(ndx_con+1)

def read_resname(lig_fname):
    mol2 = np.genfromtxt(lig_fname, delimiter="\n", dtype='str')
    for i in range(len(mol2)):
        if "@<TRIPOS>ATOM" in mol2[i]:
            resname = mol2[i+1].split()[7]
    return resname

def write_leap(VAR, fname, two_lig_func):
    msj = "source leaprc.gaff \n\n"

    msj += "loadamberparams " + "TMP/"+VAR["LIG1_FILE"][:-5]+".frcmod\n"
    msj += "loadamberparams " + VAR["DEPENDS"]+"/PARAMS.frcmod\n\n"

    msj += read_resname(VAR["LIG1_FILE"]) + " = loadmol3 " + "TMP/"+VAR["LIG1_FILE"]+"\n"
    msj += "check " + read_resname(VAR["LIG1_FILE"]) + "\n"
    msj += "saveoff " + read_resname(VAR["LIG1_FILE"]) + " " + "TMP/"+VAR["LIG1_FILE"][:-5]+".lib\n\n"
    if two_lig_func:
        msj += "loadamberparams " + "TMP/"+VAR["LIG2_FILE"][:-5]+".frcmod\n"
        msj += read_resname(VAR["LIG2_FILE"]) + " = loadmol3 " + "TMP/"+VAR["LIG2_FILE"]+"\n"
        msj += "check " + read_resname(VAR["LIG2_FILE"]) + "\n"
        msj += "saveoff " + read_resname(VAR["LIG2_FILE"]) + " " + "TMP/"+VAR["LIG2_FILE"][:-5]+".lib\n\n"

    msj += "loadamberparams " + VAR["DEPENDS"]+"/AU.frcmod\n"
    msj += "loadamberparams " + VAR["DEPENDS"]+"/ST.frcmod\n"
    msj += "AU = loadmol3 " + VAR["DEPENDS"]+"/AU.mol2\n"
    msj += "ST = loadmol3 " + VAR["DEPENDS"]+"/ST.mol2\n\n"

    msj += "loadoff " + "TMP/"+VAR["LIG1_FILE"][:-5]+".lib\n"
    if two_lig_func:
        msj += "loadoff " + "TMP/"+VAR["LIG2_FILE"][:-5]+".lib\n"

    msj += VAR["NAME"] + " = loadpdb " + "TMP/"+VAR["NAME"]+".pdb \n"
    msj += "saveamberparm " +  VAR["NAME"] + " " + "TMP/"+VAR["NAME"]+".prmtop" + " " + "TMP/"+VAR["NAME"]+".inpcrd \n"
    msj += "quit"
    out = open(fname, "w")
    out.write(msj)
    out.close()
