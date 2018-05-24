import numpy as np

def check_VAR(VAR, log):
    log += "Checking input options..."
    if VAR["LIG1_FRAC"] < 0 or VAR["LIG1_FRAC"] > 1.0:
        sys.exit("LIG1_FRAC must be between 0 and 1.")
    if VAR["MORPHOLOGY"] != "random" and VAR["MORPHOLOGY"] != "janus" and VAR["MORPHOLOGY"] != "stripe":
        sys.exit("Unsupported morphology. So far we support 'random', 'janus', and 'stripe' coatings.")
    if VAR["STRIPES"] < 1:
        sys.exit("The number of stripes must be at least one.")
    return log

def check_mol2(mol2, log):
    MOLECULE = False
    ATOM = False
    BOND = False
    for i in mol2:
        if "@<TRIPOS>MOLECULE" in i:
            MOLECULE = True
        elif "@<TRIPOS>ATOM" in i:
            ATOM = True
        elif "@<TRIPOS>BOND" in i:
            BOND = True

    if not MOLECULE:
        sys.exit("Keyword '@<TRIPOS>MOLECULE' not found in mol2 file.")
    if not ATOM:
        sys.exit("Keyword '@<TRIPOS>ATOM' not found in mol2 file.")
    if not BOND:
        sys.exit("Keyword '@<TRIPOS>BOND' not found in mol2 file.")

    N_lig_file=len(mol2)
    ATOM=False
    atoms = []
    names = []
    res_names = []
    for i in range(N_lig_file):
        if ATOM:
            if "@<TRIPOS>" in mol2[i]:
                break
            atoms.append(mol2[i].split())
            names.append(mol2[i].split()[1])
            res_names.append(mol2[i].split()[7])
        elif "@<TRIPOS>ATOM" in mol2[i]:
            ATOM = True
    log += "{} atoms were found in the mol2 file...\n".format(len(atoms))

    log += "Checking if columns 3, 4, and 5 correspond to floating numbers...\n"
    for i in range(len(atoms)):
        float(atoms[i][2]), float(atoms[i][3]), float(atoms[i][4])

    if len(atoms)!=np.unique(np.array(res_names), return_counts=True)[1][0]:
        sys.exit("There seems to be more than one residue type in the input mol2 file")
    return log

def check_frcmod(fname, log):
    frcmod = np.genfromtxt(fname, delimiter="\n", dtype='str')
    errors = []
    for i in range(len(frcmod)):
        if "ATTN, need revision" in frcmod[i]:
            errors.append(frcmod[i])
    if errors:
        log += "The following parameters in the ligand were impossible to obtain...\n"
        log += "Consider adding you own frcmod file with the missing parameters...\n"
        for i in range(len(errors)):
            log += errors[i]+"\n"
    return log

def read_resname(fname):
    mol2 = np.genfromtxt(fname, delimiter="\n", dtype='str')
    for i in range(len(mol2)):
        if "@<TRIPOS>ATOM" in mol2[i]:
            resname = mol2[i+1].split()[7]
    return resname

def write_leap(VAR, TMP, two_lig_func):
    msj = "source leaprc.gaff \n\n"

    msj += "loadamberparams {}/LIG1.frcmod\n".format(TMP)
    msj += "{} = loadmol3 {}/LIG1.mol2\n".format(read_resname(TMP+"/LIG1.mol2"), TMP)
    msj += "check {}\n".format(read_resname(TMP+"/LIG1.mol2"))
    msj += "saveoff {} {}/LIG1.lib\n\n".format(read_resname(TMP+"/LIG1.mol2"), TMP)
    if two_lig_func:
        msj += "loadamberparams {}/LIG2.frcmod\n".format(TMP)
        msj += "{} = loadmol3 {}/LIG2.mol2\n".format(read_resname(TMP+"/LIG2.mol2"), TMP)
        msj += "check {}\n".format(read_resname(TMP+"/LIG2.mol2"))
        msj += "saveoff {} {}/LIG2.lib\n\n".format(read_resname(TMP+"/LIG2.mol2"), TMP)

    msj += "loadamberparams PARAMS/AU.frcmod\n"
    msj += "loadamberparams PARAMS/AUS.frcmod\n"
    msj += "loadamberparams PARAMS/AUL.frcmod\n"
    msj += "loadamberparams PARAMS/PARAMS.frcmod\n\n"
    if VAR["FRCMOD"] != "0":
        msj += "loadamberparams {}\n".format(VAR["FRCMOD"])
    msj += "AU = loadmol3 PARAMS/AU.mol2\n"
    msj += "AUS = loadmol3 PARAMS/AUS.mol2\n"
    msj += "AUL = loadmol3 PARAMS/AUL.mol2\n"

    msj += "loadoff {}/LIG1.lib\n".format(TMP)
    if two_lig_func:
        msj += "loadoff {}/LIG2.lib\n".format(TMP)

    msj += "{} = loadpdb {}/{}.pdb\n".format(VAR["NAME"], TMP, VAR["NAME"])
    msj += "saveamberparm {} {}/{}.prmtop {}/{}.inpcrd\n".format(VAR["NAME"], TMP, VAR["NAME"], TMP, VAR["NAME"])
    msj += "quit"
    out = open(TMP+"/TLeap.in", "w")
    out.write(msj)
    out.close()
