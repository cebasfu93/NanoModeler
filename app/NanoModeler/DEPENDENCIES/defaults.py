import numpy as np
import logging

logger = logging.getLogger('nanomodeler')
logger.addHandler(logging.NullHandler())

def check_VAR(VAR, rep):
    if VAR["LIG1_FRAC"] < 0 or VAR["LIG1_FRAC"] > 1.0:
        excp_txt = "ATTENTION! LIG1_FRAC must be between 0 and 1."
        rep.write(excp_txt + "\n")
        raise Exception(excp_txt)
    if VAR["MORPHOLOGY"] != "random" and VAR["MORPHOLOGY"] != "janus" and VAR["MORPHOLOGY"] != "stripe":
        excp_txt = "ATTENTION! Unsupported morphology. So far we support 'random', 'janus', and 'stripe' coatings."
        rep.write(excp_txt + "\n")
        raise Exception(excp_txt)
    if VAR["STRIPES"] < 1:
        excp_txt = "ATTENTION! The number of stripes must be at least one."
        rep.write(excp_txt + "\n")
        raise Exception(excp_txt)

def check_mol2(mol2, rep):
    MOLECULE = False
    ATOM = False
    BOND = False
    for i in mol2:  #Looks for needed keywords in the mol2 file
        if "@<TRIPOS>MOLECULE" in i:
            MOLECULE = True
        elif "@<TRIPOS>ATOM" in i:
            ATOM = True
        elif "@<TRIPOS>BOND" in i:
            BOND = True

    if not MOLECULE:
        excp_txt = "ATTENTION! Keyword '@<TRIPOS>MOLECULE' not found in mol2 file."
        rep.write(excp_txt + "\n")
        raise Exception(excp_txt)
    if not ATOM:
        excp_txt = "ATTENTION! Keyword '@<TRIPOS>ATOM' not found in mol2 file."
        rep.write(excp_txt + "\n")
        raise Exception(excp_txt)
    if not BOND:
        excp_txt = "ATTENTION! Keyword '@<TRIPOS>BOND' not found in mol2 file."
        rep.write(excp_txt + "\n")
        raise Exception(excp_txt)

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
    logger.info("\t{} atoms were found in the mol2 file...".format(len(atoms)))

    logger.info("\tChecking if columns 3, 4, and 5 correspond to floating numbers...")
    for i in range(len(atoms)):
        float(atoms[i][2]), float(atoms[i][3]), float(atoms[i][4])

    logger.info("\tChecking if there is only one residue in the input structure...")
    if len(atoms)!=np.unique(np.array(res_names), return_counts=True)[1][0]:
        excp_txt = "ATTENTION! There seems to be more than one residue type in the input mol2 file"
        rep.write(excp_txt + "\n")
        raise Exception(excp_txt)
    return True

def check_frcmod(fname, rep):
    frcmod = np.genfromtxt(fname, delimiter="\n", dtype='str')
    errors = []
    for i in range(len(frcmod)):
        if "ATTN, need revision" in frcmod[i]:
            errors.append(frcmod[i])
    if errors:
        warn_txt = "\tATTENTION! The following parameters in the ligand were impossible to obtain..."
        rep.write(warn_txt + "\n")
        logger.warning(warn_txt)
        warn_txt = "\tConsider adding you own frcmod file with the missing parameters..."
        rep.write(warn_txt + "\n")
        logger.info(warn_txt)
        for i in range(len(errors)):
            err_txt = "\t{}".format(errors[i])
            rep.write(err_txt + "\n")
            logger.info(err_txt)
    return True

def read_resname(fname):
    mol2 = np.genfromtxt(fname, delimiter="\n", dtype='str')
    for i in range(len(mol2)):
        if "@<TRIPOS>ATOM" in mol2[i]:
            resname = mol2[i+1].split()[7]
    return resname

def write_leap(VAR, TMP, two_lig_func):
    msj = "source leaprc.gaff2 \n\n"

    #Loads parameters of the ligands
    msj += "loadamberparams {}/LIG1.frcmod\n".format(TMP)
    if two_lig_func:
        msj += "loadamberparams {}/LIG2.frcmod\n".format(TMP)
    #Loads parameters corrections of Hakkinen
    msj += "loadamberparams PARAMS/PARAMS.frcmod\n\n"
    if VAR["FRCMOD"]:
        msj += "loadamberparams {}\n".format(VAR["FRCMOD"].name)     #If the user gave an frcmod in overwrites all the previous ones

    #Loads structures of the gold atoms
    msj += "AU = loadmol3 PARAMS/AU.mol2\n"
    msj += "AUS = loadmol3 PARAMS/AUS.mol2\n"
    msj += "AUL = loadmol3 PARAMS/AUL.mol2\n"

    #Loads structures of the ligands
    msj += "{} = loadmol3 {}/LIG1.mol2\n".format(read_resname(TMP+"/LIG1.mol2"), TMP)
    msj += "check {}\n".format(read_resname(TMP+"/LIG1.mol2"))
    msj += "saveoff {} {}/LIG1.lib\n\n".format(read_resname(TMP+"/LIG1.mol2"), TMP)
    if two_lig_func:
        msj += "{} = loadmol3 {}/LIG2.mol2\n".format(read_resname(TMP+"/LIG2.mol2"), TMP)
        msj += "check {}\n".format(read_resname(TMP+"/LIG2.mol2"))
        msj += "saveoff {} {}/LIG2.lib\n\n".format(read_resname(TMP+"/LIG2.mol2"), TMP)

    #Loads libraries of the ligands
    msj += "loadoff {}/LIG1.lib\n".format(TMP)
    if two_lig_func:
        msj += "loadoff {}/LIG2.lib\n".format(TMP)

    #Loads pdb and saves parameters to amber files
    msj += "NP = loadpdb {}/NP.pdb\n".format(TMP)
    msj += "saveamberparm NP {}/NP.prmtop {}/NP.inpcrd\n".format(TMP, TMP)
    msj += "quit"
    out = open(TMP+"/TLeap.in", "w")
    out.write(msj)
    out.close()
