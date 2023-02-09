import numpy as np
import logging
import tempfile
import io

logger = logging.getLogger("nanomodeler")
logger.addHandler(logging.NullHandler())

report = logging.getLogger("nanomodeler.report")


def check_VAR(VAR):
    if VAR["LIG1_FRAC"] < 0 or VAR["LIG1_FRAC"] > 1.0:
        excp_txt = "ATTENTION! LIG1_FRAC must be between 0 and 1."
        report.error(excp_txt)
        raise Exception(excp_txt)
    if (
        VAR["MORPHOLOGY"] != "random"
        and VAR["MORPHOLOGY"] != "janus"
        and VAR["MORPHOLOGY"] != "stripe"
    ):
        excp_txt = "ATTENTION! Unsupported morphology. So far we support 'random', 'janus', and 'stripe' coatings."
        report.error(excp_txt)
        raise Exception(excp_txt)
    if VAR["STRIPES"] < 1:
        excp_txt = "ATTENTION! The number of stripes must be at least one."
        report.error(excp_txt)
        raise Exception(excp_txt)


def check_mol2(mol2):
    MOLECULE = False
    ATOM = False
    BOND = False

    for i in mol2:  # Looks for needed keywords in the mol2 file
        if "@<TRIPOS>MOLECULE" in i:
            MOLECULE = True
        elif "@<TRIPOS>ATOM" in i:
            ATOM = True
        elif "@<TRIPOS>BOND" in i:
            BOND = True

    if not MOLECULE:
        message = "Keyword '@<TRIPOS>MOLECULE' not found in mol2 file."
        report.error(message)
        raise Exception(message)
    if not ATOM:
        message = "Keyword '@<TRIPOS>ATOM' not found in mol2 file."
        report.error(message)
        raise Exception(message)
    if not BOND:
        message = "Keyword '@<TRIPOS>BOND' not found in mol2 file."
        report.error(message)
        raise Exception(message)

    N_lig_file = len(mol2)

    ATOM = False
    atoms = []
    names = []
    res_names = []

    for i in range(N_lig_file):
        if ATOM:
            if "@<TRIPOS>" in mol2[i]:
                break
            if mol2[i].split() != []:
                atoms.append(mol2[i].split())
                names.append(mol2[i].split()[1])
                res_names.append(mol2[i].split()[7])
        elif "@<TRIPOS>ATOM" in mol2[i]:
            ATOM = True

    report.info("\t{} atoms were found in the mol2 file...".format(len(atoms)))
    report.info("\tChecking if columns 3, 4, and 5 correspond to floating numbers...")

    for i in range(len(atoms)):
        float(atoms[i][2]), float(atoms[i][3]), float(atoms[i][4])

    report.info("\tChecking if there is only one residue in the input structure...")

    if len(atoms) != np.unique(np.array(res_names), return_counts=True)[1][0]:
        message = "There seems to be more than one residue type in the input mol2 file"
        report.error(message)
        raise Exception(message)

    return True


def check_frcmod(fname):
    frcmod = np.genfromtxt(fname, delimiter="\n", dtype="str")
    errors = []
    for i in range(len(frcmod)):
        if "ATTN, need revision" in frcmod[i]:
            errors.append(frcmod[i])
    if errors:
        warn_txt = "\tATTENTION! The following parameters in the ligand were impossible to obtain..."
        report.warning(warn_txt + "\n")
        logger.warning(warn_txt)
        warn_txt = (
            "\tConsider adding you own frcmod file with the missing parameters..."
        )
        report.info(warn_txt + "\n")
        logger.info(warn_txt)
        for i in range(len(errors)):
            err_txt = "\t{}".format(errors[i])
            report.info(err_txt + "\n")
            logger.info(err_txt)
    return True


def read_resname(fname):
    mol2 = np.genfromtxt(fname, delimiter="\n", dtype="str")
    for i in range(len(mol2)):
        if "@<TRIPOS>ATOM" in mol2[i]:
            resname = mol2[i + 1].split()[7]
    return resname


def write_leap(VAR, TMP, two_lig_func):
    msj = "source leaprc.gaff2 \n\n"

    # Loads parameters of the ligands
    msj += "loadamberparams {}/LIG1.frcmod\n".format(TMP)
    if two_lig_func:
        msj += "loadamberparams {}/LIG2.frcmod\n".format(TMP)
    # Loads parameters corrections of Hakkinen
    msj += "loadamberparams PARAMS/PARAMS.frcmod\n\n"

    frcmod_file = VAR["FRCMOD"]
    # If the user gave an frcmod in overwrites all the previous ones
    if frcmod_file:
        # if frcmod_file is an instance of StringIO I need to create a temporary file
        if isinstance(frcmod_file, io.StringIO):
            frcmod_tmp_file = tempfile.NamedTemporaryFile(delete=False, dir=TMP)
            report.info(
                "creating frcmod temporary file: '{}'".format(frcmod_tmp_file.name)
            )
            frcmod_tmp_file.write(frcmod_file.read().encode("utf-8"))
            frcmod_tmp_file.flush()
            frcmod_file.close()

            frcmod_file = frcmod_tmp_file

        msj += "loadamberparams {}\n".format(frcmod_file.name)

    # Loads structures of the gold atoms
    msj += "AU = loadmol3 PARAMS/AU.mol2\n"
    msj += "AUS = loadmol3 PARAMS/AUS.mol2\n"
    msj += "AUL = loadmol3 PARAMS/AUL.mol2\n"

    # Loads structures of the ligands
    msj += "{} = loadmol3 {}/LIG1.mol2\n".format(read_resname(TMP + "/LIG1.mol2"), TMP)
    msj += "check {}\n".format(read_resname(TMP + "/LIG1.mol2"))
    msj += "saveoff {} {}/LIG1.lib\n\n".format(read_resname(TMP + "/LIG1.mol2"), TMP)
    if two_lig_func:
        msj += "{} = loadmol3 {}/LIG2.mol2\n".format(
            read_resname(TMP + "/LIG2.mol2"), TMP
        )
        msj += "check {}\n".format(read_resname(TMP + "/LIG2.mol2"))
        msj += "saveoff {} {}/LIG2.lib\n\n".format(
            read_resname(TMP + "/LIG2.mol2"), TMP
        )

    # Loads libraries of the ligands
    msj += "loadoff {}/LIG1.lib\n".format(TMP)
    if two_lig_func:
        msj += "loadoff {}/LIG2.lib\n".format(TMP)

    # Loads pdb and saves parameters to amber files
    msj += "NP = loadpdb {}/NP.pdb\n".format(TMP)
    msj += "saveamberparm NP {}/NP.prmtop {}/NP.inpcrd\n".format(TMP, TMP)
    msj += "quit"

    out = open(TMP + "/TLeap.in", "w")
    out.write(msj)
    out.close()
