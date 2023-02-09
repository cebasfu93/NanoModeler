# -*- coding: utf-8 -*-

import logging

logger = logging.getLogger("nanomodeler")
logger.addHandler(logging.NullHandler())


__VERSION__ = "1.0.9"


def getVersion():
    return __VERSION__


def NanoModeler(
    LIG1_FILE=None,
    CAP1=[],
    LIG1_C=0,
    LIG1_S=0,
    LIG1_FRAC=1.0,
    MORPHOLOGY="random",
    RSEED=None,
    STRIPES=1,
    LIG2_FILE=None,
    CAP2=[],
    LIG2_C=0,
    LIG2_S=0,
    FRCMOD=None,
    CORE=None,
    ELONGATED=False,
):
    version = __VERSION__

    VAR = {
        # ByteArray of the mol2 of ligand1
        "LIG1_FILE": LIG1_FILE,
        # List of atom numbers (as in the mol2, likely to start in 1) to remove
        "CAP1": CAP1,
        # Signed integer, atom number of carbon atom in LIG1_FILE used as anchor
        "LIG1_C": LIG1_C,
        # Signed integer, atom number in the original mol2 file of ligand 1 corresponding to the anchoring S atom
        "LIG1_S": LIG1_S,
        # Float, fraction of ligand1 to place (0-1.0)
        "LIG1_FRAC": LIG1_FRAC,
        # String, morphology to distribute ligands1 and 2. random, janus, and stripe are allowed
        "MORPHOLOGY": MORPHOLOGY,
        # Integer, random seed for random morphology
        "RSEED": RSEED,
        # Signed integer, number of stripes for stripe morphology. It will start (bottom up with ligand 1)
        "STRIPES": STRIPES,
        # ByteArray of the mol2 of ligand2
        "LIG2_FILE": LIG2_FILE,
        # List of atom numbers (as in the mol2, likely to start in 1) to remove
        "CAP2": CAP2,
        # Signed integer, atom number of carbon atom in LIG2_FILE used as anchor
        "LIG2_C": LIG2_C,
        # Signed integer, atom number in the original mol2 file of ligand 1 corresponding to the anchoring S atom
        "LIG2_S": LIG2_S,
        # ByteArray of the frcmod provided by the user
        "FRCMOD": FRCMOD,
        # ByteArray of the pdb of the core
        "CORE": CORE,
        # If set to true, the first carbon of the cores will be ignored and placed 1.8A away from the S atoms along the O-S axis
        "ELONGATED": ELONGATED,
    }

    logger.info("WELCOME TO NANOMODELER v{}".format(__VERSION__))
    logger.info("Importing sys library...")
    import sys

    logger.info("Importing numpy library...")
    import numpy as np

    logger.info("Importing random library...")
    import random

    logger.info("Importing math library...")
    import math

    logger.info("Importing scipy library...")
    from scipy.spatial import distance

    logger.info("Importing sklearn library...")
    from sklearn.decomposition import PCA

    logger.info("Importing os library...")
    import os

    logger.info("Importing shutil library...")
    import shutil

    logger.info("Importing tempfile...")
    import tempfile

    logger.info("Importing zipfile...")
    import zipfile

    logger.info("Importing atexit library...")
    import atexit

    logger.info("Importing re library...")
    import re

    logger.info("Importing pickle library...")
    import pickle

    logger.info("Importing subprocess library...")
    import subprocess

    logger.info("Importing signal library...")
    import signal

    logger.info("Importing traceback library...")
    import traceback

    logger.info("Importing datetime library...")
    import datetime

    logger.info("Importing pathlib library...")
    import pathlib

    logger.info("Importing acpype.py...")
    from DEPENDENCIES.acpype import MolTopol

    logger.info("Importing NP_builder dependency...")
    from DEPENDENCIES.NP_builder import (
        assign_morph,
        coat_NP,
        get_ligand_pill,
        get_stones,
        init_core_pdb,
        init_lig_mol2,
        print_NP_pdb,
    )

    logger.info("Importing staples dependency...")
    from DEPENDENCIES.staples import (
        get_ndxs,
        load_gro,
        load_top,
        make_blocks,
        write_angles,
        write_bonds,
        write_topology,
    )

    logger.info("Importing default function dependency...")
    from DEPENDENCIES.defaults import check_frcmod, check_mol2, check_VAR, write_leap

    logger.info("Importing transformations dependency...")
    from DEPENDENCIES.transformations import (
        affine_matrix_from_points,
        quaternion_matrix,
        vector_norm,
    )

    logger.info("Importing subunits dependency...")
    import DEPENDENCIES.subunits

    logger.info("Importing rewrite_mol2 dependency...")
    from DEPENDENCIES.rewrite_mol2 import rewrite_mol2

    logger.info("Importing cleanup dependency...")
    from DEPENDENCIES.cleanup import cleanup_error

    logger.info("Creating temporary folder...")
    TMP = tempfile.mkdtemp(dir="./")

    # get singleton instance of logger (report level)
    report = logging.getLogger("nanomodeler.report")
    report.setLevel(logging.INFO)

    # clean all handler
    report.handlers = []

    # create a file handler
    reportFileHandler = logging.FileHandler(os.path.join(TMP, "report.log"), "w")
    reportFileHandler.setLevel(logging.INFO)

    # set the formatter
    formatter = logging.Formatter("%(message)s")
    reportFileHandler.setFormatter(formatter)

    # add the filehandler to the logger
    report.addHandler(reportFileHandler)

    # From this point, the cleanup_error will be excecuted after NanoModeler finishes (or crashes)
    atexit.register(cleanup_error, reportFileHandler)

    # Checks validity of values in the variables
    logger.info("Checking input options...")
    check_VAR(VAR)

    # Based on the morphology, number of stripes, and fraction of ligand1, it is decided if there are one or two ligands
    if VAR["MORPHOLOGY"] == "stripe":
        if VAR["STRIPES"] == 1:
            two_lig = False
        else:
            two_lig = True
    else:
        two_lig = VAR["LIG1_FRAC"] < 1.0

    logger.info("Imported options:")
    report.info("Imported options:")
    for i in VAR:
        if (
            i == "LIG1_FILE" or i == "LIG2_FILE" or i == "CORE" or i == "FRCMOD"
        ) and VAR[i]:
            if VAR[i]:
                op_txt = "\t{:<20}{:>20}\n".format(i, str(VAR[i].name))
                logger.info(op_txt)
                report.info(op_txt)
        else:
            op_txt = "\t{:<20}{:>20}".format(i, str(VAR[i]))
            logger.info(op_txt)
            report.info(op_txt)

    logger.info("One ligand was found...")
    logger.info("Checking ligand1 mol2 file...")

    # Read lines and eliminates the \n character from them
    LIG1_MOL2 = VAR["LIG1_FILE"].readlines()
    LIG1_MOL2 = [s.replace("\n", "") for s in LIG1_MOL2]
    check_mol2(LIG1_MOL2)

    # Writes new mol2 without capping, S atom at the end, and properly numbered
    logger.info("Rewriting ligand1 file...")
    VAR["LIG1_C"] = rewrite_mol2(
        LIG1_MOL2,
        VAR["CAP1"],
        VAR["LIG1_S"],
        VAR["LIG1_C"],
        TMP + "/LIG1.mol2",
        VAR["ELONGATED"],
    )

    if two_lig:
        logger.info("Two ligands were found...")
        logger.info("Checking ligand2 mol2 file...")

        # Read lines and eliminates the \n character from them
        LIG2_MOL2 = VAR["LIG2_FILE"].readlines()
        LIG2_MOL2 = [s.replace("\n", "") for s in LIG2_MOL2]

        check_mol2(LIG2_MOL2)

        # Writes new mol2 without capping, S atom at the end, and properly numbered
        logger.info("Rewriting ligand2 file...")
        VAR["LIG2_C"] = rewrite_mol2(
            LIG2_MOL2,
            VAR["CAP2"],
            VAR["LIG2_S"],
            VAR["LIG2_C"],
            TMP + "/LIG2.mol2",
            VAR["ELONGATED"],
        )

    ##############################NP_builder########################

    logger.info("Initializing core...")
    # Read lines and eliminates the \n character from them
    CORE_PDB = VAR["CORE"].readlines()
    CORE_PDB = [s.replace("\n", "") for s in CORE_PDB]
    # Extracts xyz coordinates, names, and residues of the core
    xyz_core, names_core, res_core = init_core_pdb(CORE_PDB, VAR["ELONGATED"])

    logger.info("Initializing ligand1...")
    # Reads rewritten mol2 file and puts the C atom in the origin
    xyz_lig1, names_lig1, res_lig1 = init_lig_mol2(
        TMP + "/LIG1.mol2", VAR["LIG1_S"], VAR["LIG1_C"]
    )
    if two_lig:
        logger.info("Initializing ligand2...")
        # Reads rewritten mol2 file and puts the C atom in the origin
        xyz_lig2, names_lig2, res_lig2 = init_lig_mol2(
            TMP + "/LIG2.mol2", VAR["LIG1_S"], VAR["LIG2_C"]
        )
        if res_lig1[0] == res_lig2[0]:
            error_msg = f"The two ligands provided have the same name ({res_lig1[0]}). Rename one of them and try again."
            logger.error(error_msg)
            raise ValueError(error_msg)
    else:
        xyz_lig2, names_lig2, res_lig2 = [], [], []

    logger.info("Running PCA for ligand1...")

    #  Gives the coordinates of the C atom and the projections of 2 pseudo random atoms into PCA1

    xyz_pillars1 = get_ligand_pill(xyz_lig1, VAR["LIG1_C"], VAR["LIG1_S"])
    if two_lig:
        logger.info("Running PCA for ligand2...")
        #  Gives the coordinates of the C atom and the projections of 2 pseudo random atoms into PCA1

        xyz_pillars2 = get_ligand_pill(xyz_lig2, VAR["LIG2_C"], VAR["LIG2_S"])
    else:
        xyz_pillars2 = []

    N_S = len(names_core[names_core == "ST"])

    logger.info("Assigning morphology...")
    # Gives the coordinates of the C atoms in the core corresponding to each ligand depending on the morphology specified
    xyz_anchors1, xyz_anchors2 = assign_morph(
        xyz_core,
        names_core,
        VAR["LIG1_FRAC"],
        VAR["RSEED"],
        VAR["MORPHOLOGY"],
        VAR["STRIPES"],
    )

    # Gives a 3D array with the xyz coordinates for all the stones of all the anchors of ligand 1
    xyz_stones1 = get_stones(
        xyz_core, names_core, xyz_anchors1, xyz_pillars1, VAR["LIG1_S"]
    )
    if two_lig:
        # Gives a 3D array with the xyz coordinates for all the stones of all the anchors of ligand 2
        xyz_stones2 = get_stones(
            xyz_core, names_core, xyz_anchors2, xyz_pillars2, VAR["LIG2_S"]
        )
    else:
        xyz_stones2 = []

    logger.info("Coating nanoparticle...")
    # Merges xyz coordinates and names of the core and the ligands into one coated NP
    xyz_coated_NP, names_coated_NP, res_coated_NP = coat_NP(
        xyz_core,
        names_core,
        xyz_lig1,
        names_lig1,
        xyz_pillars1,
        xyz_stones1,
        xyz_lig2,
        names_lig2,
        xyz_pillars2,
        xyz_stones2,
        res_lig1,
        res_lig2,
        VAR["LIG1_S"],
        VAR["LIG2_S"],
        VAR["ELONGATED"],
    )

    logger.info("Writing pdb of the coated nanoparticle...")
    # Writes pdb of the nanoparticle
    print_NP_pdb(
        xyz_coated_NP,
        names_coated_NP,
        res_coated_NP,
        xyz_anchors1,
        xyz_anchors2,
        xyz_lig1,
        xyz_lig2,
        TMP + "/NP.pdb",
    )

    ###########################Parameters (parmchk2, tleap, acpype)#####################################

    logger.info("Running parmchk2 for ligand1...")
    # Generated frcmod of the ligand (includes parameters with the S atom)
    os.system(
        "parmchk2 -i {}/LIG1.mol2 -f mol2 -o {}/LIG1.frcmod -a y".format(TMP, TMP)
    )

    logger.info("Checking parameters for ligand1...")

    check_frcmod(TMP + "/LIG1.frcmod")

    # If parameters could not be assigned, they are shown in the log file
    if two_lig:
        # Generated frcmod of the ligand (includes parameters with the S atom)
        logger.info("Running parmchk2 for ligand2...")
        os.system(
            "parmchk2 -i {}/LIG2.mol2 -f mol2 -o {}/LIG2.frcmod -a y".format(TMP, TMP)
        )
        logger.info("Checking parameters for ligand2...")
        # If parameters could not be assigned, they are shown in the log file
        check_frcmod(TMP + "/LIG2.frcmod")

    logger.info("Writing tleap input file...")

    # Writes file to be run by tleap
    write_leap(VAR, TMP, two_lig)

    logger.info("Running tleap...")
    os.system("tleap -sf {}/TLeap.in > {}/TLeap.log".format(TMP, TMP))

    logger.info("Running acpype...")
    try:
        logger.info("Running acpype...")
        acpype_system = MolTopol(
            acFileXyz="{}/NP.inpcrd".format(TMP),
            acFileTop="{}/NP.prmtop".format(TMP),
            debug=False,
            basename="{}/NP".format(TMP),
            verbose=False,
            gmx45=True,
            disam=False,
            direct=False,
            is_sorted=False,
            chiral=False,
        )
        acpype_system.writeGromacsTopolFiles(amb2gmx=True)
    except OSError as e:
        raise ValueError(
            "Could not generate NP.inpcrd. This usually happens because the atom types in your ligand mol2 file are not in AMBER format. Check NanoModeler's documentation for more information."
        )

    # to save without tmp path in the molecule name
    replace_str = "{}/".format(TMP)
    with open("{}/NP.gro".format(TMP), "r") as gro_file:
        gro_filedata = gro_file.read()

    gro_filedata = gro_filedata.replace(replace_str, "")

    with open("{}/NP.gro".format(TMP), "w") as gro_file:
        gro_file.write(gro_filedata)

    ##############################Staples########################

    logger.info("Reading gro file of the coated nanoparticle...")
    # Reads gro file without velocities
    xyz_sys, names_sys = load_gro(TMP + "/NP.gro")

    logger.info("Reading top file of the unlinked nanoparticle...")
    # Saves types and residue name of all atoms in the system
    types_sys, res_sys = load_top(TMP + "/NP.top")

    logger.info("Looking for anchoring carbon atoms for ligand1...")
    # Get indexes of all the C atoms
    ndx_C1 = get_ndxs(xyz_sys, names_sys, types_sys, res_sys, res_lig1, xyz_anchors1)
    if two_lig:
        logger.info("Looking for anchoring carbon atoms for ligand2...")
        # Get indexes of all the C atoms
        ndx_C2 = get_ndxs(
            xyz_sys, names_sys, types_sys, res_sys, res_lig2, xyz_anchors2
        )
    else:
        ndx_C2 = [], []
    # Makes list of blocks consisting of 1C, 1S and 2Au atoms
    blocks = make_blocks(xyz_sys, names_sys, names_core, res_core, ndx_C1)
    if two_lig:
        # Appends list of blocks consisting of 1C, 1S and 2Au atoms
        blocks = blocks + make_blocks(xyz_sys, names_sys, names_core, res_core, ndx_C2)

    # Manually assigns the intra- and inter-blocks parameters
    logger.info("Writing bonds parameters...")
    write_bonds(blocks, TMP + "/bonds.top", xyz_sys, names_sys)
    logger.info("Writing angles parameters...")
    write_angles(blocks, TMP + "/angles.top", xyz_sys, names_sys)
    logger.info("Writing final topology file...")
    write_topology(TMP + "/NP.top", TMP + "/bonds.top", TMP + "/angles.top")
    #############################################################

    # If NanoModeler crashes, not it wont run anything after
    atexit.unregister(cleanup_error)

    # flush the filehandler
    reportFileHandler.flush()
    # remove the filehandler from the logger
    report.handlers.remove(reportFileHandler)
    # close the filehandler
    reportFileHandler.close()

    logger.info("Compressing final files...")

    # List of files to compress to output
    copy = ["LIG1.mol2", "NP.pdb", "NP.top", "NP.gro"]
    if two_lig:
        copy.append("LIG2.mol2")

    zip_path = TMP + ".zip"
    zip_tmp = zipfile.ZipFile(zip_path, "w")

    # Saves the list of output files in zip file
    for i in copy:
        zip_tmp.write("{}/{}".format(TMP, i))

    zip_tmp.close()

    # Opens and reads the saved zip file
    zf = open(zip_path, "rb")
    zip_data = zf.read()
    zf.close()

    logger.info("Cleaning...")
    # List of files of delete at the end of the run (including the zip file)
    bye = ["ANTECHAMBER.FRCMOD", "leap.log", "md.mdp", "em.mdp", zip_path]
    for i in bye:
        os.remove(i)

    # shutil.rmtree(TMP)     #Deletes the temporary folder

    logger.info(
        '"Señoras y señores, bienvenidos al party, agarren a su pareja (de la cintura) y preparense porque lo que viene no esta facil, no esta facil no."\n\tIvy Queen.'
    )
    logger.info("NanoModeler terminated normally. Que gracias.")

    return (1, zip_data)


if __name__ == "__main__":
    NanoModeler(
        LIG1_FILE=open("user_test/pmba-protonated.mol2"),
        CAP1=[15],
        LIG1_C=6,
        LIG1_S=0,
        LIG1_FRAC=0.5,
        MORPHOLOGY="janus",
        RSEED=666,
        STRIPES=4,
        LIG2_FILE=open("user_test/pmba-deprotonated.mol2"),
        CAP2=[14],
        LIG2_C=6,
        LIG2_S=0,
        # open("LIGANDS/usr5.frcmod"),
        FRCMOD=None,
        CORE=open("CORES/au25SR18_NM.pdb"),
        ELONGATED=False,
    )
