# -*- coding: utf-8 -*-

import logging
logger = logging.getLogger('nanomodeler')
logger.addHandler(logging.NullHandler())

__VERSION__ = "1.0.0"

def getVersion():
    return __VERSION__

def NanoModeler(LIG1_FILE=None, CAP1=[], LIG1_C=0, LIG1_S=0, LIG1_FRAC=1.0, MORPHOLOGY="random", RSEED=666, STRIPES=1, LIG2_FILE=None, CAP2=[], LIG2_C=0, LIG2_S=0, FRCMOD=None, CORE=None, ELONGATED=False):

    VAR = {
    "LIG1_FILE" : LIG1_FILE,   #ByteArray of the mol2 of ligand1
    "CAP1" : CAP1,		    #List of atom numbers (as in the mol2, likely to start in 1) to remove
    "LIG1_C"    : LIG1_C,    #Signed integer, atom number of carbon atom in LIG1_FILE used as anchor
    "LIG1_S"    : LIG1_S,            #Signed integer, atom number in the original mol2 file of ligand 1 corresponding to the anchoring S atom
    "LIG1_FRAC" : LIG1_FRAC,          #Float, fraction of ligand1 to place (0-1.0)
    "MORPHOLOGY" : MORPHOLOGY,    #String, morphology to distribute ligands1 and 2. random, janus, and stripe are allowed
    "RSEED" : RSEED,              #Integer, random seed for random morphology
    "STRIPES" : STRIPES,              #Signed integer, number of stripes for stripe morphology. It will start (bottom up with ligand 1)

    "LIG2_FILE" : LIG2_FILE,   #ByteArray of the mol2 of ligand2
    "CAP2" : CAP2,               #List of atom numbers (as in the mol2, likely to start in 1) to remove
    "LIG2_C"    : LIG2_C,    #Signed integer, atom number of carbon atom in LIG2_FILE used as anchor
    "LIG2_S"    : LIG2_S,            #Signed integer, atom number in the original mol2 file of ligand 1 corresponding to the anchoring S atom

    "FRCMOD"    : FRCMOD,       #ByteArray of the frcmod provided by the user
    "CORE" : CORE,    #ByteArray of the pdb of the core
    "ELONGATED" : ELONGATED    #If set to true, the first carbon of the cores will be ignored and placed 1.8A away from the S atoms along the O-S axis
    }

    logger.info("WELCOME TO NANOMODELER")
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
    logger.info("Importing acpype.py...")
    from DEPENDENCIES.acpype import MolTopol
    logger.info("Importing NP_builder dependency...")
    from DEPENDENCIES.NP_builder import init_lig_mol2, init_core_pdb, get_ligand_pill, assign_morph, get_stones, coat_NP, print_NP_pdb
    logger.info("Importing staples dependency...")
    from DEPENDENCIES.staples import load_gro, load_top, get_ndxs, make_blocks, write_bonds, write_angles, write_topology, change_gro_mol_name
    logger.info("Importing default function dependency...")
    from DEPENDENCIES.defaults import check_VAR, check_mol2, check_frcmod, write_leap
    logger.info("Importing transformations dependency...")
    from  DEPENDENCIES.transformations import affine_matrix_from_points, vector_norm, quaternion_matrix
    logger.info("Importing subunits dependency...")
    import DEPENDENCIES.subunits
    logger.info("Importing rewrite_mol2 dependency...")
    from DEPENDENCIES.rewrite_mol2 import rewrite_mol2
    logger.info("Importing cleanup dependency...")
    from DEPENDENCIES.cleanup import cleanup_error

    logger.info("Creating temporary folder...")
    TMP = tempfile.mkdtemp(dir="./")

    rf = open(TMP + "/report.log", "w")
    report = logging.getLogger('report')
    report.setLevel(logging.INFO)
    report.addHandler(logging.StreamHandler(stream = rf))
    atexit.register(cleanup_error, TMP, rf)    #From this point, the cleanup_error will be excecuted after NanoModeler finishes (or crashes)

    logger.info("Checking input options...")
    check_VAR(VAR)      #Checks validity of values in the variables

    #Based on the morphology, number of stripes, and fraction of ligand1, it is decided if there are one or two ligands
    if (VAR["MORPHOLOGY"] == "stripe"):
        if VAR["STRIPES"] == 1:
            two_lig = False
        else:
            two_lig = True
    else:
        two_lig = (VAR["LIG1_FRAC"] < 1.0)

    logger.info("Imported options:")
    report.info("Imported options:")
    for i in VAR:
        if (i == "LIG1_FILE" or i == "LIG2_FILE"or i == "CORE" or i == "FRCMOD") and VAR[i]:
            if VAR[i]:
                op_txt = "\t{:<20}{:>20}".format(i, str(VAR[i].name))
                logger.info(op_txt)
                report.info(op_txt + "\n")
        else:
            op_txt = "\t{:<20}{:>20}".format(i, str(VAR[i]))
            logger.info(op_txt)
            report.info(op_txt)

    logger.info("One ligand was found...")
    logger.info("Checking ligand1 mol2 file...")
    LIG1_MOL2 = VAR["LIG1_FILE"].readlines()        #Read lines and eliminates the \n character from them
    LIG1_MOL2 = [s.replace("\n", "") for s in LIG1_MOL2]
    check_mol2(LIG1_MOL2)
    logger.info("Rewriting ligand1 file...")
    VAR["LIG1_C"] = rewrite_mol2(LIG1_MOL2, VAR["CAP1"], VAR["LIG1_S"], VAR["LIG1_C"], TMP+"/LIG1.mol2", VAR["ELONGATED"])        #Writes new mol2 without capping, S atom at the end, and properly numbered

    if two_lig:
        logger.info("Two ligands were found...")
        logger.info("Checking ligand2 mol2 file...")
        LIG2_MOL2 = VAR["LIG2_FILE"].readlines()        #Read lines and eliminates the \n character from them
        LIG2_MOL2 = [s.replace("\n", "") for s in LIG2_MOL2]
        check_mol2(LIG2_MOL2)
        logger.info("Rewriting ligand2 file...")
        VAR["LIG2_C"] = rewrite_mol2(LIG2_MOL2, VAR["CAP2"], VAR["LIG2_S"], VAR["LIG2_C"], TMP+"/LIG2.mol2", VAR["ELONGATED"])        #Writes new mol2 without capping, S atom at the end, and properly numbered

    ##############################NP_builder########################

    logger.info("Initializing core...")
    CORE_PDB = VAR["CORE"].readlines()          #Read lines and eliminates the \n character from them
    CORE_PDB = [s.replace("\n", "") for s in CORE_PDB]
    xyz_core, names_core, res_core = init_core_pdb(CORE_PDB, VAR["ELONGATED"])      #Extracts xyz coordinates, names, and residues of the core

    logger.info("Initializing ligand1...")
    xyz_lig1, names_lig1, res_lig1 = init_lig_mol2(TMP+"/LIG1.mol2", VAR["LIG1_S"], VAR["LIG1_C"])      #Reads rewritten mol2 file and puts the C atom in the origin
    if two_lig:
        logger.info("Initializing ligand2...")
        xyz_lig2, names_lig2, res_lig2 = init_lig_mol2(TMP+"/LIG2.mol2", VAR["LIG1_S"], VAR["LIG2_C"])      #Reads rewritten mol2 file and puts the C atom in the origin
    else:
        xyz_lig2, names_lig2, res_lig2 = [], [], []

    logger.info("Running PCA for ligand1...")
    xyz_pillars1 = get_ligand_pill(xyz_lig1, VAR["LIG1_C"], VAR["LIG1_S"])        #Gives the coordinates of the C atom and the projections of 2 pseudo random atoms into PCA1
    if two_lig:
        logger.info("Running PCA for ligand2...")
        xyz_pillars2 = get_ligand_pill(xyz_lig2, VAR["LIG2_C"], VAR["LIG2_S"])        #Gives the coordinates of the C atom and the projections of 2 pseudo random atoms into PCA1
    else:
        xyz_pillars2 = []

    N_S = len(names_core[names_core=='ST'])

    logger.info("Assigning morphology...")
    xyz_anchors1, xyz_anchors2 = assign_morph(xyz_core, names_core, VAR["LIG1_FRAC"], VAR["RSEED"], VAR["MORPHOLOGY"], VAR["STRIPES"])        #Gives the coordinates of the C atoms in the core corresponding to each ligand depending on the morphology specified

    xyz_stones1 = get_stones(xyz_core, names_core, xyz_anchors1, xyz_pillars1, VAR["LIG1_S"])       #Gives a 3D array with the xyz coordinates for all the stones of all the anchors of ligand 1
    if two_lig:
        xyz_stones2 = get_stones(xyz_core, names_core, xyz_anchors2, xyz_pillars2, VAR["LIG2_S"])       #Gives a 3D array with the xyz coordinates for all the stones of all the anchors of ligand 2
    else:
        xyz_stones2 = []

    logger.info("Coating nanoparticle...")  #Merges xyz coordinates and names of the core and the ligands into one coated NP
    xyz_coated_NP, names_coated_NP, res_coated_NP = coat_NP(xyz_core, names_core, xyz_lig1, names_lig1, xyz_pillars1, xyz_stones1, xyz_lig2, names_lig2, xyz_pillars2, xyz_stones2, res_lig1, res_lig2, VAR["LIG1_S"], VAR["LIG2_S"], VAR["ELONGATED"])

    logger.info("Writing pdb of the coated nanoparticle...")
    print_NP_pdb(xyz_coated_NP, names_coated_NP, res_coated_NP, xyz_anchors1, xyz_anchors2, xyz_lig1, xyz_lig2, TMP+"/NP.pdb")      #Writes pdb of the nanoparticle

    ################################################################

    logger.info("Running parmchk2 for ligand1...")
    os.system("parmchk2 -i {}/LIG1.mol2 -f mol2 -o {}/LIG1.frcmod -a y".format(TMP, TMP))       #Generated frcmod of the ligand (includes parameters with the S atom)
    logger.info("Checking parameters for ligand1...")
    check_frcmod(TMP+"/LIG1.frcmod")         #If parameters could not be assigned, they are shown in the log file
    if two_lig:
        logger.info("Running parmchk2 for ligand2...")
        os.system("parmchk2 -i {}/LIG2.mol2 -f mol2 -o {}/LIG2.frcmod -a y".format(TMP, TMP))       #Generated frcmod of the ligand (includes parameters with the S atom)
        logger.info("Checking parameters for ligand2...")
        check_frcmod(TMP+"/LIG2.frcmod")     #If parameters could not be assigned, they are shown in the log file

    logger.info("Writing tleap input file...")
    write_leap(VAR, TMP, two_lig)           #Writes file to be run by tleap
    logger.info("Running tleap...")
    os.system("tleap -sf {}/TLeap.in > {}/TLeap.log".format(TMP, TMP))

    logger.info("Running acpype...")
    acpype_system = MolTopol(acFileXyz = "{}/NP.inpcrd".format(TMP), acFileTop = "{}/NP.prmtop".format(TMP),
                      debug = False, basename = "{}/NP".format(TMP),
                      verbose = False,  gmx45 = True,
                      disam = False,     direct = False,
                      is_sorted = False, chiral = False)
    acpype_system.writeGromacsTopolFiles(amb2gmx = True)

    ##############################Staples########################
    logger.info("Reading gro file of the coated nanoparticle...")
    xyz_sys, names_sys = load_gro(TMP+"/NP.gro")        #Reads gro file without velocities

    logger.info("Reading top file of the unlinked nanoparticle...")
    types_sys, res_sys = load_top(TMP+"/NP.top")        #Saves types and residue name of all atoms in the system

    logger.info("Looking for anchoring carbon atoms for ligand1...")
    ndx_C1 = get_ndxs(xyz_sys, names_sys, types_sys, res_sys, res_lig1)     #Get indexes of all the C atoms
    if two_lig:
        logger.info("Looking for anchoring carbon atoms for ligand2...")
        ndx_C2 = get_ndxs(xyz_sys, names_sys, types_sys, res_sys, res_lig2)      #Get indexes of all the C atoms
    else:
        ndx_C2 = [], []

    blocks = make_blocks(xyz_sys, names_sys, names_core, res_core, ndx_C1)      #Makes list of blocks consisting of 1C, 1S and 2Au atoms
    if two_lig:
        blocks = blocks + make_blocks(xyz_sys, names_sys, names_core, res_core, ndx_C2)     #Appends list of blocks consisting of 1C, 1S and 2Au atoms

    #Manually assigns the intra- and inter-blocks parameters
    logger.info("Writing bonds parameters...")
    write_bonds(blocks, TMP+"/bonds.top", xyz_sys, names_sys)
    logger.info("Writing angles parameters...")
    write_angles(blocks, TMP+"/angles.top", xyz_sys, names_sys)
    logger.info("Writing final topology file...")
    write_topology(TMP+"/NP.top", TMP+"/bonds.top", TMP+"/angles.top")
    logger.info("Making nicer the results...")
    change_gro_mol_name(TMP+"/NP.gro")
    #############################################################

    atexit.unregister(cleanup_error)        #If NanoModeler crashes, now it wont run anything after
    report.handlers[0].flush()
    rf.close()

    logger.info("Compressing final files...")

    copy = ["LIG1.mol2", "NP.pdb", "NP.top", "NP.gro", "report.log"]      #List of files to compress to output
    if two_lig:
        copy.append("LIG2.mol2")

    zip_path = TMP + ".zip"
    zip_tmp = zipfile.ZipFile(zip_path, "w")
    for i in copy:
        zip_tmp.write("{}/{}".format(TMP, i),"result/{}".format(i))       #Saves the list of output files in zip file
    zip_tmp.close()
    zf = open(zip_path, 'rb')       #Opens and reads the saved zip file
    zip_data = zf.read()
    zf.close()

    logger.info("Cleaning...")
    bye = ["ANTECHAMBER.FRCMOD", "leap.log", "md.mdp", "em.mdp", zip_path]      #List of files of delete at the end of the run (including the zip file)
    for i in bye:
        os.remove(i)
    shutil.rmtree(TMP)     #Deletes the temporary folder

    logger.info("\"Señoras y señores, bienvenidos al party, agarren a su pareja (de la cintura) y preparense porque lo que viene no esta facil, no esta facil no.\"\n\tIvy Queen.")
    logger.info("NanoModeler terminated normally. Que gracias.")
    return (1, zip_data)

if __name__ == "__main__":
    NanoModeler(LIG1_FILE=open("LIG2.mol2"),
        CAP1=[],
        LIG1_C=1,
        LIG1_S=0,

        LIG1_FRAC=1.0,
        MORPHOLOGY="random",
        RSEED=666,
        STRIPES=1,

        LIG2_FILE=None, #open("LIG3.mol2"),
        CAP2=[],
        LIG2_C=1,
        LIG2_S=0,

        FRCMOD=None, #open("NP2.frcmod"),

        CORE=open("CORES/au144SR60_NM.pdb"),
        ELONGATED=False)
