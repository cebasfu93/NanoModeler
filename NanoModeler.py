def NanoModeler(LIG1_FILE=None, CAP1=[], LIG1_C=0, LIG1_S=0, LIG1_FRAC=1.0, MORPHOLOGY="random", RSEED=None, STRIPES=1, LIG2_FILE=None, CAP2=[], LIG2_C=0, LIG2_S=0, FRCMOD=None, CORE=None, ELONGATED=False):
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

    log = "WELCOME TO NANOMODELER\n\n"
    log += "Importing sys library...\n"
    import sys
    log += "Importing numpy library...\n"
    import numpy as np
    log += "Importing random library...\n"
    import random
    log += "Importing math library...\n"
    import math
    log += "Importing scipy library...\n"
    from scipy.spatial import distance
    log += "Importing sklearn library...\n"
    from sklearn.decomposition import PCA
    log += "Importing os library...\n"
    import os
    log += "Importing shutil library...\n"
    import shutil
    log += "Importing tempfile...\n"
    import tempfile
    log += "Importing zipfile...\n"
    import zipfile
    log += "Importing atexit library...\n"
    import atexit
    log += "Importing re library...\n"
    import re
    log += "Importing pickle library...\n"
    import pickle
    log += "Importing subprocess library...\n"
    import subprocess
    log += "Importing signal library...\n"
    import signal
    log += "Importing traceback library...\n"
    import traceback
    log += "Importing datetime library...\n"
    import datetime
    log += "Importing acpype.py...\n"
    from DEPENDENCIES.acpype import MolTopol
    log += "Importing NP_builder dependency...\n"
    from DEPENDENCIES.NP_builder import init_lig_mol2, init_core_pdb, get_ligand_pill, assign_morph, get_stones, coat_NP, print_NP_pdb
    log += "Importing staples dependency...\n"
    from DEPENDENCIES.staples import load_gro, load_top, get_ndxs, make_blocks, write_bonds, write_angles, write_topology
    log += "Importing default function dependency...\n"
    from DEPENDENCIES.defaults import check_VAR, check_mol2, check_frcmod, write_leap
    log += "Importing transformations dependency...\n"
    from  DEPENDENCIES.transformations import affine_matrix_from_points, vector_norm, quaternion_matrix
    log += "Importing subunits dependency...\n"
    import DEPENDENCIES.subunits
    log += "Importing rewrite_mol2 dependency...\n"
    from DEPENDENCIES.rewrite_mol2 import rewrite_mol2
    log += "Importing cleanup dependency...\n\n"
    from DEPENDENCIES.cleanup import cleanup_error

    log += "Creating temporary folder...\n"
    TMP = tempfile.mkdtemp(dir="./")
    atexit.register(cleanup_error, TMP, log)    #From this point, the cleanup_error will be excecuted after NanoModeler finishes (or crashes)

    log += "Checking input options...\n"
    check_VAR(VAR)      #Checks validity of values in the variables

    #Based on the morphology, number of stripes, and fraction of ligand1, it is decided if there are one or two ligands
    if (VAR["MORPHOLOGY"] == "stripe"):
        if VAR["STRIPES"] == 1:
            two_lig = False
        else:
            two_lig = True
    else:
        two_lig = (VAR["LIG1_FRAC"] < 1.0)
    log += "Imported options:\n"
    for i in VAR:
        if (i == "LIG1_FILE" or i == "LIG2_FILE" or i == "CORE" or i == "FRCMOD") and VAR[i]:
            if VAR[i]:
                log += "\t{:<20}{:>20}\n".format(i, str(VAR[i].name))
        else:
            log += "\t{:<20}{:>20}\n".format(i, str(VAR[i]))

    log += "\nOne ligand was found...\n"
    log += "Checking ligand1 mol2 file...\n"
    LIG1_MOL2 = VAR["LIG1_FILE"].readlines()        #Read lines and eliminates the \n character from them
    LIG1_MOL2 = [s.replace("\n", "") for s in LIG1_MOL2]
    log = check_mol2(LIG1_MOL2, log)
    log += "Rewriting ligand1 file...\n"
    VAR["LIG1_C"], log = rewrite_mol2(LIG1_MOL2, VAR["CAP1"], VAR["LIG1_S"], VAR["LIG1_C"], TMP+"/LIG1.mol2", VAR["ELONGATED"], log)        #Writes new mol2 without capping, S atom at the end, and properly numbered

    if two_lig:
        log += "Two ligands were found...\n"
        log += "Checking ligand2 mol2 file...\n"
        LIG2_MOL2 = VAR["LIG2_FILE"].readlines()        #Read lines and eliminates the \n character from them
        LIG2_MOL2 = [s.replace("\n", "") for s in LIG2_MOL2]
        log = check_mol2(LIG2_MOL2, log)
        log += "Rewriting ligand2 file...\n"
        VAR["LIG2_C"], log = rewrite_mol2(LIG2_MOL2, VAR["CAP2"], VAR["LIG2_S"], VAR["LIG2_C"], TMP+"/LIG2.mol2", VAR["ELONGATED"], log)        #Writes new mol2 without capping, S atom at the end, and properly numbered

    ##############################NP_builder########################

    log += "Initializing core...\n"
    CORE_PDB = VAR["CORE"].readlines()          #Read lines and eliminates the \n character from them
    CORE_PDB = [s.replace("\n", "") for s in CORE_PDB]
    xyz_core, names_core, res_core = init_core_pdb(CORE_PDB, VAR["ELONGATED"])      #Extracts xyz coordinates, names, and residues of the core

    log += "Initializing ligand1...\n"
    xyz_lig1, names_lig1, res_lig1 = init_lig_mol2(TMP+"/LIG1.mol2", VAR["LIG1_S"], VAR["LIG1_C"])      #Reads rewritten mol2 file and puts the C atom in the origin
    if two_lig:
        log += "Initializing ligand2...\n"
        xyz_lig2, names_lig2, res_lig2 = init_lig_mol2(TMP+"/LIG2.mol2", VAR["LIG1_S"], VAR["LIG2_C"])      #Reads rewritten mol2 file and puts the C atom in the origin
    else:
        xyz_lig2, names_lig2, res_lig2 = [], [], []

    log += "Running PCA for ligand1...\n"
    xyz_pillars1, log = get_ligand_pill(xyz_lig1, VAR["LIG1_C"], VAR["LIG1_S"], log)        #Gives the coordinates of the C atom and the projections of 2 pseudo random atoms into PCA1
    if two_lig:
        log += "Running PCA for ligand2...\n"
        xyz_pillars2, log = get_ligand_pill(xyz_lig2, VAR["LIG2_C"], VAR["LIG2_S"], log)        #Gives the coordinates of the C atom and the projections of 2 pseudo random atoms into PCA1
    else:
        xyz_pillars2 = []

    N_S = len(names_core[names_core=='ST'])

    log +="Assigning morphology...\n"
    xyz_anchors1, xyz_anchors2, log = assign_morph(xyz_core, names_core, VAR["LIG1_FRAC"], VAR["RSEED"], VAR["MORPHOLOGY"], VAR["STRIPES"], log)        #Gives the coordinates of the C atoms in the core corresponding to each ligand depending on the morphology specified

    xyz_stones1 = get_stones(xyz_core, names_core, xyz_anchors1, xyz_pillars1, VAR["LIG1_S"])       #Gives a 3D array with the xyz coordinates for all the stones of all the anchors of ligand 1
    if two_lig:
        xyz_stones2 = get_stones(xyz_core, names_core, xyz_anchors2, xyz_pillars2, VAR["LIG2_S"])       #Gives a 3D array with the xyz coordinates for all the stones of all the anchors of ligand 2
    else:
        xyz_stones2 = []

    log += "Coating nanoparticle...\n"  #Merges xyz coordinates and names of the core and the ligands into one coated NP
    xyz_coated_NP, names_coated_NP, res_coated_NP, log = coat_NP(xyz_core, names_core, xyz_lig1, names_lig1, xyz_pillars1, xyz_stones1, xyz_lig2, names_lig2, xyz_pillars2, xyz_stones2, res_lig1, res_lig2, VAR["LIG1_S"], VAR["LIG2_S"], VAR["ELONGATED"], log)

    log += "Writing pdb of the coated nanoparticle...\n"
    print_NP_pdb(xyz_coated_NP, names_coated_NP, res_coated_NP, xyz_anchors1, xyz_anchors2, xyz_lig1, xyz_lig2, TMP+"/NP.pdb")      #Writes pdb of the nanoparticle

    ###########################Parameters (parmchk2, tleap, acpype)#####################################

    log += "Running parmchk2 for ligand1...\n"
    os.system("parmchk2 -i {}/LIG1.mol2 -f mol2 -o {}/LIG1.frcmod -a y".format(TMP, TMP))       #Generated frcmod of the ligand (includes parameters with the S atom)
    log += "Checking parameters for ligand1...\n"
    log = check_frcmod(TMP+"/LIG1.frcmod", log)         #If parameters could not be assigned, they are shown in the log file
    if two_lig:
        log += "Running parmchk2 for ligand2...\n"
        os.system("parmchk2 -i {}/LIG2.mol2 -f mol2 -o {}/LIG2.frcmod -a y".format(TMP, TMP))       #Generated frcmod of the ligand (includes parameters with the S atom)
        log += "Checking parameters for ligand2...\n"
        log = check_frcmod(TMP+"/LIG2.frcmod", log)     #If parameters could not be assigned, they are shown in the log file

    log += "Writing tleap input file...\n"
    write_leap(VAR, TMP, two_lig)           #Writes file to be run by tleap
    log += "Running tleap...\n"
    os.system("tleap -sf {}/TLeap.in > {}/TLeap.log".format(TMP, TMP))

    log += "Running acpype...\n"
    acpype_system = MolTopol(acFileXyz = "{}/NP.inpcrd".format(TMP), acFileTop = "{}/NP.prmtop".format(TMP),
                      debug = False, basename = "{}/NP".format(TMP),
                      verbose = False,  gmx45 = True,
                      disam = False,     direct = False,
                      is_sorted = False, chiral = False)
    acpype_system.writeGromacsTopolFiles(amb2gmx = True)

    ##############################Staples########################
    log += "Reading gro file of the coated nanoparticle...\n"
    xyz_sys, names_sys = load_gro(TMP+"/NP.gro")        #Reads gro file without velocities

    log += "Reading top file of the unlinked nanoparticle...\n"
    types_sys, res_sys = load_top(TMP+"/NP.top")        #Saves types and residue name of all atoms in the system

    log += "Looking for anchoring carbon atoms for ligand1...\n"
    ndx_C1 = get_ndxs(xyz_sys, names_sys, types_sys, res_sys, res_lig1)     #Get indexes of all the C atoms
    if two_lig:
        log += "Looking for anchoring carbon atoms for ligand2...\n"
        ndx_C2 = get_ndxs(xyz_sys, names_sys, types_sys, res_sys, res_lig2)      #Get indexes of all the C atoms
    else:
        ndx_C2 = [], []

    blocks = make_blocks(xyz_sys, names_sys, names_core, res_core, ndx_C1)      #Makes list of blocks consisting of 1C, 1S and 2Au atoms
    if two_lig:
        blocks = blocks + make_blocks(xyz_sys, names_sys, names_core, res_core, ndx_C2)     #Appends list of blocks consisting of 1C, 1S and 2Au atoms

    #Manually assigns the intra- and inter-blocks parameters
    log += "Writing bonds parameters...\n"
    write_bonds(blocks, TMP+"/bonds.top", xyz_sys, names_sys)
    log += "Writing angles parameters...\n"
    write_angles(blocks, TMP+"/angles.top", xyz_sys, names_sys)
    log += "Writing final topology file...\n"
    write_topology(TMP+"/NP.top", TMP+"/bonds.top", TMP+"/angles.top")
    #############################################################

    atexit.unregister(cleanup_error)        #If NanoModeler crashes, not it wont run anything after

    log += "Compressing final files...\n"

    copy = ["LIG1.mol2", "NP.pdb", "NP.top", "NP.gro"]      #List of files to compress to output
    if two_lig:
        copy.append("LIG2.mol2")

    zip_path = TMP + ".zip"
    zip_tmp = zipfile.ZipFile(zip_path, "w")
    for i in copy:
        zip_tmp.write("{}/{}".format(TMP, i))       #Saves the list of output files in zip file
    zip_tmp.close()
    zf = open(zip_path, 'rb')       #Opens and reads the saved zip file
    zip_data = zf.read()
    zf.close()

    log += "Cleaning...\n"
    bye = ["ANTECHAMBER.FRCMOD", "leap.log", "md.mdp", "em.mdp", zip_path]      #List of files of delete at the end of the run (including the zip file)
    for i in bye:
        os.remove(i)
    #shutil.rmtree(TMP)     #Deletes the temporary folder

    log += "\"Senoras y senores, bienvenidos al party, agarren a su pareja (de la cintura) y preparense porque lo que viene no esta facil, no esta facil no.\"\n\tIvy Queen.\n"
    log += "NanoModeler terminated normally. Que gracias.\n"
    print(log)
    return (1, log, zip_data)

if __name__ == "__main__":
    NanoModeler(LIG1_FILE=open("LIGANDS/LIGCA.mol2"),
        CAP1=[],
        LIG1_C=1,
        LIG1_S=0,

        LIG1_FRAC=1.0,
        MORPHOLOGY="random",
        RSEED=666,
        STRIPES=4,

        LIG2_FILE=None, #open("Tests/lig2.mol2"),
        CAP2=[],
        LIG2_C=0,
        LIG2_S=0,

        FRCMOD=None, #open("Tests/frcMod-AP3.frcmod"),

        CORE=open("CORES/au144SR60_NM.pdb"),
        ELONGATED=False)
