def NanoModeler(NAME="test", LIG1_FILE="LIG1.mol2", CAP1="0", LIG1_C=0, LIG1_S=0, LIG1_FRAC=1.0, MORPHOLOGY="random", RSEED=666, STRIPES=1, LIG2_FILE="XXX.mol2", CAP2="0", LIG2_C=0, LIG2_S=0, FRCMOD="0", CORE="au144SR60_NM.pdb", ELONGATED=False):
    VAR = {
    "NAME": NAME,             #Name of the project
    "LIG1_FILE" : LIG1_FILE,   #Name of the mol2 of ligand1 (must be in the working directory)
    "CAP1" : CAP1,		    #Atom numbers (as in the mol2, likely to start in 1) to remove. Numbers separated by commas
    "LIG1_C"    : LIG1_C,    #Atom number of carbon atom in LIG1_FILE used as anchor
    "LIG1_S"    : LIG1_S,            #Atom number in the original mol2 file of ligand 1 corresponding to the anchoring S atom
    "LIG1_FRAC" : LIG1_FRAC,          #Fraction of ligand1 to place (0-1.0)
    "MORPHOLOGY" : MORPHOLOGY,    #Morphology to distribute ligands1 and 2. random, janus, and stripe are allowed
    "RSEED" : RSEED,              #Random seed for random morphology
    "STRIPES" : STRIPES,              #Number of stripes for stripe morphology. It will start (bottom up with ligand 1)

    "LIG2_FILE" : LIG2_FILE,   #Name of the mol2 of ligand2 (must be in the working directory)
    "CAP2" : CAP2,               #Atom numbers (as in the mol2, likely to start in 1) to remove. Numbers separated by commas
    "LIG2_C"    : LIG2_C,    #Atom number of carbon atom in LIG2_FILE used as anchor
    "LIG2_S"    : LIG2_S,            #Atom number in the original mol2 file of ligand 1 corresponding to the anchoring S atom

    "FRCMOD"    : FRCMOD,       #Path to a frcmod provided by the user
    "CORE" : CORE,    #Name of the core to coat. Found in CORES/CORE
    "ELONGATED" : ELONGATED    #If set to true, the first carbon of the cores will be ignored and placed along the O-S axis
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
    from DEPENDENCIES.rewrite_mol2 import rewrite_mol2_with_C, rewrite_mol2_with_S
    log += "Importing cleanup dependency...\n\n"
    from DEPENDENCIES.cleanup import cleanup_error
    log += "Creating temporary folder...\n"

    TMP = tempfile.mkdtemp(dir="./")
    atexit.register(cleanup_error, TMP, log)

    log = check_VAR(VAR, log)

    if (VAR["MORPHOLOGY"] == "stripe"):
        if VAR["STRIPES"] == 1:
            two_lig = False
        else:
            two_lig = True
    else:
        two_lig = (VAR["LIG1_FRAC"] < 1.0)

    log += "Imported options:\n"
    for i in VAR:
        log += (i.ljust(20) + str(VAR[i]).ljust(50)+"\n")

    log += "\n\nOne ligand was found...\n"
    log += "Checking ligand1 mol2 file...\n"
    log = check_mol2(VAR["LIG1_FILE"], log)
    log += "Rewriting ligand1 file...\n"
    if LIG1_S != 0:
        VAR["LIG1_C"], log = rewrite_mol2_with_S(VAR["LIG1_FILE"], VAR["CAP1"], VAR["LIG1_S"], VAR["LIG1_C"], TMP+"/"+VAR["LIG1_FILE"], VAR["ELONGATED"], log)
    else:
        VAR["LIG1_C"], log = rewrite_mol2_with_C(VAR["LIG1_FILE"], VAR["CAP1"], VAR["LIG1_C"], TMP+"/"+VAR["LIG1_FILE"], log)

    if two_lig:
        log += "Two ligands were found...\n"
        log += "Checking ligand2 mol2 file...\n"
        log = check_mol2(VAR["LIG2_FILE"], log)
        log += "Rewriting ligand2 file...\n"
        if LIG2_S != 0:
            VAR["LIG2_C"], log = rewrite_mol2_with_S(VAR["LIG2_FILE"], VAR["CAP2"], VAR["LIG2_S"], VAR["LIG2_C"], TMP+"/"+VAR["LIG2_FILE"], VAR["ELONGATED"], log)
        else:
            VAR["LIG2_C"], log = rewrite_mol2_with_C(VAR["LIG2_FILE"], VAR["CAP2"], VAR["LIG2_C"], TMP+"/"+VAR["LIG2_FILE"], log)

    ##############################NP_builder########################

    log += "Initializing core...\n"
    xyz_core, names_core, res_core = init_core_pdb("./CORES/"+VAR["CORE"], VAR["ELONGATED"])

    log += "Initializing ligand1...\n"
    xyz_lig1, names_lig1, res_lig1 = init_lig_mol2(TMP+"/"+VAR["LIG1_FILE"], VAR["LIG1_S"], VAR["LIG1_C"])
    if two_lig:
        log += "Initializing ligand2...\n"
        xyz_lig2, names_lig2, res_lig2 = init_lig_mol2(TMP+"/"+VAR["LIG2_FILE"], VAR["LIG1_S"], VAR["LIG2_C"])
    else:
        xyz_lig2, names_lig2, res_lig2 = [], [], []

    log += "Running PCA for ligand1...\n"
    xyz_pillars1, log = get_ligand_pill(xyz_lig1, VAR["LIG1_C"], VAR["LIG1_S"], log)
    if two_lig:
        log += "Running PCA for ligand2...\n"
        xyz_pillars2, log = get_ligand_pill(xyz_lig2, VAR["LIG2_C"], VAR["LIG2_S"], log)
    else:
        xyz_pillars2 = []

    N_S = len(names_core[names_core=='ST'])

    xyz_anchors1, xyz_anchors2, log = assign_morph(xyz_core, names_core, VAR["LIG1_FRAC"], VAR["RSEED"], VAR["MORPHOLOGY"], VAR["STRIPES"], log)

    xyz_stones1 = get_stones(xyz_core, names_core, xyz_anchors1, xyz_pillars1, VAR["LIG1_S"])
    if two_lig:
        xyz_stones2 = get_stones(xyz_core, names_core, xyz_anchors2, xyz_pillars2, VAR["LIGs_S"])
    else:
        xyz_stones2 = []

    log += "Coating nanoparticle...\n"
    xyz_coated_NP, names_coated_NP, res_coated_NP, log = coat_NP(xyz_core, names_core, xyz_lig1, names_lig1, xyz_pillars1, xyz_stones1, xyz_lig2, names_lig2, xyz_pillars2, xyz_stones2, res_lig1, res_lig2, VAR["LIG1_S"], VAR["LIG2_S"], VAR["ELONGATED"], log)

    log += "Writing pdb of the coated nanoparticle...\n"
    print_NP_pdb(xyz_coated_NP, names_coated_NP, res_coated_NP, xyz_anchors1, xyz_anchors2, xyz_lig1, xyz_lig2, TMP+"/"+VAR["NAME"]+".pdb")

    ################################################################

    log += "Running parmchk2 for ligand1...\n"
    os.system("parmchk2 -i {} -f mol2 -o {} -a y".format(TMP+"/"+VAR["LIG1_FILE"], TMP+"/"+VAR["LIG1_FILE"][:-5]+".frcmod"))
    log += "Checking parameters for ligand1...\n"
    log = check_frcmod(TMP+"/"+VAR["LIG1_FILE"][:-5]+".frcmod", log)
    if two_lig:
        log += "Running parmchk2 for ligand2...\n"
        os.system("parmchk2 -i {} -f mol2 -o {} -a y".format(TMP+"/"+VAR["LIG2_FILE"], TMP+"/"+VAR["LIG2_FILE"][:-5]+".frcmod"))
        log += "Checking parameters for ligand2...\n"
        log = check_frcmod(TMP+"/"+VAR["LIG2_FILE"][:-5]+".frcmod", log)

    log += "Writing tleap input file...\n"
    write_leap(VAR, TMP, two_lig)
    log += "Running tleap...\n"
    os.system("tleap -sf " + TMP + "/TLeap.in > " + TMP + "/TLeap.log")

    log += "Running acpype...\n"
    os.system("python DEPENDENCIES/acpype.py -p {} -x {} -r".format(TMP+"/"+VAR["NAME"]+".prmtop", TMP+"/"+VAR["NAME"]+".inpcrd -b " + VAR["NAME"] + " -c user > acpype.log"))
    os.system("mv {}_GMX.top {}.top".format(VAR["NAME"], TMP+"/"+VAR["NAME"]))
    os.system("mv {}_GMX.gro {}.gro".format(VAR["NAME"], TMP+"/"+VAR["NAME"]))

    ##############################Staples########################
    log += "Reading gro file of the coated nanoparticle...\n"
    xyz_sys, names_sys = load_gro(TMP+"/"+VAR["NAME"]+".gro")

    log += "Reading top file of the unlinked nanoparticle...\n"
    types_sys, res_sys = load_top(TMP+"/"+VAR["NAME"]+".top")

    log += "Looking for anchoring carbon atoms for ligand1...\n"
    ndx_C1 = get_ndxs(xyz_sys, names_sys, types_sys, res_sys, res_lig1)
    if two_lig:
        log += "Looking for anchoring carbon atoms for ligand2...\n"
        ndx_C2 = get_ndxs(xyz_sys, names_sys, types_sys, res_sys, res_lig2)
    else:
        ndx_C2 = [], []

    blocks = make_blocks(xyz_sys, names_sys, names_core, res_core, ndx_C1)
    if two_lig:
        blocks = blocks + make_blocks(xyz_sys, names_sys, names_core, res_core, ndx_C2)

    log += "Writing bonds parameters...\n"
    write_bonds(blocks, TMP+"/bonds.top", xyz_sys, names_sys)
    log += "Writing angles parameters...\n"
    write_angles(blocks, TMP+"/angles.top", xyz_sys, names_sys)
    log += "Writing final topology file...\n"
    write_topology(TMP+"/"+VAR["NAME"]+".top", TMP+"/bonds.top", TMP+"/angles.top")
    #############################################################

    atexit.unregister(cleanup_error)

    log += "Compressing final files...\n"
    #os.mkdir(VAR["NAME"])
    copy = [VAR["LIG1_FILE"], VAR["NAME"]+".pdb", VAR["NAME"]+".top", VAR["NAME"]+".gro"]
    if two_lig:
        copy.append(VAR["LIG2_FILE"])

    final_zip = zipfile.ZipFile(VAR["NAME"]+".zip", "w")
    for i in copy:
        final_zip.write(TMP+"/"+i)
    final_zip.close()

    log += "Cleaning...\n"
    bye = ["ANTECHAMBER.FRCMOD", "leap.log", "md.mdp", "em.mdp", "acpype.log"]
    for i in bye:
        os.remove(i)
    #shutil.rmtree(TMP)

    log += "NanoModeler terminated normally. Que gracias.\n"
    print(log)
    return (1, log, final_zip)

NanoModeler(NAME="test",

    LIG1_FILE="LIG2.mol2",
    CAP1="0",
    LIG1_C=1,
    LIG1_S=0,

    LIG1_FRAC=1.0,
    MORPHOLOGY="random",
    RSEED=666,
    STRIPES=1,

    LIG2_FILE="XXX.mol2",
    CAP2="0",
    LIG2_C=0,
    LIG2_S=0,

    FRCMOD="0",
    CORE="au25SR18_NM.pdb",
    ELONGATED=False)
