def NanoModeler(LIG1_FILE=None, CAP1=[], LIG1_C=0, LIG1_S=0, LIG1_FRAC=1.0, MORPHOLOGY="random", RSEED=666, STRIPES=1, LIG2_FILE=None, CAP2=[], LIG2_C=0, LIG2_S=0, FRCMOD=None, CORE=None, ELONGATED=False):
    VAR = {
    "LIG1_FILE" : LIG1_FILE,   #ByteArray of the mol2 of ligand1 (must be in the working directory)
    "CAP1" : CAP1,		    #Atom numbers (as in the mol2, likely to start in 1) to remove. Numbers separated by commas
    "LIG1_C"    : LIG1_C,    #Atom number of carbon atom in LIG1_FILE used as anchor
    "LIG1_S"    : LIG1_S,            #Atom number in the original mol2 file of ligand 1 corresponding to the anchoring S atom
    "LIG1_FRAC" : LIG1_FRAC,          #Fraction of ligand1 to place (0-1.0)
    "MORPHOLOGY" : MORPHOLOGY,    #Morphology to distribute ligands1 and 2. random, janus, and stripe are allowed
    "RSEED" : RSEED,              #Random seed for random morphology
    "STRIPES" : STRIPES,              #Number of stripes for stripe morphology. It will start (bottom up with ligand 1)

    "LIG2_FILE" : LIG2_FILE,   #ByteArray of the mol2 of ligand2 (must be in the working directory)
    "CAP2" : CAP2,               #Atom numbers (as in the mol2, likely to start in 1) to remove. Numbers separated by commas
    "LIG2_C"    : LIG2_C,    #Atom number of carbon atom in LIG2_FILE used as anchor
    "LIG2_S"    : LIG2_S,            #Atom number in the original mol2 file of ligand 1 corresponding to the anchoring S atom

    "FRCMOD"    : FRCMOD,       #ByteArray to a frcmod provided by the user
    "CORE" : CORE,    #ByteArray of the pdb of the core
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
    from DEPENDENCIES.rewrite_mol2 import rewrite_mol2
    log += "Importing cleanup dependency...\n\n"
    from DEPENDENCIES.cleanup import cleanup_error

    log += "Creating temporary folder...\n"
    TMP = tempfile.mkdtemp(dir="./")
    atexit.register(cleanup_error, TMP, log)

    log += "Checking input options...\n"
    check_VAR(VAR)

    if (VAR["MORPHOLOGY"] == "stripe"):
        if VAR["STRIPES"] == 1:
            two_lig = False
        else:
            two_lig = True
    else:
        two_lig = (VAR["LIG1_FRAC"] < 1.0)
    log += "Imported options:\n"
    for i in VAR:
        log += "\t{:<20}{:>20}\n".format(i, str(VAR[i]))

    log += "\nOne ligand was found...\n"
    log += "Checking ligand1 mol2 file...\n"
    LIG1_MOL2 = VAR["LIG1_FILE"].readlines()
    LIG1_MOL2 = [s.replace("\n", "") for s in LIG1_MOL2]
    log = check_mol2(LIG1_MOL2, log)
    log += "Rewriting ligand1 file...\n"
    VAR["LIG1_C"], log = rewrite_mol2(LIG1_MOL2, VAR["CAP1"], VAR["LIG1_S"], VAR["LIG1_C"], TMP+"/LIG1.mol2", VAR["ELONGATED"], log)

    if two_lig:
        log += "Two ligands were found...\n"
        log += "Checking ligand2 mol2 file...\n"
        LIG2_MOL2 = VAR["LIG2_FILE"].readlines()
        LIG2_MOL2 = [s.replace("\n", "") for s in LIG2_MOL2]
        log = check_mol2(LIG2_MOL2, log)
        log += "Rewriting ligand2 file...\n"
        VAR["LIG2_C"], log = rewrite_mol2(LIG2_MOL2, VAR["CAP2"], VAR["LIG2_S"], VAR["LIG2_C"], TMP+"/LIG2.mol2", VAR["ELONGATED"], log)

    ##############################NP_builder########################

    log += "Initializing core...\n"
    CORE_PDB = VAR["CORE"].readlines()
    CORE_PDB = [s.replace("\n", "") for s in CORE_PDB]
    xyz_core, names_core, res_core = init_core_pdb(CORE_PDB, VAR["ELONGATED"])

    log += "Initializing ligand1...\n"
    xyz_lig1, names_lig1, res_lig1 = init_lig_mol2(TMP+"/LIG1.mol2", VAR["LIG1_S"], VAR["LIG1_C"])
    if two_lig:
        log += "Initializing ligand2...\n"
        xyz_lig2, names_lig2, res_lig2 = init_lig_mol2(TMP+"/LIG2.mol2", VAR["LIG1_S"], VAR["LIG2_C"])
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

    log +="Assigning morphology...\n"
    xyz_anchors1, xyz_anchors2, log = assign_morph(xyz_core, names_core, VAR["LIG1_FRAC"], VAR["RSEED"], VAR["MORPHOLOGY"], VAR["STRIPES"], log)

    xyz_stones1 = get_stones(xyz_core, names_core, xyz_anchors1, xyz_pillars1, VAR["LIG1_S"])
    if two_lig:
        xyz_stones2 = get_stones(xyz_core, names_core, xyz_anchors2, xyz_pillars2, VAR["LIG2_S"])
    else:
        xyz_stones2 = []

    log += "Coating nanoparticle...\n"
    xyz_coated_NP, names_coated_NP, res_coated_NP, log = coat_NP(xyz_core, names_core, xyz_lig1, names_lig1, xyz_pillars1, xyz_stones1, xyz_lig2, names_lig2, xyz_pillars2, xyz_stones2, res_lig1, res_lig2, VAR["LIG1_S"], VAR["LIG2_S"], VAR["ELONGATED"], log)

    log += "Writing pdb of the coated nanoparticle...\n"
    print_NP_pdb(xyz_coated_NP, names_coated_NP, res_coated_NP, xyz_anchors1, xyz_anchors2, xyz_lig1, xyz_lig2, TMP+"/NP.pdb")

    ################################################################

    log += "Running parmchk2 for ligand1...\n"
    os.system("parmchk2 -i {}/LIG1.mol2 -f mol2 -o {}/LIG1.frcmod -a y".format(TMP, TMP))
    log += "Checking parameters for ligand1...\n"
    log = check_frcmod(TMP+"/LIG1.frcmod", log)
    if two_lig:
        log += "Running parmchk2 for ligand2...\n"
        os.system("parmchk2 -i {}/LIG2.mol2 -f mol2 -o {}/LIG2.frcmod -a y".format(TMP, TMP))
        log += "Checking parameters for ligand2...\n"
        log = check_frcmod(TMP+"/LIG2.frcmod", log)

    log += "Writing tleap input file...\n"
    write_leap(VAR, TMP, two_lig)
    log += "Running tleap...\n"
    os.system("tleap -sf {}/TLeap.in > {}/TLeap.log".format(TMP, TMP))

    log += "Running acpype...\n"
    os.system("python DEPENDENCIES/acpype.py -b {}/NP -c user -p {}/NP.prmtop -x {}/NP.inpcrd > {}/acpype.log".format(TMP, TMP, TMP, TMP))
    os.system("mv {}/NP_GMX.top {}/NP.top".format(TMP, TMP))
    os.system("mv {}/NP_GMX.gro {}/NP.gro".format(TMP, TMP))

    ##############################Staples########################
    log += "Reading gro file of the coated nanoparticle...\n"
    xyz_sys, names_sys = load_gro(TMP+"/NP.gro")

    log += "Reading top file of the unlinked nanoparticle...\n"
    types_sys, res_sys = load_top(TMP+"/NP.top")

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
    write_topology(TMP+"/NP.top", TMP+"/bonds.top", TMP+"/angles.top")
    #############################################################

    atexit.unregister(cleanup_error)

    log += "Compressing final files...\n"

    copy = ["LIG1.mol2", "NP.pdb", "NP.top", "NP.gro"]
    if two_lig:
        copy.append("LIG2.mol2")

    final_zip = zipfile.ZipFile(TMP+".zip", "w")
    for i in copy:
        final_zip.write("{}/{}".format(TMP, i))
    final_zip.close()

    log += "Cleaning...\n"
    bye = ["ANTECHAMBER.FRCMOD", "leap.log", "md.mdp", "em.mdp"]
    for i in bye:
        os.remove(i)
    #shutil.rmtree(TMP)

    log += "\"Señoras y señores, bienvenidos al party, agarren a su pareja (de la cintura) y preparense porque lo que viene no esta facil, no esta facil no.\"\n\tIvy Queen.\n"
    log += "NanoModeler terminated normally. Que gracias.\n"
    print(log)
    return (1, log, final_zip)

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

        FRCMOD=None,

        CORE=open("CORES/au25SR18_NM.pdb"),
        ELONGATED=False)
