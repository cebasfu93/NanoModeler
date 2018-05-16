def NanoModeler(NAME="test", LIG1_FILE="LIG1.mol2", CAP1="N", LIG1_FRAC="1.0", MORPHOLOGY="random", RSEED="666", STRIPES="1", LIG2_FILE="XXX.mol2", CAP2="N", CORE="au144SR60_NM.pdb", COREDIR="/DATA/SoftwareSFU/IN-HOUSE/NanoModeler/CORES", DEPENDS="/DATA/SoftwareSFU/IN-HOUSE/NanoModeler/DEPENDENCIES"):
    VAR = {
    "NAME": NAME,             #Name of the project
    "LIG1_FILE" : LIG1_FILE,   #Name of the mol2 of ligand1 (must be in the working directory)
    "CAP1" : CAP1,		    #Atom numbers (as in the mol2, likely to start in 1) to remove. Numbers separated by commas
    "LIG1_FRAC" : LIG1_FRAC,          #Fraction of ligand1 to place (0-1.0)
    "MORPHOLOGY" : MORPHOLOGY,    #Morphology to distribute ligands1 and 2. random, janus, and stripe are allowed
    "RSEED" : RSEED,              #Random seed for random morphology
    "STRIPES" : STRIPES,              #Number of stripes for stripe morphology. It will start (bottom up with ligand 1)

    "LIG2_FILE" : LIG2_FILE,   #Name of the mol2 of ligand2 (must be in the working directory)
    "CAP2" : CAP2,               #Atom numbers (as in the mol2, likely to start in 1) to remove. Numbers separated by commas

    "CORE" : CORE,    #Name of the core to coat. Found in COREDIR/CORE

    "COREDIR" : COREDIR,   #Path to folder containing all the available cores
    "DEPENDS" : DEPENDS #Path to the folder containing all the dependencies that come with NanoModeler
    }

    log = open("NanoModeler.log", "w")
    log.write("WELCOME TO NANOMODELER\n")
    log.write("Importing sys library...\n")
    import sys
    sys.stdout = log
    sys.stderr = log
    print("Importing numpy library...")
    import numpy as np
    print("Importing random library...")
    import random
    print("Importing math library...")
    import math
    print("Importing scipy library...")
    from scipy.spatial import distance
    print("Importing sklearn library...")
    from sklearn.decomposition import PCA
    print("Importing os library...")
    import os
    print("Importing shutil library...")
    import shutil
    print("Importing atexit library...")
    import atexit
    print("Looking for folder with dependencies...")
    sys.path.append(VAR["DEPENDS"])
    print("Importing default function dependency...")
    from defaults import check_VAR, check_mol2, write_leap
    print("Importing transformations...")
    from  transformations import affine_matrix_from_points, vector_norm, quaternion_matrix
    print("Importing subunits dependency...")
    import subunits
    print("Importing rewrite_mol2 dependency...")
    from rewrite_mol2 import rewrite_mol2
    print("Importing NP_builder dependency...")
    from NP_builder import init_lig_mol2, init_core_pdb, get_ligand_pill, assign_morph, get_stones, coat_NP, print_NP_pdb
    print("Importing staples dependency...")
    from staples import load_gro, load_top, get_ndxs, make_blocks, write_bonds, write_angles, write_topology
    print("Importing cleanup dependency...")
    from cleanup import cleanup_error, cleanup_normal
    atexit.register(cleanup_error)
    check_VAR(VAR)

    if (VAR["MORPHOLOGY"] == "stripe"):
        if int(VAR["STRIPES"]) == 1:
            two_lig = False
        else:
            two_lig = True
    else:
        two_lig = (float(VAR["LIG1_FRAC"]) < 1.0)

    print("Imported options:")
    for i in VAR:
        print(i.ljust(20) + VAR[i].ljust(50))

    print("\n\nCreating folder...")
    os.mkdir("TMP")

    print("One ligand was found...")
    print("Checking ligand1 mol2 file...")
    check_mol2(VAR["LIG1_FILE"])
    print("Rewriting ligand1 file...")
    rewrite_mol2(VAR["LIG1_FILE"], VAR["CAP1"], "TMP/"+VAR["LIG1_FILE"])

    if two_lig:
        print("Two ligands were found...")
        print("Checking ligand2 mol2 file...")
        check_mol2(VAR["LIG2_FILE"])
        print("Rewriting ligand2 file...")
        rewrite_mol2(VAR["LIG2_FILE"], VAR["CAP2"], "TMP/"+VAR["LIG2_FILE"])

    ##############################NP_builder########################

    print("Initializing ligand1...")
    xyz_lig1, names_lig1, anchor_ndx1, name_anchor1, res_anchor1, res_lig1 = init_lig_mol2("TMP/"+VAR["LIG1_FILE"], VAR["CAP1"])
    if two_lig:
        print("Initializing ligand2...")
        xyz_lig2, names_lig2, anchor_ndx2, name_anchor2, res_anchor2, res_lig2 = init_lig_mol2("TMP/"+VAR["LIG2_FILE"], VAR["CAP2"])
    else:
        xyz_lig2, names_lig2, anchor_ndx2, name_anchor2, res_anchor2, res_lig2 = [], [], [], [], [], []

    print("Initializing core...")
    xyz_core, names_core, res_core = init_core_pdb(VAR["COREDIR"]+"/"+VAR["CORE"])

    print("Running PCA for ligand1...")
    xyz_pillars1 = get_ligand_pill(xyz_lig1, anchor_ndx1)
    if two_lig:
        print("Running PCA for ligand2...")
        xyz_pillars2 = get_ligand_pill(xyz_lig2, anchor_ndx2)
    else:
        xyz_pillars2 = []

    N_S = len(names_core[names_core=='ST'])

    xyz_anchors1, xyz_anchors2 = assign_morph(xyz_core, names_core, float(VAR["LIG1_FRAC"]), int(VAR["RSEED"]), VAR["MORPHOLOGY"], int(VAR["STRIPES"]))

    xyz_stones1 = get_stones(xyz_anchors1, xyz_pillars1)
    if two_lig:
        xyz_stones2 = get_stones(xyz_anchors2, xyz_pillars2)
    else:
        xyz_stones2 = []

    print("Coating nanoparticle...")
    xyz_coated_NP, names_coated_NP, res_coated_NP = coat_NP(xyz_core, names_core, xyz_lig1, names_lig1, xyz_pillars1, xyz_stones1, xyz_lig2, names_lig2, xyz_pillars2, xyz_stones2, res_lig1, res_lig2)

    print("Writing pdb of the coated nanoparticle...")
    print_NP_pdb(xyz_coated_NP, names_coated_NP, res_coated_NP, xyz_anchors1, xyz_anchors2, xyz_lig1, xyz_lig2, "TMP/"+VAR["NAME"]+".pdb")

    ################################################################

    print("Running parmchk2 for ligand1...")
    os.system("parmchk2 -i {} -f mol2 -o {} -a y".format("TMP/"+VAR["LIG1_FILE"], "TMP/"+VAR["LIG1_FILE"][:-5]+".frcmod"))
    if two_lig:
        print("Running parmchk2 for ligand2...")
        os.system("parmchk2 -i {} -f mol2 -o {} -a y".format("TMP/"+VAR["LIG2_FILE"], "TMP/"+VAR["LIG2_FILE"][:-5]+".frcmod"))

    print("Writing tleap input file...")
    write_leap(VAR, "TMP/TLeap.in", two_lig)
    print("Running tleap...")
    os.system("tleap -sf TMP/TLeap.in > TMP/TLeap.log")

    print("Running acpype...")
    os.system("python {}/acpype.py -p {} -x {} -r".format(VAR["DEPENDS"], "TMP/"+VAR["NAME"]+".prmtop", "TMP/"+VAR["NAME"]+".inpcrd -b " + VAR["NAME"] + " -c user > acpype.log"))
    os.system("mv {}_GMX.top {}.top".format(VAR["NAME"], "TMP/"+VAR["NAME"]))
    os.system("mv {}_GMX.gro {}.gro".format(VAR["NAME"], "TMP/"+VAR["NAME"]))

    ##############################Staples########################
    print("Reading gro file of the coated nanoparticle...")
    xyz_sys, names_sys = load_gro("TMP/"+VAR["NAME"]+".gro")

    print("Reading top file of the unlinked nanoparticle...")
    types_sys, res_sys = load_top("TMP/"+VAR["NAME"]+".top")

    print("Looking for anchoring carbon atoms for ligand1 and their hydrogens if applicable...")
    ndx_C1, ndx_H1 = get_ndxs(xyz_sys, types_sys, names_sys, res_sys, name_anchor1, res_anchor1)
    if two_lig:
        print("Looking for anchoring carbon atoms for ligand2 and their hydrogens if applicable...")
        ndx_C2, ndx_H2 = get_ndxs(xyz_sys, types_sys, names_sys, res_sys, name_anchor2, res_anchor2)
    else:
        ndx_C2, ndx_H2 = [], []

    blocks = make_blocks(xyz_core, names_core, xyz_sys, ndx_C1, ndx_H1)
    if two_lig:
        blocks = blocks + make_blocks(xyz_core, names_core, xyz_sys, ndx_C2, ndx_H2)

    print("Writing bonds parameters...")
    write_bonds(blocks, "TMP/bonds.top", xyz_sys, names_sys)
    print("Writing angles parameters...")
    write_angles(blocks, "TMP/angles.top", xyz_sys, names_sys, res_core)
    print("Writing final topology file...")
    write_topology("TMP/"+VAR["NAME"]+".top", "TMP/bonds.top", "TMP/angles.top")
    #############################################################

    atexit.unregister(cleanup_error)

    print("Copying final files...")
    os.mkdir(VAR["NAME"])
    copy = [VAR["LIG1_FILE"], VAR["NAME"]+".pdb", VAR["NAME"]+".top", VAR["NAME"]+".gro"]
    if two_lig:
        copy.append(VAR["LIG2_FILE"])
    for i in copy:
        shutil.copyfile("TMP/"+i, VAR["NAME"]+"/"+i)

    atexit.register(cleanup_normal, VAR, log)

NanoModeler(NAME="test",
    LIG1_FILE="LIG2.mol2",
    CAP1="N",
    LIG1_FRAC="1.0",
    MORPHOLOGY="random",
    RSEED="666",
    STRIPES="1",
    LIG2_FILE="XXX.mol2",
    CAP2="N",
    CORE="au144SR60_NM.pdb",
    COREDIR="./CORES",
    DEPENDS="./DEPENDENCIES")
