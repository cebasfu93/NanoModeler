log = open("NanoModeler.log", "w")
log.write("WELCOME TO NANOMODELER\n")
log.write("Importing sys library...\n")
import sys
#sys.stdout = log
#sys.stderr = log
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
print("Importing argparse library...")
import os
print("Importing shutil library...")
import shutil
print("Importing default variables...")
from defaults import VAR, write_leap
print("Looking for folder with dependencies...")
sys.path.append(VAR["DEPENDS"])
print("Importing transformations...")
from  transformations import *
print("Importing subunits...")
import subunits
print("Importing collections...")
import collections
print("Importing rewrite_mol2...")
from rewrite_mol2 import rewrite_mol2
print("Importing NP_builder...")
from NP_builder import init_lig_mol2, init_core_pdb, get_ligand_pill, assign_morph, get_stones, coat_NP, print_NP_pdb
print("Importing staples...")
from staples import load_gro, load_top, get_ndxs, make_blocks, write_bonds, write_angles, write_topology
print("Importing checking functions...")
from check import check_mol2, check_VAR
check_VAR(VAR)

print("Importing input options...\n \n")
inp=np.genfromtxt(sys.argv[1], dtype="str")
for i in range(len(inp)):
    VAR[inp[i][0]] = inp[i][1]

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
print("Rewriting ligand1 file...")
rewrite_mol2(VAR["LIG1_FILE"], VAR["CAP1"], "TMP/"+VAR["LIG1_FILE"])
print("Checking ligand1 mol2 file...")
check_mol2("TMP/"+VAR["LIG1_FILE"])
if two_lig:
    print("Two ligands were found...")
    print("Rewriting ligand2 file...")
    rewrite_mol2(VAR["LIG2_FILE"], VAR["CAP2"], "TMP/"+VAR["LIG2_FILE"])
    print("Checking ligand2 mol2 file...")
    check_mol2("TMP/"+VAR["LIG2_FILE"])

##############################NP_builder########################

print("Initializing ligand1...")
xyz_lig1, names_lig1, anchor_ndx1, name_anchor1, res_anchor1, res_lig1 = init_lig_mol2(VAR["LIG1_FILE"], VAR["CAP1"])
if two_lig:
    print("Initializing ligand2...")
    xyz_lig2, names_lig2, anchor_ndx2, name_anchor2, res_anchor2, res_lig2 = init_lig_mol2(VAR["LIG2_FILE"], VAR["CAP2"])
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
write_leap("TMP/TLeap.in", two_lig)
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

print("Copying final files...")
os.mkdir(VAR["NAME"])
copy = [VAR["LIG1_FILE"], VAR["NAME"]+".pdb", VAR["NAME"]+".top", VAR["NAME"]+".gro"]
if two_lig:
    copy.append(VAR["LIG2_FILE"])
for i in copy:
    shutil.copyfile("TMP/"+i, VAR["NAME"]+"/"+i)

print("Cleaning...")
bye = ["ANTECHAMBER.FRCMOD", "leap.log", "md.mdp", "em.mdp", "acpype.log"]
for i in bye:
    os.remove(i)
shutil.rmtree("TMP")

print("Compressing files to output...")
print("NANOMODELER terminated normally. Que gracias.")
log.close()
shutil.copyfile("NanoModeler.log", VAR["NAME"]+"/NanoModeler.log")
os.remove("NanoModeler.log")

os.system("tar -zcvf {}.tar.gz {}".format(VAR["NAME"], VAR["NAME"]))
#shutil.rmtree(VAR["NAME"])
