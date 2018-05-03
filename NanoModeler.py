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
print("Importing NP_builder...")
from NP_builder import init_lig_mol2, init_core_pdb, get_ligand_pill, assign_morph, get_stones, coat_NP, print_NP_pdb
print("Importing staples...")
from staples import load_gro, load_top, get_ndxs, write_bonds, write_angles, write_topology
print("Importing checking functions...")
from check import check_mol2, check_VAR
check_VAR(VAR)

print("Importing input options...\n \n")
inp=np.genfromtxt(sys.argv[1], dtype="str")
for i in range(len(inp)):
    VAR[inp[i][0]] = inp[i][1]
two_lig = float(VAR["LIG1_FRAC"]) < 1.0
print("Imported options:")
for i in VAR:
    print(i.ljust(20) + VAR[i].ljust(50))

print("\n\nCreating folder...")
os.mkdir("TMP")

print("Copying ligand1 file...")
shutil.copyfile(VAR["LIG1_FILE"], "TMP/"+VAR["LIG1_FILE"])
print("Checking ligand1 mol2 file...")
check_mol2("TMP/"+VAR["LIG1_FILE"])
if not VAR["LIG2_FILE"] == "XXX.mol2":
    print("Copying ligand2 file...")
    shutil.copyfile(VAR["LIG2_FILE"], "TMP/"+VAR["LIG2_FILE"])
    print("Checking ligand2 mol2 file...")
    check_mol2("TMP/"+VAR["LIG2_FILE"])
##############################NP_builder########################

print("Initializing ligand1...")
xyz_lig1, names_lig1, anchor_ndx1, name_anchor1, res_anchor1, res_lig1 = init_lig_mol2(VAR["LIG1_FILE"])
if two_lig:
    print("Initializing ligand2...")
    xyz_lig2, names_lig2, anchor_ndx2, name_anchor2, res_anchor2, res_lig2 = init_lig_mol2("LIG2_FILE")
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

xyz_anchors1, xyz_anchors2 = assign_morph(xyz_core, names_core, float(VAR["LIG1_FRAC"]), int(VAR["RSEED"]), VAR["MORPHOLOGY"])

xyz_stones1 = get_stones(xyz_anchors1, xyz_pillars1)
if two_lig:
    xyz_stones2 = get_stones(xyz_anchors2, xyz_pillars2)
else:
    xyz_stones2 = []

print("Coating nanoparticle...")
xyz_coated_NP, names_coated_NP, res_coated_NP = coat_NP(xyz_core, names_core, float(VAR["LIG1_FRAC"]), xyz_lig1, names_lig1, xyz_pillars1, xyz_stones1, xyz_lig2, names_lig2, xyz_pillars2, xyz_stones2, res_lig1, res_lig2)

print("Writing pdb of the coated nanoparticle...")
print_NP_pdb(xyz_coated_NP, names_coated_NP, res_coated_NP, xyz_anchors1, xyz_anchors2, xyz_lig1, xyz_lig2, float(VAR["LIG1_FRAC"]), "TMP/"+VAR["NAME"]+".pdb")

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

ndx_C1, ndx_H1 = get_ndxs(xyz_sys, types_sys, names_sys, res_sys, name_anchor1, res_anchor1)
if two_lig:
    ndx_C2, ndx_H2 = get_ndxs(xyz_sys, types_sys, names_sys, res_sys, name_anchor2, res_anchor2)
else:
    ndx_C2, ndx_H2 = [], []



print("Writing final file...")
write_bonds(staples, "TMP/bonds.top", xyz_sys, names_sys, types_sys)
write_angles(staples, "TMP/angles.top", xyz_sys, names_sys, types_sys)
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
