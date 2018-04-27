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
print("Importing scipy library...")
from scipy.spatial import distance
print("Importing sklearn library...")
from sklearn.decomposition import PCA
print("Importing argparse library...")
import argparse
print("Importing os library...")
import os
print("Importing shutil library...")
import shutil
print("Importing subprocess library...")
import subprocess
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
from staples import load_gro, load_top, get_gro_ndx, get_lig_info, make_blocks, make_staples, classify_staples, write_bonds, write_angles, write_topology, staple_to_residues, print_pdb

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

shutil.copyfile(VAR["CORE"], "TMP/"+VAR["CORE"])
print("Copying ligand1 file...")
shutil.copyfile(VAR["LIG1_FILE"], "TMP/"+VAR["LIG1_FILE"])
if not VAR["LIG2_FILE"] == "XXX.mol2":
    print("Copying ligand2 file...")
    shutil.copyfile(VAR["LIG2_FILE"], "TMP/"+VAR["LIG2_FILE"])

##############################NP_builder########################

print("Initializing ligand1...")
xyz_lig1, names_lig1, anchor_ndx1 = init_lig_mol2(VAR["LIG1_FILE"])
if two_lig:
    print("Initializing ligand2...")
    xyz_lig2, names_lig2, anchor_ndx2 = init_lig_mol2("LIG2_FILE")
else:
    xyz_lig2, names_lig2, anchor_ndx2 = [], [], []

print("Initializing metallic core...")
xyz_core, names_core = init_core_pdb(VAR["CORE"])

print("Running PCA for ligand1...")
xyz_pillars1 = get_ligand_pill(xyz_lig1, anchor_ndx1)
if two_lig:
    print("Running PCA for ligand2...")
    xyz_pillars2 = get_ligand_pill(xyz_lig2, anchor_ndx2)
else:
    xyz_pillars2 = []

N_S = len(names_core[names_core=='ST'])

xyz_anchors1, xyz_anchors2 = assign_morph(xyz_core, names_core, VAR["COREANCHOR"], float(VAR["LIG1_FRAC"]), int(VAR["RSEED"]), VAR["MORPHOLOGY"])

xyz_stones1 = get_stones(xyz_anchors1, xyz_pillars1)
if two_lig:
    xyz_stones2 = get_stones(xyz_anchors2, xyz_pillars2)
else:
    xyz_stones2 = []

print("Coating nanoparticle...")
xyz_coated_NP, names_coated_NP = coat_NP(xyz_core, names_core, float(VAR["LIG1_FRAC"]), xyz_lig1, names_lig1, xyz_pillars1, xyz_stones1, xyz_lig2, names_lig2, xyz_pillars2, xyz_stones2)

print("Writing pdb of the coated nanoparticle...")
print_NP_pdb(xyz_coated_NP, names_coated_NP, xyz_anchors1, xyz_anchors2, xyz_lig1, xyz_lig2, float(VAR["LIG1_FRAC"]), VAR["LIG1_NAME"], VAR["LIG2_NAME"], "TMP/"+VAR["NAME"]+".pdb")

################################################################

print("Running parmchk2 for ligand1...")
os.system("parmchk2 -i {} -f mol2 -o {} -a y".format("TMP/"+VAR["LIG1_FILE"], "TMP/"+VAR["LIG1_FILE"][:-5]+".frcmod"))
if two_lig:
    print("Running parmchk2 for ligand2...")
    os.system("parmchk2 -i {} -f mol2 -o {} -a y".format("TMP/"+VAR["LIG2_FILE"], "TMP/"+VAR["LIG2_FILE"][:-5]+".frcmod"))

print("Writing tleap input file...")
write_leap("TMP/"+VAR["LEAPFILE"]+".in", two_lig)
print("Running tleap...")
os.system("tleap -sf {}.in > {}.log".format("TMP/"+VAR["LEAPFILE"], "TMP/"+VAR["LEAPFILE"]))

print("Running acpype...")
os.system("python {}/acpype.py -p {} -x {} -r".format(VAR["DEPENDS"], "TMP/"+VAR["NAME"]+".prmtop", "TMP/"+VAR["NAME"]+".inpcrd -b " + VAR["NAME"] + " -c user > acpype.log"))
os.system("mv {}_GMX.top {}.top".format(VAR["NAME"], "TMP/"+VAR["NAME"]))
os.system("mv {}_GMX.gro {}.gro".format(VAR["NAME"], "TMP/"+VAR["NAME"]))

##############################Staples########################
print("Reading gro file of the coated nanoparticle...")
xyz_sys, names_sys = load_gro("TMP/"+VAR["NAME"]+".gro")
print("Reading top file of the unlinked nanoparticle...")
types_sys = load_top("TMP/"+VAR["NAME"]+".top")

print("Looking for gold and sulfur atoms...")
ndx_AU = get_gro_ndx(names_sys, 'AU')
ndx_ST = get_gro_ndx(names_sys, 'ST')

print("Looking for the connecting atom of ligand1...")
name_C1, type_C1, name_H1, type_H1 = get_lig_info("TMP/"+VAR["LIG1_FILE"])
ndx_C1 = get_gro_ndx(names_sys, name_C1)
print("Looking for hydrogen atoms at the interface with ligand1...")
ndx_H1a = get_gro_ndx(names_sys, name_H1[0])
ndx_H1b = get_gro_ndx(names_sys, name_H1[1])

if two_lig:
    print("Looking for the connecting atom of ligand2...")
    name_C2, type_C2, name_H2, type_H2 = get_lig_info("TMP/"+VAR["LIG2_FILE"])
    print("Looking for hydrogen atoms at the interface with ligand2...")
    ndx_C2 = get_gro_ndx(names_sys, name_C2)
    ndx_H2a = get_gro_ndx(names_sys, name_H2[0])
    ndx_H2b = get_gro_ndx(names_sys, name_H2[1])

else:
    name_C2, type_C2, name_H2, type_H2, ndx_C2, ndx_H2a, ndx_H2b = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])

ndx_C = np.unique(np.concatenate((ndx_C1, ndx_C2)).astype('int'))
ndx_H = np.unique(np.concatenate((ndx_H1a, ndx_H1b, ndx_H2a, ndx_H2b)).astype('int'))
print("Building Au-S-Au blocks...")
blocks = make_blocks(ndx_AU, ndx_ST, ndx_C, ndx_H, xyz_sys, ndx_ST)
print("Identifying staples...")
staples = make_staples(blocks)
print("Classifying staples...")
staples = classify_staples(staples)

#residues = staple_to_residues(staples)
#print_pdb(residues, "TMP/staples.pdb')

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

print("NANOMODELER terminated normally. Que gracias.")

log.close()
shutil.copyfile("NanoModeler.log", VAR["NAME"]+"/NanoModeler.log")
os.remove("NanoModeler.log")
