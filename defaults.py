VAR = {
"NAME": "test",             #Name of the project
"LIG1_FILE" : "LF1.mol2",   #Name of the mol2 of ligand1 (must be in the working directory)
"CAP1" : "N",		    #Atom numbers (as in the mol2, likely to start in 1) to remove. Numbers separated by commas	

"LIG2_FILE" : "LF2.mol2",   #Name of the mol2 of ligand2 (must be in the working directory)
"CAP2" : "N",
"MORPHOLOGY" : "random",    #Morphology to distribute ligands1 and 2. random, janus, and stripe are allowed
"LIG1_FRAC" : 1.0,          #Fraction of ligand1 to place
"RSEED" : 666,              #Random seed for random morphology
"STRIPES" : 1,              #Number of stripes for stripe morphology. It will start (bottom up with ligand 1)

"CORE" : "au144SR60_NM.pdb",

"COREDIR" : "/DATA/SoftwareSFU/IN-HOUSE/NanoModeler/CORES",
"DEPENDS" : "/DATA/SoftwareSFU/IN-HOUSE/NanoModeler/DEPENDENCIES"
}

import numpy as np
def read_resname(lig_fname):
    mol2 = np.genfromtxt(lig_fname, delimiter="\n", dtype='str')
    for i in range(len(mol2)):
        if "@<TRIPOS>ATOM" in mol2[i]:
            resname = mol2[i+1].split()[7]
    return resname

def write_leap(fname, two_lig_func):
    msj = "source leaprc.gaff \n\n"
    msj += "loadamberparams " + "TMP/"+VAR["LIG1_FILE"][:-5]+".frcmod\n"

    msj += read_resname(VAR["LIG1_FILE"]) + " = loadmol3 " + "TMP/"+VAR["LIG1_FILE"]+"\n"
    msj += "check " + read_resname(VAR["LIG1_FILE"]) + "\n"
    msj += "saveoff " + read_resname(VAR["LIG1_FILE"]) + " " + "TMP/"+VAR["LIG1_FILE"][:-5]+".lib\n\n"
    if two_lig_func:
        msj += "loadamberparams " + "TMP/"+VAR["LIG2_FILE"][:-5]+".frcmod\n"
        msj += read_resname(VAR["LIG2_FILE"]) + " = loadmol3 " + "TMP/"+VAR["LIG2_FILE"]+"\n"
        msj += "check " + read_resname(VAR["LIG2_FILE"]) + "\n"
        msj += "saveoff " + read_resname(VAR["LIG2_FILE"]) + " " + "TMP/"+VAR["LIG2_FILE"][:-5]+".lib\n\n"

    msj += "loadamberparams " + VAR["DEPENDS"]+"/AU.frcmod\n"
    msj += "loadamberparams " + VAR["DEPENDS"]+"/ST.frcmod\n"
    msj += "AU = loadmol3 " + VAR["DEPENDS"]+"/AU.mol2\n"
    msj += "ST = loadmol3 " + VAR["DEPENDS"]+"/ST.mol2\n\n"

    msj += "loadoff " + "TMP/"+VAR["LIG1_FILE"][:-5]+".lib\n"
    if two_lig_func:
        msj += "loadoff " + "TMP/"+VAR["LIG2_FILE"][:-5]+".lib\n"

    msj += VAR["NAME"] + " = loadpdb " + "TMP/"+VAR["NAME"]+".pdb \n"
    msj += "saveamberparm " +  VAR["NAME"] + " " + "TMP/"+VAR["NAME"]+".prmtop" + " " + "TMP/"+VAR["NAME"]+".inpcrd \n"
    msj += "quit"
    out = open(fname, "w")
    out.write(msj)
    out.close()
