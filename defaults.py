VAR = {
"NAME": "test",             #Name of the project
"LIG1_FILE" : "LF1.mol2",   #Name of the mol2 of ligand1 (must be in the working directory)
"CAP1" : []
"LIG1_FRAC" : 1.0,          #Fraction of ligand1 to place

"LIG2_FILE" : "LF2.mol2",   #Name of the mol2 of ligand2 (must be in the working directory)
"CAP2" : []
"MORPHOLOGY" : "random",    #Morphology to distribute ligands1 and 2. random and janus are allowed
"RSEED" : 666,              #Random seed for random morphology

"CORE" : "au144SR60_NM.pdb",

"COREDIR" : "/DATA/SoftwareSFU/IN-HOUSE/NanoModeler/CORES",
"DEPENDS" : "/DATA/SoftwareSFU/IN-HOUSE/NanoModeler/DEPENDENCIES"
}

def write_leap(fname, two_lig_func):
    msj = "source leaprc.gaff \n"
    msj += "loadamberparams " + "TMP/"+VAR["LIG1_FILE"][:-5]+".frcmod\n"
    msj += VAR["LIG1_FILE"][:-5] + " = loadmol3 " + "TMP/"+VAR["LIG1_FILE"]+"\n"
    msj += "check " + VAR["LIG1_FILE"][:-5] + "\n"
    msj += "saveoff " + VAR["LIG1_FILE"][:-5] + " " + "TMP/"+VAR["LIG1_FILE"][:-5]+".lib\n"
    if two_lig_func:
        msj += "loadamberparams " + "TMP/"+VAR["LIG2_FILE"][:-5]+".frcmod\n"
        msj += VAR["LIG2_FILE"][:-5] + " = loadmol3 " + "TMP/"+VAR["LIG2_FILE"]+"\n"
        msj += "check " + VAR["LIG2_FILE"][:-5] + "\n"
        msj += "saveoff " + VAR["LIG2_FILE"][:-5] + " " + "TMP/"+VAR["LIG2_FILE"][:-5]+".lib\n"

    msj += "loadamberparams " + VAR["DEPENDS"]+"/AU.frcmod\n"
    msj += "loadamberparams " + VAR["DEPENDS"]+"/ST.frcmod\n"
    msj += "AU = loadmol3 " + VAR["DEPENDS"]+"/AU.mol2\n"
    msj += "ST = loadmol3 " + VAR["DEPENDS"]+"/ST.mol2\n"
    msj += "loadoff " + "TMP/"+VAR["LIG1_FILE"][:-5]+".lib\n"
    msj += "loadamberparams "+ "TMP/"+VAR["LIG1_FILE"][:-5]+".frcmod\n"
    msj += VAR["LIG1_FILE"][:-5] + " = loadmol3 " + "TMP/"+VAR["LIG1_FILE"]+"\n"
    if two_lig_func:
        msj += "loadoff " + "TMP/"+VAR["LIG2_FILE"][:-5]+".lib\n"
        msj += "loadamberparams "+ "TMP/"+VAR["LIG2_FILE"][:-5]+".frcmod\n"
        msj += VAR["LIG1_FILE"][:-5] + " = loadmol3 " + "TMP/"+VAR["LIG2_FILE"]+"\n"
    msj += VAR["NAME"] + " = loadpdb " + "TMP/"+VAR["NAME"]+".pdb \n"
    msj += "saveamberparm " +  VAR["NAME"] + " " + "TMP/"+VAR["NAME"]+".prmtop" + " " + "TMP/"+VAR["NAME"]+".inpcrd \n"
    msj += "quit"
    out = open(fname, "w")
    out.write(msj)
    out.close()
