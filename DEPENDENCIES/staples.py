import numpy as np
import sys
from scipy.spatial import distance
import subunits

signature = " NanoModeler"

def load_gro(gro_fname):
    #Loads the gro file written by NP_builder.py and return the coordinates and names of the atoms in the system
    gro_file = np.genfromtxt(gro_fname, dtype='str', skip_header=2, skip_footer=1, delimiter="\n")
    xyz = []
    names = []
    resids = []
    for line in gro_file:
        xyz.append([line[-24:-16],line[-16:-8], line[-8:]])
        names.append(line[-34:-29].strip())
        resids.append(line[-39:-34].strip())
    return np.array(xyz).astype('float'), np.array(names)

def load_top(top_fname):
    #Loads the topology of the system writen by acpype.py and returns the types of the atoms
    top_file = np.genfromtxt(top_fname, dtype='str', delimiter = "\n")
    types = []
    residues = []
    for i in range(len(top_file)):
        if "[ atoms ]" in top_file[i]:
            ini = i + 2
        elif "[ bonds ]" in top_file[i]:
            fin = i
    for i in range(ini, fin):
        types.append(top_file[i].split()[1])
        residues.append(top_file[i].split()[3])
    return np.array(types), np.array(residues)

def get_ndxs(xyz_sys_func, types_sys_func, names_sys_func, res_sys_func, name_anchor_func, res_anchor_func):
    ndx_C = np.where(np.logical_and(names_sys_func==name_anchor_func, res_sys_func==res_anchor_func))[0]
    type_anchor_func = types_sys_func[ndx_C]

    print("Checking if the assigned atom type for the anchors is supported...")
    if not (type_anchor_func[0]=="CT" or type_anchor_func[0]=="CA"):
        sys.exit("One of the anchors was assigned an unsupported atom type. Those supported are CT and CA.")
    N_anch = len(ndx_C)

    D_C_all = distance.cdist(xyz_sys_func[ndx_C], xyz_sys_func)
    if type_anchor_func[0]=="CT":
        print("Looking for closest hydrogen atoms to anchors...")
        ndx_H = np.argsort(D_C_all)[:,1:3]
        if not np.all(types_sys_func[ndx_H.flatten().astype("int")]=="HC"):
            sys.exit("There are no parameters for the hydrogen atoms next to the anchor, or the atoms next to the anchor are not hydrogen atoms. The hydrogens next to CT anchor must be HC...")
    elif type_anchor_func[0]=="CA":
        ndx_H = np.array([])
    return ndx_C, ndx_H

def make_blocks(xyz_core_func, names_core_func, xyz_sys_func, ndx_C_func, ndx_H_func):
    all_Au = np.append(np.where(names_core_func=="AUS")[0], np.where(names_core_func=="AUL")[0])
    all_S = np.where(names_core_func=="ST")[0]
    blocks = []
    D_C_CORE = distance.cdist(xyz_sys_func[ndx_C_func], xyz_core_func)
    for i in range(len(ndx_C_func)):
        ndx_S = all_S[np.argsort(D_C_CORE[i, all_S])[0]]
        D_S_Au = distance.cdist([xyz_sys_func[ndx_S]], xyz_sys_func[all_Au])[0]
        ndx_Au = all_Au[np.argsort(D_S_Au)[0:2]]
        tipos_Au = names_core_func[ndx_Au]
        if i ==0:
            print(ndx_C_func[i], ndx_S, ndx_Au, ndx_H_func[i])
        if np.any(np.logical_and(tipos_Au != "AUS", tipos_Au != "AUL")):
            sys.exit("There was a problem recognizing if some gold atoms where type AUL or AUS.")
        if len(ndx_H_func)!=0:
            blocks.append(subunits.Block(ndx_S=ndx_S, ndx_Au=ndx_Au, ndx_C=ndx_C_func[i], ndx_H=ndx_H_func[i], types_Au=tipos_Au))
        else:
            blocks.append(subunits.Block(ndx_S=ndx_S, ndx_Au=ndx_Au, ndx_C=ndx_C_func[i], ndx_H=[], types_Au=tipos_Au))
    return blocks

def write_bonds(blocks_list, fname, xyz_sys_func, names_sys_func):
    #Goes through the S atoms of every staple, looks for the closest Au and C atoms, and assign bond parameters
    bonds = open(fname, 'w')
    func_type = str(1)
    for i in range(len(blocks_list)):
        b = blocks_list[i]
        #S - Au bonds
        cons = 62730
        for j in range(2):
            if b.typesAu[j]=="AUL":
                zero = 0.233
            elif b.typesAu[j]=="AUS":
                zero = 0.241
            bonds.write(str(b.S+1).rjust(6)+str(b.Au[j]+1).rjust(7)+str(func_type).rjust(4)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ;\t"+names_sys_func[b.S]+" - "+names_sys_func[b.Au[j]]+signature+"\n")

        #S - C bonds
        if len(b.H) == 2:
            cons = 99113.0
            zero = 0.184
        elif len(b.H) == 0:
            cons = 198321.6
            zero = 0.175
        else:
            sys.exit("There is something wrong with the anchors' hydrogen indexing.")
        bonds.write(str(b.S+1).rjust(6)+str(b.C+1).rjust(7)+str(func_type).rjust(4)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ;\t"+names_sys_func[b.S]+" - "+names_sys_func[b.C]+signature+"\n")
    bonds.close()

def write_angles(blocks_list, fname, xyz_sys_func, names_sys_func, res_core_func):
    #Goes through every staple and wirtes the parameters for the angles involving S atoms. Then a particular case is used for the S-Aul-S bond
    angles = open(fname, 'w')
    func_type = str(1)
    Au_taken = []
    for i in range(len(blocks_list)):
        b = blocks_list[i]

        #AuL - S - AuL
        if np.all(b.typesAu == "AUL"):
            if np.all(res_core_func[b.Au]=="STV"):
                cons = 1460.24
                zero = 119.2
            elif np.all(res_core_func[b.Au]=="STC"):
                cons = 460.24
                zero = 100.0
            else:
                sys.exit("There was a problem recognizing the staple type when trying to write an Au - S - Au angle.")
            angles.write(str(b.Au[0]+1).rjust(6)+str(b.S+1).rjust(7)+str(b.Au[1]+1).rjust(7)+str(func_type).rjust(7)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ;\t"+names_sys_func[b.Au[0]]+" - "+names_sys_func[b.S]+" - "+names_sys_func[b.Au[1]]+signature+"\n")

        #AuL - S - AuS
        if np.any(b.typesAu == "AUS"):
            cons = 460.240
            zero = 91.3
            angles.write(str(b.Au[0]+1).rjust(6)+str(b.S+1).rjust(7)+str(b.Au[1]+1).rjust(7)+str(func_type).rjust(7)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ;\t"+names_sys_func[b.Au[0]]+" - "+names_sys_func[b.S]+" - "+names_sys_func[b.Au[1]]+signature+"\n")

        #Au - S -  C
        for j in range(2):
            cons = 146.370
            if b.typesAu[j] == "AUL":
                zero = 106.8
            elif b.typesAu[j] == "AUS":
                zero = 111.6
            angles.write(str(b.Au[j]+1).rjust(6)+str(b.S+1).rjust(7)+str(b.C+1).rjust(7)+str(func_type).rjust(7)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ;\t"+names_sys_func[b.Au[j]]+" - "+names_sys_func[b.S]+" - "+names_sys_func[b.C]+signature+"\n")

        #S - C - H
        if len(b.H) == 2:
            cons = 418.40
            zero = 107.0
            angles.write(str(b.S+1).rjust(6)+str(b.C+1).rjust(7)+str(b.H[0]+1).rjust(7)+str(func_type).rjust(7)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ;\t"+names_sys_func[b.S]+" - "+names_sys_func[b.C]+" - "+names_sys_func[b.H[0]]+signature+"\n")
            angles.write(str(b.S+1).rjust(6)+str(b.C+1).rjust(7)+str(b.H[1]+1).rjust(7)+str(func_type).rjust(7)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ;\t"+names_sys_func[b.S]+" - "+names_sys_func[b.C]+" - "+names_sys_func[b.H[1]]+signature+"\n")

        #S - AuL - S
        all_S = np.where(names_sys_func=="ST")[0]
        for j in range(2):
            if b.typesAu[j] == "AUL" and b.Au[j] not in Au_taken:
                Au_taken.append(b.Au[j])
                cons = 460.240
                zero = 172.4

                D_AuL_S = distance.cdist([xyz_sys_func[b.Au[j]]], xyz_sys_func[all_S])
                near_S = all_S[np.argsort(D_AuL_S[0])[0:2]]
                angles.write(str(near_S[0]+1).rjust(6)+str(b.Au[j]+1).rjust(7)+str(near_S[1]+1).rjust(7)+str(func_type).rjust(7)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ;\t"+names_sys_func[near_S[0]]+" - "+names_sys_func[b.Au[j]]+" - "+names_sys_func[near_S[1]]+signature+"\n")

    angles.close()

def write_topology(fname, bonds, angles):
    #Copies the previous topology file writen by acpype.py and inserts the new bond and angles at the beggining of their respective sections
    top_file = open(fname, "r")
    cont_top_file = top_file.readlines()
    top_file.close()

    final_top=open(fname, "w")
    final_top.close()
    final_top=open(fname, "a")

    for i in range(len(cont_top_file)):
        final_top.writelines(cont_top_file[i])
        if ";   ai     aj funct   r             k" in cont_top_file[i]:
            bonds_file=open(bonds,"r")
            bonds_contents=bonds_file.readlines()
            bonds_file.close()
            for j in range(len(bonds_contents)):
                final_top.writelines(bonds_contents[j])
        if ";   ai     aj     ak    funct   theta         cth" in cont_top_file[i]:
            angles_file=open(angles,"r")
            angles_contents=angles_file.readlines()
            angles_file.close()
            for j in range(len(angles_contents)):
                final_top.writelines(angles_contents[j])
