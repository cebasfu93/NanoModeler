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

def get_ndxs(xyz_sys_func, types_sys_func, names_sys_func, res_sys_func, name_anchor_func, res_anchor_func, log):
    ndx_C = np.where(np.logical_and(names_sys_func==name_anchor_func, res_sys_func==res_anchor_func))[0]
    type_anchor_func = types_sys_func[ndx_C]

    log += "Checking if the assigned atom type for the anchors is supported...\n"
    if not (type_anchor_func[0]=="CT" or type_anchor_func[0]=="CA"):
        sys.exit("One of the anchors was assigned an unsupported atom type. Those supported are CT and CA.")
    N_anch = len(ndx_C)
    return ndx_C

def make_blocks(xyz_sys_func, names_sys_func, names_core_func, res_core_func, ndx_C_func):
    core_Au = np.where(np.logical_or(names_core_func=="AUS", names_core_func=="AUL"))[0]
    ndx_S = np.where(names_sys_func=="ST")[0]
    ndx_Au = np.where(np.logical_or(names_sys_func=="AUS", names_sys_func=="AUL"))[0]
    D_S_C = distance.cdist(xyz_sys_func[ndx_S], xyz_sys_func[ndx_C_func])
    D_S_Au = distance.cdist(xyz_sys_func[ndx_S], xyz_sys_func[ndx_Au])

    blocks = []
    for i in range(len(ndx_S)):
        now_S = ndx_S[i]
        now_C = ndx_C_func[np.argsort(D_S_C[i])[0]]
        sort_S_Au = np.argsort(D_S_Au[i])[0:2]
        old_Au = core_Au[sort_S_Au]
        now_Au = ndx_Au[sort_S_Au]
        tipos_Au = names_sys_func[now_Au]
        blocks.append(subunits.Block(ndx_S=now_S, ndx_Au=now_Au, ndx_C=now_C, types_Au=tipos_Au, staple=res_core_func[old_Au]))
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
    bonds.close()

def write_angles(blocks_list, fname, xyz_sys_func, names_sys_func):
    #Goes through every staple and wirtes the parameters for the angles involving S atoms. Then a particular case is used for the S-Aul-S bond
    angles = open(fname, 'w')
    func_type = str(1)
    Au_taken = []
    for i in range(len(blocks_list)):
        b = blocks_list[i]

        #AuL - S - AuL
        if np.all(b.typesAu == "AUL"):
            if np.all(b.staple=="STV"):
                cons = 1460.24
                zero = 119.2
            elif np.all(b.staple=="STC"):
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
