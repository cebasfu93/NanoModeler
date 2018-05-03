import numpy as np
import sys
from scipy.spatial import distance
import subunits
import collections

signature = " NanoModeler"

def load_gro(gro_fname):
    #Loads the gro file written by NP_builder.py and return the coordinates and names of the atoms in the system
    gro_file = np.genfromtxt(gro_fname, dtype='str', skip_header=2, skip_footer=1)
    xyz = []
    names = []
    resids = []
    for line in gro_file:
        xyz.append(line[4:7])
        names.append(line[2])
        resids.append(line[1])
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

def get_gro_ndx(names_array, search_object):
    #Returns he indexes of the atoms with a given name
    ndx = []
    for i in range(len(names_array)):
        if search_object == str(names_array[i]):
            ndx.append(i)
    return np.array(ndx)

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
        ndx_H = np.zeros((N_anch, 2))
        for i in range(N_anch):
            D_sort = np.argsort(D_C_all[i,:])
            ndx_H[i,0] = D_sort[1]
            ndx_H[i,1] = D_sort[2]
        if not np.all(types_sys_func[ndx_H.flatten()]=="HC"):
            sys.exit("There are no parameters for the hydrogen atoms next to the anchor, or the atoms next to the anchor are not hydrogen atoms. The hydrogens next to CT anchor must be HC...")
    elif type_anchor_func[0]=="CA":
        ndx_H = np.array([])
    return ndx_C, ndx_H

def write_bonds(staples_list, fname, xyz_sys_func, names_sys_func, types_sys_func):
    #Goes through the S atoms of every staple, looks for the closest Au and C atoms, and assign bond parameters
    bonds = open(fname, 'w')
    func_type = str(1)
    for i in range(len(staples_list)):
        s = staples_list[i]
        D_S_Au = distance.cdist(xyz_sys_func[s.S], xyz_sys_func[s.Au])
        D_S_C = distance.cdist(xyz_sys_func[s.S], xyz_sys_func[s.C])
        for j in range(len(s.S)):
            #S - Au bonds
            near_Au = s.Au[D_S_Au[j].argsort()[0:2]]
            cons = 62730
            for k in range(len(near_Au)):
                if near_Au[k] in s.Au_l:
                    zero = 0.233
                else:
                    zero = 0.241
                bonds.write(str(s.S[j]+1).rjust(6)+str(near_Au[k]+1).rjust(7)+str(func_type).rjust(4)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ;\t"+names_sys_func[s.S[j]]+" - "+names_sys_func[near_Au[k]]+signature+"\n")

            #S - C bonds
            near_C = s.C[D_S_C[j].argsort()[0]]
            if types_sys_func[near_C] == 'CT':
                cons = 99113.0
                zero = 0.184
                bonds.write(str(s.S[j]+1).rjust(6)+str(near_C+1).rjust(7)+str(func_type).rjust(4)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ;\t"+names_sys_func[s.S[j]]+" - "+names_sys_func[near_C]+signature+"\n")
            elif types_sys_func[near_C] == 'CA':
                cons = 198321.6
                zero = 0.175
                bonds.write(str(s.S[j]+1).rjust(6)+str(near_C+1).rjust(7)+str(func_type).rjust(4)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ;\t"+names_sys_func[s.S[j]]+" - "+names_sys_func[near_C]+" m\n")
            else:
                print("Unrecognized bond type")

def write_angles(staples_list, fname, xyz_sys_func, names_sys_func, types_sys_func):
    #Goes through every staple and wirtes the parameters for the angles involving S atoms. Then a particular case is used for the S-Aul-S bond
    angles = open(fname, 'w')
    func_type = str(1)
    for i in range(len(staples_list)):
        s = staples_list[i]
        D_S_Au = distance.cdist(xyz_sys_func[s.S], xyz_sys_func[s.Au])
        D_S_C = distance.cdist(xyz_sys_func[s.S], xyz_sys_func[s.C])
        for j in range(len(s.S)):
            #Au - S - Au angles
            near_Au = s.Au[D_S_Au[j].argsort()[0:2]]
            cons = 460.24
            if ((near_Au[0] in s.Au_l) and (near_Au[1] not in s.Au_l) or (near_Au[0] not in s.Au_l) and (near_Au[1] in s.Au_l)):
                zero = 91.3
            elif ((near_Au[0] in s.Au_l) and (near_Au[1] in s.Au_l)):
                if s.tipo == "STV":
                    cons = 1460.24
                    zero = 119.2
                elif s.tipo == "STC":
                    zero = 100.0
            else:
                print("There is an unsupported Au-S-Au bond")
            angles.write(str(near_Au[0]+1).rjust(6)+str(s.S[j]+1).rjust(7)+str(near_Au[1]+1).rjust(7)+str(func_type).rjust(7)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ;\t"+names_sys_func[near_Au[0]]+" - "+names_sys_func[s.S[j]]+" - "+names_sys_func[near_Au[1]]+signature+"\n")

            #Au - S - C angles
            near_C = s.C[D_S_C[j].argsort()[0]]
            if types_sys_func[near_C] == "CT":
                cons = 146.37
                for k in range(len(near_Au)):
                    if near_Au[k] in s.Au_l:
                        zero = 106.8
                    else:
                        zero = 111.6
                    angles.write(str(near_Au[k]+1).rjust(6)+str(s.S[j]+1).rjust(7)+str(near_C+1).rjust(7)+str(func_type).rjust(7)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ;\t"+names_sys_func[near_Au[k]]+" - "+names_sys_func[s.S[j]]+" - "+names_sys_func[near_C]+signature+"\n")

            #S - C - H angles
            D_C_H = distance.cdist([xyz_sys_func[near_C]], xyz_sys_func[s.H])
            near_H = s.H[D_C_H[0].argsort()[0:2]]
            cons = 418.40
            zero = 107.0
            for k in range(len(near_H)):
                if types_sys_func[near_H[k]] == "HC":
                    angles.write(str(s.S[j]+1).rjust(6)+str(near_C+1).rjust(7)+str(near_H[k]+1).rjust(7)+str(func_type).rjust(7)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ;\t"+names_sys_func[s.S[j]]+" - "+names_sys_func[near_C]+" - "+names_sys_func[near_H[k]]+signature+"\n")
                else:
                    print("There are no parameters for angles involving this kind of hydrogen atoms")

        #Aul - S - Aul angles
        D_Aul_S = distance.cdist(xyz_sys_func[s.Au_l], xyz_sys_func[s.S])
        cons = 460.240
        zero = 172.4
        for j in range(len(s.Au_l)):
            near_S = s.S[D_Aul_S[j].argsort()[0:2]]
            angles.write(str(near_S[0]+1).rjust(6)+str(s.Au_l[j]+1).rjust(7)+str(near_S[1]+1).rjust(7)+str(func_type).rjust(7)+"{:.4e}".format(zero).rjust(14)+"{:.4e}".format(cons).rjust(14)+" ;\t"+names_sys_func[near_S[0]]+" - "+names_sys_func[s.Au_l[j]]+" - "+names_sys_func[near_S[1]]+signature+"\n")

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
