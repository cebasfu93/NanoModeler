import numpy as np
from scipy import optimize
from scipy.spatial import distance
import math
import subunits

def calc_angle(a, b, c):
    #Calculates the angle formed by a-b-c
    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)*180/math.pi
    return angle

def center(coords):
    coords = np.array(coords, dtype='float')
    COM = np.mean(coords, axis=0)
    coords = coords - COM
    return coords

def penalty(p, S_xyz, Au_xyz):
    D_S = distance.cdist([p], [S_xyz])
    angles=[]
    for i in range(2):
        angles.append(angle(p, S_xyz, Au_xyz[i]))
    angles = np.array(angles)
    #print(angles)
    angles = np.abs((angles-109.2)/109.2)
    D_S = np.abs((D_S-1.8)/1.8)
    #print(angles, np.sum(angles), D_S)
    return float(np.sum(angles)+D_S)

def read_Au25SR18(fname):
    gro = np.genfromtxt(fname, dtype='str', delimiter="\n", skip_header=2, skip_footer=1)
    xyz = []
    names= []
    for i in range(len(gro)):
        if "Au" in gro[i]:
            xyz.append([gro[i][-24:-16],gro[i][-16:-8], gro[i][-8:]])
            names.append('AU')
    for i in range(len(gro)):
        if "S" in gro[i]:
            xyz.append([gro[i][-24:-16],gro[i][-16:-8], gro[i][-8:]])
            names.append('ST')
    for i in range(len(gro)):
        if "CT1" in gro[i]:
            xyz.append([gro[i][-24:-16],gro[i][-16:-8], gro[i][-8:]])
            names.append('C')
    xyz = center(xyz)*10
    names = np.array(names, dtype='str')
    return xyz, names

def read_Au68SR34(fname):
    fxyz = np.genfromtxt(fname, skip_header=1, dtype="str")
    xyz = []
    names = []
    N_xyz = len(fxyz)
    for i in range(N_xyz):
        if fxyz[i,0]=="Au":
            xyz.append(fxyz[i,1:])
            names.append("AU")
        elif fxyz[i,0]=="S":
            xyz.append(fxyz[i,1:])
            names.append("ST")
        elif fxyz[i,0]=="H":
            xyz.append(fxyz[i,1:])
            names.append("C")

    xyz = center(xyz)
    names = np.array(names)
    return xyz, names

def read_Au102SR44(fname):
    gro = np.genfromtxt(fname, dtype='str', delimiter="\n", skip_header=2, skip_footer=1)
    xyz = []
    names= []
    for i in range(len(gro)):
        if "Au" in gro[i]:
            xyz.append([gro[i][-24:-16],gro[i][-16:-8], gro[i][-8:]])
            names.append('AU')
    for i in range(len(gro)):
        if "S" in gro[i]:
            xyz.append([gro[i][-24:-16],gro[i][-16:-8], gro[i][-8:]])
            names.append('ST')
    for i in range(len(gro)):
        if "CZ" in gro[i]:
            xyz.append([gro[i][-24:-16],gro[i][-16:-8], gro[i][-8:]])
            names.append('C')
    xyz = center(xyz)*10
    names = np.array(names, dtype='str')
    return xyz, names

def read_Au144SR60(fname):
    pdb=np.genfromtxt(fname, delimiter='\n', dtype=str)
    N_pdb = len(pdb)
    xyz = []
    names = []
    for i in range(N_pdb):
        if "ATOM" in pdb[i]:
            at_act = pdb[i].split()
            if at_act[2] == 'Au':
                xyz.append(at_act[5:8])
                names.append('AU')
            elif at_act[2] == 'S':
                xyz.append(at_act[5:8])
                names.append('ST')
            elif at_act[2] == 'C':
                xyz.append(at_act[5:8])
                names.append('C')

    xyz = center(xyz)
    names = np.array(names, dtype='str')
    return xyz, names

def read_Au314SH96(fname):
    fxyz = np.genfromtxt(fname, skip_header=1, dtype=str)
    xyz = []
    names = []
    N_xyz = len(fxyz)
    for i in range(N_xyz):
        if fxyz[i,0] == "Au":
            xyz.append(fxyz[i,1:])
            names.append("AU")
        elif fxyz[i,0] == "S":
            xyz.append(fxyz[i,1:])
            names.append("ST")
        #elif fxyz[i,0] == "H":
            #xyz.append(fxyz[i,1:])
            #names.append("C")
    xyz = center(xyz)
    for i in range(len(xyz)):
        if names[i] == "ST":
            norm = np.linalg.norm(xyz[i])
            xyz = np.vstack((xyz, xyz[i]*(norm+1.83)/norm))
            names.append("C")

    xyz = center(xyz)
    names = np.array(names)
    return xyz, names

def classify_staples(xyz_sys, names_sys):

    ndx_AU = np.where(names_sys=='AU')[0]
    ndx_ST = np.where(names_sys=='ST')[0]
    ndx_C = np.where(names_sys=='C')[0]

    N_ST = len(ndx_ST)

    blocks = []
    D_ST_AU = distance.cdist(xyz_sys[ndx_ST],xyz_sys[ndx_AU])
    D_ST_C = distance.cdist(xyz_sys[ndx_ST],xyz_sys[ndx_C])
    for i in range(N_ST):
        ndx_au = np.argsort(D_ST_AU[i])[0:2]
        ndx_c = np.argsort(D_ST_C[i])[0]
        blocks.append(subunits.Block(ndx_S=ndx_ST[i], ndx_Au=ndx_AU[ndx_au], ndx_C=ndx_C[ndx_c]))


    N_blocks = len(blocks)
    ganchos = []
    for i in range(N_blocks):
        taken = False
        for j in range(len(ganchos)):
            if any(x in ganchos[j].Au for x in blocks[i].Au):
                ganchos[j].add(blocks[i].S, blocks[i].Au, blocks[i].C)
                taken = True
        if not taken:
            ganchos.append(subunits.Gancho(ndx_S=blocks[i].S, ndx_Au=blocks[i].Au, ndx_C=blocks[i].C))

    N_ganchos = len(ganchos)

    new_ganchos = []

    for i in range(N_ganchos):
        taken = False
        for j in range(len(new_ganchos)):
            for k in range(len(new_ganchos[j].Au)):
                if new_ganchos[j].Au[k] in ganchos[i].Au:
                    taken = True
                    old = j
        print(len(new_ganchos), len(ganchos))
        if not taken:
            new_ganchos.append(ganchos[i])
        elif taken:
            new_ganchos[old].add(ganchos[i].S, ganchos[i].Au, ganchos[i].C)
    N_ganchos = len(new_ganchos)

    staples = []
    for i in range(N_ganchos):
        staples.append(subunits.Staple(ndx_S=new_ganchos[i].S, ndx_Au=new_ganchos[i].Au, ndx_C=new_ganchos[i].C))
        #staples.append(subunits.Staple(ndx_S=ganchos[i].S, ndx_Au=ganchos[i].Au, ndx_C=ganchos[i].C))

    #Depending in the number of sulphur and gold atoms in each staple, it clssifies it. For STC and STV it calculates the angle Aul-S-Aul and with an tolerance of +/-9 degrees, the staple is classified
    for i in range(len(staples)):
        staple_act = staples[i]
        N_S = len(staple_act.S)
        N_Au = len(staple_act.Au)
        #print(N_S, N_Au)
        if N_S == 1 and N_Au == 2:
            staple_act.change_tipo('STP')
        elif N_S == 2 and N_Au == 3:
            staple_act.change_tipo('STR')
        elif N_S == 3 and N_Au == 4:
            D_S_Au = distance.cdist(xyz_sys[staple_act.S], xyz_sys[staple_act.Au])
            for j in range(len(staple_act.S)):
                near_Au = staple_act.Au[D_S_Au[j].argsort()[0:2]]
                if np.all(np.in1d(near_Au, staple_act.Au_l)):
                    angle = calc_angle(xyz_sys[near_Au[0]], xyz_sys[staple_act.S[j]], xyz_sys[near_Au[1]])
                    print(angle)
                    #if angle <= 109.0 and angle >= 91.0:
                    if angle <= 109.0 and angle >= 87.0:
                        staple_act.change_tipo('STV')
                    elif angle <= 128.2 and angle >= 110.2:
                        staple_act.change_tipo('STC')
                    else:
                        print("One of the staples should be STC or STV but the AuL-S-AuL angle displays an odd value. Unrecognized staple")
        else:
            print("Unrecognized staple")
    print(len(staples))
    return staples

def write_pdb(xyz_sys, names_sys, staples, fname):
    pdb=open(fname, 'w')
    pdb.close()
    all_busy_Au = np.array([], dtype='int')
    N_staples = len(staples)

    N_AU = len(np.where(names_sys=="AU")[0])

    for i in range(N_staples):
        all_busy_Au = np.append(all_busy_Au, staples[i].Au)

    res = 1
    at = 1

    for i in range(len(xyz_sys)):
        if "AU" in names_sys[i] and i not in all_busy_Au:
            write_pdb_block("AU", "AU", xyz_sys[i], res, at, fname)
            at+=1
            res+=1

    for i in range(N_staples):
        st_act = staples[i]
        for j in range(len(st_act.Au_s)):
            write_pdb_block("AUS", st_act.tipo, xyz_sys[st_act.Au_s[j]], res, at, fname)
            at+=1
        for j in range(len(st_act.Au_l)):
            write_pdb_block("AUL", st_act.tipo, xyz_sys[st_act.Au_l[j]], res, at, fname)
            at+=1
        for j in range(len(st_act.S)):
            write_pdb_block("ST", st_act.tipo, xyz_sys[st_act.S[j]], res, at, fname)
            at+=1
        res+=1
    for i in range(N_staples):
        st_act = staples[i]
        for j in range(len(st_act.C)):
            write_pdb_block("C", st_act.tipo, xyz_sys[st_act.C[j]], res, at, fname)
            at+=1
        res+=1

def write_pdb_block(atname_func, res_name_func, xyz_func, resnum, atnum, out_filename):
    #Writes one line of a generic pdb file
    xyz_func=np.round(xyz_func, decimals=4)
    coords=open(out_filename, 'a')
    coords.write('ATOM'.ljust(6))
    coords.write(str(atnum).rjust(5))
    coords.write(' '+str(atname_func).ljust(4))
    coords.write(' '+str(res_name_func).ljust(3))
    coords.write('  '+str(resnum).rjust(4))
    coords.write('    '+ str(xyz_func[0]).rjust(8))
    coords.write(str(xyz_func[1]).rjust(8))
    coords.write(str(xyz_func[2]).rjust(8)+"\n")
    coords.close()

def print_xyz(coords, nombres, fname):
    fxyz = open(fname + ".xyz", "w")
    fxyz.write("{} \n\n".format(len(nombres)))
    for i in range(len(nombres)):
        fxyz.write("{}\t\t{:.3f}   {:.3f}   {:.3f}\n".format(nombres[i], coords[i,0], coords[i,1], coords[i,2]))
    fxyz.close()
"""
#AU25SR18 (61 atoms)
xyz_Au25SR18, names_Au25SR18 = read_Au25SR18("AU25SR18/au25_pet18.gro")
staples_Au25SR18 = classify_staples(xyz_Au25SR18, names_Au25SR18)
write_pdb(xyz_Au25SR18, names_Au25SR18, staples_Au25SR18, "au25SR18_NM.pdb")

#AU38SR24 (86 atoms)
xyz_Au38SR24, names_Au38SR24 = read_Au25SR18("AU38SR24/au38_pet24.gro")
staples_Au38SR24 = classify_staples(xyz_Au38SR24, names_Au38SR24)
write_pdb(xyz_Au38SR24, names_Au38SR24, staples_Au38SR24, "au38SR24_NM.pdb")

#AU68SR34 (136 atoms)
xyz_Au68SR34, names_Au68SR34 = read_Au314SH96("AU68SR34/Au68SH34-I1.xyz")
staples_Au68SR34 = classify_staples(xyz_Au68SR34, names_Au68SR34)
write_pdb(xyz_Au68SR34, names_Au68SR34, staples_Au68SR34, "au68SR34-I1_NM.pdb")

xyz_Au68SR34, names_Au68SR34 = read_Au314SH96("AU68SR34/Au68SH34-I2.xyz")
staples_Au68SR34 = classify_staples(xyz_Au68SR34, names_Au68SR34)
write_pdb(xyz_Au68SR34, names_Au68SR34, staples_Au68SR34, "au68SR34-I2_NM.pdb")

xyz_Au68SR34, names_Au68SR34 = read_Au314SH96("AU68SR34/Au68SH34-I3.xyz")
staples_Au68SR34 = classify_staples(xyz_Au68SR34, names_Au68SR34)
write_pdb(xyz_Au68SR34, names_Au68SR34, staples_Au68SR34, "au68SR34-I3_NM.pdb")
"""

xyz_Au68SR34, names_Au68SR34 = read_Au314SH96("AU68SR34/Au68SH34-I4.xyz")
staples_Au68SR34 = classify_staples(xyz_Au68SR34, names_Au68SR34)
write_pdb(xyz_Au68SR34, names_Au68SR34, staples_Au68SR34, "au68SR34-I4_NM.pdb")

#AU102SR44 (190 atoms)
#xyz_Au102SR44, names_Au102SR44 = read_Au102SR44("AU102SR44/au102_pmba44.gro")
#staples_Au102SR44 = classify_staples(xyz_Au102SR44, names_Au102SR44)
#write_pdb(xyz_Au102SR44, names_Au102SR44, staples_Au102SR44, "au102SR44_NM.pdb")
"""
#AU144SR60 (264 atoms)
xyz_Au144SR60, names_Au144SR60 = read_Au144SR60("AU144SR60/au144SR60.pdb")
staples_Au144SR60 = classify_staples(xyz_Au144SR60, names_Au144SR60)
write_pdb(xyz_Au144SR60, names_Au144SR60, staples_Au144SR60, "au144SR60_NM.pdb")
#print_xyz(xyz_Au144SR60, names_Au144SR60, "au144SR60_NM")


#AU314SR96 (506 atoms)
xyz_Au314SR96, names_Au314SR96 = read_Au314SH96("AU314SR96/au314SH96.xyz")
staples_Au314SR96 = classify_staples(xyz_Au314SR96, names_Au314SR96)
write_pdb(xyz_Au314SR96, names_Au314SR96, staples_Au314SR96, "au314SR96_NM.pdb")
"""
