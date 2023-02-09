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
    if np.abs(cosine_angle-1.0) < 0.001:
        angle = 0
    else:
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
    angles = np.abs((angles-109.2)/109.2)
    D_S = np.abs((D_S-1.8)/1.8)
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

def read_Au36CP24(fname):
    fxyz = np.genfromtxt(fname, skip_header=1, dtype=str)
    xyz = []
    xyz_C = []
    names = []
    N_xyz = len(fxyz)
    for i in range(N_xyz):
        if fxyz[i,0] == "Au":
            xyz.append(fxyz[i,1:])
            names.append("AU")
        elif fxyz[i,0] == "S":
            xyz.append(fxyz[i,1:])
            names.append("ST")
        elif fxyz[i,0] == "C":
            xyz_C.append(fxyz[i,1:])
    xyz, xyz_C = np.array(xyz, dtype='float'), np.array(xyz_C, dtype='float')
    names = np.array(names)
    D_ST_C = distance.cdist(xyz[np.where(names=="ST")[0]], xyz_C)
    for i in range(len(D_ST_C)):
        xyz = np.vstack((xyz, xyz_C[np.argsort(D_ST_C[i,:])[0]]))
        names = np.append(names, "C")
    xyz = center(xyz)
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

def read_A133L52(fname):
    txt = np.genfromtxt(fname, dtype='str', delimiter="\n")
    letters = []
    xyz = []
    names = []
    xs, ys, zs = [], [], []
    ini = 0

    for i in range(len(txt)):
        if "A17" in txt[i]:
            ini = i
            break
        if "A"in txt[i]:
            xyz.append([float(txt[i+1]), float(txt[i+2]), float(txt[i+3])])
            names.append('AU')

    X, Y, Z = False, False, False
    for i in range(ini, len(txt)):
        if "A" in txt[i]:
            names.append("AU")
        elif "S" in txt[i]:
            names.append("ST")
        elif "C" in txt[i]:
            names.append("C")
        elif "x" not in txt[i] and "y" not in txt[i] and "z" not in txt[i]:
            if X:
                xs.append(txt[i])
            elif Y:
                ys.append(txt[i])
            elif Z:
                zs.append(txt[i])

        elif "x" in txt[i]:
            X, Y, Z = True, False, False
        elif "y" in txt[i]:
            X, Y, Z = False, True, False
        elif "z" in txt[i]:
            X, Y, Z = False, False, True

    for i in range(len(xs)):
        new_coord = [xs[i], ys[i], zs[i]]
        xyz = np.vstack((xyz, new_coord))
    xyz = np.array(xyz, dtype='float')
    xyz[:,0]=xyz[:,0]*30.14
    xyz[:,1]=xyz[:,1]*30.44
    xyz[:,2]=xyz[:,2]*43.69
    names = np.array(names)

    for i in range(len(xyz)):
        if xyz[i,1]<4:
            xyz[i,1] = xyz[i,1] + 30.44

    where_S = np.where(names=="ST")[0]
    where_C = np.where(names=="C")[0]
    new_xyz = xyz[np.where(names=="AU")[0]]
    new_names = []
    N_AU = len(new_xyz)
    for i in range(N_AU):
        new_names.append("AU")
    new_xyz = np.vstack((new_xyz, xyz[where_S]))
    for i in range(len(where_S)):
        new_names.append("ST")
    D_S_C = distance.cdist(xyz[where_S], xyz[where_C])
    for i in range(len(D_S_C)):
        new_xyz = np.vstack((new_xyz, xyz[where_C[np.argsort(D_S_C[i])[0]]]))
        new_names.append("C")


    new_xyz = center(new_xyz)
    new_names = np.array(new_names)
    return new_xyz, new_names

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

def read_Au68SR34(fname):
    fxyz = np.genfromtxt(fname, skip_header=2, dtype=str)
    xyz = []
    names = []
    N_xyz = len(fxyz)
    for i in range(N_xyz):
        if fxyz[i,0] == "Xx":
            xyz.append(fxyz[i,1:])
            names.append("AU")
        elif fxyz[i,0] == "S":
            xyz.append(fxyz[i,1:])
            names.append("ST")
        elif fxyz[i,0] == "C":
            xyz.append(fxyz[i,1:])
            names.append("C")

    xyz = center(xyz)
    names = np.array(names)
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

def read_Au44SR28(fname):
    fxyz = np.genfromtxt(fname, skip_header=1, dtype='str')
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
        elif fxyz[i,0] == "C":
            xyz.append(fxyz[i,1:])
            names.append("C")

    xyz = center(xyz)
    names = np.array(names)
    return xyz, names

def classify_staples(xyz_sys, names_sys):

    ndx_AU = np.where(names_sys=='AU')[0]
    ndx_ST = np.where(names_sys=='ST')[0]
    ndx_C = np.where(names_sys=='C')[0]

    N_AU = len(ndx_AU)
    N_ST = len(ndx_ST)

    blocks = []
    D_ST_AU = distance.cdist(xyz_sys[ndx_ST],xyz_sys[ndx_AU])
    D_ST_C = distance.cdist(xyz_sys[ndx_ST],xyz_sys[ndx_C])
    for i in range(N_ST):
        ndx_c = np.argsort(D_ST_C[i])[0]
        ndx_au = np.argsort(D_ST_AU[i])[0:2]
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
        if not taken:
            new_ganchos.append(ganchos[i])
        elif taken:
            new_ganchos[old].add(ganchos[i].S, ganchos[i].Au, ganchos[i].C)
    N_ganchos = len(new_ganchos)

    staples = []
    for i in range(N_ganchos):
        staples.append(subunits.Staple(ndx_S=new_ganchos[i].S, ndx_Au=new_ganchos[i].Au, ndx_C=new_ganchos[i].C))

    #Depending in the number of sulphur and gold atoms in each staple, it clssifies it. For STC and STV it calculates the angle Aul-S-Aul and with an tolerance of +/-9 degrees, the staple is classified
    print("Au{} SR{}".format(len(ndx_AU), len(ndx_ST)))

    new_staples = []
    for i in range(len(staples)):
        staple_act = staples[i]
        N_S = len(staple_act.S)
        N_Au = len(staple_act.Au)
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
                    if angle <= 109.0 and angle >= 83.0:
                        staple_act.change_tipo('STC')
                    elif angle <= 128.2 and angle >= 110.2:
                        staple_act.change_tipo('STV')
                    else:
                        print("One of the staples should be STC or STV but the AuL-S-AuL angle displays an odd value. Unrecognized staple")
        if not (N_S == 4 and N_Au == 5 and N_AU == 102):
            new_staples.append(staple_act)

        else:
            print("Unrecognized staple")

    if N_AU == 102:
        new_staples.append(subunits.Staple(ndx_S=[116, 117], ndx_Au=[29, 88, 88, 39], ndx_C=[160, 161], tipo="STR"))
        new_staples.append(subunits.Staple(ndx_S=[122, 123], ndx_Au=[38, 85, 85], ndx_C = [166, 167], tipo="STR"))
        new_staples.append(subunits.Staple(ndx_S=[134, 135], ndx_Au=[68, 97, 97, 78], ndx_C=[178, 179], tipo="STR"))
        new_staples.append(subunits.Staple(ndx_S=[138, 139], ndx_Au=[77, 94, 94], ndx_C=[182, 183], tipo="STR"))

    return new_staples

def write_pdb(xyz_sys, names_sys, staples, fname):
    if len(xyz_sys) == 192:
        xyz_sys = np.delete(xyz_sys, 79, axis = 0)
        xyz_sys = np.delete(xyz_sys, 39, axis = 0)
        names_sys = np.delete(names_sys, 79, axis = 0)
        names_sys = np.delete(names_sys, 39, axis = 0)

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
        #print(vars(staples[i]))
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
#AU25SR18 (61 atoms, Pohjolainen/Hakkinen, JCTC, 2016)
xyz_Au25SR18, names_Au25SR18 = read_Au25SR18("AU25SR18/au25_pet18.gro")
staples_Au25SR18 = classify_staples(xyz_Au25SR18, names_Au25SR18)
write_pdb(xyz_Au25SR18, names_Au25SR18, staples_Au25SR18, "au25SR18_NM.pdb")

#AU36SR24 (84 atoms, Das/Jin, J. Phys. Chem. A, 2015)
xyz_Au36SR24, names_Au36SR24 = read_Au36CP24("AU36SR24/Au36CP24.xyz")
staples_Au36SR24 = classify_staples(xyz_Au36SR24, names_Au36SR24)
write_pdb(xyz_Au36SR24, names_Au36SR24, staples_Au36SR24, "au36SR24_NM.pdb")

#AU38SR24 (86 atoms, Pohjolainen/Hakkinen, JCTC, 2016... Qian/Jin, JACS, 2010)
xyz_Au38SR24, names_Au38SR24 = read_Au25SR18("AU38SR24/au38_pet24.gro")
staples_Au38SR24 = classify_staples(xyz_Au38SR24, names_Au38SR24)
write_pdb(xyz_Au38SR24, names_Au38SR24, staples_Au38SR24, "au38SR24_NM.pdb")

#AU44SR28 (100 atoms, Pei/Liu, JACS, 2013)
xyz_Au44SR28, names_Au44SR28 = read_Au44SR28("AU44SR28/Au44SR28.xyz")
staples_Au44SR28 = classify_staples(xyz_Au44SR28, names_Au44SR28)
write_pdb(xyz_Au44SR28, names_Au44SR28, staples_Au44SR28, "au44SR28_NM.pdb")

#AU68SR32 (132 atoms, Wu/Zeng, Sci. Adv, 2015)
xyz_Au68SR32, names_Au68SR32 = read_Au314SH96("AU68SR32/Au68SH32-I1.xyz")
staples_Au68SR32 = classify_staples(xyz_Au68SR32, names_Au68SR32)
write_pdb(xyz_Au68SR32, names_Au68SR32, staples_Au68SR32, "au68SR32-I1_NM.pdb")

xyz_Au68SR32, names_Au68SR32 = read_Au314SH96("AU68SR32/Au68SH32-I2.xyz")
staples_Au68SR32 = classify_staples(xyz_Au68SR32, names_Au68SR32)
write_pdb(xyz_Au68SR32, names_Au68SR32, staples_Au68SR32, "au68SR32-I2_NM.pdb")

xyz_Au68SR32, names_Au68SR32 = read_Au314SH96("AU68SR32/Au68SH32-I3.xyz")
staples_Au68SR32 = classify_staples(xyz_Au68SR32, names_Au68SR32)
write_pdb(xyz_Au68SR32, names_Au68SR32, staples_Au68SR32, "au68SR32-I3_NM.pdb")

xyz_Au68SR32, names_Au68SR32 = read_Au314SH96("AU68SR32/Au68SH32-I4.xyz")
staples_Au68SR32 = classify_staples(xyz_Au68SR32, names_Au68SR32)
write_pdb(xyz_Au68SR32, names_Au68SR32, staples_Au68SR32, "au68SR32-I4_NM.pdb")

#AU68SR34 (136 atoms, Xu/Gao, J. Phys. Chem C, 2015)
xyz_Au68SR34, names_Au68SR34 = read_Au314SH96("AU68SR34/Au68SH34-I1.xyz")
staples_Au68SR34 = classify_staples(xyz_Au68SR34, names_Au68SR34)
write_pdb(xyz_Au68SR34, names_Au68SR34, staples_Au68SR34, "au68SR34-I1_NM.pdb")

xyz_Au68SR34, names_Au68SR34 = read_Au314SH96("AU68SR34/Au68SH34-I2.xyz")
staples_Au68SR34 = classify_staples(xyz_Au68SR34, names_Au68SR34)
write_pdb(xyz_Au68SR34, names_Au68SR34, staples_Au68SR34, "au68SR34-I2_NM.pdb")

xyz_Au68SR34, names_Au68SR34 = read_Au314SH96("AU68SR34/Au68SH34-I3.xyz")
staples_Au68SR34 = classify_staples(xyz_Au68SR34, names_Au68SR34)
write_pdb(xyz_Au68SR34, names_Au68SR34, staples_Au68SR34, "au68SR34-I3_NM.pdb")

xyz_Au68SR34, names_Au68SR34 = read_Au314SH96("AU68SR34/Au68SH34-I4.xyz")
staples_Au68SR34 = classify_staples(xyz_Au68SR34, names_Au68SR34)
write_pdb(xyz_Au68SR34, names_Au68SR34, staples_Au68SR34, "au68SR34-I4_NM.pdb")

#AU102SR44 (190 atoms, Pohjolainen/Hakkinen, JCTC, 2016)
#xyz_Au102SR44, names_Au102SR44 = read_Au102SR44("AU102SR44/au102_pmba44.gro")
#staples_Au102SR44 = classify_staples(xyz_Au102SR44, names_Au102SR44)
#write_pdb(xyz_Au102SR44, names_Au102SR44, staples_Au102SR44, "au102SR44_NM.pdb")

#AU133SR52 (237 atoms, Zeng/Jin, Sci. Adv, 2015)
xyz_Au133SR52, names_Au133SR52 = read_A133L52("AU133SR52/Au133SR52.txt")
staples_Au133SR52 = classify_staples(xyz_Au133SR52, names_Au133SR52)
write_pdb(xyz_Au133SR52, names_Au133SR52, staples_Au133SR52,"au133SR52_NM.pdb")

#AU144SR60 (264 atoms, Lopez-Acevedo/Hakkinen, J. Phys. Chem. C. Lett, 2009)
xyz_Au144SR60, names_Au144SR60 = read_Au144SR60("AU144SR60/au144SR60.pdb")
staples_Au144SR60 = classify_staples(xyz_Au144SR60, names_Au144SR60)
write_pdb(xyz_Au144SR60, names_Au144SR60, staples_Au144SR60, "au144SR60_NM.pdb")
#print_xyz(xyz_Au144SR60, names_Au144SR60, "au144SR60_NM")

#AU314SR96 (506 atoms, Malola/Hakkinen, ACS Nano, 2013)
xyz_Au314SR96, names_Au314SR96 = read_Au314SH96("AU314SR96/au314SH96.xyz")
staples_Au314SR96 = classify_staples(xyz_Au314SR96, names_Au314SR96)
write_pdb(xyz_Au314SR96, names_Au314SR96, staples_Au314SR96, "au314SR96_NM.pdb")
"""

####Prepared with minimized staples
core_dic = {"AU25SR18/test-25/test-25_CORE.xyz": "au25SR18_NM.pdb", \
            "AU36SR24/test-36/test-36_CORE.xyz": "au36SR24_NM.pdb", \
            "AU38SR24/test-38/test-38_CORE.xyz": "au38SR24_NM.pdb", \
            "AU44SR28/test-44/test-44_CORE.xyz": "au44SR28_NM.pdb", \
            "AU68SR32/test-I1/test-I1_CORE.xyz": "au68SR32-I1_NM.pdb", \
            "AU68SR32/test-I2/test-I2_CORE.xyz": "au68SR32-I2_NM.pdb", \
            "AU68SR32/test-I3/test-I3_CORE.xyz": "au68SR32-I3_NM.pdb", \
            "AU68SR32/test-I4/test-I4_CORE.xyz": "au68SR32-I4_NM.pdb", \
            "AU68SR34/test-I1/test-I1_CORE.xyz": "au68SR34-I1_NM.pdb", \
            "AU68SR34/test-I2/test-I2_CORE.xyz": "au68SR34-I2_NM.pdb", \
            "AU68SR34/test-I3/test-I3_CORE.xyz": "au68SR34-I3_NM.pdb", \
            "AU68SR34/test-I4/test-I4_CORE.xyz": "au68SR34-I4_NM.pdb", \
            "AU133SR52/test-133/test-133_CORE.xyz": "au133SR52_NM.pdb", \
            "AU144SR60/test-144/test-144_CORE.xyz": "au144SR60_NM.pdb", \
            "AU314SR96/test-314/test-314_CORE.xyz": "au314SR96_NM.pdb", \
}

for i in core_dic:
    xyz_core, names_core = read_Au68SR34(i)
    staples_core = classify_staples(xyz_core, names_core)
    write_pdb(xyz_core, names_core, staples_core, core_dic[i])

"""
#AU102SR44 (190 atoms, Pohjolainen/Hakkinen, JCTC, 2016)
xyz_Au102SR44, names_Au102SR44 = read_Au102SR44("AU102SR44/au102_pmba44.gro")
staples_Au102SR44 = classify_staples(xyz_Au102SR44, names_Au102SR44)
write_pdb(xyz_Au102SR44, names_Au102SR44, staples_Au102SR44, "au102SR44_NM.pdb")
"""
