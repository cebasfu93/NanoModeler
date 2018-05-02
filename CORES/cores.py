import numpy as np
from scipy import optimize
from scipy.spatial import distance
import math
import subunits

def angle(a, b, c):
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
    fxyz = np.genfromtxt(fname, delimiter='\n', dtype=str, skip_header=2)
    N_fxyz = len(fxyz)
    xyz = []
    names = []
    N_S = 0
    for i in range(N_fxyz):
        at_act = fxyz[i].split()
        if "Au" in at_act:
            xyz.append(at_act[1:4])
            names.append('AU')
        elif "S" in at_act:
            N_S+=1
            xyz.append(at_act[1:4])
            names.append("ST")

    xyz = center(xyz)
    xyz = xyz[np.argsort(names)]
    names = np.sort(names)
    S_xyz = xyz[np.where(names=="ST")]
    Au_xyz = xyz[np.where(names=="AU")]
    D_S_Au = distance.cdist(S_xyz, Au_xyz)
    ndx = np.argsort(D_S_Au, axis=1)
    for i in range(N_S):
        s_xyz = S_xyz[i]
        au_xyz = Au_xyz[ndx[i,:2]]
        s_norm = np.linalg.norm(s_xyz)
        x0 = (s_norm + 1.8)/s_norm * s_xyz
        #sol = optimize.minimize(penalty, x0, (s_xyz, au_xyz)).x
        dic={"args":(s_xyz, au_xyz)}
        sol = optimize.basinhopping(penalty, x0, minimizer_kwargs=dic, stepsize=0.2, niter=10).x
        xyz = np.vstack((xyz, sol))
        names = np.append(names, "C")
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
            if any(x in blocks[i].Au for x in ganchos[j].Au):
                ganchos[j].add(blocks[i].S, blocks[i].Au, blocks[i].C)
                taken = True
        if not taken:
            ganchos.append(subunits.Gancho(ndx_S=blocks[i].S, ndx_Au=blocks[i].Au, ndx_C=blocks[i].C))

    N_ganchos = len(ganchos)
    staples = []
    for i in range(N_ganchos):
        staples.append(subunits.Staple(ndx_S=ganchos[i].S, ndx_Au=ganchos[i].Au, ndx_C=ganchos[i].C))

    #Depending in the number of sulphur and gold atoms in each staple, it clssifies it. For STC and STV it calculates the angle Aul-S-Aul and with an tolerance of +/-9 degrees, the staple is classified
    for i in range(len(staples)):
        staple_act = staples[i]
        N_S = len(staple_act.S)
        N_Au = len(staple_act.Au)
        if N_S == 1 and N_Au == 2:
            staple_act.change_tipo('STP')
        elif N_S == 2 and N_Au == 3:
            staple_act.change_tipo('STR')
        elif N_S == 3 and N_Au == 4:
            D_S_Aul = distance.cdist(xyz_sys[staple_act.S], xyz_sys[staple_act.Au_l])
            for j in range(len(staple_act.S)):
                near_Au = staple_act.S[D_S_Aul[j].argsort()[0:2]]
                if np.all(np.in1d(near_Au, staple_act.Au_l)):
                    angle = angle(xyz_sys[near_Au[0]], xyz_sys[staple_act.S[j]], xyz_sys[near_Au[1]])
                    if angle <= 109.0 and angle >= 91.0:
                        staple_act.change_tipo('STC')
                    elif angle <= 128.2 and angle >= 110.2:
                        staple_act.change_tipo('STV')
                    else:
                        print("One of the staples should be STC or STV but the AuL-S-AuL angle displays an odd value. Unrecognized staple")
        else:
            print("Unrecognized staple")
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
    for i in range(N_AU):
        if i not in all_busy_Au:
            write_pdb_block("AU", "AU", xyz_sys[i], res, at, fname)
            at+=1
            res+=1

    for i in range(N_staples):
        st_act = staples[i]
        for j in range(len(st_act.Au_l)):
            write_pdb_block("AUL", st_act.tipo, xyz_sys[st_act.Au_l[j]], res, at, fname)
            at+=1
        for j in range(len(st_act.Au_s)):
            write_pdb_block("AUS", st_act.tipo, xyz_sys[st_act.Au_s[j]], res, at, fname)
            at+=1
        for j in range(len(st_act.S)):
            write_pdb_block("ST", st_act.tipo, xyz_sys[st_act.S[j]], res, at, fname)
            at+=1
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

xyz_Au144SR60, names_Au144SR60 = read_Au144SR60("au144SR60.pdb")
staples_Au144SR60 = classify_staples(xyz_Au144SR60, names_Au144SR60)
write_pdb(xyz_Au144SR60, names_Au144SR60, staples_Au144SR60, "au144SR60_NM.pdb")
#print_xyz(xyz_Au144SR60, names_Au144SR60, "au144SR60_NM")

#xyz_Au314SR96, names_Au314SR96 = fix_Au314SH96("au314SH96.xyz")
#print_xyz(xyz_Au314SR96, names_Au314SR96, "au314SR96_NM")
