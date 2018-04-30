import numpy as np
from scipy import optimize
from scipy.spatial import distance
import math

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

def fix_Au144SR60(fname):
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
    return xyz, names

def fix_Au314SH96(fname):
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


def print_xyz(coords, nombres, fname):
    fxyz = open(fname + ".xyz", "w")
    fxyz.write("{} \n\n".format(len(nombres)))
    for i in range(len(nombres)):
        fxyz.write("{}\t\t{:.3f}   {:.3f}   {:.3f}\n".format(nombres[i], coords[i,0], coords[i,1], coords[i,2]))
    fxyz.close()

xyz_Au144SR60, names_Au144SR60 = fix_Au144SR60("au144SR60.pdb")
print_xyz(xyz_Au144SR60, names_Au144SR60, "au144SR60_NM")

xyz_Au314SR96, names_Au314SR96 = fix_Au314SH96("au314SH96.xyz")
print_xyz(xyz_Au314SR96, names_Au314SR96, "au314SR96_NM")
