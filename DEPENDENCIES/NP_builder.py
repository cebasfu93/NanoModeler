import numpy as np
from  DEPENDENCIES.transformations import *
import random
from scipy.spatial import distance
from sklearn.decomposition import PCA
import math

def rot_mat(p, u, t):
    ct = math.cos(t)
    st = math.sin(t)
    x = u[0]
    y = u[1]
    z = u[2]
    rot = np.array([[ct + x**2*(1-ct), x*y*(1-ct)-z*st, x*z*(1-ct)+y*st], \
    [x*y*(1-ct)+z*st, ct+y**2*(1-ct), y*z*(1-ct)-x*st],\
    [x*z*(1-ct)-y*st, y*z*(1-ct)+x*st, ct+z**2*(1-ct)]])
    return np.dot(rot, p.T).T

def phi(xyz):
    return math.acos(xyz[2]/np.linalg.norm(xyz))

def init_lig_mol2(fname, lig_s, lig_c):
    #Imports ligand mol2 file
    mol2=np.genfromtxt(fname, delimiter='\n', dtype='str')
    N_lig_file=len(mol2)
    found_ATOM=0
    names_lig_func=[]
    xyz_lig_func=[]
    res_lig_func=[]
    for i in range(N_lig_file):
        if found_ATOM:
            if "@<TRIPOS>BOND" in mol2[i]:
                break
            at_file = mol2[i].split()
            names_lig_func.append(at_file[1])
            xyz_lig_func.append(at_file[2:5])
            res_lig_func.append(at_file[7])
        elif "@<TRIPOS>ATOM" in mol2[i]:
            found_ATOM = True

    xyz_lig_func, names_lig_func, res_lig_func = np.array(xyz_lig_func, dtype='float'), np.array(names_lig_func), np.array(res_lig_func)

    origin = xyz_lig_func[lig_c-1,:]

    #Moves the ligand so that the S is in (0,0,0)
    xyz_lig_func = xyz_lig_func - origin
    return xyz_lig_func/10., names_lig_func, res_lig_func

def init_core_pdb(fname, elong):
    #Imports core pdb file. Centers the core in (0,0,0) and returns xyz coordinates and names
    pdb=np.genfromtxt(fname, delimiter='\n', dtype=str)
    names_core_func = []
    res_core_func = []
    xyz_core_func = []
    for i in range(len(pdb)):
        at_act=pdb[i]
        names_core_func.append(at_act[12:16].strip())
        res_core_func.append(at_act[17:20].strip())
        xyz_core_func.append([at_act[30:38], at_act[38:46], at_act[46:54]])
    xyz_core_func = np.array(xyz_core_func, dtype="float")
    names_core_func = np.array(names_core_func)
    res_core_func = np.array(res_core_func)
    if elong:
        xyz_C_new = []
        xyz_all_S = xyz_core_func[names_core_func=="ST"]
        norm_all_S = np.linalg.norm(xyz_all_S, axis=1)
        for j in range(len(xyz_all_S)):
            xyz_C_new.append(xyz_all_S[j]*(1.8+norm_all_S[j])/norm_all_S[j])
        xyz_C_new = np.array(xyz_C_new)
        xyz_core_func[names_core_func=="C"] = xyz_C_new
    return xyz_core_func/10., names_core_func, res_core_func

def get_ligand_pill(xyz_lig_func, lig_c, lig_s, log):
    #Runs a PCA and takes the first eigenvector as the best fitting line.
    N_at = len(xyz_lig_func)
    pca = PCA(n_components=3)
    if lig_s == 0:
        pca.fit(xyz_lig_func[:-1])
    else:
        pca.fit(xyz_lig_func)
    pca1 = pca.components_[0]
    if np.sum(np.mean(xyz_lig_func, axis=0)>=0)<2:
        pca1=-pca1
    var1 = pca.explained_variance_[0]/np.sum(pca.explained_variance_)*100

    log += "PCA1 explains: {:.1f}% of the points' variance...\n".format(var1)
    log += "Consider this is a measure on how linear the input ligand is. The lower this value, the more likely it will be to get clashes in the final structure.\n"

    #Randomly takes 2 other atoms in the ligand and project their positions in PCA1
    random.seed(666)
    rango = list(range(N_at))
    rango.remove(N_at-1)
    rango.remove(lig_c-1)
    pillars_ndx = random.sample(rango, 2)
    pillars_func = [xyz_lig_func[lig_c-1]]
    for p in pillars_ndx:
        pillars_func = np.append(pillars_func, np.array([np.dot(xyz_lig_func[p], pca1) * pca1]), axis=0)
    return pillars_func, log

def assign_morph(xyz_core_func, names_core_func, frac_lig1_func, rseed_func, morph_func, stripes_func, log):
    #Distributes all the anchors in lig1 and lig2 dependending in the specified morphology
    xyz_anchors_func = xyz_core_func[names_core_func=='C',:]
    N_anchors = len(xyz_anchors_func)
    for_lig1 = round(N_anchors*frac_lig1_func)
    indexes = list(range(N_anchors))
    if morph_func == "random":
        log += "Assigning a random distribution of the ligands...\n"
        random.seed(rseed_func)
        random.shuffle(indexes)
        lig1_ndx = indexes[:for_lig1]
        lig2_ndx = indexes[for_lig1:]
    elif morph_func == "janus":
        log += "Assigning a janus distribution for the ligands...\n"
        bottom = xyz_anchors_func[np.argsort(xyz_anchors_func[:,2])[0]]
        D_bottom_anch = distance.cdist([bottom], xyz_anchors_func)
        lig1_ndx = D_bottom_anch[0].argsort()[:for_lig1]
        lig2_ndx = list(set(indexes) - set(lig1_ndx))
    elif morph_func == "stripe":
        log += "Assigning a striped distribution for the ligands...\n"
        phis = np.arccos(np.divide(xyz_anchors_func[:,2], np.linalg.norm(xyz_anchors_func, axis=1)))
        dphi = (math.pi+0.00001)/stripes_func
        lig1_ndx = []
        lig2_ndx = []
        for i in range(N_anchors):
            if phi(xyz_anchors_func[i])//dphi%2 == 0:
                lig1_ndx.append(i)
            elif phi(xyz_anchors_func[i])//dphi%2 == 1:
                lig2_ndx.append(i)

    log += "The nanoparticle will have {} of ligand 1...\n".format(len(lig1_ndx))
    log += "The nanoparticle will have {} of ligand 2...\n".format(len(lig2_ndx))
    xyz_anchors1_func=xyz_anchors_func[lig1_ndx]
    xyz_anchors2_func=xyz_anchors_func[lig2_ndx]
    return xyz_anchors1_func, xyz_anchors2_func, log

def get_stones(xyz_core_func, names_core_func, xyz_anchorsi_func, xyz_pillarsi_func, lig_s):
    #Return a 3D array with the xyz coordinates for all the stones of all the anchors
    n_stones_lig = len(xyz_pillarsi_func)
    n_anchors = len(xyz_anchorsi_func)
    xyz_stones = np.zeros((n_anchors, n_stones_lig+1, 3))
    ndx_core_ST = np.where(names_core_func=="ST")[0]
    xyz_ST = xyz_core_func[ndx_core_ST]
    D_anch_ST = distance.cdist(xyz_anchorsi_func, xyz_ST)
    sort_D_anch_ST = np.argsort(D_anch_ST, axis=1)
    #Takes the COM-C vectors and scale them to match the distance between staples in the ligand's file
    for i in range(n_anchors):
        xyz_stones[i,0,:] = xyz_anchorsi_func[i,:]

        mag_C = np.linalg.norm(xyz_anchorsi_func[i,:])

        f1 = (mag_C + np.linalg.norm(xyz_pillarsi_func[1,:]))/mag_C
        xyz_stones[i,1,:] = xyz_anchorsi_func[i,:]*f1
        f2 = (mag_C + np.linalg.norm(xyz_pillarsi_func[2,:]))/mag_C
        xyz_stones[i,2,:] = xyz_anchorsi_func[i,:]*f2

        xyz_stones[i,3,:] = xyz_core_func[ndx_core_ST[sort_D_anch_ST[i,0]]]

    return xyz_stones

def solve_clashes(xyz_coated_tmp, trans_lig_tmp, xyz_stone_act, resnum, log):
    n_clash_iter = 100
    thresh = 0.1
    D_clash = distance.cdist(trans_lig_tmp, xyz_coated_tmp)
    clash_dis = np.min(D_clash)
    theta = 0
    trans_lig_best = trans_lig_tmp
    CLASH = False
    if clash_dis < thresh:
        CLASH = True
        log += "Clashes were found while placing residue {}...\n".format(resnum)
        log += "Trying to solve the clashes...\n"

    while theta < 6.28:
        theta += 2*math.pi/n_clash_iter
        trans_lig_try = trans_lig_tmp
        unit_u = xyz_stone_act/np.linalg.norm(xyz_stone_act)
        trans_lig_try = rot_mat(trans_lig_tmp, unit_u, theta)
        D_clash = distance.cdist(trans_lig_try, xyz_coated_tmp)
        if np.min(D_clash) > clash_dis:
            clash_dis = np.min(D_clash)
            trans_lig_best = trans_lig_try
        if theta >= 6.28 and clash_dis < thresh:
            log += "It was not possible to solve all the clashes. Residue {} has a close contact of {:.2f} nm...\n".format(resnum, clash_dis)
            log += "Revise the final geometry...\n"
        if theta >= 6.28 and clash_dis > thresh and CLASH:
            log += "The clash was solved, i.e., the minimum distance between atoms is at least {} nm...\n".format(thresh)

    return trans_lig_best, log

def coat_NP(xyz_core_func, names_core_func, xyz_lig1_func, names_lig1_func, xyz_pillars1_func, xyz_stones1_func, xyz_lig2_func, names_lig2_func, xyz_pillars2_func, xyz_stones2_func, res_lig1_func, res_lig2_func, lig1_s, lig2_s, elong, log):
    #Merges xyz coordinates and names of the core and the ligands into one coated NP
    keep_rows=[]
    for i in range(len(names_core_func)):
        if names_core_func[i]!='ST' and names_core_func[i]!='C':
            keep_rows.append(i)

    xyz_coated_func=xyz_core_func[keep_rows,:]
    names_coated_func=names_core_func[keep_rows]
    res_coated_func=names_core_func[keep_rows]

    #Transforms and appends rototranslated ligand 1
    xyz_lig1_func_conv=np.insert(xyz_lig1_func, 3, 1, axis=1).T
    for i in range(len(xyz_stones1_func[:,0,0])):
        xyz_stones_now = xyz_stones1_func[i,:-1,:]
        trans_matrix=affine_matrix_from_points(xyz_pillars1_func.T, xyz_stones_now.T, shear=False, scale=False, usesvd=True)
        trans_lig=np.dot(trans_matrix, xyz_lig1_func_conv).T[:,:3]
        trans_lig = trans_lig[:-1]

        if lig1_s == 0 or elong:
            trans_lig, log = solve_clashes(xyz_coated_func, trans_lig, xyz_stones_now[0,:], len(keep_rows)+i+1, log)

        else:
            D_clash = distance.cdist(trans_lig, xyz_coated_func)
            clash_dis = np.min(D_clash)
            if clash_dis < 0.1:
                log += "The sulphur atom was given in the mol2 file of ligand 1...\n"
                log += "There are no degrees of freedom available to prevent clashes...\n"
                log += "Clashes were found while placing residue {}...\n".format(len(keep_rows)+i+1)
                log += "Consider parametrizing ligand 1 without the thiol sulphur atom, then NanoModeler will try to find a conformation without clashes...\n"

        trans_lig = np.append(trans_lig, [xyz_stones1_func[i,-1,:]], axis=0)
        xyz_coated_func=np.append(xyz_coated_func, trans_lig, axis=0)
        names_coated_func=np.append(names_coated_func, names_lig1_func, axis=0)
        res_coated_func=np.append(res_coated_func, res_lig1_func, axis=0)

    #Transforms and appends rototranslated ligand 2
    if len(xyz_lig2_func)!=0:
        xyz_lig2_func_conv=np.insert(xyz_lig2_func, 3, 1, axis=1).T
        for i in range(len(xyz_stones2_func[:,0,0])):
            xyz_stones_now = xyz_stones2_func[i,:-1,:]
            trans_matrix=affine_matrix_from_points(xyz_pillars2_func.T, xyz_stones_now.T, shear=False, scale=False, usesvd=True)
            trans_lig=np.dot(trans_matrix, xyz_lig2_func_conv).T[:,:3]
            trans_lig = trans_lig[:-1]

            if lig2_s == 0 or elong:
                trans_lig, log = solve_clashes(xyz_coated_func, trans_lig, xyz_stones_now[0,:], len(keep_rows)+len(xyz_stones1_func[:,0,0])+i+1, log)

            else:
                D_clash = distance.cdist(trans_lig, xyz_coated_func)
                clash_dis = np.min(D_clash)
                if clash_dis < 0.1:
                    log += "The sulphur atom was given in the mol2 file of ligand 2...\n"
                    log += "There are no degrees of freedom available to prevent clashes...\n"
                    log += "Clashes were found while placing residue {}...\n".format(len(keep_rows)+len(xyz_stones1_func[:,0,0])+i+1)
                    log += "Consider parametrizing ligand 2 without the thiol sulphur atom, then NanoModeler will try to find a conformation without clashes...\n"

            trans_lig = np.append(trans_lig, [xyz_stones2_func[i,-1,:]], axis=0)
            xyz_coated_func=np.append(xyz_coated_func, trans_lig, axis=0)
            names_coated_func=np.append(names_coated_func, names_lig2_func, axis=0)
            res_coated_func=np.append(res_coated_func, res_lig2_func, axis=0)
    return xyz_coated_func, names_coated_func, res_coated_func, log

def print_NP_pdb(xyz_coated_func, names_coated_func, res_coated_func, xyz_anchors1_func, xyz_anchors2_func, xyz_lig1_func, xyz_lig2_func, out_fname):
    xyz_coated_func = xyz_coated_func*10
    xyz_anchors1_func = xyz_anchors1_func*10
    xyz_anchors2_func = xyz_anchors2_func*10
    xyz_lig1_func = xyz_lig1_func*10
    xyz_lig2_func = xyz_lig2_func*10
    N_at_lig1 = len(xyz_lig1_func[:,0])
    if len(xyz_lig2_func)!=0:
        N_at_lig2 = len(xyz_lig2_func[:,0])
    else:
        N_at_lig2 = 0

    #Writes the pdb of the core and the placed stones
    N_lig1 = len(xyz_anchors1_func)
    N_tot_lig1 = N_lig1 * N_at_lig1
    N_lig2 = len(xyz_anchors2_func)
    N_tot_lig2 = N_lig2 * N_at_lig2
    N_core = len(xyz_coated_func) - N_lig1*N_at_lig1 - N_lig2*N_at_lig2

    at=0
    res=0
    output=open(out_fname, "w")
    #Writes the core
    for i in range(N_core):
        at+=1
        res+=1
        at_name_act=names_coated_func[i]
        r_name_act=res_coated_func[i]
        write_pdb_block(at_name_act, r_name_act, xyz_coated_func[i,:], res, at, out_fname)

    #Writes ligand 1
    lig_atoms=0
    for i in range(N_tot_lig1):
        at+=1
        at_name_act=names_coated_func[i+N_core]
        r_name_act=res_coated_func[i+N_core]
        if(lig_atoms%N_at_lig1==0):
            res+=1
        lig_atoms+=1
        write_pdb_block(at_name_act, r_name_act, xyz_coated_func[i+N_core,:], res, at, out_fname)

    #Writes ligand 2
    lig_atoms=0
    for i in range(N_tot_lig2):
        at+=1
        at_name_act=names_coated_func[i+N_core+N_tot_lig1]
        r_name_act=res_coated_func[i+N_core+N_tot_lig1]
        if(lig_atoms%N_at_lig2==0):
            res+=1
        lig_atoms+=1
        write_pdb_block(at_name_act, r_name_act, xyz_coated_func[i+N_core+N_tot_lig1,:], res, at, out_fname)

    output.close()

def write_pdb_block(atname_func, res_name_func, xyz_func, resnum, atnum, out_filename):
    #Writes the information of one atom in a pdb file
    xyz_func=np.round(xyz_func, decimals=4)
    coords=open(out_filename, 'a')
    coords.write('ATOM'.ljust(6))
    coords.write(str(atnum).rjust(5))
    coords.write(' ' + str(atname_func).ljust(4))
    coords.write(' '+str(res_name_func).ljust(3))
    coords.write('  '+str(resnum).rjust(4))
    coords.write('    ' + str(round(xyz_func[0],3)).rjust(8))
    coords.write(str(round(xyz_func[1],3)).rjust(8))
    coords.write(str(round(xyz_func[2],3)).rjust(8)+"\n")
    coords.close()
