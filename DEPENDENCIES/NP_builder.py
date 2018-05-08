import numpy as np
from  transformations import *
import random
from scipy.spatial import distance
from sklearn.decomposition import PCA

def rot_mat(p, u, t):
    ct = math.cos(t)
    st = math.sin(t)
    x = u[0]
    y = u[1]
    z = u[2]
    rot = np.array([[ct + x**2*(1-ct), x*y*(1-ct)-z*st, x*z*(1-ct)+y*st], \
    [x*y*(1-ct)+z*st, ct+y**2*(1-ct), y*z*(1-ct)-x*st],\
    [x*z*(1-ct)-y*st, y*z*(1-ct)+x*st, ct+z**2*(1-ct)]])
    return np.dot(rot, p)

def init_lig_mol2(fname, cap):
    #Imports ligand mol2 file. Returns xyz coordinates, names, and index corresponding to the anchor
    mol2=np.genfromtxt(fname, delimiter='\n', dtype='str')
    N_lig_file=len(mol2)
    found_ATOM=0
    names_lig_func=[]
    xyz_lig_func=[]
    resID_func = []
    res_lig_func=[]
    for i in range(N_lig_file):
        if found_ATOM:
            if "@<TRIPOS>" in mol2[i]:
                break
            at_file = mol2[i].split()
            names_lig_func.append(at_file[1])
            xyz_lig_func.append(at_file[2:5])
            resID_func.append(at_file[6])
            res_lig_func.append(at_file[7])
        elif "@<TRIPOS>ATOM" in mol2[i]:
            found_ATOM = True

    xyz_lig_func, names_lig_func, res_lig_func, resID_func = np.array(xyz_lig_func, dtype='float'), np.array(names_lig_func), np.array(res_lig_func), np.array(resID_func, dtype="int")

    if cap=="N":
        print("There are no capping atoms in the structure...")
    else:
        cap = np.array(cap.split(","), dtype="int")-1
        print(cap)
        print("Removing capping atoms...")
        xyz_lig_func, names_lig_func, res_lig_func, resID_func = np.delete(xyz_lig_func, cap, axis=0), np.delete(names_lig_func, cap), np.delete(res_lig_func, cap), np.delete(resID_func, cap)

    for i in range(N_lig_file):
        if "@<TRIPOS>RESIDUECONNECT" in mol2[i]:
            name_anchor_func = mol2[i+1].split()[1]
            anchor_ndx_func = np.where(np.logical_and(names_lig_func==name_anchor_func, resID_func==int(mol2[i+1].split()[0])))[0][0]
            res_anchor_func = res_lig_func[anchor_ndx_func]

    anchor_pos = np.copy(xyz_lig_func)[anchor_ndx_func,:]
    #Moves the ligand so that the anchor is in (0,0,0)
    for i in range(len(xyz_lig_func[:,0])):
        xyz_lig_func[i,:] = xyz_lig_func[i,:] - anchor_pos
    return xyz_lig_func, names_lig_func, anchor_ndx_func, name_anchor_func, res_anchor_func, res_lig_func

def init_core_pdb(fname):
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
    return xyz_core_func, names_core_func, res_core_func

def get_ligand_pill(xyz_lig_func, anchor_ndx_func):
    #Runs a PCA and takes the first eigenvector as the best fitting line.
    pca = PCA(n_components=3)
    pca.fit(xyz_lig_func)
    pca1 = pca.components_[0]
    var1 = pca.explained_variance_[0]/np.sum(pca.explained_variance_)*100

    print("PCA1 explains: {:.1f}% of the points' variance...".format(var1))
    print("Consider this is a measure on how linear the input ligand is. The lower this value, the more likely it will be to get clashes in the final structure.")

    #Randomly takes 2 other atoms in the ligand and project their positions in PCA1
    random.seed(666)
    rango = list(range(len(xyz_lig_func[:,0])))
    rango.remove(anchor_ndx_func)
    pillars_ndx = random.sample(rango, 2)
    pillars_func = np.array([0.0, 0.0, 0.0]) #This corresponds to the first stone (i.e. the anchor) which will always be in (0,0,0)
    for i in pillars_ndx:
        pillars_func = np.vstack((pillars_func, np.dot(xyz_lig_func[i], pca1) * pca1))
    return pillars_func

def assign_morph(xyz_core_func, names_core_func, frac_lig1_func, rseed_func, morph_func, stripes_func, first_func):
    #Distributes all the anchors in lig1 and lig2 dependending in the specified morphology
    xyz_anchors_func = xyz_core_func[names_core_func=='C',:]
    N_anchors = len(xyz_anchors_func)
    for_lig1 = round(N_anchors*frac_lig1_func)
    indexes = list(range(N_anchors))
    if morph_func == "random":
        print("Assigning a random distribution of the ligands...")
        random.seed(rseed_func)
        random.shuffle(indexes)
        lig1_ndx = indexes[:for_lig1]
        lig2_ndx = indexes[for_lig1:]
    elif morph_func == "janus":
        print("Assigning a janus distribution for the ligands...")
        bottom = xyz_anchors_func[np.argsort(xyz_anchors_func[:,2])[0]]
        D_bottom_anch = distance.cdist([bottom], xyz_anchors_func)
        lig1_ndx = D_bottom_anch[0].argsort()[:for_lig1]
        lig2_ndx = list(set(indexes) - set(lig1_ndx))
    elif morph_func == "stripe":
        print("Assigning a striped distribution for the ligands...")
        max_Z = np.max(xyz_anchors_func[:,2])+0.00001
        min_Z = np.min(xyz_anchors_func[:,2])
        dZ = (max_Z - min_Z)/stripes_func
        lig1_ndx = []
        lig2_ndx = []
        if first_func == 1:
            print("The stripes will be generated from the bottom up starting with ligand 1...")
            firstOne = True
        elif first_func == 2:
            print("The stripes will be generated from the bottom up starting with ligand 2...")
            firstOne = False
        for i in range(N_anchors):
            if firstOne:
                if (xyz_anchors_func[i,2]-min_Z)//dZ%2 == 0:
                    lig1_ndx.append(i)
                elif (xyz_anchors_func[i,2]-min_Z)//dZ%2 == 1:
                    lig2_ndx.append(i)
            else:
                if (xyz_anchors_func[i,2]-min_Z)//dZ%2 == 0:
                    lig2_ndx.append(i)
                elif (xyz_anchors_func[i,2]-min_Z)//dZ%2 == 1:
                    lig1_ndx.append(i)        

    xyz_anchors1_func=xyz_anchors_func[lig1_ndx]
    xyz_anchors2_func=xyz_anchors_func[lig2_ndx]
    return xyz_anchors1_func, xyz_anchors2_func

def get_stones(xyz_anchorsi_func, xyz_pillarsi_func):
    #Return a 3D array with the xyz coordinates for all the stones of all the anchors
    n_stones_lig = len(xyz_pillarsi_func)
    n_anchors = len(xyz_anchorsi_func)
    xyz_stones = np.zeros((n_anchors, n_stones_lig, 3))

    #Takes the COM-C vectors and scale them to match the distance between staples in the ligand's file
    for i in range(n_anchors):
        mag_C = np.linalg.norm(xyz_anchorsi_func[i,:])
        for j in range(n_stones_lig):
            scaling = (mag_C + np.linalg.norm(xyz_pillarsi_func[j,:]))/mag_C
            xyz_stones[i,j,:]=xyz_anchorsi_func[i,:]*scaling
    return xyz_stones

def solve_clashes(xyz_coated_tmp, trans_lig_tmp, xyz_stone_act, resnum):
    n_clash_iter = 100
    thresh = 1.0
    D_clash = distance.cdist(trans_lig_tmp, xyz_coated_tmp)
    clash_dis = np.min(D_clash)
    theta = 0
    trans_lig_best = trans_lig_tmp
    if clash_dis < thresh:
        print("Clashes were found while placing residue {}...".format(resnum))
        print("Trying to solve the clashes...")

    while theta < 6.28:
        theta += 2*math.pi/n_clash_iter
        trans_lig_try = trans_lig_tmp
        unit_u = xyz_stone_act/np.linalg.norm(xyz_stone_act)
        for k in range(len(trans_lig_tmp)):
            trans_lig_try[k,:] = rot_mat(trans_lig_tmp[k,:], unit_u, theta)
        D_clash = distance.cdist(trans_lig_try, xyz_coated_tmp)
        if np.min(D_clash) > clash_dis:
            clash_dis = np.min(D_clash)
            trans_lig_best = trans_lig_try
        if theta >= 6.28 and clash_dis < thresh:
            print("It was not possible to solve all the clashes. Residue {} has a close contact of {:.2f} nm...".format(resnum, clash_dis/10))
            print("Revise the final geometry...")
            break

    return trans_lig_best

def coat_NP(xyz_core_func, names_core_func, frac_lig1_func, xyz_lig1_func, names_lig1_func, xyz_pillars1_func, xyz_stones1_func, xyz_lig2_func, names_lig2_func, xyz_pillars2_func, xyz_stones2_func, res_lig1_func, res_lig2_func):
    #Merges xyz coordinates and names of the core and the ligands into one coated NP
    keep_rows=[]
    for i in range(len(names_core_func)):
        if names_core_func[i]!='C':
            keep_rows.append(i)

    xyz_coated_func=xyz_core_func[keep_rows,:]
    names_coated_func=names_core_func[keep_rows]
    res_coated_func=names_core_func[keep_rows]
    res_coated_func[np.where(names_coated_func=="AUL")[0]]="AU"
    res_coated_func[np.where(names_coated_func=="AUS")[0]]="AU"
    names_coated_func[np.where(names_coated_func=="AUL")[0]]="AU"
    names_coated_func[np.where(names_coated_func=="AUS")[0]]="AU"


    #Transforms and appends rototranslated ligand 1
    xyz_lig1_func_conv=np.insert(xyz_lig1_func, 3, 1, axis=1).T
    for i in range(len(xyz_stones1_func[:,0,0])):
        xyz_stones_now = xyz_stones1_func[i,:,:]
        trans_matrix=affine_matrix_from_points(xyz_pillars1_func.T, xyz_stones_now.T, shear=False, scale=False, usesvd=True)
        trans_lig=np.dot(trans_matrix, xyz_lig1_func_conv).T[:,:3]

        trans_lig = solve_clashes(xyz_coated_func, trans_lig, xyz_stones1_func[i,0,:], len(keep_rows)+i+1)

        xyz_coated_func=np.append(xyz_coated_func, trans_lig, axis=0)
        names_coated_func=np.append(names_coated_func, names_lig1_func, axis=0)
        res_coated_func=np.append(res_coated_func, res_lig1_func, axis=0)

    #Transforms and appends rototranslated ligand 2
    if frac_lig1_func < 1.0:
        xyz_lig2_func_conv=np.insert(xyz_lig2_func, 3, 1, axis=1).T
        for i in range(len(xyz_stones2_func[:,0,0])):
            xyz_stones_now = xyz_stones2_func[i,:,:]
            trans_matrix=affine_matrix_from_points(xyz_pillars2_func.T, xyz_stones_now.T, shear=False, scale=False, usesvd=True)
            trans_lig=np.dot(trans_matrix, xyz_lig2_func_conv).T[:,:3]

            trans_lig = solve_clashes(xyz_coated_func, trans_lig, xyz_stones2_func[i,0,:], len(keep_rows)+len(xyz_stones1_func[:,0,0])+i+1)

            xyz_coated_func=np.append(xyz_coated_func, trans_lig, axis=0)
            names_coated_func=np.append(names_coated_func, names_lig2_func, axis=0)
            res_coated_func=np.append(res_coated_func, res_lig2_func, axis=0)
    return xyz_coated_func, names_coated_func, res_coated_func

def print_NP_pdb(xyz_coated_func, names_coated_func, res_coated_func, xyz_anchors1_func, xyz_anchors2_func, xyz_lig1_func, xyz_lig2_func, frac_lig1_func, out_fname):
    N_at_lig1 = len(xyz_lig1_func[:,0])
    if frac_lig1_func < 1.0:
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
