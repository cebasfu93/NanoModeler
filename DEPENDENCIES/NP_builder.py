import numpy as np
from  transformations import *
import random
from scipy.spatial import distance
from sklearn.decomposition import PCA

def init_lig_mol2(fname):
    #Imports ligand mol2 file. Returns xyz coordinates, names, and index corresponding to the anchor
    mol2=np.genfromtxt(fname, delimiter='\n', dtype='str')
    N_lig_file=len(mol2)
    found_ATOM=0
    names_lig_func=[]
    xyz_lig_func=[]
    res_lig_func=[]
    for i in range(N_lig_file):
        if found_ATOM:
            if "@<TRIPOS>" in mol2[i]:
                break
            at_file = mol2[i].split()
            names_lig_func.append(at_file[1])
            xyz_lig_func.append(at_file[2:5])
            res_lig_func.append(at_file[7])
        elif "@<TRIPOS>ATOM" in mol2[i]:
            found_ATOM = True

    xyz_lig_func, names_lig_func, res_lig_func = np.array(xyz_lig_func, dtype='float'), np.array(names_lig_func), np.array(res_lig_func)
    for i in range(N_lig_file):
        if "@<TRIPOS>RESIDUECONNECT" in mol2[i]:
            anchor_ndx_func = np.where(names_lig_func==mol2[i+1].split()[1])[0][0]

    anchor_pos = np.copy(xyz_lig_func)[anchor_ndx_func,:]

    #Moves the ligand so that the anchor is in (0,0,0)
    for i in range(len(xyz_lig_func[:,0])):
        xyz_lig_func[i,:] = xyz_lig_func[i,:] - anchor_pos
    return xyz_lig_func, names_lig_func, anchor_ndx_func, res_lig_func

def init_core_xyz(fname):
    #Imports core pdb file. Centers the core in (0,0,0) and returns xyz coordinates and names
    fxyz=np.genfromtxt(fname, delimiter='\n', dtype=str, skip_header=1)
    names_core_func = []
    xyz_core_func = []
    for i in range(len(fxyz)):
        at_act=fxyz[i].split()
        names_core_func.append(at_act[0])
        xyz_core_func.append(at_act[1:4])
    xyz_core_func=np.array(xyz_core_func, dtype='float')
    names_core_func=np.array(names_core_func)
    return xyz_core_func, names_core_func

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

def assign_morph(xyz_core_func, names_core_func, frac_lig1_func, rseed_func, morph_func):
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
        seed = xyz_anchors_func[0]
        D_seed_anch = distance.cdist([seed], xyz_anchors_func)
        lig1_ndx = D_seed_anch[0].argsort()[:for_lig1]
        lig2_ndx = list(set(indexes) - set(lig1_ndx))
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

def coat_NP(xyz_core_func, names_core_func, frac_lig1_func, xyz_lig1_func, names_lig1_func, xyz_pillars1_func, xyz_stones1_func, xyz_lig2_func, names_lig2_func, xyz_pillars2_func, xyz_stones2_func, res_lig1_func, res_lig2_func):
    #Merges xyz coordinates and names of the core and the ligands into one coated NP
    keep_rows=[]
    for i in range(len(names_core_func)):
        if names_core_func[i]=='AU' or names_core_func[i]=='ST':
            keep_rows.append(i)

    xyz_coated_func=xyz_core_func[keep_rows,:]
    names_coated_func=names_core_func[keep_rows]
    res_coated_func=names_core_func[keep_rows]

    #Transforms and appends rototranslated ligand 1
    xyz_lig1_func_conv=np.insert(xyz_lig1_func, 3, 1, axis=1).T
    for i in range(len(xyz_stones1_func[:,0,0])):
        xyz_stones_now = xyz_stones1_func[i,:,:]
        trans_matrix=affine_matrix_from_points(xyz_pillars1_func.T, xyz_stones_now.T, shear=False, scale=False, usesvd=True)
        trans_lig=np.dot(trans_matrix, xyz_lig1_func_conv).T[:,:3]

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

    res=0
    at=0
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
