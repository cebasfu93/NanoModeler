import numpy as np
import sys
from sklearn.decomposition import PCA

def rewrite_mol2(fname, cap, oname, log):
    mol2 = np.genfromtxt(fname, delimiter='\n', dtype='str')
    out = open(oname, "w")
    if cap=="N":
        log += "There are no capping atoms...\n"
        cap="0"
    else:
        log += "Capping atoms {} will be removed...\n".format(cap)
    ATOM=False
    BOND=False
    atoms = []
    old_at_num = []
    res_names = []
    bonds = []
    charge_cap = []
    cap = np.array(cap.split(","), dtype="int")
    for i in range(len(mol2)):
        if "@<TRIPOS>MOLECULE" in mol2[i]:
            log += "Reading @<TRIPOS>MOLECULE section...\n"
            mol_name = mol2[i+1].split()[0]
        elif "@<TRIPOS>RESIDUECONNECT"in mol2[i]:
            log += "Reading @<TRIPOS>RESIDUECONNECT section...\n"
            connect = mol2[i+1]
        elif "@<TRIPOS>BOND" in mol2[i]:
            log += "Reading @<TRIPOS>BOND section...\n"
            BOND=True
            ATOM=False
        elif BOND:
            at1 = int(mol2[i].split()[1])
            at2 = int(mol2[i].split()[2])
            if at1 not in cap and at2 not in cap:
                bonds.append(np.array(mol2[i].split()))
            if "@<TRIPOS>" in mol2[i+1]:
                BOND=False
        elif "@<TRIPOS>ATOM" in mol2[i]:
            log += "Reading @<TRIPOS>ATOM section...\n"
            ATOM=True
        elif ATOM:
            at = int(mol2[i].split()[0])
            if at not in cap:
                old_at_num.append(at)
                res_names.append(mol2[i].split()[7])
                atoms.append(np.array(mol2[i].split()))
            else:
                charge_cap.append(float(mol2[i].split()[8]))

    if len(atoms)!=np.unique(np.array(res_names), return_counts=True)[1][0]:
        sys.exit("There seems to be more than one residue type in the input mol2 file")

    xyz = []
    names = np.array([])
    for atom in atoms:
        names = np.append(names, atom[1])
        xyz.append(np.array([float(atom[2]), float(atom[3]), float(atom[4])]))
    xyz = np.array(xyz)
    anch_name = connect.split()[1]
    anch_ndx = np.where(names == anch_name)[0][0]
    xyz_anch = xyz[anch_ndx]
    xyz = np.subtract(xyz, xyz_anch)
    pca = PCA(n_components=3)
    pca.fit(xyz)
    pca1 = pca.components_[0]
    if np.sum(np.mean(xyz, axis=0)>=0)<2:
        pca1=-pca1
    new_pt = pca1/np.linalg.norm(pca1)*1.8
    xyz = np.append(xyz, np.array([new_pt]), axis=0)

    atoms.append(['0', 'ST', str(new_pt[0]), str(new_pt[1]), str(new_pt[2]), "S", atoms[anch_ndx][6], atoms[anch_ndx][7], "0.0"])
    old_at_num.append(len(atoms))
    bonds.append([str(len(bonds)+1), str(len(atoms)), str(anch_ndx+1), "1"])

    N_at = len(atoms)
    N_bo = len(bonds)
    charge_per_atom = np.sum(charge_cap)/(N_at-1)

    log += "The capping group has a total charge of {:.3f}...\n".format(np.sum(charge_cap))
    log += "A charge of {:.3f} will be added to all atoms in the ligand...\n".format(charge_per_atom)
    log += "Writing @<TRIPOS>MOLECULE section..."
    out.write("@<TRIPOS>MOLECULE\n{}\n\t{}\t{}\t1\nSMALL\nUSER_CHARGES\n".format(mol_name, N_at, N_bo))

    log += "Writing @<TRIPOS>ATOM section...\n"
    out.write("@<TRIPOS>ATOM\n")
    at = 0
    for atom in atoms:
        at+=1
        if at != len(atoms)-1:
            out.write("{0:>4} {1:>4} {2:>13.4f} {3:>9.4f} {4:>9.4f} {5:>4} {6} {7} {8:>7.4f}\n".format(\
            at, atom[1], xyz[at-1,0], xyz[at-1,1], xyz[at-1,2], atom[5], atom[6], atom[7], float(atom[8])+charge_per_atom))
        else:
            out.write("{0:>4} {1:>4} {2:>13.4f} {3:>9.4f} {4:>9.4f} {5:>4} {6} {7} {8:>7.4f}\n".format(\
            at, atom[1], xyz[at-1,0], xyz[at-1,1], xyz[at-1,2], atom[5], atom[6], atom[7], float(atom[8])+charge_per_atom))

    log += "Writing @<TRIPOS>BOND section...\n"
    out.write("@<TRIPOS>BOND\n")
    bo = 0
    old_at_num = np.array(old_at_num)
    for bond in bonds:
        bo+=1
        new_at1 = np.where(old_at_num==int(bond[1]))[0][0]+1
        new_at2 = np.where(old_at_num==int(bond[2]))[0][0]+1
        out.write("{0:>5} {1:>5} {2:>5} {3:>2}\n".format(bo, str(new_at1), str(new_at2), bond[3]))

    log += "Writing @<TRIPOS>SUBSTRUCTURE section...\n"
    out.write("@<TRIPOS>SUBSTRUCTURE\n")
    out.write("\t1 {}\t\t\t1 ****\t\t\t0 ****  ****\n".format(atoms[0][7]))
    log += "Writing @<TRIPOS>RESIDUECONNECT section...\n"
    out.write("@<TRIPOS>RESIDUECONNECT\n")
    out.write(connect+"\n")
    out.close()
