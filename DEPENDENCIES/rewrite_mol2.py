import numpy as np
import sys
import logging
from sklearn.decomposition import PCA

logger = logging.getLogger('nanomodeler')
logger.addHandler(logging.NullHandler())


def rewrite_mol2(mol2, cap, lig_s, lig_c, oname, elong):
    # Saves new mol2 without capping atoms, with an S atom, renumbered and S atom bonded to C atom
    out = open(oname, "w")

    if not cap:
        logger.info("\tNo capping atoms were found...")
    else:
        logger.info("\tCapping atoms {} will be removed...".format(cap))

    if lig_s != 0:
        # Reads bonds section and stays with the atom bonded to the S atom that is not in the capping group
        lig_c = find_C(mol2, cap, lig_s)

    ATOM=False
    BOND=False
    atoms = []
    old_at_num = []
    res_names = []
    bonds = []
    charge_cap = []

    for i in range(len(mol2)):
        if "@<TRIPOS>MOLECULE" in mol2[i]:
            logger.info("\tReading @<TRIPOS>MOLECULE section...")
            mol_name = mol2[i+1].split()[0]
        elif "@<TRIPOS>BOND" in mol2[i]:
            logger.info("\tReading @<TRIPOS>BOND section...")
            BOND=True
            # ATOM=False
        elif BOND:
            if mol2[i].split() != []:
                at1 = int(mol2[i].split()[1])
                at2 = int(mol2[i].split()[2])
                if at1 not in cap and at2 not in cap and at1 != lig_s and at2 != lig_s:
                    # Saves all the atoms that dont involve capping atoms or the S atom
                    bonds.append(np.array(mol2[i].split()))
            if "@<TRIPOS>" in mol2[i+1]:
                BOND=False
        elif "@<TRIPOS>ATOM" in mol2[i]:
            logger.info("\tReading @<TRIPOS>ATOM section...")
            ATOM=True
        elif ATOM:
            if mol2[i].split() != []:
                at = int(mol2[i].split()[0])
                if at == lig_c:
                    # Saves the C atom
                    anch_name = mol2[i].split()[1]
                    c_atom = np.array(mol2[i].split())
                if at == lig_s:
                    # Saves the S atom
                    s_atom = np.array(mol2[i].split())
                elif at not in cap:
                    # Saves atoms not in the capping group
                    old_at_num.append(at)
                    res_names.append(mol2[i].split()[7])
                    atoms.append(np.array(mol2[i].split()))
                else:
                    # Saves the charge of the atoms in the capping group
                    charge_cap.append(float(mol2[i].split()[8]))
            if "@<TRIPOS>" in mol2[i+1]:
                ATOM=False
    if lig_s != 0 and "S" != s_atom[1][0].upper():
        raise ValueError(f"The chosen Sulfur atom (atom number {lig_s}) does not start with 'S' ({s_atom[1]})")
    if "C" != c_atom[1][0].upper():
        raise ValueError(f"The C atom found (atom number {lig_c}) does not start with 'C' ({c_atom[1]})")


    xyz = []
    names = np.array([])
    for atom in atoms:
        # If the xyz columns are not float, python will exit
        names = np.append(names, atom[1])
        xyz.append(np.array([float(atom[2]), float(atom[3]), float(atom[4])]))
    xyz = np.array(xyz)

    # If there is no S or if elong is set to true, the S atom is placed 1.8A away from the C atom in the opposite direction of PCA1
    if elong or not lig_s:
        anch_ndx = np.where(names == anch_name)[0][0]
        xyz_anch = xyz[anch_ndx]
        xyz = np.subtract(xyz, xyz_anch)
        pca = PCA(n_components=3)
        pca.fit(xyz)
        pca1 = pca.components_[0]
        # Sklearn always outputs the PCA with two positive components, in this way, the PCA is aligned to the coordinates
        if np.sum(np.mean(xyz, axis=0)>=0)<2:
            pca1=-pca1
        new_pt = pca1/np.linalg.norm(pca1)*1.8
    else:
        new_pt = np.array(s_atom[2:5], dtype='float')

    xyz = np.append(xyz, np.array([new_pt]), axis=0)
    old_at_num.append(len(atoms))
    old_at_num = np.array(old_at_num, dtype='int')

    # +1 includes the S atom that is placed in the following if statement
    N_at = len(atoms) + 1
    # +1 includes the S-C bond that is placed after the following if statement
    N_bo = len(bonds) + 1
    # The S atom is added to the atoms list
    if not lig_s:
        # The charge of the capping atom is divided in the rest of the atoms (excluding the S atom)
        charge_per_atom = np.sum(charge_cap)/(N_at-1)
        for atom in atoms:
            atom[8] = str(float(atom[8])+charge_per_atom)
        atoms.append(['0', 'ST', str(new_pt[0]), str(new_pt[1]), str(new_pt[2]), 'S', atoms[0][6], atoms[0][7], "0.0"])
    else:
        # The charge of the capping atom is divided in the rest of the atoms (including the S atom)
        charge_per_atom = np.sum(charge_cap)/(N_at)
        for atom in atoms:
            atom[8] = str(float(atom[8])+charge_per_atom)
        atoms.append(['0', 'ST', str(new_pt[0]), str(new_pt[1]), str(new_pt[2]), 'S', s_atom[6], s_atom[7], float(s_atom[8])+charge_per_atom])
    # The S atom is bonded to the C atom
    bonds.append(['0', str(len(atoms)), str(np.where(old_at_num==lig_c)[0][0]+1), 1])

    logger.info("\tThe capping group has a total charge of {:.3f}...".format(np.sum(charge_cap)))
    logger.info("\tA charge of {:.3f} will be added to all atoms in the ligand...".format(charge_per_atom))
    logger.info("\tWriting @<TRIPOS>MOLECULE section...")

    out.write("@<TRIPOS>MOLECULE\n{}\n\t{}\t{}\t1\nSMALL\nUSER_CHARGES\n".format(mol_name, N_at, N_bo))

    logger.info("\tWriting @<TRIPOS>ATOM section...")

    out.write("@<TRIPOS>ATOM\n")
    at = 0
    # The atoms are renumbered and written in the format of mol2 file (format from internet)
    for atom in atoms:
        at+=1
        out.write("{0:>4} {1:>4} {2:>13.4f} {3:>9.4f} {4:>9.4f} {5:>4} {6} {7} {8:>7.4f}\n".format(\
        at, atom[1], xyz[at-1,0], xyz[at-1,1], xyz[at-1,2], atom[5], atom[6], atom[7], float(atom[8])))

    logger.info("\tWriting @<TRIPOS>BOND section...")

    out.write("@<TRIPOS>BOND\n")
    bo = 0

    # The bonds are renumbered and written in the format of mol2 file (format from internet)
    for bond in bonds:
        bo+=1
        if bo != len(bonds):
            new_at1 = np.where(old_at_num==int(bond[1]))[0][0]+1
            new_at2 = np.where(old_at_num==int(bond[2]))[0][0]+1
            out.write("{0:>5} {1:>5} {2:>5} {3:>2}\n".format(bo, str(new_at1), str(new_at2), bond[3]))
        else:
            out.write("{0:>5} {1:>5} {2:>5} {3:>2}\n".format(bo, str(bond[1]), str(bond[2]), bond[3]))

    logger.info("\tWriting @<TRIPOS>SUBSTRUCTURE section...")
    out.write("@<TRIPOS>SUBSTRUCTURE\n")
    out.write("\t1 {}\t\t\t1 ****\t\t\t0 ****  ****\n".format(atoms[0][7]))
    out.close()

    # The C atom numbered is updated according to the new numbering
    new_lig_c = np.where(old_at_num==lig_c)[0][0]+1
    return new_lig_c

def find_C(mol2, cap, lig_s):
    BOND = False
    for i in range(len(mol2)):
        if BOND:
            if mol2[i].split() != []:
                at1 = int(mol2[i].split()[1])
                at2 = int(mol2[i].split()[2])
                if at1 not in cap and at2 not in cap:
                    if at1 == lig_s:
                        found_c = at2
                    elif at2 == lig_s:
                        found_c = at1
            if "@<TRIPOS>" in mol2[i+1]:
                break
        if "@<TRIPOS>BOND" in mol2[i]:
            BOND = True
    return found_c
