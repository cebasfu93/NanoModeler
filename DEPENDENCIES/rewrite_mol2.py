import numpy as np

def rewrite_mol2(fname, cap, oname):
    mol2 = np.genfromtxt(fname, delimiter='\n', dtype='str')
    out = open(oname, "w")
    if cap=="N":
        for i in range(len(mol2)):
            out.write(mol2[i]+"\n")
        out.close()
    else:
        ATOM=False
        BOND=False
        atoms = []
        old_at_num = []
        bonds = []
        charge_cap = []
        cap = np.array(cap.split(","), dtype="int")
        for i in range(len(mol2)):
            print(mol2[i])
            if "@<TRIPOS>MOLECULE" in mol2[i]:
                mol_name = mol2[i+1].split()[0]
            elif "@<TRIPOS>RESIDUECONNECT"in mol2[i]:
                connect = mol2[i+1]
            elif "@<TRIPOS>BOND" in mol2[i]:
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
                ATOM=True
            elif ATOM:
                at = int(mol2[i].split()[0])
                if at not in cap:
                    old_at_num.append(at)
                    atoms.append(np.array(mol2[i].split()))
                else:
                    charge_cap.append(float(mol2[i].split()[8]))

        N_at = len(atoms)
        N_bo = len(bonds)
        charge_per_atom = np.sum(charge_cap)/N_at

        out.write("@<TRIPOS>MOLECULE\n{}\n\t{}\t{}\nSMALL\nUSER_CHARGES\n".format(mol_name, N_at, N_bo))

        out.write("@<TRIPOS>ATOM\n")
        at = 0
        for atom in atoms:
            at+=1
            out.write("{0:>4} {1:>4} {2:>13.4f} {3:>9.4f} {4:>9.4f} {5:>4} {6} {7} {8:>7.4f}\n".format(\
            at, atom[1], float(atom[2]), float(atom[3]), float(atom[4]), atom[5], atom[6], atom[7], float(atom[8])+charge_per_atom))

        out.write("@<TRIPOS>BOND\n")
        bo = 0
        old_at_num = np.array(old_at_num)
        for bond in bonds:
            bo+=1
            new_at1 = np.where(old_at_num==int(bond[1]))[0][0]+1
            new_at2 = np.where(old_at_num==int(bond[2]))[0][0]+1
            out.write("{0:>5} {1:>5} {2:>5} {3:>2}\n".format(bo, str(new_at1), str(new_at2), bond[3]))

        out.write("@<TRIPOS>RESIDUECONNECT\n")
        out.write(connect+"\n")
        out.close()

#rewrite_mol2("../LIG2_CAP.mol2", "1,2,3,4", "LIG2_rewriten.mol2")
