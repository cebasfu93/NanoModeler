import numpy as np

def get_size(fname):
    pdb = np.genfromtxt(fname, delimiter='\n', dtype=str)
    xyz = []
    for line in pdb:
        if "AU" in line:
            xyz.append([line[-8:], line[-16:-8], line[-24:-16]])
    xyz = np.array(xyz, dtype='float')
    sizes = np.subtract(np.max(xyz, axis=0), np.min(xyz, axis=0))
    diam = np.mean(sizes)/10
    print("{} has a diameter of {:.2f} nm.".format(fname, diam))

files = ["au25SR18_NM.pdb", \
        "au36SR24_NM.pdb", \
        "au38SR24_NM.pdb", \
        "au44SR28_NM.pdb", \
        "au68SR32-I1_NM.pdb", \
        "au68SR32-I2_NM.pdb", \
        "au68SR32-I3_NM.pdb", \
        "au68SR32-I4_NM.pdb", \
        "au68SR34-I1_NM.pdb", \
        "au68SR34-I2_NM.pdb", \
        "au68SR34-I3_NM.pdb", \
        "au68SR34-I4_NM.pdb", \
        "au102SR44_NM.pdb", \
        "au133SR52_NM.pdb", \
        "au144SR60_NM.pdb", \
        "au314SR96_NM.pdb"]

for f in files:
    get_size(f)
