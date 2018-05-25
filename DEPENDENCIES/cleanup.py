import os

def cleanup_error(TMP, log):
    log += "Cleaning up...\n"
    bye = ["ANTECHAMBER.FRCMOD", "leap.log", "md.mdp", "em.mdp"]
    for i in bye:
        if os.path.isfile(i):
            os.remove(i)
    #os.system("rm -r " + TMP)
    log += "NanoModeler terminated with errors.\n"
    print(log)
    return (0, log, 0)
