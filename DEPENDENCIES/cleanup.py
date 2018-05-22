import os

def cleanup_error(TMP, log):
    log += "Cleaning up...\n"
    bye = ["ANTECHAMBER.FRCMOD", "leap.log", "md.mdp", "em.mdp", "acpype.log"]
    for i in bye:
        os.system("rm " + str(i))
    #os.system("rm -r " + TMP)
    log += "NanoModeler terminated with errors.\n"
    print(log)
    return (0, log, 0)

def cleanup_normal(log):
    print(log)
    return (1, log)
