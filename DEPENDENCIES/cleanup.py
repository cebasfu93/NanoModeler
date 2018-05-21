import os
import shutil

def cleanup_error(TMP, log):
    log += "Cleaning up...\n"
    bye = ["ANTECHAMBER.FRCMOD", "leap.log", "md.mdp", "em.mdp", "acpype.log"]
    for i in bye:
        os.system("rm " + str(i))
    #os.system("rm -r TMP")
    log += "NanoModeler terminated with errors.\n"
    print(log)
