import os
import shutil

def cleanup_error(log):
    log += "Cleaning up...\n"
    bye = ["ANTECHAMBER.FRCMOD", "leap.log", "md.mdp", "em.mdp", "acpype.log"]
    for i in bye:
        os.system("rm " + str(i))
    #os.system("rm -r TMP")
    log += "NanoModeler terminated with errors.\n"
    print(log)


def cleanup_normal(VAR, log):
    log += "Cleaning...\n"
    bye = ["ANTECHAMBER.FRCMOD", "leap.log", "md.mdp", "em.mdp", "acpype.log"]
    for i in bye:
        os.remove(i)

    #shutil.rmtree("TMP")
    log += "Compressing files to output...\n"
    log += "NanoModeler terminated normally. Que gracias.\n"

    print(log)
    #os.system("tar -zcvf {}.tar.gz {}".format(VAR["NAME"], VAR["NAME"]))
    #shutil.rmtree(VAR["NAME"])
