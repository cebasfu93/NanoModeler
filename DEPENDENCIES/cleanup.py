import os
import shutil

def cleanup_error():
    print("Cleaning up...")
    bye = ["ANTECHAMBER.FRCMOD", "leap.log", "md.mdp", "em.mdp", "acpype.log"]
    for i in bye:
        os.system("rm " + str(i))
    os.system("rm -r TMP")
    print("NanoModeler terminated with errors.")


def cleanup_normal(VAR, log):
    print("Cleaning...")
    bye = ["ANTECHAMBER.FRCMOD", "leap.log", "md.mdp", "em.mdp", "acpype.log"]
    for i in bye:
        os.remove(i)
    shutil.rmtree("TMP")

    print("Compressing files to output...")
    print("NanoModeler terminated normally. Que gracias.")
    log.close()
    shutil.copyfile("NanoModeler.log", VAR["NAME"]+"/NanoModeler.log")
    os.remove("NanoModeler.log")

    os.system("tar -zcvf {}.tar.gz {}".format(VAR["NAME"], VAR["NAME"]))
    shutil.rmtree(VAR["NAME"])
