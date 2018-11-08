import os
import logging

logger = logging.getLogger('nanomodeler')
logger.addHandler(logging.NullHandler())

def cleanup_error(TMP, rep):
    logger.warning('Cleaning up...')
    bye = ['ANTECHAMBER.FRCMOD', 'leap.log', 'md.mdp', 'em.mdp']
    for i in bye:
        if os.path.isfile(i):
            try:
                os.remove(i)
            except OSError as e:
                rm_txt = "ATTENTION! An error occurred while removing file: {} (details: {})".format(i,str(e))
                logger.error(rm_txt)
                rep.write(rm_txt + "\n")
    err_txt = "NanoModeler terminated with errors."
    logger.warning(err_txt)
    rep.write(err_txt + "\n")
    rep.close()
    return (0,0)
