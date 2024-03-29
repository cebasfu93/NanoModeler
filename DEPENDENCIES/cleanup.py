import os
import logging

logger = logging.getLogger('nanomodeler')
logger.addHandler(logging.NullHandler())

report = logging.getLogger('nanomodeler.report')

def cleanup_error(logFileHandler):
    logger.warning('Cleaning up...')
    bye = ['ANTECHAMBER.FRCMOD', 'leap.log', 'md.mdp', 'em.mdp']
    for i in bye:
        if os.path.isfile(i):
            try:
                os.remove(i)
            except OSError as e:
                rm_txt = "ATTENTION! An error occurred while removing file: {} (details: {})".format(i,str(e))
                logger.error(rm_txt)
                report.error(rm_txt)
                
    err_txt = "NanoModeler terminated with errors."
    logger.warning(err_txt)
    report.warning(err_txt)

    logFileHandler.flush()
    logFileHandler.close()

    report.handlers.remove(logFileHandler);

    return (0,0)
