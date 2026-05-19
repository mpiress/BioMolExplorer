from kernel.header_builder import HeaderBuilder

__doc__ = HeaderBuilder.build(

    module_title="CMolD: Collection of Molecular Descriptors",

    module_description=(
        "Module responsible for defining framework logger patterns."
        "Used by mapping warnings, errors, and execution performance of BioMolExplorer."
    ),

    module_version="1.0.3"
)
#----------------------------------------------------------------------------------------------
import warnings
warnings.filterwarnings("ignore")
#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
import logging
import os
from pathlib import Path
#----------------------------------------------------------------------------------------------


class LoggerManager:
    _loggers = {}

    @classmethod
    def get_logger(cls, name, log_file=None, level=logging.ERROR):
        path    = str(Path.cwd())
        logpath = path + '/logs/'
        
        if not os.path.exists(logpath):
            os.makedirs(logpath, exist_ok=True)

        if name not in cls._loggers:
            logger = logging.getLogger(name)
            logger.setLevel(level)

            if log_file:
                file_handler = logging.FileHandler(log_file)
                file_handler.setLevel(level)
                formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
                file_handler.setFormatter(formatter)
                logger.addHandler(file_handler)

            cls._loggers[name] = logger

        return cls._loggers[name]
