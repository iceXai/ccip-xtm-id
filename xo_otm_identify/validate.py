#!/usr/bin/env python3
#coding: utf-8
"""
@author: Dr. Stephan Paul (AWI/iceXai; stephan.paul@awi.de)
"""

# In[]


# In[]


class Validator(object):
    def __init__(self, cfg: Configuration) -> None:
        self.cfg = cfg
        
    def validate(self) -> None:
        VALID = []
        VALID.append(self.validate_sensor_carrier_combination())
        VALID.append(self.validate_version())
        VALID.append(self.validate_aois())
        if not all(VALID):
            sys.exit()
    
    def validate_sensor_carrier_combination(self) -> bool:
        #get configuration settings
        SENSOR = self.cfg.sensor  
        CARRIER = self.cfg.carrier
        COMBO = f'{CARRIER}/{SENSOR}'
        #valid combinations
        VALID_COMBOS = ['terra/modis','aqua/modis','s3a/olci','s3a/slstr',
                        's3b/olci','s3b/slstr','snpp/viirs','jpss1/viirs']
        #validate choices
        if COMBO not in VALID_COMBOS:
            msg = f'Combination of {CARRIER}/{SENSOR} is invalid/unsupported!'
            logger.critical(msg)
            return False
        else:
            return True
    
    def validate_version(self) -> bool:
        VERSION = self.cfg.version
        SENSOR = self.cfg.sensor  
        META_FILE = f'{SENSOR}_{VERSION}.yaml'
        META_PATH = os.path.join(os.getcwd(), 'meta', META_FILE)
        #validate existence
        if not os.path.isfile(META_PATH):
            msg = f'Version {VERSION} is incorrect, no corresponding meta '+\
                f'could file found!'
            logger.critical(msg)
            return False
        else:
            return True

    def validate_aois(self) -> bool:
        USER_AOIS = self.cfg.user_aois
        with open(os.path.join(os.getcwd(), 'aoi', 'list_of_aois.yaml')) as f:
            AOIS = yaml.safe_load(f)
        INVALID = False
        #validate choice
        for aoi in USER_AOIS:
            if aoi not in AOIS['aois'].keys():
                msg = f'Choice of {aoi} is not supported!'
                logger.critical(msg)
                INVALID = True
        if INVALID:
            return False
        else:
            return True