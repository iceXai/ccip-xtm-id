#!/usr/bin/env python3
#coding: utf-8
"""
@author: Dr. Stephan Paul (AWI/iceXai; stephan.paul@awi.de)
"""

# In[]

import os

from .cfg import Configuration

# In[]


class Validator(object):
    def __init__(self, cfg: Configuration) -> None:
        self.cfg = cfg
        
    def validate(self) -> None:
        VALID = []
        VALID.append(self.validate_carrier_combination())
        VALID.append(self.validate_l1p_pathing())
        VALID.append(self.validate_l2i_pathing())
        VALID.append(self.validate_aoi())
        if not all(VALID):
            sys.exit()
    
    def validate_carrier_combination(self) -> bool:
        #get configuration settings
        CARRIER1 = self.cfg.carrier1
        CARRIER2 = self.cfg.carrier2
        COMBO = f'{CARRIER1}/{CARRIER2}'
        #valid combinations
        VALID_COMBOS = ['cryosat2/envisat','envisat/ers2','ers2/ers1']
        #validate choices
        if COMBO not in VALID_COMBOS:
            msg = f'{CARRIER1}/{CARRIER2} combination is invalid/unsupported!'
            logger.critical(msg)
            return False
        else:
            return True
    
    def validate_l1p_pathing(self) -> bool:
        INVALID = False
        for yy, mm in zip(self.cfg.yy, self.cfg.mm):
            PATH = self.cfg.path_to_carrier1_l1p(yy, mm)
            if not os.path.isdir(PATH):
                msg = f'L1p pathing for {CARRIER1} is incorrect for {yy}/{mm}'
                logger.critical(msg)
                INVALID = True
            PATH = self.cfg.path_to_carrier2_l1p(yy, mm)
            if not os.path.isdir(PATH):
                msg = f'L1p pathing for {CARRIER2} is incorrect for {yy}/{mm}'
                logger.critical(msg)
                INVALID = True
        if INVALID:
            return False            
        else:
            return True
        
    def validate_l2i_pathing(self) -> bool:
        INVALID = False
        for yy, mm in zip(self.cfg.yy, self.cfg.mm):
            PATH = self.cfg.path_to_carrier1_l2i(yy, mm)
            if not os.path.isdir(PATH):
                msg = f'L2i pathing for {CARRIER1} is incorrect for {yy}/{mm}'
                logger.critical(msg)
                INVALID = True
            PATH = self.cfg.path_to_carrier2_l2i(yy, mm)
            if not os.path.isdir(PATH):
                msg = f'L2i pathing for {CARRIER2} is incorrect for {yy}/{mm}'
                logger.critical(msg)
                INVALID = True
        if INVALID:
            return False            
        else:
            return True

    def validate_aoi(self) -> bool:
        USER_AOI = self.cfg.aoi        
        #validate choice
        if USER_AOI != 'arc' and USER_AOI != 'ant':
            return False
        else:
            return True