# -*- coding: utf-8 -*-
"""
@author: Dr. Stephan Paul (AWI/iceXai; stephan.paul@awi.de)
"""

# In[]

import os
import yaml

import numpy as np

from typing import Dict, List
from loguru import logger


# In[]

class Parameters:
    def __init__(self, user_parameters: list):
        #status
        logger.info('Load supported parameters...')
        #loads the yaml file and reads it content
        PARFILE = 'list_of_supported_parameters.yaml'
        path = os.path.join(os.getcwd(), 'cfg', PARFILE)
        with open(path) as f:
            self.par = yaml.safe_load(f)
        #add source to l1p 
        self.par['l1p']['src'] = ['', 'mission_data_source']
        #store user parameters
        self.usr_par = user_parameters
        
    @property
    def l1p_groups(self) -> np.array:
        return np.unique([self.l1p_parameters[key][0] 
                          for key in self.l1p_parameters.keys()])
    
    def parameters_by_group(self, group: str) -> list:
        return [key for key in self.l1p_parameters.keys() 
                if self.l1p_parameters[key][0] == group]

    def validate_cfg_parameters(self, parameters: list) -> bool:
        #store for later use
        self.cfg_parameters = parameters
        #validate input arguments for support
        all_valid = True
        if parameters is None:
            #status
            logger.critical(f'No parameters for CSV output specified!')
            return False
        else:
            for par in parameters: 
                NOT_L1P = par not in self.par['l1p']
                NOT_L2I = par not in self.par['l2i']
                if NOT_L1P and NOT_L2I:
                    #status
                    logger.critical(f'Specified parameter {par} currently '+\
                                    f'not supported!')
                    all_valid = False
            #return status
            return all_valid
        
    @property
    def l1p_parameters(self) -> Dict[str, str]:
        #set mandatory file source variable/path
        pdict = {}
        pdict['src'] = self.par['l1p']['src']
        #loop over params and add them
        for par in self.usr_par:
            l1pdict = self.par['l1p']
            if par in l1pdict.keys():
                pdict[par] = l1pdict[par]
        #return to caller
        return pdict
    
    @property
    def l2i_parameters(self) -> Dict[str, str]:
        #set mandatory file source variable/path
        pdict = {}
        #loop over params and add them
        for par in self.usr_par:
            l2idict = self.par['l2i']
            if par in l2idict.keys():
                pdict[par] = l2idict[par]
        #return to caller
        return pdict

