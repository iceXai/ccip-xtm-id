# -*- coding: utf-8 -*-
"""
@author: Dr. Stephan Paul (AWI/iceXai; stephan.paul@awi.de)
"""

# In[]
from loguru import logger

import yaml
import datetime
import os
import sys
import importlib


# In[]
"""
Configuration
"""

class Configuration(object):
    def __init__(self, cfg_file_name: str):
        #status
        logger.info('Load configuration file...')
        #loads the yaml file and reads it content
        path = os.path.join(os.getcwd(), 'cfg', cfg_file_name)
        with open(path) as f:
            self.config = yaml.safe_load(f)
        #configure the logger
        self.configure_logger()
            
    """ Logger Setup """
    def configure_logger(self) -> None:
        path = self.output_path
        LOG_PATH = os.path.join(path,'log')
        if not os.path.isdir(LOG_PATH):
            os.makedirs(LOG_PATH)
        CARRIER1 = self.carrier1
        CARRIER2 = self.carrier2
        AOI = self.aoi
        MATCH = self.matchtype
        BUFFER = self.buffersize
        DT = self.delta_t
        NAME = f'{CARRIER1}_{CARRIER2}_{MATCH}_{DT}h_{BUFFER}m_{AOI}.log'
        logger.add(f'{LOG_PATH}/{NAME}')
    
    """ General """
    @property
    def carrier1(self):
        return self.config['carrier']['reference'].lower()
    
    @property
    def carrier2(self):
        return self.config['carrier']['match'].lower()
    
    @property
    def aoi(self):
        return self.config['aoi'].lower()
    
    @property
    def parameter(self):
        return self.config['parameter']       
    
    """ I/O """
    @property
    def output_path(self):
        return self.config['output']
    
    @property
    def input_l1p(self):
        return self.config['input']['l1p']
    
    @property
    def input_l2i(self):
        return self.config['input']['l2i']
    
    @property
    def override(self):
        return self.config['override']['l2i']
    
    @property
    def l1p_version_carrier1(self):
        return self.config['version']['reference']['l1p']
    
    @property
    def l2i_version_carrier1(self):
        return self.config['version']['reference']['l2i']
    
    @property
    def l1p_version_carrier2(self):
        return self.config['version']['match']['l1p']
    
    @property
    def l2i_version_carrier2(self):
        return self.config['version']['match']['l2i']
    
    """ Date """
    @property
    def yy(self):
        return self.config['date']['year']
    
    @property
    def mm(self):
        return self.config['date']['month']
    
    """ Matching """
    @property
    def matchtype(self):
        return self.config['matching']['type'].lower()
    
    @property
    def buffersize(self):
        return self.config['matching']['buffer']
    
    @property
    def delta_t(self):
        return self.config['matching']['dt']    