# -*- coding: utf-8 -*-
"""
@author: Dr. Stephan Paul (AWI/iceXai; stephan.paul@awi.de)
"""

# In[]
from loguru import logger
from typing import Dict, List

from .par import Parameters

import yaml
import os
import sys
import shutil

# In[]

class Configuration(object):
    def __init__(self, cfg_file_name: str):
        #configure the logger
        self.configure_logger()
        #status
        logger.info('Load configuration file...')
        #loads the yaml file and reads it content
        path = os.path.join(os.getcwd(), 'cfg', cfg_file_name)
        with open(path) as f:
            self.config = yaml.safe_load(f)
        #initiate parameters
        self.par = Parameters(self.user_parameter)
        #store configuration file for later reference
        self._copy_cfg_file(cfg_file_name)
            
    """ Logger Setup """
    def configure_logger(self) -> None:
        path = self._output_path
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
    def carrier1_tag(self):
        return self.carrier_tag(self.carrier1)
    
    @property
    def carrier2_tag(self):
        return self.carrier_tag(self.carrier2)
    
    def carrier_tag(self, carrier: str) -> str:
        if carrier == 'cryosat2':
            return 'cs2'
        if carrier == 'envisat':
            return 'env'
        if carrier == 'ers2' or carrier == 'ers1':
            return carrier
    
    @property
    def aoi(self):
        return self.config['aoi'].lower()
    
    @property
    def user_parameter(self):
        return self.config['parameter'] 

    @property
    def hemisphere_long(self):
        if self.aoi=='arc':
            return 'north'
        else:
            return 'south'
    
    @property
    def hemisphere_short(self):
        if self.aoi=='arc':
            return 'nh'
        else:
            return 'sh'
    
    @property
    def epsg(self):
        if self.aoi=='arc':
            return 6931
        else:
            return 6932
            
    """ I/O """
    @property
    def _output_path(self):
        return self.config['output']
    
    @property
    def _input_l1p(self):
        return self.config['input']['l1p']
    
    @property
    def _input_l2i(self):
        return self.config['input']['l2i']
    
    @property
    def override(self):
        return self.config['override']
    
    @property
    def _l1p_version_carrier1(self):
        return self.config['version']['reference']['l1p']
    
    @property
    def _l2i_version_carrier1(self):
        return self.config['version']['reference']['l2i']
    
    @property
    def _l1p_version_carrier2(self):
        return self.config['version']['match']['l1p']
    
    @property
    def _l2i_version_carrier2(self):
        return self.config['version']['match']['l2i']
    
    """ Pathing """
    def path_to_carrier1_l1p(self, year: str, month: str) -> str:
        CARRIER = self.carrier1
        VERSION = self._l1p_version_carrier1
        DATE = os.path.join(year, month)
        return self._path_to_l1p(CARRIER, VERSION, DATE)
    
    def path_to_carrier2_l1p(self, year: str, month: str) -> str:
        CARRIER = self.carrier2
        VERSION = self._l1p_version_carrier2
        DATE = os.path.join(year, month)
        return self._path_to_l1p(CARRIER, VERSION, DATE)
    
    def _path_to_l1p(self, carrier: str, version: str, date: str) -> str:
        INPATH = self._input_l1p
        if carrier == 'cryosat2':
            PRODUCT = 'ipf1-e/rep/l1p'
        else:
            PRODUCT = 'l1p'
        HEMISPHERE = self.hemisphere_long
        return os.path.join(INPATH, carrier, PRODUCT,
                            version, HEMISPHERE, date)
    
    def path_to_carrier1_l2i(self, year: str, month: str) -> str:
        CARRIER = self.carrier1
        VERSION = self._l2i_version_carrier1
        DATE = os.path.join(year, month)
        return self._path_to_l2i(CARRIER, VERSION, DATE)
    
    def path_to_carrier2_l2i(self, year: str, month: str) -> str:
        CARRIER = self.carrier2
        VERSION = self._l2i_version_carrier2
        DATE = os.path.join(year, month)
        return self._path_to_l2i(CARRIER, VERSION, DATE)
    
    def _path_to_l2i(self, carrier: str, version: str, date: str) -> str:
        INPATH = self._input_l2i
        HEMISPHERE = self.hemisphere_short
        PRODUCT = 'l2i' 
        return os.path.join(INPATH, carrier, version, 
                            HEMISPHERE, PRODUCT, date)
        
    """ Output files """
    def _copy_cfg_file(self, cfg_file: str) -> None:
        INPATH = os.path.join(os.getcwd(), 'cfg', cfg_file)
        OUTPATH = self._output_path
        ARCHIVE = 'archive'
        DESTPATH = os.path.join(OUTPATH, ARCHIVE)
        if not os.path.isdir(DESTPATH):
            os.makedirs(DESTPATH)
        shutil.copy(INPATH, DESTPATH)
    
    def _output_shp_name(self, year: str, month: str) -> str:
        MATCHTYPE = self.matchtype
        CARRIER1 = self.carrier1
        CARRIER2 = self.carrier2
        AOI = self.aoi
        BUFFER = f'{str(self.buffersize).zfill(5)}m'
        return f'{MATCHTYPE}_{CARRIER1}_{CARRIER2}_'+\
            f'{year}_{month}_{AOI}_{BUFFER}.shp'
        
    def output_to_shp(self, year: str, month: str) -> str:
        OUTPATH = self._output_path
        FILENAME = self._output_shp_name(year, month)
        return os.path.join(OUTPATH, FILENAME)

    def _output_carrier_csv_name(self, subtype: str, year: str, month: str, 
                                 product: str, dt: int) -> str:
        MATCHTYPE = self.matchtype
        AOI = self.aoi
        BUFFER = f'{str(self.buffersize).zfill(5)}m'
        return f'{MATCHTYPE}_{year}_{month}_{AOI}_{BUFFER}_{subtype}'+\
            f'_{product}_dt{dt+1}.csv'
            
    def _output_meta_csv_name(self, year: str, month: str, dt: int) -> str:
        MATCHTYPE = self.matchtype
        AOI = self.aoi
        BUFFER = f'{str(self.buffersize).zfill(5)}m'
        return f'{MATCHTYPE}_{year}_{month}_{AOI}_{BUFFER}_meta_dt{dt+1}.csv'
        
    def output_to_csv(self, year: str, month: str, dt: int) -> Dict[str, str]:
        OUTPATH = self._output_path
        CSVDICT = {}
        #carrier1 data
        TYPE = self.carrier1
        CSVDICT[TYPE] = {}
        PRODUCT = 'l1p'
        NAME = self._output_carrier_csv_name(TYPE, year, month, PRODUCT, dt)
        CSVDICT[TYPE][PRODUCT] = os.path.join(OUTPATH, NAME)
        PRODUCT = 'l2i'
        NAME = self._output_carrier_csv_name(TYPE, year, month, PRODUCT, dt)
        CSVDICT[TYPE][PRODUCT] = os.path.join(OUTPATH, NAME)
        #carrier2 data
        TYPE = self.carrier2
        CSVDICT[TYPE] = {}
        PRODUCT = 'l1p'
        NAME = self._output_carrier_csv_name(TYPE, year, month, PRODUCT, dt)
        CSVDICT[TYPE][PRODUCT] = os.path.join(OUTPATH, NAME)
        PRODUCT = 'l2i'
        NAME = self._output_carrier_csv_name(TYPE, year, month, PRODUCT, dt)
        CSVDICT[TYPE][PRODUCT] = os.path.join(OUTPATH, NAME)
        #metadata
        NAME = self._output_meta_csv_name(year, month, dt)
        CSVDICT['meta'] = os.path.join(OUTPATH, NAME)
        #return to caller
        return CSVDICT
    
    def shp_exists(self, year: str, month: str) -> bool:
        FILENAME = self.output_to_shp(year, month)
        return os.path.isfile(FILENAME)
        
        
    """ Date """
    @property
    def yy(self):
        return [str(y) for y in self.config['date']['year']]
    
    @property
    def mm(self):
        return [str(m).zfill(2) for m in self.config['date']['month']]
    
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
    
