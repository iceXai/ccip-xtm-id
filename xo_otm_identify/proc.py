#!/usr/bin/env python3
#coding: utf-8
"""
@author: Dr. Stephan Paul (AWI/iceXai; stephan.paul@awi.de)
"""

# In[]

import os
import pathlib

import pandas as pd
import geopandas as gpd
import xarray as xr
import numpy as np
import datetime as dt

from typing import List, Dict
from loguru import logger

from .cfg import Configuration
from .par import Parameters



# In[]

class Processor:
    """
    Class to handle the processing of the given sensor/carrier combo, i.e.,
    the compilation of parameters and waveforms from XOs/OTMs to CSV
    """
    def __init__(self, cfg: Configuration, year: str, month: str):
        self.cfg = cfg
        self.year = year
        self.month = month
        self.path_carrier1_l2i = cfg.path_to_carrier1_l2i(year, month)
        self.path_carrier2_l2i = cfg.path_to_carrier2_l2i(year, month)
        self.files_carrier1_l2i = os.listdir(self.path_carrier1_l2i)
        self.files_carrier2_l2i = os.listdir(self.path_carrier2_l2i)
    
    def run(self) -> None:
        #import geodata
        gdf = self._import_shp()
        #compile csv data
        csvdict = self._compile_csv(gdf)
        #export to csv
        self._export(csvdict)
    
    def _import_shp(self) -> gpd.GeoDataFrame:
        INNAME = self.cfg.output_to_shp(self.year, self.month)
        return gpd.read_file(INNAME)

    def _compile_csv(self, gdf: gpd.GeoDataFrame) -> Dict[str, pd.DataFrame]:
        #status
        logger.info(f'Compiling data for CSV output for '+\
                    f'{self.year}/{self.month}...')
        #allocate output
        csvdict = {}
        TYPE = self.cfg.carrier1
        csvdict[TYPE] = {'l1p': pd.DataFrame(),
                         'l2i': pd.DataFrame(),
                         }
        TYPE = self.cfg.carrier2
        csvdict[TYPE] = {'l1p': pd.DataFrame(),
                         'l2i': pd.DataFrame(),
                         }
        csvdict['meta'] = pd.DataFrame()
        #number of total matches
        N_MATCHES = len(gdf)
        #loop over identified matches
        for idx, f1, f2, f1r0, f1r1, f2r0, f2r1, dt, geom in gdf.itertuples():
            #check time difference mismatch and skip
            if dt > self.cfg.delta_t:
                continue
            
            #otherwise continue processing with importing L1p data
            path_to_f1 = os.path.join(self.cfg._input_l1p, f1)
            l1b_c1, c1_df = self._import_l1p(path_to_f1, f1r0, f1r1)
            path_to_f2 = os.path.join(self.cfg._input_l1p, f2)
            l1b_c2, c2_df = self._import_l1p(path_to_f2, f2r0, f2r1)
            
            #insert xo/otm id
            c1_df.insert(0, self.cfg.matchtype+'_idx', idx)
            c2_df.insert(0, self.cfg.matchtype+'_idx', idx)
            
            #identify l2i matches for l1p files
            paths_l2i_c1 = self.files_carrier1_l2i
            f1_l2i = self._identify_l2i_by_l1p(path_to_f1, paths_l2i_c1)
            paths_l2i_c2 = self.files_carrier2_l2i
            f2_l2i = self._identify_l2i_by_l1p(path_to_f2, paths_l2i_c2)
            
            #import L2i data
            path_to_f1_l2i = os.path.join(self.cfg._input_l2i, f1_l2i)
            c1_l2i_df = self._import_l2i(path_to_f1_l2i, f1r0, f1r1, True)
            path_to_f2_l2i = os.path.join(self.cfg._input_l2i, f2_l2i)
            c2_l2i_df = self._import_l2i(path_to_f2_l2i, f2r0, f2r1, False)
            
            #insert xo/otm id
            c1_l2i_df.insert(0, self.cfg.matchtype+'_idx', idx)
            c2_l2i_df.insert(0, self.cfg.matchtype+'_idx', idx)
            
            #append l1p/l2i to output
            TYPE = self.cfg.carrier1
            csvdict[TYPE]['l1p'] = pd.concat([csvdict[TYPE]['l1p'], c1_df], 
                                             ignore_index=True)
            csvdict[TYPE]['l2i'] = pd.concat([csvdict[TYPE]['l2i'], c1_l2i_df], 
                                             ignore_index=True)
            TYPE = self.cfg.carrier2
            csvdict[TYPE]['l1p'] = pd.concat([csvdict[TYPE]['l1p'], c2_df], 
                                             ignore_index=True)
            csvdict[TYPE]['l2i'] = pd.concat([csvdict[TYPE]['l2i'], c2_l2i_df], 
                                             ignore_index=True)
            
            #create meta data entries
            meta = {self.cfg.matchtype+'_idx': idx,
                    'l1p'+self.cfg.carrier1: f1,
                    'l1p'+self.cfg.carrier2: f2,
                    'l1b'+self.cfg.carrier1: l1b_c1,
                    'l1b'+self.cfg.carrier2: l1b_c2,
                    'l2i'+self.cfg.carrier1: f1_l2i,
                    'l2i'+self.cfg.carrier2: f2_l2i,
                    self.cfg.matchtype+'_dt': dt,
                    }
            
            #convert to a dataframe
            meta_df = pd.DataFrame([meta])
            
            #append it
            csvdict['meta'] = pd.concat([csvdict['meta'], meta_df],
                                        ignore_index=True)
        #status
        TYPE = self.cfg.carrier1
        N_L1P = len(csvdict[TYPE]['l1p'])
        N_L2I = len(csvdict[TYPE]['l2i'])
        logger.info(f'Processed L1p/L2i waveforms ({TYPE}): {N_L1P}/{N_L2I}')
        TYPE = self.cfg.carrier2
        N_L1P = len(csvdict[TYPE]['l1p'])
        N_L2I = len(csvdict[TYPE]['l2i'])
        logger.info(f'Processed L1p/L2i waveforms ({TYPE}): {N_L1P}/{N_L2I}')
        META = csvdict['meta']
        N_META = len(META)
        MATCH = self.cfg.matchtype.upper()
        logger.info(f'Used {N_META} out of {N_MATCHES} potential {MATCH}s')
        logger.info(f'Processing CSV output complete!')
        #return to caller
        return csvdict

    def _import_l1p(self, path: str, r0: int, r1: int) -> str, pd.DataFrame:
        #compile product-level specific parameter list
        parameters = self.cfg.par.l1p_parameters()
        #retrieve parameters groups for L1p
        l1p_groups = self.cfg.par.l1p_groups
        #create empty DataFrame 
        df = pd.DataFrame()
        #loop over groups in NC file
        for group in l1p_groups:
            #open file connection
            nc = xr.open_dataset(path, group = group)
            if group  == '':
                #load mandatory file source information
                l1b_source = nc.attrs[parameters['src'][1]]
                #close file connection and continue with next group
                nc.close()
                continue
            #get parameters within that group
            pars_in_group = self.cfg.par.parameters_by_group(group)
            #load parameters
            for par in pars_in_group:
                #get in-file variable name
                VAR = parameters[par][1]
                if par == 'pwr' or par == 'rng':
                    DATA = nc[VAR].values[r0:r1,:]
                    #convert to DataFrame
                    BINS = [par+'b'+str(rbin+1).zfill(3) 
                            for rbin in range(DATA.shape[1])]
                    par_df = pd.DataFrame(DATA, columns = BINS)
                    #add it to the existing one
                    df = pd.concat([df, par_df],axis=1)
                else:
                    #try to get variable
                    try:
                        DATA = {}
                        DATA[par] = nc[VAR].values[r0:r1]
                        #convert it to DataFrame
                        par_df = pd.DataFrame(DATA)
                        #add it to the existing one
                        df = pd.concat([df, par_df],axis=1)
                    except KeyError:
                        logger.warning(f'Parameter {VAR} not found in L1p'+\
                                       f' file')
                        continue
            #close file connection for next group
            nc.close()
        #return to caller
        return l1b_source, df

    def _identify_l2i_by_l1p(self, l1p_path: str, l2i_paths: str) -> str:
        #get l1p file name
        l1p = pathlib.Path(l1p_path).name
        #retrieve time/date tag from file name
        tag = l1p.split('-')[6:8]
        tag = [dt.datetime.strptime(tag[0],'%Y%m%dT%H%M%S'),
               dt.datetime.strptime(tag[1],'%Y%m%dT%H%M%S')]
        #get l2i files
        l2i_tags = [f.split('-')[6:8] for f in l2i_paths]
        #id correct l2i match
        l2i_source = None
        for idx, l2i_tag in enumerate(l2i_tags):
            #convert to datetime
            l2i_tag = [dt.datetime.strptime(l2i_tag[0],'%Y%m%dT%H%M%S'),
                       dt.datetime.strptime(l2i_tag[1][:-5],'%Y%m%dT%H%M%S')]
            #calculate time difference in total seconds
            dt_start = l2i_tag[0]-tag[0]
            total_dt_start = dt_start.total_seconds()
            dt_end = l2i_tag[1]-tag[1]
            total_dt_end = dt_end.total_seconds()
            #only keep in case dt is exact match
            if total_dt_start == 0 and total_dt_end == 0:
                #path/file name
                l2i = l2i_paths[idx]
                #load mandatory file source information
                l2i_source = pathlib.Path(l2i).name
                #break loop in case one is found
                break
        if l2i_source is None:
            logger.warning(f'No L2i match found for L1p file: {l1p}')
        #return to caller
        return l2i_source

    def _import_l2i(self, path: str, r0: int, r1: int, 
                    ref: bool) -> pd.DataFrame:
        #compile product-level specific parameter list
        parameters = self.cfg.par.l2i_parameters()
        #open file connection
        nc = xr.open_dataset(path)
        #load the data
        DATA = {}
        #load common parameters
        DATA['sft'] = nc['surface_type'].values[r0:r1]
        for par in parameters:
            #get in-file variable name
            VAR = parameters[par]
            #try to get variable
            try:
                DATA[par] = nc[VAR].values[r0:r1]
            except KeyError:
                logger.warning(f'Parameter {VAR} not found in L2i file')
                continue    
        if ref:
            #load reference sensor specific things
            DATA['frb'] = nc['sea_ice_freeboard'].values[r0:r1]
        #convert it to DataFrame
        df = pd.DataFrame(DATA)
        if no ref:
            #load match sensor specific things
            frb = nc['threshold_freeboards'].values[r0:r1,:]
            ths = nc['tfmra_thresholds'].values
            #convert it to DataFrame w/ correct column names
            frb_columns = ['frb_th'+str(format(np.round(th,2),'.2f'))[2:]
                           for th in ths]
            frb_df = pd.DataFrame(frb, columns = frb_columns)
            #prepend multi-threshold freeboard data
            par_df = pd.concat([frb_df, df],axis=1)
        #close file connection
        nc.close()
        #return to caller
        return df
        








    def _export(self) -> None:
        pass
    