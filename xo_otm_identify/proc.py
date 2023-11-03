#!/usr/bin/env python3
#coding: utf-8
"""
@author: Dr. Stephan Paul (AWI/iceXai; stephan.paul@awi.de)
"""

# In[]

import os

import pandas as pd
import geopandas as gpd
import xarray as xr
import numpy as np

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
            #check time difference mismatch
            if dt > self.cfg.delta_t:
                continue
            #otherwise continue processing



    def _import_l1p(self, carrier: str, r0: int, r1: int) -> pd.DataFrame:
        #compile carrier specific parameter list
        parameters = self.cfg.par.compile_carrier_data(carrier)
        #create empty DataFrame 
        df = pd.DataFrame()
        #load mandatory file source
        #TODO check for xarray implementation
        # self.src[carrier] = nc_handle.getncattr(pars['src'])
        # del pars['src']


    def load_l1p_nc2pd(self, carrier: str, r0: int, r1: int) -> pd.DataFrame:
        #compile carrier specific parameter list
        pars = self.ncdict.compile_carrier_data(carrier)
        #create empty dict to fill
        data_df = pd.DataFrame()
        #load mandatory file source
        self.src[carrier] = nc_handle.getncattr(pars['src'])
        del pars['src']
        #loop over remaining
        for key in pars['l1']:
            if key == 'pwr' or key == 'rng':
                #load data
                data = nc_handle[pars['l1'][key]][r0:r1,:]
                #convert to DataFrame
                df = pd.DataFrame(data,
                                  columns=[key+'b'+str(rbin+1).zfill(3) 
                                           for rbin in range(data.shape[1])])
                #add it to the existing one
                data_df = pd.concat([data_df,df],axis=1)
            else:
                #load data to dict
                data = {}
                data[key] = nc_handle[pars['l1'][key]][r0:r1]
                #convert it to DataFrame
                df = pd.DataFrame(data)
                #add it to the existing one
                data_df = pd.concat([data_df,df],axis=1)
        #return it
        return data_df



    def _import_l2i(self, carrier: str, r0: int, r1: int) -> pd.DataFrame:
        pass


    ##xo data stack function
    def compile4csv(self) -> None:
        #status
        print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
                  '] > Wrap-up crossover data to csv-file')
        #allocate output
        csv_dict = {'c1':pd.DataFrame(),'c2':pd.DataFrame(),
                      'meta':pd.DataFrame()}
        #prepare for l2i output as well
        if self.return_l2i_status():
            csv_dict['c1_l2i'] = pd.DataFrame()
            csv_dict['c2_l2i'] = pd.DataFrame()
        #number ot total xo/otm
        n_xo = len(self.out)
        #loop over identified crossover entries
        for gdfid,f1,f2,f1r0,f1r1,f2r0,f2r1,dt_xo,xogeo in self.out.itertuples():
            #check delta_t
            if dt_xo<self.delta_t:
                #load l1b file contents
                nc1 = Dataset(os.path.join(self.path_c1,f1))
                nc2 = Dataset(os.path.join(self.path_c2,f2))
                #get l1p data
                c1_data = self.load_l1p_nc2pd(self.c1,f1r0,f1r1,nc1)
                c2_data = self.load_l1p_nc2pd(self.c2,f2r0,f2r1,nc2)
                #close file link
                nc1.close()
                nc2.close()
                #insert xo id
                c1_data.insert(0,self.matchtype+'_idx',gdfid)
                c2_data.insert(0,self.matchtype+'_idx',gdfid)
                #append to carrier stack
                csv_dict['c1'] = pd.concat([csv_dict['c1'],
                                              c1_data],
                                             ignore_index=True)
                csv_dict['c2'] = pd.concat([csv_dict['c2'],
                                              c2_data],
                                             ignore_index=True)
                
                #load l2i files and get l2i data
                if self.return_l2i_status():
                    #retrieve data as pandas DataFrame
                    c1_l2i_data = self.load_l2i_nc2pd(f1,self.c1,f1r0,f1r1,True)
                    c2_l2i_data = self.load_l2i_nc2pd(f2,self.c2,f2r0,f2r1,False)
                    
                    #insert xo id
                    c1_l2i_data.insert(0,self.matchtype+'_idx',gdfid)
                    c2_l2i_data.insert(0,self.matchtype+'_idx',gdfid)
                    #append to carrier stack
                    csv_dict['c1_l2i'] = pd.concat([csv_dict['c1_l2i'],
                                                      c1_l2i_data],
                                                     ignore_index=True)
                    #append to carrier stack
                    csv_dict['c2_l2i'] = pd.concat([csv_dict['c2_l2i'],
                                                      c2_l2i_data],
                                                     ignore_index=True)
                    #add meta data w/ l2i
                    meta_entry = {self.matchtype+'_idx': gdfid,
                                  'l1p_'+self.c1: f1,
                                  'l1p_'+self.c2: f2,
                                  'l1b_'+self.c1: self.src[self.c1],
                                  'l1b_'+self.c2: self.src[self.c2],
                                  'l2i_'+self.c1: self.srcl2i[self.c1],
                                  'l2i_'+self.c2: self.srcl2i[self.c2],
                                  self.matchtype+'_dt': dt_xo}
                    csv_dict['meta'] = pd.concat([csv_dict['meta'],
                                                  pd.DataFrame([meta_entry])],
                                                  ignore_index=True)
                else:
                    #add meta data w/o l2i
                    meta_entry = {self.matchtype+'_idx': gdfid, 
                                  'l1p_'+self.c1: f1,
                                  'l1p_'+self.c2: f2,
                                  'l1b_'+self.c1: self.src[self.c1],
                                  'l1b_'+self.c2: self.src[self.c2],
                                  self.matchtype+'_dt': dt_xo}
                    csv_dict['meta'] = pd.concat([csv_dict['meta'],
                                                  pd.DataFrame([meta_entry])],
                                                  ignore_index=True)
        #save output
        if len(csv_dict['meta'])>0:
            #status
            if self.return_l2i_status():
                print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
                      '] - Processed L1p/L2i waveforms ('+self.c1+'):'+\
                          str(len(csv_dict['c1']))+'/'+\
                          str(len(csv_dict['c1_l2i'])))
                print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
                      '] - Processed L1p/L2i waveforms ('+self.c2+'):'+\
                          str(len(csv_dict['c2']))+'/'+\
                          str(len(csv_dict['c2_l2i'])))
            else:
                print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
                      '] - Processed L1p waveforms ('+self.c1+'):'+\
                          str(len(csv_dict['c1'])))
                print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
                      '] - Processed L1p waveforms ('+self.c2+'):'+\
                          str(len(csv_dict['c2'])))
            print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
                  '] - Processed '+str(len(csv_dict['meta']))+\
                  ' out of '+str(n_xo)+' potential '+self.matchtype.upper()+'s')        
            print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
                  '] - Processing complete - compiling CSV output...')
            self.data4csv = csv_dict
        else:
            #status
            print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
                  '] - Processing complete - nothing to compile...')










    def _export(self) -> None:
        pass
    