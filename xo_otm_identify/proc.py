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
        self.path_carrier1 = cfg.path_to_carrier1_l2i(year, month)
        self.path_carrier2 = cfg.path_to_carrier2_l2i(year, month)
        self.files_carrier1 = os.listdir(self.path_carrier1)
        self.files_carrier2 = os.listdir(self.path_carrier2)
    
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
            c1_df = self._import_l1p(path_to_f1, f1r0, f1r1)
            path_to_f2 = os.path.join(self.cfg._input_l1p, f2)
            c2_df = self._import_l1p(path_to_f2, f2r0, f2r1)
            #insert xo/otm id
            c1_df.insert(0, self.cfg.matchtype+'_idx', idx)
            c2_df.insert(0, self.cfg.matchtype+'_idx', idx)
            #append to output
            TYPE = self.cfg.carrier1
            csvdict[TYPE] = pd.concat([csvdict[TYPE],c1_df],ignore_index=True)
            TYPE = self.cfg.carrier2
            csvdict[TYPE] = pd.concat([csvdict[TYPE],c2_df],ignore_index=True)



    def _import_l1p(self, path: str, r0: int, r1: int) -> pd.DataFrame:
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
                #TODO has yet to be returned!!!
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
        return df

    def _import_l2i_ref(self, path: str, r0: int, r1: int) -> pd.DataFrame:
        #compile product-level specific parameter list
        parameters = self.cfg.par.l2i_parameters()
        #get l1p file name
        l1p = pathlib.Path(path).name
        #retrieve time/date tag from file name
        tag = l1p.split('-')[6:8]
        tag = [dt.datetime.strptime(tag[0],'%Y%m%dT%H%M%S'),
               dt.datetime.strptime(tag[1],'%Y%m%dT%H%M%S')]
        #get l2i files
        l2i = self.files_carrier1
        l2i_tag = [f.split('-')[6:8] for f in l2i]




    def _import_l2i_match(self, path: str, r0: int, r1: int) -> pd.DataFrame:
        #compile product-level specific parameter list
        parameters = self.cfg.par.l2i_parameters()



    def _identify_l2i_by_l1p(self, path: str):
        pass




    def load_l2i_nc2pd(self,sensorfile: str, carrier: str,
                       r0: int, r1: int, is_ref: bool) -> pd.DataFrame:
        #compile carrier specific parameter list
        pars = self.ncdict.compile_carrier_data(carrier)
            
        #retrieve time/date tag from file name
        tag = sensorfile.split('-')[6:8]
        tag = [dt.datetime.strptime(tag[0],'%Y%m%dT%H%M%S'),
               dt.datetime.strptime(tag[1],'%Y%m%dT%H%M%S')]
        #create search
        #search_list = self.l2imeta['files'][carrier]
        
        #id correct l2i match
        tmp = [l2i.split('-')[6:8] for l2i in self.l2imeta['files'][carrier]]
        for i,t in enumerate(tmp):
            #convert to datetime
            t_dt = [dt.datetime.strptime(t[0],'%Y%m%dT%H%M%S'),
                    dt.datetime.strptime(t[1][:-5],'%Y%m%dT%H%M%S')]
            #calculate time difference in total seconds
            dt_start = t_dt[0]-tag[0]
            dt_end = t_dt[1]-tag[1]
            #only keep in case dt is exact match
            if dt_start.total_seconds() == 0 and \
                dt_end.total_seconds() == 0:
                #file name
                l2i = self.l2imeta['files'][carrier][i]
                #path
                pathl2i = self.l2imeta['paths'][carrier]
                #compile full file path
                fp = os.path.join(pathl2i,l2i)
                #store file name for later use in meta data
                self.srcl2i[carrier] = l2i
                #break loop in case one is found
                break
        
        #open file connection
        nc_handle = Dataset(fp)
        #create empty dict to fill
        data_df = pd.DataFrame()
        #load the data
        if is_ref:
            data = {}
            data['frb'] = nc_handle['sea_ice_freeboard'][r0:r1]
            data['sft'] = nc_handle['surface_type'][r0:r1]
            if self.ncdict.has_l2_pars:
                for key in pars['l2']:
                    data[key] = nc_handle[pars['l2'][key]][r0:r1]
            nc_handle.close()
            #convert it to DataFrame
            df = pd.DataFrame(data)
            #add it to the existing one
            data_df = pd.concat([data_df,df],axis=1)
        else:
            data = {}
            frb = nc_handle['threshold_freeboards'][r0:r1,:]
            ths = nc_handle['tfmra_thresholds'][:]
            data['sft'] = nc_handle['surface_type'][r0:r1]
            if self.ncdict.has_l2_pars:
                for key in pars['l2']:
                    data[key] = nc_handle[pars['l2'][key]][r0:r1]
            nc_handle.close()
            #convert it to DataFrame w/ correct column names
            df = pd.DataFrame(frb,
                              columns=['frb_th'+\
                                       str(format(np.round(th,2),'.2f'))[2:]
                                       for th in ths])
            #append surface-type data
            df = pd.concat([df,pd.DataFrame(data)],axis=1)
            #add it to the existing one
            data_df = pd.concat([data_df,df],axis=1)
        #return it
        return data_df




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
    