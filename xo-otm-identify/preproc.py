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

from typing import List, Dict
from loguru import logger



# In[]

class PreProcessor:
    """
    Class to handle the preprocessing of the given sensor/carrier combo, i.e.,
    the compilation of lat/lon and file indices from all available L1p files
    """
    def __init__(self, path_to_carrier1: str, path_to_carrier2: str):
        self.path_carrier1 = path_to_carrier1
        self.path_carrier2 = path_to_carrier2
        self.files_carrier1 = os.listdir(self.path_carrier1)
        self.files_carrier2 = os.listdir(self.path_carrier2)
        
    def run(self) -> None:
        #preprocess
        self.carrier1 = self._preprocess(self.path_carrier1, 
                                         self.files_carrier1)
        self.carrier2 = self._preprocess(self.path_carrier2, 
                                         self.files_carrier2)
        #identify xo's/otm's
        
        
        #export to shapefiles
        
            
    def _preprocess(self, 
             paths: List[str],
             files: List[str]) -> Dict[str, gpd.GeoDataFrame]:
        #number of files
        N_FILES = len(files)
    
        #compile the pd.df from the raw files
        xo_preproc = []
        for ncidx, ncfile in enumerate(files):
            #status
            logger.info(f'Processing: {ncfile} - {ncidx+1} out of {N_FILES}')
            #open file
            nc = xr.open_dataset(os.path.join(path_to_files,ncfile))
            #get data
            lat = nc['time_orbit/latitude'].values
            lon = nc['time_orbit/longitude'].values
            tmp = nc['time_orbit/timestamp'].values
            sft = nc['surface_type/flag'].values
            #dist2coast only available for ENVISAT(?)
            try:
                d2c = nc['classifier/dist_coast'].values
            except:
                d2c_available = False
            else:
                d2c_available = True
            #close file link
            nc.close()
            
            #limit data to predefined AOI's
            if self.aoi=='arc':
                west1 = (lat>=65.0) & (lon<=(-85.0)) & (lon>=(-180.0))
                west2 = (lat>=65.0) & (lon<=(275.0)) & (lon>=(180.0))
                east = (lat>=70.0) & (lon>=70.0) & (lon<=180.0)
                aoi_idx = np.where(west1 | west2 | east)
            elif self.aoi=='ant':
                aoi_idx = np.where(lat <= -55.0)
            else:
                logger.critical(f'{self.aoi} undefined!')

                
            #limit data to surface type 'ocean'
            idx2pick = np.where(sft[aoi_idx]==1)
            
            #limit to minimum distant to coast according to the user-selected 
            #buffer range
            if d2c_available:
                idx2pick = np.where(d2c[idx2pick] >= self.buffer)
                
            #append to pd.df
            if len(aoi_idx)>0:
                xo_preproc.append(pd.DataFrame({'id': ncfile,
                                                'time': tmp[idx2pick],
                                                'lat':  lat[idx2pick],
                                                'lon':  lon[idx2pick],
                                                'fidx': idx2pick[0]}))
        #status
        logger.info(f'Converting data to pandas df...')             
        #concatenate and convert to a single pd.df  
        xo_preproc = pd.concat(xo_preproc,ignore_index=True)
        #status
        logger.info(f'Complete!')
        #status
        logger.info(f'Converting data into ESAE2-based geopandas gdf...')    
        #create gpd.gdf and set/change crs
        xo_gdf_dict = {'pts':[],'lin':[]}
        lon = xo_preproc.lon
        lat = xo_preproc.lat
        point_geometry = gpd.points_from_xy(lon,lat)
        xo_gdf_dict['pts'] = gpd.GeoDataFrame(xo_preproc,
                                              geometry=point_geometry)
        xo_gdf_dict['pts'] = xo_gdf_dict['pts'].set_crs(epsg=4326)
        #transform crs
        xo_gdf_dict['pts'] = xo_gdf_dict['pts'].to_crs(epsg=self.epsg)
        #status
        logger.info(f'Complete!')
        
        #status
        logger.info(f'Compile orbit tracks...') 
        #convert it to spatial lines
        pts = xo_gdf_dict['pts']
        pts_grp = pts.groupby(['id'],as_index=False)['geometry']
        lin = pts_grp.apply(lambda x: shapely.geometry.LineString(x.tolist())\
                            if len(x)>=2 else None)   
        lin_id = lin.id
        xo_gdf_dict['lin'] = gpd.GeoDataFrame(lin_id,geometry=lin.geometry,
                                              crs="EPSG:"+str(self.epsg))
        #status
        logger.info(f'Complete!')
        #assign preproc gdf
        return xo_gdf_dict#
    
    
    def _identify(self):
        pass
    
    
    def _export(self):
        pass
    