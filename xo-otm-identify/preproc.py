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
        #compile geodata
        carrier1_gdfdict = self._compile(self.path_carrier1, 
                                         self.files_carrier1)
        carrier2_gdfdict = self._compile(self.path_carrier2, 
                                         self.files_carrier2)
        #identify xo's/otm's
        id_gdf = self.identify(carrier1_gdfdict, carrier2_gdfdict)
        
        #export to shapefiles
        
            
    def _compile(self, 
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
    
    def _identify_intersection(self, line1, line2):
        if self.matchtype == 'xo':
            return lin1.intersection(lin2)
        if self.matchtype == 'otm':
            return lin1.intersection(lin2.buffer(self.buffer))
        
    def _identify_XOs(self, match, line1, line2):
        #point subset by file id
        pts1_sub = pts1[pts1['id']==xo_id1]
        pts2_sub = pts2[pts2['id']==xo_id2]
        pts1_sub_geo = gpd.GeoSeries(pts1_sub.geometry)
        pts2_sub_geo = gpd.GeoSeries(pts2_sub.geometry)
        
        #buffer intersect
        xo_buff = xo.buffer(self.buffer)
        
        #retrieve indices of points
        pts1_xo_idx = []
        for pt1_idx, pt1 in enumerate(pts1_sub_geo):
            point_intersects = pt1.intersects(xo_buff)
            if point_intersects:
                pts1_xo_idx.append(pt1_idx)
        pts2_xo_idx = []
        for pt2_idx, pt2 in enumerate(pts2_sub_geo):
            point_intersects = pt2.intersects(xo_buff)
            if point_intersects:
                pts2_xo_idx.append(pt2_idx)
                
        #check whether points exist in buffer
        if len(pts1_xo_idx)>0 and len(pts2_xo_idx)>0:
            #check time difference
            dt_idx1 = int(len(pts1_xo_idx)/2)
            dt_idx2 = int(len(pts2_xo_idx)/2)
            t1 = pts1_sub['time'].iloc[pts1_xo_idx[dt_idx1]]
            t2 = pts2_sub['time'].iloc[pts2_xo_idx[dt_idx2]]
            
            #get correct indices in files
            fidx1 = pts1_sub['fidx'].iloc[pts1_xo_idx]
            fidx2 = pts2_sub['fidx'].iloc[pts2_xo_idx]
            
            #convert to datetime obj
            dt_pts1 = dt.datetime.fromtimestamp(t1,dt.timezone.utc)
            dt_pts2 = dt.datetime.fromtimestamp(t2,dt.timezone.utc)
            
            #calculate average time difference [in hours]
            delta_dt_xo = dt_pts1-dt_pts2
            abs_delta_hh = abs(delta_dt_xo.total_seconds()/(3600))
            #store attributes/geometry to attributes and 
            #geometry list
            #[idx_files1_list, idx_files2_list,
            # rng1_start, rng1_end, 
            # rng2_start, rng2_end, dt]
            attr.append((xo_id1,xo_id2,
                         fidx1.iloc[0],
                         fidx1.iloc[-1],
                         fidx2.iloc[0],
                         fidx2.iloc[-1],
                         np.round(abs_delta_hh,3)))
            geom.append(xo)
    
    def identify_OTMs(self, match, line1, line2):
        #point subset by file id
        pts1_sub = pts1[pts1['id']==xo_id1]
        pts2_sub = pts2[pts2['id']==xo_id2]
        pts1_sub_geo = gpd.GeoSeries(pts1_sub.geometry)
        pts2_sub_geo = gpd.GeoSeries(pts2_sub.geometry)
        
        #buffer intersect
        xo_buff = xo.buffer(5000)
        
        #retrieve indices of points
        pts1_xo_idx = []
        for pt1_idx, pt1 in enumerate(pts1_sub_geo):
            point_intersects = xo_buff.contains(pt1)
            if point_intersects:
                pts1_xo_idx.append(pt1_idx)
        pts2_xo_idx = []
        for pt2_idx, pt2 in enumerate(pts2_sub_geo):
            point_intersects = xo_buff.contains(pt2)
            if point_intersects:
                pts2_xo_idx.append(pt2_idx)
    
        #check whether points exist in buffer
        if len(pts1_xo_idx)>50 and len(pts2_xo_idx)>50:
            #get correct indices in files
            fidx1 = pts1_sub['fidx'].iloc[pts1_xo_idx]
            fidx2 = pts2_sub['fidx'].iloc[pts2_xo_idx]
            
            #check and calculate time difference [in hours]
            t1 = pts1_sub['time'].iloc[pts1_xo_idx[0]]
            t2 = pts2_sub['time'].iloc[pts2_xo_idx[0]]
            dt_pts1 = dt.datetime.fromtimestamp(t1,dt.timezone.utc)
            dt_pts2 = dt.datetime.fromtimestamp(t2,dt.timezone.utc)
            delta_dt_xo = dt_pts1-dt_pts2
            abs_delta_hh = abs(delta_dt_xo.total_seconds()/(3600))
        
            #store attributes/geometry to attributes and 
            #geometry list
            #[idx_files1_list, idx_files2_list,
            # rng1_start, rng1_end, 
            # rng2_start, rng2_end, dt]
            attr.append((xo_id1,xo_id2,
                         fidx1.iloc[0],
                         fidx1.iloc[-1],
                         fidx2.iloc[0],
                         fidx2.iloc[-1],
                         np.round(abs_delta_hh,3)))
            geom.append(xo)
    
    def _identify(self, 
                  gdfdict_c1: Dict[str, gpd.GeoDataFrame], 
                  gdfdict_c2: Dict[str, gpd.GeoDataFrame]) -> gpd.GeoDataFrame:
        #allocate sum-up
        attr = []
        geom = []
        
        #load point data from sensor 1/2
        POINTS1 = gdfdict_c1['pts']
        POINTS2 = gdfdict_c2['pts']
        
        #load line data from sensor 1/2
        LINES1 = gdfdict_c1['lin']
        LINES2 = gdfdict_c2['lin']
        
        #convert to datetime to speed up the id process
        avg_t1 = POINTS1[['id','time']].groupby(['id'], as_index=False).mean()
        dt1_df = pd.to_datetime(avg_t1.time,unit='s')
        avg_t2 = POINTS2[['id','time']].groupby(['id'], as_index=False).mean()
        dt2_df = pd.to_datetime(avg_t2.time,unit='s')
        
        #loop over all sensor 1 orbits to identify crossovers and 
        #orbit-trajectory matches with sensor 2 orbits
        N_LINES_CARRIER1 = len(lns1)
        for idx1, xo_id1, lin1 in LINES1.itertuples():
            #status
            logger.info(f'Processing {idx1+1} out of {N_LINES_CARRIER1}')
            logger.info(f'Identifying orbit intersections...')

            #create time-dependent subset of sensor2 orbits to save time
            #by calculating absolute date time difference [in hours]
            abs_dtd_hh = (dt2_df - dt1_df.iloc[idx1]).abs()
            abs_dtd_hh = abs_dtd_hh / pd.Timedelta(hours=1)
            #subset indices to use by specified maximum time difference
            LINES_IDX = abs_dtd_hh[abs_dtd_hh <= self.delta_t].index

            #loop over the sensor2 orbit subset and collect information
            for idx2, xo_id2, lin2 in LINES2.iloc[LINES_IDX,:].itertuples():
                #summarize data
                if lin1 is None or lin2 is None:
                    continue

                #find intersection
                MATCH = self._identify_intersection(lin1, lin2)
                #identify XOs and only use point-like intersects 
                #(no close "fly-by's")
                if type(MATCH) is shapely.geometry.point.Point:
                    self.identify_XOs()                        
                #identify OTMs
                if type(MATCH) is shapely.geometry.linestring.LineString:
                    self.identify_OTMs()
                    
        #store data into a GeoDataFrame
        print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
              '] - Saving results...')
        self.out = gpd.GeoDataFrame(attr,
                                   geometry=geom,
                                   crs="EPSG:"+str(self.epsg),
                                   columns=[self.matchtype+'_f1',
                                            self.matchtype+'_f2',
                                            'f1_rng_0',
                                            'f1_rng_1',
                                            'f2_rng_0',
                                            'f2_rng_1',
                                            'dt_'+self.matchtype])
        #return output
        print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
              '] < Process complete. :)') 
    
    
    def _export(self):
        pass
    