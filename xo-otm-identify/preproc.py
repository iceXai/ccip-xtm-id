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

from cfg import Configuration



# In[]

class PreProcessor:
    """
    Class to handle the preprocessing of the given sensor/carrier combo, i.e.,
    the compilation of lat/lon and file indices from all available L1p files
    """
    def __init__(self, cfg: Configuration):
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
        self._export(id_gdf)
        
        #import already existing hsape files
        self._import()
            
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
    
    def _identify_intersection(self, line1, line2) -> None:
        if self.matchtype == 'xo':
            return lin1.intersection(lin2)
        if self.matchtype == 'otm':
            return lin1.intersection(lin2.buffer(self.buffer))
        
    def _identify_points_in_match(self, line_points) -> None:
        #get the geometrie
        line_point_geometries = gpd.GeoSeries(line_points.geometry)
        #buffer intersect
        match_buffered = match.buffer(self.buffer)
        #retrieve indices of points
        line_points_in_match = []
        for pt_idx, pt in enumerate(line_point_geometries):
            point_intersects = pt.intersects(match_buffered)
            if point_intersects:
                line_points_in_match.append(pt_idx)
        #return to caller
        return line_points_in_match
    
    def _identify_XOs(self, 
                      match,
                      line1_points, 
                      line2_points) -> None:
        #identify line points in match intersect
        line2_points_in_match = self._identify_points_in_match(line1_points)
        line2_points_in_match = self._identify_points_in_match(line2_points)
        #get the geometrie
        line1_point_geometries = gpd.GeoSeries(line1_points.geometry)
        line2_point_geometries = gpd.GeoSeries(line2_points.geometry)
                
        #check whether points exist in buffer
        if len(line1_points_in_match)>0 and len(line2_points_in_match)>0:
            #get correct indices in files
            fidx1 = line1_points['fidx'].iloc[line1_points_in_match]
            fidx2 = line2_points['fidx'].iloc[line2_points_in_match]
            
            #get time of all points in the match intersect
            dt_idx1 = int(len(line1_points_in_match)/2)
            dt_idx2 = int(len(line2_points_in_match)/2)
            t1 = line1_points['time'].iloc[line1_points_in_match[dt_idx1]]
            t2 = line2_points['time'].iloc[line2_points_in_match[dt_idx2]]
            
            #convert to datetime obj
            dt_pts1 = dt.datetime.fromtimestamp(t1,dt.timezone.utc)
            dt_pts2 = dt.datetime.fromtimestamp(t2,dt.timezone.utc)
            
            #calculate average time difference [in hours]
            delta_dt_match = dt_pts1-dt_pts2
            abs_delta_hh = abs(delta_dt_match.total_seconds()/(3600))
            #return attributes of the match with file index ranges and average
            #time difference
            ATTRS = (line1_points['id'].iloc[0],
                     line2_points['id'].iloc[0],
                     fidx1.iloc[0],
                     fidx1.iloc[-1],
                     fidx2.iloc[0],
                     fidx2.iloc[-1],
                     np.round(abs_delta_hh,3))
            return ATTRS
        else: 
            return None
    
    def _identify_OTMs(self, 
                       match,
                       line1_points, 
                       line2_points) -> None:
        #identify line points in match intersect
        line2_points_in_match = self._identify_points_in_match(line1_points)
        line2_points_in_match = self._identify_points_in_match(line2_points)
        #get the geometrie
        line1_point_geometries = gpd.GeoSeries(line1_points.geometry)
        line2_point_geometries = gpd.GeoSeries(line2_points.geometry)
    
        #check whether enough points exist in buffer
        if len(line1_points_in_match)>50 and len(line2_points_in_match)>50:
            #get correct indices in files
            fidx1 = line1_points['fidx'].iloc[line1_points_in_match]
            fidx2 = line2_points['fidx'].iloc[line2_points_in_match]
            
            #get time of all points in the match intersect
            dt_idx1 = int(len(line1_points_in_match)/2)
            dt_idx2 = int(len(line2_points_in_match)/2)
            t1 = line1_points['time'].iloc[line1_points_in_match[dt_idx1]]
            t2 = line2_points['time'].iloc[line2_points_in_match[dt_idx2]]
            
            #convert to datetime obj
            dt_pts1 = dt.datetime.fromtimestamp(t1,dt.timezone.utc)
            dt_pts2 = dt.datetime.fromtimestamp(t2,dt.timezone.utc)
            
            #calculate average time difference [in hours]
            delta_dt_match = dt_pts1-dt_pts2
            abs_delta_hh = abs(delta_dt_match.total_seconds()/(3600))
            #return attributes of the match with file index ranges and average
            #time difference
            ATTRS = (line1_points['id'].iloc[0],
                     line2_points['id'].iloc[0],
                     fidx1.iloc[0],
                     fidx1.iloc[-1],
                     fidx2.iloc[0],
                     fidx2.iloc[-1],
                     np.round(abs_delta_hh,3))
            return ATTRS
        else: 
            return None
    
    def _identify(self, 
                  gdfdict_c1: Dict[str, gpd.GeoDataFrame], 
                  gdfdict_c2: Dict[str, gpd.GeoDataFrame]) -> gpd.GeoDataFrame:
        #allocate sum-up
        attributes = []
        geometries = []
        
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
                #get matching points
                points_on_line1 = POINTS1[POINTS1['id']==xo_id1]
                points_on_line2 = POINTS2[POINTS2['id']==xo_id2]
                #identify XOs and only use point-like intersects 
                #(no close "fly-by's")
                if type(MATCH) is shapely.geometry.point.Point:
                    ATTRS = self._identify_XOs(MATCH, 
                                                points_on_line1, 
                                                points_on_line2)                           
                #identify OTMs
                elif type(MATCH) is shapely.geometry.linestring.LineString:
                    ATTRS = self._identify_OTMs(MATCH, 
                                                points_on_line1, 
                                                points_on_line2)
                else:
                    logger.critical(f'Unsupported intersect type: '+\
                                    f'{type(MATCH)}!')
                    sys.exit()
                if ATTRS is not None:
                    attributes.append(ATTRS)
                    geometries.append(MATCH)
                    
        #store data into a GeoDataFrame
        logger.info(f'Saving results to GeoDataFrame...')
        COLUMNS = [self.matchtype+'_f1',
                   self.matchtype+'_f2',
                   'f1_rng_0',
                   'f1_rng_1',
                   'f2_rng_0',
                   'f2_rng_1',
                   'dt_'+self.matchtype,
                   ]
        return gpd.GeoDataFrame(attributes,
                                geometry = geometries,
                                crs = "EPSG:"+str(self.epsg),
                                columns = COLUMNS)
        logger.info(f'Process complete!')    
    
    def _export(self, id_gdf: gpd.GeoDataFrame) -> None:
        
        
        
        
        #set-up output base name
        if self.matchtype == 'otm':
            self.out2shp = os.path.join(self.outpath,
                                    '_'.join([self.matchtype,self.c1,self.c2,
                                              self.year,self.month,
                                              self.aoi]))
            self.out2csv = os.path.join(self.outpath,
                                        '_'.join([self.matchtype,self.year,self.month,
                                                  self.aoi]))
        else:
            buff_tag = str(self.buffer).zfill(5)+'m'
            self.out2shp = os.path.join(self.outpath,
                                        '_'.join([self.matchtype,self.c1,self.c2,
                                                  self.year,self.month,
                                                  self.aoi,buff_tag]))
            self.out2csv = os.path.join(self.outpath,
                                        '_'.join([self.matchtype,self.year,self.month,
                                                  self.aoi,buff_tag]))
            
        #save identified crossovers to shp file
        def data2shp(self) -> None:
            #make sure crs is set
            self.out = self.out.set_crs(epsg=self.epsg)
            #save
            self.out.to_file(''.join([self.out2shp,'.shp']))
            
        #load identified crossovers from shp file
        def shp2data(self) -> None:
            if os.path.isfile(self.out2shp+'.shp'):
                self.out = gpd.read_file(self.out2shp+'.shp')
            else:
                print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
                      '] ! No corresponding .shp file found')
                sys.exit()
                
                
    def _import(self):
        
            
            
    