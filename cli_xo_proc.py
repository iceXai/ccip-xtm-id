#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# In[]
##libs
import argparse
import numpy as np
#import xarray as xr
from netCDF4 import Dataset
import shapely
import os
import sys
import geopandas as gpd
import pandas as pd
import datetime as dt
from typing import List, Dict


# In[]
# In[]
def xo_wrapper():
    ##commandline stuff
    parser = argparse.ArgumentParser(description='CommandLine Interface for '+\
                                     'Dual Mission Crossovers (XO)')
    #required flags/arguments
    parser.add_argument('--carrier',required=True,metavar=('CARRIER'),nargs=2,
                        default=None,
                        help='Specify carrier satellites [Two out of '+\
                            '\'ERS1/2\'/\'ENVISAT\' or \'CRYOSAT2\'] '+\
                            'NOTE: Reference Carrier should go FIRST!')
    parser.add_argument('--l1p-input',required=True,metavar='L1P-INPUT',
                        nargs=1,dest='l1p_input',
                        help='Path to L1p inut data structure as provided/used '+\
                            'by pySIRAL [up to carrier level]')
    parser.add_argument('--output',required=True,metavar='OUTPUT DIR',nargs=1,
                        help='Specify output directory')
    parser.add_argument('--aoi',required=True,metavar='AOI',nargs=1,
                        help='Specify an area of interest to be used for the '+\
                            'search of crossovers '+\
                            '[either \'ARC\' or \'ANT\']')  
    parser.add_argument('--year',required=True,metavar=('YEAR'),nargs=1,
                        help='Specify year')
    parser.add_argument('--month',required=True,metavar=('MONTH'),nargs=1,
                        help='Specify month')
    parser.add_argument('--match-type',required=True,metavar='matchtype',nargs=1,
                        dest='matchtype',
                        help='Specify wether to use actual point-like '+\
                            'crossovers or orbit track matches '+\
                            '[either \'XO\' or \'OTM\']')     
    
    #additional non-required flags
    parser.add_argument('--buffer-size',required=False,metavar='buffer',
                        nargs=1,dest='buffer',default=12500,type=int,
                        help='Specify a buffer size [in meters] around each '+\
                            'crossover to collect waveforms')
    parser.add_argument('--max-delta-t',required=False,metavar='DT',nargs=1,
                        dest='delta',default=[1],type=int,
                        help='Specify a maximum time difference [in hours] '+\
                            'allowed for a crossover '+\
                            'to be considered')
    parser.add_argument('--csv',action='store_true',required=False,
                        dest='csv',default=False,
                        help='Specify to create csv output of XO/OTM '+\
                            'waveform parameters for t <= --max-delta-t')
    parser.add_argument('--include-l2i',action='store_true',required=False,
                        dest='l2i',default=False,
                        help='Specify to also include parameters/information '+\
                            'from corresponding L2i data')
    parser.add_argument('--l2i-input',required=False,metavar='L2I-INPUT',
                        nargs=1,dest='l2i_input',default=None,
                        help='Path to L2i inut data structure as provided/used '+\
                            'by pySIRAL [up to carrier level]')
    parser.add_argument('--params',required=False,metavar='PARAMS',nargs='+',
                        default=None,dest='par',
                        help='Specify parameters to be output to CSV'+\
                            ' (requires --csv flag). These include:'+\
                            ' \'pwr\', \'rng\', \'rdm\', \'alt\', \'flg\','+\
                            ' \'fmi\', \'ppk\', \'lew\', \'lep\', \'leq\','+\
                            ' \'nsp\', \'sig\', \'eps\', \'sku\', \'spk\','+\
                            ' \'ssd\', \'ssk\' among others.\n'+\
                            'L2i output is fixed and unaffected by choices '+\
                            'made with --params.') 
    parser.add_argument('--c1-l1p-vers',required=False,metavar='c1_l1p_v',
                        nargs=1,dest='c1_l1p_v',default='v1p1',
                        help='Specify a Carrier1 L1p data version'+\
                            ' (first entry with --carrier); e.g., \'v1p3\'')
    parser.add_argument('--c2-l1p-vers',required=False,metavar='c2_l1p_v',
                        nargs=1,dest='c2_l1p_v',default='v1p1',
                        help='Specify a Carrier2 L1p data version'+\
                            ' (second entry with --carrier); e.g., \'v1p3\'')
    parser.add_argument('--c1-l2i-vers',required=False,metavar='c1_l2i_v',
                        nargs=1,dest='c1_l2i_v',default='v3p0-rc2',
                        help='Specify a Carrier1 L2i data version'+\
                            ' (first entry with --carrier); e.g., \'v3p0-rc2\'')
    parser.add_argument('--c2-l2i-vers',required=False,metavar='c2_l2i_v',
                        nargs=1,dest='c2_l2i_v',default='v1p1',
                        help='Specify a Carrier2 L2i data version with'+\
                            ' multi-thresholds'+\
                            ' (second entry with --carrier); e.g., '+\
                                '\'v3p0-preview5-mt-v1p3\'') 
    
    
    #additional flags to skip preceeding processing steps
    parser.add_argument('--skip-preproc',action='store_true',required=False,
                        dest='skip',default=False,
                        help='Specify in case the preprocessing, i.e., '+\
                            'the XO/OTM identification was already done')   
    

    #splash
    print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
          '] > CommandLine Interface for Dual Mission Crossovers (XO)')    
    #parse arguments
    args = parser.parse_args()
    
    #initialize l1/l2 dict container
    ncdict = par_dict_container()
        
    #validate args
    xo_args_val(args,ncdict)
    
    #init job
    xo_job = xo_processor(args,ncdict)
    
    #check process status
    do_full_process = xo_job.return_process_status()
    
    #prepoc
    if do_full_process:
        #preproc
        xo_job.preproc()
        #id
        xo_job.identify()
        #save data set for later re-use
        xo_job.data2shp()
    else:
        #load preprocessed file
        xo_job.shp2data()
        
    #check csv output status
    do_output2csv = xo_job.return_csv_status()

    #collect waveform parameters per file from xo data and store as csv
    if do_output2csv:
        xo_job.compile4csv()
        xo_job.data2csv()
    
    ##end :)
    print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
              '] < Done. :)') 
    


# In[]
##xo process handler class
class xo_processor(object):
    #initial setup and parsing
    def __init__(self,args,ncdict):
        ##allocate placeholders
        self.out       = None
        self.data4csv  = None
        self.src       = {}
        self.srcl2i    = {}
        
        ##convert and store user input for later use
        
        #parse the nc file dictionary
        self.ncdict = ncdict
        
        #convert date and zero-pad months
        self.year  = str(args.year[0])
        self.month = str(args.month[0].zfill(2))
        
        #convert dt
        self.delta_t = float(args.delta[0])
            
        #parse csv
        self.csv = args.csv
        
        #parse l2i
        self.l2i = args.l2i
        
        #parse buffer
        self.buffer = args.buffer
        
        #parse skip/aoi and pathing info
        self.skip        = args.skip
        self.aoi         = args.aoi[0].lower()
        self.l1p_inpath  = args.l1p_input[0]
        self.l2i_inpath  = args.l2i_input[0]
        self.outpath     = args.output[0]
        
        #parse carrier info
        self.c1 = args.carrier[0].lower()
        self.c2 = args.carrier[1].lower()
        
        #parse xo type
        self.matchtype = args.matchtype[0].lower()
        
        #parse L1p versioning
        self.l1pvers = {self.c1: args.c1_l1p_v[0],
                        self.c2: args.c2_l1p_v[0]}
        
        #parse L1p versioning
        self.l2ivers = {self.c1: args.c1_l2i_v[0],
                        self.c2: args.c2_l2i_v[0]}
        
        #epsg
        if self.aoi=='arc':
            self.epsg = 6931
            self.hemis_long = 'north'
            self.hemis_short = 'nh'
        else:
            self.epsg = 6932
            self.hemis_long = 'south'
            self.hemis.short = 'sh'

        ##pathing
        #l1p base paths
        if self.c1 == 'cryosat2':
            c1_l1p_intm = 'ipf1-d/rep/l1p'
        else:
            c1_l1p_intm = 'l1p'
        if self.c2 == 'cryosat2':
            c2_l1p_intm = 'ipf1-d/rep/l1p'
        else:
            c2_l1p_intm = 'l1p'
        #compile paths
        self.path_c1 = os.path.join(self.l1p_inpath,
                                    self.c1,
                                    c1_l1p_intm,
                                    self.l1pvers[self.c1],
                                    self.hemis_long,
                                    self.year,
                                    self.month)
        self.path_c2 = os.path.join(self.l1p_inpath,
                                    self.c2,
                                    c2_l1p_intm,
                                    self.l1pvers[self.c2],
                                    self.hemis_long,
                                    self.year,
                                    self.month)
        #compile files
        self.files_c1 = os.listdir(self.path_c1)
        self.files_c2 = os.listdir(self.path_c2)
        
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
        
        #l2i pathing and files
        if self.l2i:
            #l2i base paths
            self.l2i_base_paths = {self.c1: os.path.join(self.l2i_inpath,
                                                         self.c1,
                                                         self.l2ivers[self.c1],
                                                         self.hemis_short,
                                                         'l2i'),
                                   self.c2: os.path.join(self.l2i_inpath,
                                                         self.c2,
                                                         self.l2ivers[self.c2],
                                                         self.hemis_short,
                                                         'l2i')}
            
            #compile current paths
            paths = {}
            paths[self.c1] = os.path.join(self.l2i_base_paths[self.c1],
                                          self.year,
                                          self.month)
            paths[self.c2] = os.path.join(self.l2i_base_paths[self.c2],
                                          self.year,
                                          self.month)
            
            #compile files
            files = {}
            files[self.c1] = os.listdir(paths[self.c1])
            files[self.c2] = os.listdir(paths[self.c2])
            
            #sum-up
            self.l2imeta = {'paths': paths,
                            'files': files}
    
    #return status of skipping or not to the wrapper
    def return_process_status(self) -> bool:
        #return value
        return not self.skip

    #return status whether csv shall be output or not to the wrapper
    def return_csv_status(self) -> bool:
        #return value
        return self.csv
        
    #return status whether l2i shall be output or not to the wrapper
    def return_l2i_status(self) -> bool:
        #return value
        return self.l2i
    
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
            
    #save identified crossovers to csv file
    def data2csv(self) -> None:
        if self.data4csv is not None:
            #data shortcuts
            meta = self.data4csv['meta'][[self.matchtype+'_idx',self.matchtype+'_dt']]
            c1_xo_id = self.data4csv['c1'][self.matchtype+'_idx']
            c2_xo_id = self.data4csv['c2'][self.matchtype+'_idx']
            if self.return_l2i_status():
                c1l2i_xo_id = self.data4csv['c1_l2i'][self.matchtype+'_idx']
                c2l2i_xo_id = self.data4csv['c2_l2i'][self.matchtype+'_idx']
            #loop over time differences in one-hourly steps
            for i, dt in enumerate(list(range(0,int(self.delta_t)))):
                #find matching XO's in meta data
                #in the first round include both sides [0:1]
                if i == 0:
                    valid_xo_idx = meta[self.matchtype+'_idx'][meta[self.matchtype+'_dt'].between(dt,(dt+1))].tolist()
                
                #in all others exclude left side ]dt:dt+1], i.e. ]1,2]!
                else:
                    valid_xo_idx = meta[self.matchtype+'_idx'][meta[self.matchtype+'_dt'].between(dt,(dt+1),
                                                                  inclusive='right')].tolist()
                
                #pick correct entries/indices for each carrier
                c1_idx_in_csv = []
                for current_idx in valid_xo_idx:
                    c1_idx_in_csv.extend(c1_xo_id.index[c1_xo_id == current_idx])
                c2_idx_in_csv = []
                for current_idx in valid_xo_idx:
                    c2_idx_in_csv.extend(c2_xo_id.index[c2_xo_id == current_idx])
                meta_idx_in_csv = []
                for current_idx in valid_xo_idx: 
                    meta_idx_in_csv.extend(meta.index[meta[self.matchtype+'_idx'] == current_idx])
                
                #with l2i to output as well
                if self.return_l2i_status():
                    c1l2i_idx_in_csv = []
                    for current_idx in valid_xo_idx:
                        c1l2i_idx_in_csv.extend(c1l2i_xo_id.index[c1l2i_xo_id == current_idx])
                    c2l2i_idx_in_csv = []
                    for current_idx in valid_xo_idx:
                        c2l2i_idx_in_csv.extend(c2l2i_xo_id.index[c2l2i_xo_id == current_idx])
                
                #output to csv
                self.data4csv['c1'].iloc[c1_idx_in_csv,:].to_csv(''.join([self.out2csv,'_',
                                                                       self.c1,
                                                                       '_l1p_dt',
                                                                       str(dt+1),
                                                                       '.csv']),
                                                              index=False)
                self.data4csv['c2'].iloc[c2_idx_in_csv,:].to_csv(''.join([self.out2csv,'_',
                                                                       self.c2,
                                                                       '_l1p_dt',
                                                                       str(dt+1),
                                                                       '.csv']),
                                                              index=False)
                self.data4csv['meta'].iloc[meta_idx_in_csv,:].to_csv(''.join([self.out2csv,
                                                                           '_meta',
                                                                           '_dt',
                                                                           str(dt+1),
                                                                           '.csv']),
                                                                  index=False)
                
                #output l2i if specified
                if self.return_l2i_status():
                    self.data4csv['c1_l2i'].iloc[c1l2i_idx_in_csv,:].to_csv(''.join([self.out2csv,
                                                                                  '_',
                                                                                  self.c1,
                                                                                  '_l2i_dt',
                                                                                  str(dt+1),
                                                                                  '.csv']),
                                                                         index=False)
                    self.data4csv['c2_l2i'].iloc[c2l2i_idx_in_csv,:].to_csv(''.join([self.out2csv,
                                                                                  '_',
                                                                                  self.c2,
                                                                                  '_l2i_dt',
                                                                                  str(dt+1),
                                                                                  '.csv']),
                                                                         index=False)

    
    
    # In[]
    ##xo preproc function
    def preproc(self) -> None:
        #preproc carrier 1
        self.preproc_c1 = self.preproc_carrier(self.path_c1,self.files_c1)
        #preproc carrier 2
        self.preproc_c2 = self.preproc_carrier(self.path_c2,self.files_c2)
        
    #preproc sub function doing the actual preproc'ing    
    def preproc_carrier(self,path_to_files: str,list_of_files: str)\
        -> Dict[str, gpd.GeoDataFrame]:
    
        #compile the pd.df from the raw files
        xo_preproc = []
        for ncidx, ncfile in enumerate(list_of_files):
            print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
                  '] > Processing file: '+ncfile+' - '+str(ncidx+1)+\
                      ' out of '+str(len(list_of_files)))
                
            #open file
            nc = Dataset(os.path.join(path_to_files,ncfile))
            
            #get data
            lat = nc['time_orbit/latitude'][:]
            lon = nc['time_orbit/longitude'][:]
            tmp = nc['time_orbit/timestamp'][:]
            sft = nc['surface_type/flag'][:]
            
            #dist2coast only available for ENVISAT(?)
            try:
                d2c = nc['classifier/dist_coast'][:]
            except:
                d2c_available = False
            else:
                d2c_available = True
                
            #close file link
            nc.close()
            
            #limit data to predefined AOI's
            if self.aoi=='arc':
                aoi_idx = np.where((lat>=65.0) & (lon<=(-85.0)) & (lon>=(-180.0)) | \
                                   (lat>=65.0) & (lon<=(275.0)) & (lon>=(180.0)) | \
                                   (lat>=70.0) & (lon>=70.0) & (lon<=180.0))
            elif self.aoi=='ant':
                aoi_idx = np.where(lat <= -55.0)
            else:
                pass #throw error soon
                
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
                
        #concatenate and convert to a single pd.df
        print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
              '] - Converting data to pandas df...')      
        xo_preproc = pd.concat(xo_preproc,ignore_index=True)
        
        #create gpd.gdf and set/change crs
        print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
              '] - Converting data to geopandas gdf and projecting to EASE2...')
        xo_gdf_dict = {'pts':[],'lin':[]}
        lon = xo_preproc.lon
        lat = xo_preproc.lat
        xo_gdf_dict['pts'] = gpd.GeoDataFrame(xo_preproc,
                                              geometry=gpd.points_from_xy(lon,lat))
        xo_gdf_dict['pts'] = xo_gdf_dict['pts'].set_crs(epsg=4326)
        
        #transform crs
        xo_gdf_dict['pts'] = xo_gdf_dict['pts'].to_crs(epsg=self.epsg)
        
        #convert it to spatial lines
        print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
              '] - Compile orbit tracks...')
        pts = xo_gdf_dict['pts']
        pts_grp = pts.groupby(['id'],as_index=False)['geometry']
        lin = pts_grp.apply(lambda x: shapely.geometry.LineString(x.tolist())\
                            if len(x)>=2 else None)   
        lin_id = lin.id
        xo_gdf_dict['lin'] = gpd.GeoDataFrame(lin_id,geometry=lin.geometry,
                                              crs="EPSG:"+str(self.epsg))
        
        #assign preproc gdf
        return xo_gdf_dict
    
    
    # In[]
    ##xo identifier function
    def identify(self) -> None:
        #allocate sum-up
        attr = []
        geom = []
        
        #load point data from sensor 1/2
        pts1 = self.preproc_c1['pts']
        pts2 = self.preproc_c2['pts']
        
        #load line data from sensor 1/2
        lns1 = self.preproc_c1['lin']
        lns2 = self.preproc_c2['lin']
        
        #convert to datetime to speed up the id process
        avg_t1 = pts1[['id','time']].groupby(['id'], as_index=False).mean()
        dt1_df = pd.to_datetime(avg_t1.time,unit='s')
        avg_t2 = pts2[['id','time']].groupby(['id'], as_index=False).mean()
        dt2_df = pd.to_datetime(avg_t2.time,unit='s')
        
        #loop over all sensor 1 orbits to identify crossovers with 
        #sensor 2 tracks
        for idx1, xo_id1, lin1 in lns1.itertuples():
            print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
                  '] > Processing '+str(idx1+1)+' out of '+\
                      str(len(lns1)))
            print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
                  '] - Identifying intersections with other orbit tracks...')

            #create time-dependent subset of sensor2 orbits to save time
            delta_dt_df = dt1_df.iloc[idx1]-dt2_df
            #calculate absolute date time difference [in hours]
            abs_dtd_hh = (dt2_df-dt1_df.iloc[idx1]).abs()/pd.Timedelta(hours=1)
            #indices to use in subset
            if self.matchtype == 'xo':    
                lns_idx = abs_dtd_hh[abs_dtd_hh <= self.delta_t].index
            else:
                lns_idx = abs_dtd_hh[abs_dtd_hh <= 1.0].index
                
            #loop over the sensor2 orbit subset and collect information
            for idx2, xo_id2, lin2 in lns2.iloc[lns_idx,:].itertuples():
                #summarize data
                if lin1 is not None and lin2 is not None:
                    #find intersection
                    if self.matchtype == 'xo':
                        xo = lin1.intersection(lin2)
                    else:
                        xo = lin1.intersection(lin2.buffer(5000))
                    
                    #only use point-like intersects (no close "fly-by's")
                    if type(xo) is shapely.geometry.point.Point:
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
                    
                    #line features
                    if type(xo) is shapely.geometry.linestring.LineString:
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
                    else:
                        #others stay unused
                        pass
        #store data into a GeoDataFrame
        print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
              '] - Saving results...')
        self.out = gpd.GeoDataFrame(attr,
                                   geometry=geom,
                                   crs="EPSG:"+str(self.epsg),
                                   columns=[self.matchtype+'_f1',self.matchtype+'_f2',\
                                            'f1_rng_0','f1_rng_1',\
                                            'f2_rng_0','f2_rng_1',\
                                            'dt_'+self.matchtype])
        #return output
        print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
              '] < Process complete. :)')               
    
    
    # In[]    
    ##data load functions to reduce code in the dtlim2csv function call
    def load_l1p_nc2pd(self, carrier: str, r0: int, r1: int,
                   nc_handle: Dataset) -> pd.DataFrame:
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

    
    # In[]
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



# In[]

##xo arguments validator
class xo_args_val(object):
    def __init__(self, args: list, ncdict: dict):
        #status
        print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
              '] < Parse user input...')
        print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
              '] > Validate user input...')
        #parse arguments
        self.args = args
        self.ncdict = ncdict
        #check validity
        self.check()

    def check(self) -> None:
        ##perform checks
        #xo-type
        if self.args.matchtype[0].lower() != 'xo' \
            and self.args.matchtype[0].lower() != 'otm':
            msg = 'Error with specified XO-type - Either \'XO\' or \'OTM\' must be chosen!'
            self.abort(msg)  
        #aoi
        if self.args.aoi[0].lower() != 'arc' \
            and self.args.aoi[0].lower() != 'ant':
            msg = 'Error with specified AOI - Either \'arc\' or \'ant\' must be chosen!'
            self.abort(msg)   
        #paths l1p/l2i
        if self.args.aoi[0].lower() == 'arc':
            hemis_long = 'north'
            hemis_short = 'nh'
        else:
            hemis_long = 'south'
            hemis_short = 'sh'
        #L1p
        if self.args.carrier[0].lower() == 'cryosat2': 
            c1_folder = os.path.join(self.args.l1p_input[0],
                                     self.args.carrier[0].lower(),
                                     'ipf1-d/rep/l1p',self.args.c1_l1p_v[0],
                                     hemis_long,
                                     str(self.args.year[0]),
                                     str(self.args.month[0].zfill(2)))
        else:
            c1_folder = os.path.join(self.args.l1p_input[0],
                                     self.args.carrier[0].lower(),
                                     'l1p',self.args.c1_l1p_v[0],
                                     hemis_long,
                                     str(self.args.year[0]),
                                     str(self.args.month[0].zfill(2)))
        if self.args.carrier[1].lower() == 'cryosat2':
            c2_folder = os.path.join(self.args.l1p_input[0],
                                     self.args.carrier[1].lower(),
                                     'ipf1-d/rep/l1p',self.args.c2_l1p_v[0],
                                     hemis_long,
                                     str(self.args.year[0]),
                                     str(self.args.month[0].zfill(2)))
        else:
            c2_folder = os.path.join(self.args.l1p_input[0],
                                     self.args.carrier[1].lower(),
                                     'l1p',self.args.c2_l1p_v[0],
                                     hemis_long,
                                     str(self.args.year[0]),
                                     str(self.args.month[0].zfill(2)))
        if len(os.listdir(c1_folder)) == 0 \
            or len(os.listdir(c2_folder)) == 0:
            msg = 'Error with specified L1p input path/versioning - '+\
                'No input files have been found!'
            self.abort(msg)
        
        #L2i
        if self.args.l2i and self.args.l2i_input != None:
            c1_folder = os.path.join(self.args.l2i_input[0],
                                     self.args.carrier[0].lower(),
                                     self.args.c1_l2i_v[0],
                                     hemis_short,
                                     'l2i',
                                     str(self.args.year[0]),
                                     str(self.args.month[0].zfill(2)))
            c2_folder = os.path.join(self.args.l2i_input[0],
                                     self.args.carrier[1].lower(),
                                     self.args.c2_l2i_v[0],
                                     hemis_short,
                                     'l2i',
                                     str(self.args.year[0]),
                                     str(self.args.month[0].zfill(2)))
            if len(os.listdir(c1_folder)) == 0 \
                or len(os.listdir(c2_folder)) == 0:
                msg = 'Error with specified L2i input path/versioning - '+\
                    'No input files have been found!'
                self.abort(msg)
        else:
            msg = 'Missing paths to L2I data but \'--include-l2i\' was given!'
            self.abort(msg)
            
        #output dir
        if not os.path.isdir(self.args.output[0]):
            #create output path
            os.mkdir(self.args.output[0])
        #output parameter choices
        if not self.ncdict.validate_usr_pars(self.args.par):
            msg = 'Unrecognized CSV output parameter(s) or none at all '+\
                'found in user input!'
            self.abort(msg)
                
    def abort(self,msg) -> None:
        print('['+str(dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))+\
              '] ! '+msg)
        sys.exit()


# In[]

#container class for the NC file entries depending on sensor and 
#production level
class par_dict_container(object):
    def __init__(self):
        #specify parameter dictionaries
        self.l1_dict = {'pwr': 'waveform/power',
                        'rng': 'waveform/range',
                        'rdm': 'waveform/radar_mode',
                        'alt': 'time_orbit/altitude',
                        'lat': 'time_orbit/latitude',
                        'lon': 'time_orbit/longitude',
                        'flg': 'surface_type/flag',
                        'fmi': 'classifier/first_maximum_index',
                        'ppk': 'classifier/peakiness',
                        'lew': 'classifier/leading_edge_width',
                        'lep': 'classifier/leading_edge_peakiness',
                        'leq': 'classifier/leading_edge_quality',
                        'nsp': 'classifier/noise_power',
                        'sig': 'classifier/sigma0',
                        'eps': 'classifier/epsilon_sec',
                        'sku': 'classifier/stack_kurtosis',
                        'spk': 'classifier/stack_peakiness',
                        'ssd': 'classifier/stack_standard_deviation',
                        'ssk': 'classifier/stack_skewness',
                        'src': 'mission_data_source',
                        }
        self.l2_dict = {'sla': 'sea_level_anomaly',
                        'mss': 'mean_sea_surface',
                        'elv': 'elevation',
                        'pdc': 'pulse_deblurring_correction',
                        'd2o': 'distance_to_ocean',
                        'd2i': 'distance_to_low_ice_concentration',
                        'miz': 'flag_miz',
                        }
        #parameter status flags
        self.has_l1_pars = False
        self.has_l2_pars = False
        
        
    def validate_usr_pars(self,pars: list) -> bool:
        #store for later use
        self.pars = pars
        #validate input arguments for presence in pre-defined dict
        all_valid = True
        if pars is not None:
            for key in pars: 
                if key not in self.l1_dict and \
                    key not in self.l2_dict: 
                    all_valid = False
            #return status
            return all_valid
        else:
            #nothing specified by user!
            return False
        
    def compile_carrier_data(self,carrier: str) -> dict:
        #set mandatory file source variable/path
        cdict = {'src': self.l1_dict['src'],
                 'l1': {},
                 'l2': {},
                 }
        #loop over params and add them
        for key in self.pars:
            #ers specific
            if key == 'eps':
                if carrier == 'ers1' or carrier == 'ers2':
                    cdict['l1'][key] = self.l1_dict[key]
            elif key == 'nsp':
                if carrier != 'ers1' or carrier != 'ers2':
                    cdict['l1'][key] = self.l1_dict[key]
            #cryosat specific
            elif key == 'sku' or key == 'spk' or key == 'ssd' or key == 'ssk':
                if carrier == 'cryosat2':
                    cdict['l1'][key] = self.l1_dict[key]
            #level2 product keys
            elif key == 'sla' or key == 'mss' or key == 'elv' or \
                key == 'pdc' or key == 'd2i' or key == 'd2o' or key == 'miz':
                cdict['l2'][key] = self.l2_dict[key]
            #all others
            else:
                cdict['l1'][key] = self.l1_dict[key]
        #set status flags
        self.set_level_flags(cdict)
        #return compiled carrier-specific dict
        return cdict
    
    def set_level_flags(self, cdict: dict) -> None:
        #set status flags according to the carreir dictionary
        if len(cdict['l1'])>0:
            self.has_l1_pars = True
        if len(cdict['l2'])>0:
            self.has_l2_pars = True
        



# In[]
if __name__ == "__main__":
    xo_wrapper()   


