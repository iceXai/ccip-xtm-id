#!/usr/bin/env python3
#coding: utf-8
"""
@author: Dr. Stephan Paul (AWI/iceXai; stephan.paul@awi.de)
"""


# In[]

from xo_otm_identify import CommandLineArgParser
from xo_otm_identify import PreProcessor
from xo_otm_identify import Processor
from xo_otm_identify import Configuration
from xo_otm_identify import Validator


# In[]

#initiate the argument parsser
clarg = CommandLineArgParser()
#clarg.parse_command_line_arguments()
clarg.set_test_cfg_arg('otm_envisat-ers2_full_cfg.yaml')
#get the arguments
args = clarg.get_args()

#initiate the configuration module
cfg = Configuration(args['cfg'])

#initiate validator module and validate configuration input
validator = Validator(cfg)
validator.validate()

#check override status
OVERRIDE_PREPROC = cfg.override

for yy, mm in zip(cfg.yy, cfg.mm):
    #check whether a shapefile already exists
    SHAPEFILE_EXISTS = cfg.shp_exists(yy, mm)
    
    #in case no cossresponding shapefiles exist or the override argument is set
    #initiate the PreProcessor() class and run it
    if not SHAPEFILE_EXISTS or SHAPEFILE_EXISTS and OVERRIDE_PREPROC:
        preproc = PreProcessor(cfg, yy, mm)
        #exectutes the preprocssor with compiling the rawdata, identifying 
        #XOs/OTMs, and exporting the geoinfo to shapefiles
        preproc.run()
    
    #otherwise just run the Processor() class taking care of loading and 
    #compiling/dumping the CSV output
    proc = Processor(cfg, yy, mm)
    #exectutes the procesor by importing the geoinfo shapefiles, compiling
    #the CSV output from identified l1p/l2i files, and exporting it
    proc.run()
        
