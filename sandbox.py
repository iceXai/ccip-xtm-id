#!/usr/bin/env python3
#coding: utf-8
"""
@author: Dr. Stephan Paul (AWI/iceXai; stephan.paul@awi.de)
"""


# In[]

from xo_otm_identify import CommandLineArgParser
from xo_otm_identify import PreProcessor
from xo_otm_identify import Configuration


# In[]

clarg = CommandLineArgParser()
#clarg.parse_command_line_arguments()
clarg.set_test_cfg_arg('otm_envisat-ers2_full_cfg.yaml')

args = clarg.get_args()

cfg = Configuration(args['cfg'])