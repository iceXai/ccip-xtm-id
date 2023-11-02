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



# In[]

class Processor:
    """
    Class to handle the processing of the given sensor/carrier combo, i.e.,
    the compilation of parameters and waveforms from XOs/OTMs to CSV
    """
    def __init__(self, cfg: Configuration, year: str, month: str):
        pass