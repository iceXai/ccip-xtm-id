#!/usr/bin/env python3
#coding: utf-8
"""
@author: Dr. Stephan Paul (AWI/iceXai; stephan.paul@awi.de)
"""


# In[] 
import argparse


# In[]
"""
Commandline Argument Parser
"""

class CommandLineArgParser(object):
    """
    Commandline argument parser class handling all input for the XO/OTM 
    identification
    """
    def __init__(self):
        #argument dict with all necessary entries for an argparse argument
        self._argsdict = {
                "config-file": {
                    "action": 'store',
                    "dest": 'cfg',
                    "required": True,
                    "type": str,
                    "help": "Config file to be used and stored in ./cfg/"},
                }
        #options list storing the argname used by the parser and a reference 
        #to the argument dict
        self._options = [
                ('--config-file', 'config-file'),
                ]
        #set parser
        self.parser = self._initiate_parser()
        
    def _initiate_parser(self) -> object:
        # create the parser
        parser = argparse.ArgumentParser()
        for option in self._options:
            ARGNAME, ARGDICT_KEY = option
            parser.add_argument(ARGNAME, **self._argsdict[ARGDICT_KEY])
        #return parser
        return parser
    
    def parse_command_line_arguments(self) -> None:
        self.args = self.parser.parse_args()
        
    def set_test_cfg_arg(self, cfg: str) -> None:
        self.args = self.parser.parse_args(["--config-file", cfg])
        
    def get_args(self) -> dict:
        return self.args.__dict__