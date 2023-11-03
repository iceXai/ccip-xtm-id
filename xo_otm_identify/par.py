# -*- coding: utf-8 -*-
"""
@author: Dr. Stephan Paul (AWI/iceXai; stephan.paul@awi.de)
"""

# In[]


# In[]

class Parameters:
    def __init__(self):
        #specify supported parameter dictionaries
        l1p_parameters = {'pwr': 'waveform/power',
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
        l2i_parameetrs = {'sla': 'sea_level_anomaly',
                          'mss': 'mean_sea_surface',
                          'elv': 'elevation',
                          'pdc': 'pulse_deblurring_correction',
                          'd2o': 'distance_to_ocean',
                          'd2i': 'distance_to_low_ice_concentration',
                          'miz': 'flag_miz',
                          }
        self.pardict = {'l1p': l1p_parameters,
                        'l2i': l2i_parameters
                        }
        #parameter status flags
        self.has_l1p_parameters = False
        self.has_l2i_parameters = False
        
        
    def validate_cfg_parameters(self, parameters: list) -> bool:
        #store for later use
        self.cfg_parameters = parameters
        #validate input arguments for support
        all_valid = True
        if parameters is None:
            #status
            logger.critical(f'No parameters for CSV output specified!')
            return False
        else:
            for par in parameters: 
                NOT_L1P = par not in self.pardict['l1p']
                NOT_L2I = par not in self.pardict['l2i']
                if NOT_L1P and NOT_L2I:
                    #status
                    logger.critical(f'Specified parameter {par} currently '+\
                                    f'not supported!')
                    all_valid = False
            #return status
            return all_valid
        
    def compile_carrier_data(self, carrier: str) -> Dict[str, str]:
        #set mandatory file source variable/path
        cdict = {'src': self.pardict['l1p']['src'],
                 'l1p': {},
                 'l2i': {},
                 }
        #loop over params and add them
        for par in self.cfg_parameters:
            l1pdict = self.pardict['l1p']
            l2idict = self.pardict['l2i']
            if par in l1pdict.keys():
                cdict['l1p'][par] = l1pdict[par]
            else:
                cdict['l2i'][par] = l2idict[par]
        #set prouct level flags
        self._set_product_flags(cdict)
        #return to caller
        return cdict
    
    def _set_product_flags(self, cdict: Dict[str, str]) -> None:
        #set status flags according to the carreir dictionary
        if len(cdict['l1p'])>0:
            self.has_l1p_parameters = True
        if len(cdict['l2i'])>0:
            self.has_l2i_parameters = True
