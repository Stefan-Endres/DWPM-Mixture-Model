#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

#%% Imports
from __future__ import division
from scipy.interpolate import interp1d
import data_handling, Van_der_Waals, numpy
VdW = Van_der_Waals.VdW()

import decimal
decimal.getcontext().prec = 100
from decimal import Decimal

try: #DEBUGGING; DELETE
    del s
    del p
    del I
except NameError:
    pass

#%% Inputs (will be called if no input container "I" is defined before exec)
def inputs():
    I = {# Model inputs
         'Compounds'    : ['carbon_dioxide','ethane'], # Compound to simulate.
         'Mixture model': 'DWMP', #  'VdW standard' 'DWMP'
         'Model'       : 'Adachi-Lu',   # Model used in the simulation, 
                                     # options:
                                      # 'Soave'                   
                                      # 'Adachi-Lu'   
         
         # Optional inputs
         'T'           : 281.15,  # 281.15
         'P'           : 6278.150329,   # 
         
                    # 6492.799322	x1 = 0.4535 y1= 0.4605
                    # Psat 1 = 6071.500625 @ 283.15
                    # Psat 2 = 6342.145033
         
         # NOTE: Using Psat_V_roots, @ T = 283.15:
          # Psat 1 = 6174.7757144388343
          # Psat 2 = 4759.7154146201701
         
         # NOTE: Using Psat_V_roots, @ T = 287.15:
          # Psat 1 = 7524.4499990006861
          # Psat 2 = 5842.1051091283525
                  
# datapoint = array([0.4535,0.4605]) # @ P = 6492.799322 T = 281.15  

                         #comp. 1  Psat = 5462.722077  
                         #comp. 2  Psat = 5714.421355
                         # HIgh -ish 6492.799322
         'Phase split' : True, # Find phase split at specified T, P.
         'Save results': True,
         'Plot binary' : True,
         'Fig. number' : 2,
         'Plot pure'   : False,       
         'Plot options': {'text.usetex' : True, # Options for all plots
                          'font.size' : 11,            
                          'font.family' : 'lmodern',
                          'text.latex.unicode': True

                          },
         }

    return I