#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
#%% Imports
from __future__ import division
from scipy.interpolate import interp1d
import logging
import data_handling
import numpy
#from common import parameters
import Van_der_Waals
VdW = Van_der_Waals.VdW()
#'''
try:
    del I # DEBUGGING; DELETE
    del s
except(NameError):
    pass
#'''
#%% Inputs (will be called if no input container "I" is defined before exec)
def inputs():
    I = {# Model inputs
         'Compound'    : ['phenol'], # Compound to simulate.
         'Model'       : 'Adachi-Lu' ,   # Model used in the simulation, 
                                     # options:
                                      # 'Soave'                   
                                      # 'Adachi-Lu'   
         
         # Optional inputs (Set to False for off)
         'T'           : 350.0,  # Specify to find P_sat and volume roots at T
                                   # set to False to skip
         'P'           : False,   # ToDO find T_sat and Volume roots at P
                                   # set to False to skip    
         'Save results': True,
         'Plot pure'   : True,       
         'Plot options': {'text.usetex' : True, # Options for all plots
                          'font.size' : 12,            
                          'font.family' : 'lmodern', 
                          'text.latex.unicode': True

                          },
         # If set to True this will force a ew optimization for the 'm' 
          # parameter for the selected 'Model', to be used if new vapour data 
          # is added.                 
         'Force paramater update' : True, 
         'Test'        : False
         }

    return I
#%% Functions
def optim_a_m(p):
    """
    Optimize the paramter 'm' for a pure component for the specified Van der 
    Waals 'a' dependancy model.
    
    Parameters
    ----------
    p : dictionary
        Dictionary containing the data and parameters to fit the parameter 'm'
        to the specified model p['Model'] over p['T'] and p['P'] with specified
        paramters p['a_c'], p['b_c'], p['P_c'], p['T_c'] and p['R'].

    Returns
    -------
    mlsq : float
        Least square fit of 'm' paramater for p['Model']
    
    Dependencies
    ------------
    Van_der_Waals.VdW, scipy, itertools
    
    """
    #%% Initialize
    from scipy.optimize import leastsq
    import itertools
    
    s = {}
    s['b'] = p['b_c'] # b = b_c
    #s['a'] = p['a_c'] # First estimate
    
    p['m'] = 1e-10 # Initial estimate for 'm'.
    #s = VdW.a_T(s, p) ## Initial 'a 'from 'm' estimate
    
    if round(p['P'][len(p['P'])-1],4) == round(p['P_c'],4): 
        p['P'] = p['P'][:len(p['P'])-1] # Trim critical data points.
        p['T'] = p['T'][:len(p['T'])-1]
    # pass
    #%% find a solution for 'a', 'V_l' 'V_v' and 'm' at every T-Psat data point
    def solve_maxwell(i, s, p, tol = 1e-10):
        # Set state 's' = to data point in 'p'
        tol = 1e-10
        s['P'], s['T'] = p['P'][i], p['T'][i]
        s = VdW.a_T(s, p) # Solve 'a' for current temperature.
        er  = 1           # Defines initial error
        while er > tol:
            a_old = s['a'] # Used to track iteration error
            s = VdW.V_root(s, p) # Solve 'V_l' 'V_v' state
            s = VdW.a_maxwell(s, p) # Solve 'a' explicitly 
            er = abs(s['a'] - a_old)
        p = VdW.a_m_sol(s, p) # Update 'm'
        return s['a'],p['m']#s, p
    #%% Compile and solve function at every data point.
    soln = map(solve_maxwell,range(len(p['P'])),
               itertools.repeat(s, len(p['P'])),
               itertools.repeat(p, len(p['P']))
               )       
    # Solution of 'a' at every data point
    a = [soln[i][0] for i in range(len(p['P'])) ]
    m = [soln[i][1] for i in range(len(p['P'])) ] 
    #%%
    def residuals(m, a, T): # (Parameter(s), Output Data, Input Data) # p['T']
        p['m']  = m
        s['T']  = T
        return a - VdW.a_T(s, p)['a'] # Error function to correlate m with T
        
    m0 = sum(m)/len(m)                    # Initial estimate for m
    mlsq = leastsq(residuals, m0, args=(numpy.array(a), numpy.array(p['T'])))
    print mlsq[0][0]
    return mlsq[0]

#%% Plot Functions
def plot_Psat(s,p,Options,fignumber=1): # I['Plot options']
    import matplotlib.pyplot as plot
    plot.figure(fignumber)
    plot.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
    plot.rcParams.update(Options)
    p['P'] = [Pa for Pa in p['P']] # Pa -> kPa
    s['P_sat store'] = [Pa for Pa in s['P_sat store']] # Pa -> kPa
    plot.plot(p['T'], p['P'], 'xr', label='Data points')
    plot.plot(s['T_sat store'],s['P_sat store'], '--r', 
              label='Van der Waals EoS %s m = %s'% (p['Model'], p['m']))
    plot.xlabel("Temperature / K")
    plot.ylabel("Pressure$^{sat}$ / Pa")
    plot.title("Van der Waals EoS correlation for $%s$" \
                % (data.c[0]['name'][0]))
    plot.legend(loc=Options['legend.loc'])
    return

#plot.plot([287.15], [7739.658686], 'ob', label='Data points')
#plot.plot([287.15], [7456.72003], 'or', label='Data points')
# plot.legend(loc=2)

#%% Main
if __name__ == '__main__':
    pass