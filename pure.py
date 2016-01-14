#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
#%% Imports
from __future__ import division
from scipy.interpolate import interp1d
import logging
import data_handling, Van_der_Waals, numpy
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
def parameters(Data,I):
    """
    Move data container to parameter output dictionary and find the critical 
    Van der Waals contants if not defined.
    
    Parameters
    ----------
    Data : Dictionary containing data loaded from the stored .csv file.
    
    I : Input dictionary, must contain key 'Model' defining the a parameter
        dependancy model.
    """
    p = {'T'   : Data['T (K)'],
         'P'   : Data['P (Pa)'],
         'T_c' : Data['T_c (K)'][0],
         'P_c' : Data['P_c (Pa)'][0],
         'V_c' : Data['V_c (m3 mol-1)'][0],
         'Z_c' : Data['Z_c'][0],
         'R'   : Data['R (m3 Pa K-1 mol-1)'][0],
         'w'   : Data['w'][0],
         'a_c' : Data['a_c (Pa m6 mol-2)'][0],
         'b_c' : Data['b_c (m3 mol-1)'][0],
         'vT'  : Data['virialT'],
         'vB'  : Data['virialB'],
         'Model': I['Model']
         }

    if I['Model'] == 'Adachi-Lu': # Find model params if not defined
        p['m'] = Data['m (Adachi-Lu)'][0]
    elif I['Model'] == 'Soave':
        p['m'] = Data['m (Soave)'][0]

    if p['a_c'] == '' or p['b_c'] == '': 
        p['b_c'] = p['R']*p['T_c']/(8*p['P_c']) 
        p['a_c'] = 27*(p['R']**2)*(p['T_c']**2)/(64.0*p['P_c'])       
    else:
        pass

    for key, value in p.iteritems(): # Filter out '' values
        if not value.__class__ == float:
            p[key] = filter(lambda a: a != '', value)
    
    return p
#%%
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
    #%% Load Data
    try: # Dectect input, use local script if not defined then import data
        Compounds = I['Compound'] # Variable to draw data in data_handling.py
        data = data_handling.ImportData()  
        data.load_pure_data(I['Compound'])
    except(KeyError,NameError): # Define local inputs if I is not found.
        I = inputs() # 
        Compounds = I['Compound'] # Variable to draw data in data_handling.py
        data = data_handling.ImportData()  
        data.load_pure_data(I['Compound'])
        #execfile('data_handling.py') # Load data

    #%% Find model parameters if not defined
    p = parameters(data.c[0],I) # Load parameters as p dictionary
 
     # Find a_c, b_c  parameters if not available and test dimensional values
    if data.c[0]['a_c (Pa m6 mol-2)'][0] == '' \
    or data.c[0]['b_c (m3 mol-1)'][0] == '': 
        try: # Check if all params are available
            data.c[0]['b_c (m3 mol-1)'] = p['R']*p['T_c']/(8.0*p['P_c']) 
            data.c[0]['a_c (Pa m6 mol-2)'] = 27*(p['R']**2)*(p['T_c']**2) \
                                             /(64*p['P_c'])
            if not numpy.round(p['b_c'],8) == \
                   numpy.round(data.c[0]['b_c (m3 mol-1)'],8):
                logging.warn('Calculated parameter for \'b_c\' does'
                              + 'not match stored data value')
                logging.warn('Changing data value for \'b_c\' from'
                             + ' \'b_c\' = {}'.format(p['b_c']) 
                             + ' to \'b_c\' = {}'.format(data.c[0]['b_c'])
                             )
                 
                p['b_c'] = data.c[0]['b_c']
            
            if not numpy.round(p['a_c'],8) == \
                   numpy.round(data.c[0]['a_c (Pa m6 mol-2)'],8):
                logging.warn(' Calculated parameter for \'a_c\' does'
                              + 'not match stored data value.')
                logging.warn(' Changing data value for \'a_c\' from'
                             + ' \'a_c\' = {}'.format(p['a_c'])
                             + ' to \'a_c\' = {}'.format(data.c[0]['a_c'])
                             )
                
                p['a_c'] = data.c[0]['a_c']
                
        except(NameError,KeyError):
            raise IOError('Missing \'P_c\' and/or \'T_c\' paramter')

     # Find 'm' paramter in a dependancy model parameter if not available
    if data.c[0]['m ({})'.format(I['Model'] )][0] == '' \
    or I['Force paramater update']:  # Detect model params
                                     # ore Force optimization if true
        try: # Check if vapour pressure data is available
            data.c[0]['P (Pa)']
            data.c[0]['T (K)']
        except(NameError,KeyError):
            raise IOError('No parameters or vapour pressure data detected.')
            
        #find_a_m() # Find params if data is available
        p['m'] = optim_a_m(p)
        exec 'data.c[0][\'m ({})\'] = p[\'m\']'.format(I['Model'])
        print 'New parameter found for %s model, m = %f'%(I['Model'],p['m'])

    #%% Find phase equilibrium at specified Temperature point (T, V_v and V_l)
    if I['T']: # Note that if I['T'] is > 0 then the boolean is 'True'
        s      = {}
        s['b'] = p['b_c'] # b = b_c
        s['a'] = p['a_c'] # First estimate
        s['T'] = I['T']
        try:
            s['P'] = interp1d(p['T'],p['P'])(s['T'])
        except(ValueError):
            raise IOError('Specified temperature {} K is larger than the criti'
                          'cal temperature {} K.'.format(s['T'],p['T_c']))

        s = VdW.Psat_V_roots(s,p)
        print('VLE at {T} K: P_sat = {P} kPa, V_v = {Vv} m3 mol-1, V_l = {Vl}'
              ' m3 mol-1'.format(T  = I['T'],
                                 P  = s['P_sat']/1000.0,
                                 Vv = s['V_v'],
                                 Vl = s['V_l']))
                                  
    #%% Find phase equilibrium at specified Pressure point (P, V_v and V_l)
    try: #TODO:
        if I['P']: # Note that if I['P'] is > 0 then the boolean is 'True'
            pass#VdW.Tsat_V_roots(s,p) # NOTE TODO!
    except(KeyError):
        pass
    
    #%% Save if True
    """ NOTE!!!: Save the data container and define first data entry [0]
    DO NOT save the dictionary container (or test first if save_to_dict.py is 
    robust))
    """
    if I['Save results']:
        from csvDict import save_dict_as_csv
        # Order of headings to save in .csv
        Order = ['T (K)', 'P (Pa)', 'T_c (K)', 'P_c (Pa)', 'V_c (m3 mol-1)', 
                 'Z_c', 'R (m3 Pa K-1 mol-1)' ,'w' ,'a_c (Pa m6 mol-2)', 
                 'b_c (m3 mol-1)', 'm (Adachi-Lu)', 'm (Soave)', 'virialT', 
                 'virialB','name']
        # Save path string     
        sstr = 'Data/Pure_Component/{}.csv'.format(data.c[0]['name'][0])
        print 'Saving new results to {}'.format(sstr)
        save_dict_as_csv(data.c[0],sstr,Order)
        
    #%% Plotting if True
    try: 
        if I['Plot pure']:
            from numpy import linspace
            s['T_sat store'] = linspace(p['T'][0],p['T'][len(p['T'])-1])
            s['P_sat store'] = []
            s['P_est'] = interp1d(p['T'],p['P'])(s['T_sat store'])
            i = 0
            for T, P in zip(s['T_sat store'][:len(s['T_sat store'])-1],  
                            s['P_est'][:len(s['T_sat store'])-1]): # Trim crit.
                s['T'] = T # Solve P_sat at this Temperature
                s['P'] = P # P est
                s = VdW.Psat_V_roots(s,p,tol=1e-1)
                s['P_sat store'].append(s['P_sat']) 
                
            s['P_sat store'].append(p['P_c']) # Append Critical point  
            plot_Psat(s,p,I['Plot options'])
    except(KeyError):
        pass
    
    #%%
    


    # TESTING AND DEBUGGING
    #s = {}
    #s['P'] = data.c[0]['P (Pa)'][10]/1000
    #s['T'] = data.c[0]['T (K)'][10]
    #s['m'] = data.c[0]['m (Soave)'][0]
    #s['a'] = 27*(p['R']**2)*(p['T_c']**2)/(64*p['P_c'])
    #s['b'] = p['R']*p['T_c']/(8*p['P_c']) 
    #VdW.Psat_V_roots(s,p)
    



    #a = VdW.a_m_sol(s,p)
    #a = VdW.a_m_sol(s,p)
    #VdW.Psat_V_roots(s,p)
