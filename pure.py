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
import scipy
#from common import parameters
import van_der_waals
VdW = van_der_waals.VdW()

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

#%% Define pure simulation function
def pure_sim(data, i=0):
    """data = data class, i compound to simulate """
    # TODO: replace all c[0] with c[i]

    # Find model parameters if not defined
    p = data_handling.parameter_build(data.c[i])
    #p = parameters(data.c[0],I)
    ## (data.c[0], data.comps)
    s = {}
     # Find a_c, b_c  parameters if not available and test dimensional values
    if (data.c[i]['a_c (Pa m6 mol-2)'][0] == ''
        or data.c[i]['b_c (m3 mol-1)'][0] == ''):

        try: # Check if all params are available
            data.c[i]['b_c (m3 mol-1)'] = p['R']*p['T_c']/(8.0*p['P_c'])
            data.c[i]['a_c (Pa m6 mol-2)'] = 27*(p['R']**2)*(p['T_c']**2) \
                                             /(64*p['P_c'])
            if not numpy.round(p['b_c'], 8) == \
                   numpy.round(data.c[i]['b_c (m3 mol-1)'], 8):
                logging.warn('Calculated parameter for \'b_c\' does'
                              + 'not match stored data value')
                logging.warn('Changing data value for \'b_c\' from'
                             + ' \'b_c\' = {}'.format(p['b_c'])
                             + ' to \'b_c\' = {}'.format(data.c[0]['b_c'])
                             )

                p['b_c'] = data.c[0]['b_c']

            if not numpy.round(p['a_c'], 8) == \
                   numpy.round(data.c[i]['a_c (Pa m6 mol-2)'], 8):

                logging.warn(' Calculated parameter for \'a_c\' does'
                              + 'not match stored data value.')
                logging.warn(' Changing data value for \'a_c\' from'
                             + ' \'a_c\' = {}'.format(p['a_c'])
                             + ' to \'a_c\' = {}'.format(data.c[i]['a_c'])
                             )

                p['a_c'] = data.c[i]['a_c']

        except(NameError,KeyError):
            raise IOError('Missing \'P_c\' and/or \'T_c\' paramter')

     # Find 'm' parameter in a dependency model parameter if not available
    if data.c[i]['m ({})'.format(data.model)][0] == '' \
    or data.force_pure_update:  # Detect model params
                                # or force optimization if true
        try: # Check if vapour pressure data is available
            data.c[i]['P (Pa)']
            data.c[i]['T (K)']
        except(NameError, KeyError):
            raise IOError('No parameters or vapour pressure data detected.')

        #find_a_m() # Find params if data is available
        p['m'] = optim_a_m(p)
        exec 'data.c[0][\'m ({})\'] = p[\'m\']'.format(data.model)
        data.c[i]['m ({})'.format(data.model)] = p['m']
        logging.info('New parameter found for'
                     ' {} model, m = {}'.format(data.model, p['m']))

    #%% Find phase equilibrium at specified Temperature point (T, V_v and V_l)
    if data.T: # Note that if data.T is > 0 then the boolean is 'True'
        s      = {}
        s['b'] = p['b_c'] # b = b_c
        s['a'] = p['a_c'] # First estimate
        s['T'] = data.T
        try:
            s['P'] = scipy.interpolate.interp1d(p['T'],p['P'])(s['T'])
        except(ValueError):
            raise IOError('Specified temperature {} K is larger than the criti'
                          'cal temperature {} K.'.format(s['T'],p['T_c']))

        s = VdW.Psat_V_roots(s, p)
        out_str =('VLE at {T} K: P_sat = {P} kPa, V_v = {Vv} m3 mol-1,'
               ' V_l = {Vl} m3 mol-1'.format(T  = data.T,
                                  P  = s['P_sat']/1000.0,
                                  Vv = s['V_v'],
                                  Vl = s['V_l']))
        print out_str
        logging.info(out_str)

    #%% Find phase equilibrium at specified Pressure point (P, V_v and V_l)
    try: #TODO:
        if data.P: # Note that if I['P'] is > 0 then the boolean is 'True'
            pass#VdW.Tsat_V_roots(s,p) # NOTE TODO!
    except(KeyError):
        pass

    #%% Save if True
    """ NOTE!!!: Save the data container and define first data entry [0]
    DO NOT save the dictionary container (or test first if save_to_dict.py is
    robust))
    """
    if data.save_pure:
        from csvDict import save_dict_as_csv
        # Order of headings to save in .csv
        Order = ['T (K)', 'P (Pa)', 'T_c (K)', 'P_c (Pa)', 'V_c (m3 mol-1)',
                 'Z_c', 'R (m3 Pa K-1 mol-1)' ,'w' ,'a_c (Pa m6 mol-2)',
                 'b_c (m3 mol-1)', 'm (Adachi-Lu)', 'm (Soave)', 'virialT',
                 'virialB','name']
        # Save path string
        sstr = os.path.join(data.datadir,
                            'Pure_Component',
                            '{}.csv'.format(data.c[i]['name'][0]))

        #sstr = 'Data/Pure_Component/{}.csv'.format(data.c[i]['name'][0])
        print 'Saving new results to {}'.format(sstr)
        save_dict_as_csv(data.c[i],sstr,Order)

    return s, p



#plot.plot([287.15], [7739.658686], 'ob', label='Data points')
#plot.plot([287.15], [7456.72003], 'or', label='Data points')
# plot.legend(loc=2)

#%% Main
if __name__ == '__main__':
    pass