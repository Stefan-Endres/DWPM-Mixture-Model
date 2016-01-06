#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to simulate phase equilibria of multicomponent systems.
"""

# %% Imports
from __future__ import division
import data_handling
import Van_der_Waals
from nComp import *
VdW = Van_der_Waals.VdW()

try:  # DEBUGGING; DELETE
    del s
    del p
    del I
except NameError:
    pass


# %% Inputs (will be called if no input container "I" is defined before exec)
def inputs():
    I = {
         # Model inputs
          # Compounds to simulate.
         #'Compounds'    : ['acetone', 'water', 'phenol'], 
         'Compounds'    : ['carbon_dioxide','ethane'], 
         #'Compounds'    : ['acetone','water'], 
         'Mixture model': 'DWPM',  # Removed 'VdW standard', set r = s = 1
         'Model'        : 'Adachi-Lu',  # Activity coefficient Model used in 
                                        # the simulation, 
                                         # options:
                                          # 'Adachi-Lu' 
                                          # 'Soave'

         'Valid phases' : ['x', 'y'], # List of valid phases in equilibrium
                                       # ex. for VLE use ['x', 'y']
                                       # Speciification does not preclude
                                       # LLE detection and calculation.

         # Optional inputs
         'T'           : 281.15,  # 281.15
         'P'           : 6278.150329,   # 
         'Phase split' : True,  # Find phase split at specified T, P.
         'Save results': True,
         'Fig. number' : 2,
         'Plot pure'   : False,
         'Plot options': {'text.usetex' : True,  # Options for all plots
                          'font.size' : 11,
                          'font.family' : 'lmodern',
                          'text.latex.unicode': True
                          },
         }

    return I


#%% TEST FUNCTION  Binary NRTL 
def g_x_test_func(s, p):
    """
    This is the test function of a binary NRTL Model from Misos et. al. (2007)
    using the parameters referenced in the paper.
    
    x_1^0 = 0.5 is an unstable point.
    """
    from math import log, e 
    t_12 = 3.00498# tau paramters
    t_21 = 4.69071
    a_12 = 0.391965 # Alpha paramter
    a_21 = 0.391965#**(-1.0) # Checked. Should be a_12 in Mitsos, see SvA p 448
    
    for i in range(1, p.m['n']+1): 
        if s.c[i]['x'] == 0.0:  # Prevent math errors from zero log call.
            s.m['g_mix'] = {}
            s.m['g_mix']['t'] = 0.0
            s.m['g_mix']['x'] = s.m['g_mix']['t']
            return s  # should be = 0 as s2['y']*log(s2['y']) = 1*log(1) = 0
            
    s.m['g_mix'] = {}
    s.m['g_mix']['t'] = ( s.c[1]['x'] * log(s.c[1]['x']) 
                        + s.c[2]['x'] * log(s.c[2]['x']) 
                        + s.c[1]['x'] * (s.c[2]['x']) 
                        * ((t_12 * e**(-a_12 * t_12)) 
                        / (s.c[2]['x'] + s.c[1]['x'] * e**(- a_12 * t_12))  
                        + (t_21 * e**(-a_21 * t_21)) 
                        / (s.c[1]['x'] + s.c[2]['x'] * e**(- a_21 * t_21))))
                        
    s.m['g_mix']['x'] = s.m['g_mix']['t']
    return s
                                                
             
# %%
if __name__ == '__main__':
    # %% Load Data
    try:  # Dectect input, use local script if not defined then import data
        data = data_handling.ImportData()
        data.load_pure_data(I['Compounds'])
        data.load_E(I['Compounds'])
    except(KeyError, NameError):  # Define local inputs if I is not found.
        I = inputs()
        data = data_handling.ImportData()
        data.load_pure_data(I['Compounds'])
        data.load_E(I['Compounds'])
#        execfile('data_handling.py') # Load data

    # %% Find pure component model parameters if not defined
    for compound in I['Compounds']:
        if False:  # TO DO ADD EXCEPTION HANDLING TO DETECT NEEDED PARAMS
            I['Compound'] = [compound]
            execfile('pure.py')

    # %% Initialize binary and single component paramters
    p = MixParameters()
    p.mixture_parameters(data.VLE, I)
    p.m['n'] = len(I['Compounds'])  # Define system size

    for i in range(p.m['n']):  # Set params for all compounds
        p.parameters(data.c[i], I)  # Defines p.c[i]

    p.m['R'] = p.c[1]['R']  # Use a component Universal gas constant

    # %% Initialize state variables
    s = state()
    s.mixed()  # Define mix state variable, call using s.m['key']
    # Define three component state variables (use index 1 and 2 for clarity)
    for i in range(1, p.m['n']+1):
        s.pure(p, i)  # Call using ex. s.c[1]['key']

    #%% TESTS
    # Test State variable
    p.m['r'], p.m['s'] = 1.0, 1.0
#    s = s.update_state(s, p, P=I['P'], T=I['T'], X=[0.1])
#    s = s.update_state(s, p, P=I['P'], T=I['T'], X=[[0.1,0.2],[0.2,0.2]])
#    s = s.update_state(s, p, P=I['P'], T=I['T'], phase=['x','y'],X=[[0.5,0.2],[0.2,0.2]])
#    s = g_mix(s, p)
        
    
    #%   
#    if False: # Test Gibbs curves
#        p.m['r'], p.m['s'] = 1.0, 1.0
#        p.m['k'][1][2] = 0.124
#        p.m['k'][2][1] = p.m['k'][1][2]
#        s = s.update_state(s, p, P=24e5, T=263.1)
#        x_r = 500
#        g_range(s, p, x_r)
#        plot_dg_mix(s,p)
#        
    #%    
    from numpy import array
#    if False: # Test Duality funcs
#        Lambda = array([1, 2])
#        X_d = array([0.4, 0.3]) 
#        Z_0 = array([0.5, 0.5]) 
#        ubd(Lambda, g_mix, X_d, Z_0, s, p)
#        lbd(X_d, g_mix, Lambda, Z_0, s, p)
#        
#    if False: # Test Binary NRTL Mitsos et al 
#        g_range_test(s, p, x_r=1000)
#        plot_dg_mix_test(s,p)
    
    #%% Equilibrium Optimization tests   
    if True: # Equilibrium Optimization tests           
        if False: #%% TEST CURVE Mitsos et al. (2007)  ##  True: Validated 
            # Single phase test
            X_d = array([0.4]) 
            Z_0 = array([0.5]) 
            # Approx. Solution  Lambda = array([-0.03244]) X_d = array([0.004557])      
            s = s.update_state(s, p, P=101e3, T=300.0,  X = X_d) 
            s = dual_equal(s, p, g_x_test_func, Z_0 , tol=1e-9)
            X_I = s.m['Z_eq']
            
            Args= (g_x_test_func, s.m['Lambda_d'], X_I, s, p, ['All'])
            from tgo import tgo
            X_II = tgo(Eq_sol, [(1e-5, 0.99999)], args=Args, n=100, k_t = 5)
            print 'EQUILIBRIUM SOLUTIONS I: {}'.format(X_I)
            print '                     II: {}'.format(X_II)  

            Z_0 = s.m['Z_eq']

            # Plot error func
            from scipy import linspace
            X_r = linspace(1e-5, 0.9999, 1000)
            #Args= (g_x_test_func, s.m['Lambda_d'], Z_0, s, p, ['x'])

            print 'Feeding Lamda_d = {} to ep. func.'.format(s.m['Lambda_d'])

            plot_g_mix(s, p, g_x_test_func,
                        Tie =[[X_II, X_I]])     
            plot_ep(Eq_sol, X_r, s, p, args=Args)
            
        #% CO2-Ethane test 
        if True: 
            X_d = [array([0.4, 0.2]), array([0.4, 0.2])]
            Z_0 = array([0.4, 0.2]) # Must be vector
            
            p.m['r'], p.m['s'] = 1.0, 1.0
            p.m['k'][1][2] = 0.124
            p.m['k'][2][1] = p.m['k'][1][2]
            s = s.update_state(s, p, P=24e5, T=263.1,  X = [[0.128],[0.229]]) 
            #Z_0 = [0.3]
            Z_0 = [0.25]
            s = dual_equal(s, p, g_mix, Z_0, tol=1e-3)
            X_I = s.m['Z_eq']
            
            Args= (g_mix, s.m['Lambda_d'], X_I, s, p, ['All'])
            from tgo import tgo
            X_II = tgo(Eq_sol, [(1e-5, 0.99999)], args=Args, n=10)
            print 'EQUILIBRIUM SOLUTIONS I: {}'.format(X_I)
            print '                     II: {}'.format(X_II)  

            # Plot error func
            from scipy import linspace
            X_r = linspace(0, 1, 1000)
            print 'Feeding Lamda_d = {} to ep. func.'.format(s.m['Lambda_d'])

            plot_g_mix(s, p, g_mix,
                        Tie =[[X_II, X_I]])     
            plot_ep(Eq_sol, X_r, s, p, args=Args)
    #%% Isotherm tests
        if True:
            # Acetone-Water
            p.m['k'][1][2] = 0.623384831942
            p.m['k'][2][1] = 0.0528180441074
            p.m['r'] = 7.53330786789
            p.m['s'] = 0.107160035705
            
#            T_isos =  p.m['T']   
#            from  more_itertools import unique_everseen
#            T_isos = list(unique_everseen(T_isos))
#            for T_i in T_isos:
#                plot_isotherm(s, p, T_plot = T_i, added_res= 50) 
    #%% Jacobian and Hessian tests.
    if False: # Binary
        # b_mix

        X_d = [[0.6], [0.6]]#array([0.6, 0.6]) 
        s = s.update_state(s, p, P=24e5, T=263.1,  X = X_d) 
        

                
        print 'd_b_mix_d_x = {}'.format(d_b_mix_d_x(s, p))
        print 'FDM est. = {}'.format( d1(b_mix, s, p, dx=1e-1))
    
        print 'd2_b_mix_d_x2 = {}'.format(d_b_mix_d_x(s, p, d=2))
        print 'FDM est. = {}'.format( d1(b_mix, s, p, d=2, dx=1e-1))
    
    
    