#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to simulate phase equilibria of multicomponent systems.
"""

# %% Imports
from __future__ import division
import data_handling
import Van_der_Waals
import numpy
from nComp import *
VdW = Van_der_Waals.VdW()

#%% TEST FUNCTION  Binary NRTL 
def g_x_test_func(s, p, k=None, ref='x'):
    """
    This is the test function of a binary NRTL Model of the water-butyl-acetate
    system from Misos et. al. (2007) using the parameters referenced in the 
    paper.
    
    x_1^0 = 0.5 is an unstable point.
    """
    from math import log, e 
    t_12 = 3.00498# tau paramters
    t_21 = 4.69071
    a_12 = 0.391965 # Alpha paramter
    a_21 = 0.391965#**(-1.0) # Checked. Should be a_12 in Mitsos, see SvA p 448
    
    for i in range(1, p.m['n']+1): 
        if s.c[i]['x'] <= 1e-20:  # Prevent math errors from zero log call.
            s.m['g_mix'] = {}
            s.m['g_mix']['t'] = 0.0
            s.m['g_mix']['x'] = s.m['g_mix']['t']
            s.m['g_mix']['ph min'] = 'x'
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
    s.m['g_mix']['ph min'] = 'x'
    return s
                                                
def g_x_test_func2(s, p, k=None, ref='x'):
    """
    This is the test function of a trenary NRTL Model of the toluene_water_
    aniline system from Misos et. al. (2007) using the parameters referenced
    in the paper.
    
    x_1^0 = [0.3, 0.2] is an unstable point.
    """
    from math import log, e 
    t_12 = 4.93035  # tau paramters
    t_21 = 7.77063
    t_13 = 1.59806
    t_31 = 0.03509
    t_23 = 4.18462
    t_32 = 1.27932
    
    a_12 = 0.2485  # Alpha paramters
    a_21 = a_12
    a_13 = 0.3000
    a_31 = a_13 
    a_23 = 0.3412
    a_32 = a_23
    
    for i in range(1, p.m['n']+1): 
        if s.c[i]['x'] <= 1e-20:  # Prevent math errors from zero log call.
            s.m['g_mix'] = {}
            s.m['g_mix']['t'] = 0.0
            s.m['g_mix']['x'] = s.m['g_mix']['t']
            s.m['g_mix']['ph min'] = 'x'
            return s  # should be = 0 as s2['y']*log(s2['y']) = 1*log(1) = 0


    s.m['g_mix'] = {}
    s.m['g_mix']['t'] = ( s.c[1]['x'] * log(s.c[1]['x']) 
                        + s.c[2]['x'] * log(s.c[2]['x']) 
                        + s.c[3]['x'] * log(s.c[3]['x']) 
                        
                        + s.c[1]['x'] * 
                        (t_21 * e**(-a_21 * t_21) * s.c[2]['x'] 
                         + t_31 * e**(-a_31 * t_31) * s.c[3]['x'])
                        / (s.c[1]['x']  
                           + e**(-a_21 * t_21) * s.c[2]['x']  
                           + e**(-a_31 * t_31) * s.c[3]['x']
                           )
                        
                        + s.c[2]['x'] * 
                        (t_12 * e**(-a_12 * t_12) * s.c[1]['x'] 
                         + t_32 * e**(-a_32 * t_32) * s.c[3]['x'])
                        / (e**(-a_12 * t_12) * s.c[1]['x'] 
                           + s.c[2]['x']  
                           + e**(-a_32 * t_32) * s.c[3]['x']
                           )
                                                
                        + s.c[3]['x'] * 
                        (t_13 * e**(-a_13 * t_13) * s.c[1]['x'] 
                         + t_23 * e**(-a_23 * t_23) * s.c[2]['x'])
                          / (e**(-a_13 * t_13) * s.c[1]['x'] 
                           + e**(-a_23 * t_23) * s.c[2]['x'] 
                           + s.c[3]['x']
                           )
                        )
                        
    s.m['g_mix']['x'] = s.m['g_mix']['t']
    s.m['g_mix']['ph min'] = 'x'
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
        if False:  # TODO ADD EXCEPTION HANDLING TO DETECT NEEDED PARAMS
            I['Compound'] = [compound]
            execfile('pure.py')

    # %% Initialize binary and single component paramters
    p = MixParameters()
    p.mixture_parameters(data.VLE, I, data)
    p.m['n'] = len(I['Compounds'])  # Define system size

    for i in range(p.m['n']):  # Set params for all compounds
        p.parameter_build(data.c[i], I)  # Defines p.c[i]

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

    #%% Gibbs surface tests
    if False: # Trenary test function
        s = s.update_state(s, p, P=101e3, T=293.15, # irrelevant in NRTL func
                           X=[0.1,0.1], Force_Update=True)  
                           
        Tie = [[-0.407050993421,# -0.3247905329, # G_P
                0.34759,        # x_1
                0.0917013,      # lambda_1
                0.07562,        # x_2
                0.46985]        # lambda_2
                ]       
                           
        plot_g_mix(s, p, g_x_test_func2, Tie=Tie, x_r=100)
        
        s = s.update_state(s, p, # irrelevant in NRTL func
                           X=[0.34759, 0.07562], Force_Update=True)  
        print g_x_test_func2(s, p, k=None, ref='x').m['g_mix']['t']
        
    #%% Equilibrium Optimization tests   
    if True: # Equilibrium Optimization tests   
        if False: #%% TEST CURVE 1 Mitsos et al. (2007)  ##  Validated 
            Z_0 = numpy.array([0.13])
            Z_0 = numpy.array([0.5])
            s = phase_equilibrium_calculation(s, p, g_x_test_func, Z_0, k=None,
                                              P=101e3, T=300.0, 
                                              tol=1e-9, 
                                              Print_Results=True, 
                                              Plot_Results=True)

        if True: #%% TEST CURVE 2 Mitsos et al. (2007)  ##  Validated 
            Z_0 = numpy.array([0.3, 0.2])
            s = phase_equilibrium_calculation(s, p, g_x_test_func2, Z_0, 
                                              k=None,
                                              P=101e3, T=300.0, 
                                              tol=1e-9, 
                                              Print_Results=True, 
                                              Plot_Results=True) 
        
        #% CO2-Ethane test 
        if False:  ##  Validated 
            p.m['r'], p.m['s'] = 1.0, 1.0
            p.m['k'][1][2] = 0.124
            p.m['k'][2][1] = p.m['k'][1][2]
            Z_0 = numpy.array([0.25])
            s = phase_equilibrium_calculation(s, p, g_mix, Z_0, k=None,
                                              P=24e5, T=263.1, 
                                              tol=1e-9, 
                                              Print_Results=True, 
                                              Plot_Results=True) 


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
    if False: # Generic tests.  ## (Successfully validated against anal. sol.
        def f_x(s, p):  # Simple test func
            return s.c[1]['x']**3.0 - s.c[1]['x']**2.0 + s.c[2]['x']**4.0
            
        X_d = [2.0, 3.0]
        s.update_state(s, p, P=24e5, T=263.1,  X = X_d, Force_Update=True) 
        print '='*25
        print 'Point  X_d = [2.0, 3.0] ' 
        print '='*25
        print 'd=1, z=1, m=1 = {}'.format(FD(f_x, s, p, d=1, z=1, m=1))
        print 'd=1, z=2, m=1 = {}'.format(FD(f_x, s, p, d=1, z=2, m=1))
        print 'd=2, z=1, m=2 = {}'.format(FD(f_x, s, p, d=2, z=1, m=2))
        print 'd=2, z=1, m=1 = {}'.format(FD(f_x, s, p, d=2, z=1, m=1))
        print 'd=2, z=2, m=2 = {}'.format(FD(f_x, s, p, d=2, z=2, m=2))
        print 'Hessian = '
        print hessian(f_x, s, p, dx=1e-6, gmix=False)
        print '='*25
        X_d = [-1.0, -6.0]
        s.update_state(s, p, P=24e5, T=263.1,  X = X_d, Force_Update=True) 
        print 'Point  X_d = [-1.0, -6.0] ' 
        print '='*25
        print 'd=1, z=1, m=1 = {}'.format(FD(f_x, s, p, d=1, z=1, m=1))
        print 'd=1, z=2, m=1 = {}'.format(FD(f_x, s, p, d=1, z=2, m=1))
        print 'd=2, z=1, m=2 = {}'.format(FD(f_x, s, p, d=2, z=1, m=2))
        print 'd=2, z=1, m=1 = {}'.format(FD(f_x, s, p, d=2, z=1, m=1))
        print 'd=2, z=2, m=2 = {}'.format(FD(f_x, s, p, d=2, z=2, m=2))
        print 'Hessian = '
        print hessian(f_x, s, p, dx=1e-6, gmix=False)
        print '='*25
        
    #%% Jacobian and Hessian tests.
    if False: # Binary
        if False: # b_mix tests
            #X_d = [[0.6, 0.4], [0.6, 0.4]]#array([0.6, 0.6]) 
            X_d = [0.5]
            s.update_state(s, p, P=24e5, T=263.1,  X = X_d) 
            print 'd_b_mix_d_x = {}'.format(d_b_mix_d_x(s, p))
            print 'FDM est. = {}'.format( FD(b_mix, s, p, dx=1e-1))
    
            print 'd2_b_mix_d_x2 = {}'.format(d_b_mix_d_x(s, p, d=2))
            print 'FDM est. = {}'.format( FD(b_mix, s, p, d=2, dx=1e-1))
            H = hessian(b_mix, s, p, gmix=False)
            import numpy
            numpy.linalg.det(H)
            H2 = numpy.linalg.eig(H)[0]
            
        
        if False: # Test func stability tests. Test func Mitsos et al. 1
            X_d = [0.13]
            s.update_state(s, p, P=24e5, T=263.1,  X = X_d) 
            H = hessian(g_x_test_func, s, p, dx=1e-6, gmix=True)


    #%% phase_seperation_detection tests.
    if False: # Binary test func Mitsos et al. 1
        if False: # Test an unstable and stable point
            X_d = [0.13]
            s.update_state(s, p, P=24e5, T=263.1,  X = X_d) 
            H = hessian(g_x_test_func, s, p, dx=1e-6, gmix=True)

            Stable = stability(X_d, g_x_test_func, s, p, k=['x'])
            print Stable
            
            X_d = [0.57]
            H = hessian(g_x_test_func, s, p, dx=1e-6, gmix=True)

            Stable = stability(X_d, g_x_test_func, s, p, k=['x'])
            print Stable
            
        if True: # Test detection 
            if False: # Test func 1 LLE
                s = phase_seperation_detection(g_x_test_func, s, p, 
                                               P=101e3, T=300.0,
                                               n=100,
                                               LLE_only=True)
                                               
                plot_g_mix(s, p, g_x_test_func,  
                           Tie= s.m['ph equil']['x'], x_r=1000)
        
            if True: # Test CO2-ethane VLE
                p.m['r'], p.m['s'] = 1.0, 1.0
                p.m['k'][1][2] = 0.124
                p.m['k'][2][1] = p.m['k'][1][2]
                Z_0 = numpy.array([0.25])
                Fd = phase_seperation_detection(g_mix, s, p, 
                                               P=24e5, T=263.1,
                                               n=100,
                                               VLE_only=True)
        
    #%% data range tests
    if False: # CO2_Ethane range
        p.m['r'], p.m['s'] = 1.0, 1.0
        p.m['k'][1][2] = 0.124
        p.m['k'][2][1] = p.m['k'][1][2]
        equilibrium_range(g_mix, s, p, n=100, VLE_only=True)














        
        