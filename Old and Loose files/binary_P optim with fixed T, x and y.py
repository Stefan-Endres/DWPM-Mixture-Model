#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

 #%% Partial Fugacities   
def ln_fug_coeff_1_partial_l(s, p):     ### NOTE ADD EXCEPTION HANDLING FOR log(Z - Beta)
    """ 
    Natural log of 
    Liquid partial fugacity coefficient of liquid phase for i = component 1
    ln(phi^ _i)
    """ 
    from math import log
    amix = a_mix(s,p) # Note already calculated as s.m['a']
    bmix = b_mix(s,p) # Note already calculated as s.m['b']
    #apartial1 =s.m['a_mix_partial']
    #b1 = s.c[1]['b']
    Beta = bmix * s.s['P'] / (p.m['R'] * s.s['T'])
    Z = s.s['P'] * s.m['V_l'] / (p.m['R'] * s.s['T'])
    I = Beta/Z
    q = amix / (bmix * p.m['R'] * s.s['T'])
    #bar_q_i = q * (1 + s.m['a_mix_partial'] / amix - s.c[1]['b'] / bmix)
    bar_q_i = q * (1 + s.m['a_mix_partial'] / amix - s.c[1]['b'] / bmix)
    return (s.c[1]['b'] / bmix) * (Z - 1) - log(Z - Beta) - bar_q_i * I


    
def ln_fug_coeff_1_partial_v(s, p):     ### NOTE ADD EXCEPTION HANDLING FOR log(Z - Beta)
    """ 
    Natural log of 
    Liquid partial fugacity coefficient of liquid phase for i = component 1
    ln(phi^ _i)
    """ 
    from math import log
    amix = a_mix(s,p) 
    bmix = b_mix(s,p)
    #apartial1 =s.m['a_mix_partial']
    #b1 = s.c[1]['b']
    Beta = bmix * s.s['P'] / (p.m['R'] * s.s['T'])
    Z = s.s['P'] * s.m['V_v'] / (p.m['R'] * s.s['T'])
    I = Beta/Z
    q = amix / (bmix * p.m['R'] * s.s['T'])
    #bar_q_i = q * (1 + s.m['a_mix_partial'] / amix - s.c[1]['b'] / bmix)
    bar_q_i = q * (1 + s.m['a_mix_partial'] / amix - s.c[1]['b'] / bmix)
    return (s.c[1]['b'] / bmix) * (Z - 1) - log(Z - Beta) - bar_q_i * I 
    
    
#%% 
def fugacity_error(P, s, p, x_1, y_1):
    """ 
    Used in P_VdW_Multicomp
    
    Returns the error in the goal function:
    x_i * phi_i^l - y_i * phi_i^v == 0
    
    To optimize for P
    
    """
    import math 
    if P < p.m['P'][0]*1e-2 or P > p.m['P'][len(p.m['P'])-1] * 1.1: ### TEEEEST 
        P = s.m['P']#p.m['P']
    # Update P and update pure comp. + mixture 
    s.s['P'], s.m['P'], s.c[1]['P'], s.c[2]['P'] = (P,)*4       
    #%%  Find g_l^R (P,V_l,T,x)
    #Update new x_1 if specified or save x1 and y1 vals.
    s.c[1]['x'], s.c[2]['x'] = x_1, (1.0-x_1)
    #s.c[1]['y'], s.c[2]['y'] = x_1, (1.0-x_1) # (find liquid values first)    
    #%% Find params at current composition and P,T conditions. Note : P CHANGES
    s.c[1]['a'] = VdW.a_T(s.c[1],p.c[1])['a'] # Used to find a_mix
    s.c[2]['a'] = VdW.a_T(s.c[2],p.c[2])['a']

    try: # Note: Highly non-linear models
        s.m['a12'] = a12(s,p)
        s.m['a21'] = a21(s,p)
        s.m['a']   = a_mix(s,p)
        #s.m['b12'] = b12(s,p)  # Not currently in use
        s.m['b']   = b_mix(s,p)
    except (ValueError, ZeroDivisionError): # DO NOT RAISE, SET PENALTY
        s.s['Math Error'] = True
        s.m['a12'], s.m['a21'], s.m['a'], s.m['b12'], s.m['b'] = (0.0,)*5     
        print 'WARNING: Math error in fugacity_error. Failed to calculated '+\
        'a12, a21, a_mix or b_mix for g_l^R, setting to 0.0'
    #%% Find Volume Roots ('V_v' and 'V_l') at P, T, a, b for both components.
    s.m = VdW.V_root(s.m, p.m) # 'V_v' and 'V_l' mixture volumes at x1, x2
    #%% Find Change in gibbs energies of mixing
    try:
        # Find the phi_i^l g_mix_x^R l'
        s.m['a_mix_partial'] = a_mix_partial(s,p)
        #s.m['g_mix_x^R l'] = g_R_l(s.m,p.m)  # Uses V_l at x_1
        s.m['g_mix_x_1 l'] = ln_fug_coeff_1_partial_l(s, p)
#        print s.m['g_mix_x_1 l'] 
        s.m['phi_i^l'] = math.e**s.m['g_mix_x_1 l'] # Uses V_l at x_1
    except (ValueError, ZeroDivisionError):
        s.s['Math Error'] = True
        print 'WARNING: Math error in fugacity_error. Failure to calculate '+\
        'g_mix_x^R l or phi_i^l for g_l^R'
        #print s.m['g_mix_x^R l']
        print s.m['g_mix_x_1 l']
    #%%  Find g_v^R (P,V_v,T,x)
    #%%
    #Update new x_1 if specified or save x1 and y1 vals.
    s.c[1]['x'],s.c[2]['x'] = y_1, (1.0-y_1) # Because y = x in a_mix 
   # s.c[1]['y'],s.c[2]['y'] = y_1, (1.0-y_1) # (find vapour values)      
    #%% Find params at current composition and P,T conditions
    s.c[1]['a'] = VdW.a_T(s.c[1],p.c[1])['a'] # Used to find a_mix 
    s.c[2]['a'] = VdW.a_T(s.c[2],p.c[2])['a'] # TO DO NO NEED TO UPDATE!!!?
    try: # Note: Highly non-linear models
        s.m['a12'] = a12(s,p)
        s.m['a21'] = a21(s,p)
        s.m['a']   = a_mix(s,p)
        #s.m['b12'] = b12(s,p)  # Not currently in use
        s.m['b']   = b_mix(s,p)
    except (ValueError, ZeroDivisionError): # DO NOT RAISE, SET PENALTY
        s.s['Math Error'] = True
        s.m['a12'], s.m['a21'], s.m['a'], s.m['b12'], s.m['b'] = (0.0,)*5
        print 'WARNING: Math error in fugacity_error. Failed to calculated '+\
        'a12, a21, a_mix or b_mix for g_v^R, setting to 0.0'
    #%% Find Volume Roots ('V_v' and 'V_l') at P, T, a, b for both components.
    s.m = VdW.V_root(s.m, p.m) # 'V_v' and 'V_l' mixture volumes at y1, y2
    #%% Find Change in gibbs energies of mixing
    try:
        # Find the phi_i^v g_mix_x^R v'
        s.m['a_mix_partial'] = a_mix_partial(s,p)
        #s.m['g_mix_x^R v'] = g_R_v(s.m,p.m)  # Uses V_v at y_1
        s.m['g_mix_x_1 v'] = ln_fug_coeff_1_partial_v(s, p)
        s.m['phi_i^v'] = math.e**s.m['g_mix_x_1 v'] # Uses V_v at y_1
    except (ValueError, ZeroDivisionError):
        s.s['Math Error'] = True
        print 'WARNING: Math error in fugacity_error. Failure to calculate '+\
        'g_mix_x^R l or phi_i^l for g_v^R'

    #%% Find error from x_i * phi_i^l - y_i * phi_i^v == 0
    #return (x_1*s.m['phi_i^l']*1e10 - y_1*s.m['phi_i^v'] *1e10 )**5
    return x_1*s.m['phi_i^l']- y_1*s.m['phi_i^v']
      


def P_VdW_Multicomp(s, p, P_guess=101.3e3, T=None, x_1=None, y_1=None, \
                    update_pure=None):
    """
    Returns the Pressure of a (binary) phase muxture at s.m['T'], s.c[1]['x'],
    s.c[1]['y']
    
    Initial P_guess=101.3 unless otherwise specified.    
    Specify T=None, x1=None, x2=None to use other values.
    
    """
    from scipy.optimize import fsolve
    #%% Save x1 and y1 if needed
    #if x_1 is None and 'x1' not in s.c[1]: # Save x
    if x_1 is None: # Vals giver
        if 'x' in s.c[1] and 'x' in s.c[2]:
            x_1, x_2 = s.c[1]['x'], s.c[2]['x'] 
        else:
            raise IOError('Error in P_VdW_Multicomp: No x_1 specified')
    #if y_1 is None and 'y1' not in s.c[1]:
    if y_1 is None:
        if 'y' in s.c[1] and 'y' in s.c[2]:
            y_1, y_2 = s.c[1]['y'], s.c[2]['y'] 
        else:
            raise IOError('Error in P_VdW_Multicomp: No y_1 specified')
    #%% Update Temperature if needed
    if T is not None:
        s.m['T'], s.s['T'], s.c[1]['T'], s.c[2]['T'] = (T,)*4
    elif 'T' not in s.m or 'T' not in s.c[1] or 'T' not in s.c[2]:
        raise IOError('Error in P_VdW_Multicomp: No T specified')
    #%% Update Pressure
    s.m['P'], P = (P_guess,)*2   
    # Find optimized Pressure from x_i * phi_i^l(P) - y_i * phi_i^v(P) == 0
    ans = fsolve(fugacity_error,P,args=(s,p,x_1,y_1),xtol=1e-8)
    
    #%% TEST CONVERGENCE:
    testPerr = fugacity_error(ans, s, p, x_1, y_1)
    if abs(testPerr) > 1e-5:
        print 'WARNING: Poor convergence in P_VdW_Multicomp detected, '\
              'x_i * phi_i^l(P) - y_i * phi_i^v(P) == {}.'.format(testPerr) +\
              ' Attempting lower intitial P guess at 0.9*P.m[\'P\']'
        #ans = fsolve(fugacity_error,P-0.9*P,args=(s,p,x_1,y_1),xtol=1e-8)
        if abs(testPerr) > 1e-30:
              print 'WARNING: Poor convergence in P_VdW_Multicomp detected, '\
              'x_i * phi_i^l(P) - y_i * phi_i^v(P) == {}'.format(testPerr) +\
              ''
    return ans 
    
#%% Optimize parameters via pressure minimization.
def P_error(p_set, s, p):
    """
    Returns the error function at specified paramter set p_set.
    
    TO DO: ADD PUNISHMENT FOR MATH ERRORS (Every for loop?)
    
    (Previously inside optim_mixture_parameter(s, p) )
    
    For 'VdW standard' model: p_set is scalar = p.m['k12'] 
    For 'DWMP' model: p_set is a vector [k12, k21, r, s]
    """
    if p.m['Model'] == 'VdW standard':
        p.m['k12'] = p_set
    if p.m['Model'] == 'DWMP':
        p.m['k12'], p.m['k21'], p.m['r'], p.m['s'] = \
        p_set[0],   p_set[1],   p_set[2], p_set[3]
        
    err = 0.0
    for i in range(len(p.m['x1'])-1): # Find error over data range    #   -133
        # Uptate T, P INCLUDE THIS IN RELEVANT FUNCTION
        s.s['P'], s.s['T'] = p.m['P'][i], p.m['T'][i]
        s.m['P'], s.c[1]['P'], s.c[2]['P'] = (s.s['P'],)*3
        s.m['T'], s.c[1]['T'], s.c[2]['T'] = (s.s['T'],)*3
        s.c[1]['x'], s.c[1]['y'] = p.m['x1'][i], p.m['y1'][i]
        s.c[2]['x'], s.c[2]['y'] = 1.0-s.c[1]['x'], 1.0-s.c[1]['y']
        try:
            P = Py_VdW_Multicomp(s, p, P_guess = s.m['P'], T=s.m['T'], \
                                x_1=s.c[1]['x'], y_1=s.c[1]['y'], \
                                update_pure=None)
        except (ValueError, ZeroDivisionError):
            err += 10.0*max(err)
            print 'WARNING: Math Error in P_error, raising goalfunc error.'
        
        err += (p.m['P'][i] - P)**2
    
    
    #if p.m['Model'] == 'DWMP':
    #    return [err, 0, 0, 0]
        
    return err


#%% Error func (over all data points)
def optim_param_set(s, p, guess=None, constraints=None):
    """
    Optimizes the parameter set 
    """
    if p.m['Model'] == 'VdW standard':
        from scipy.optimize import fsolve
        p.m['k12'] = fsolve(P_error,guess,args=(s,p)) #Scalar 
    if p.m['Model'] == 'DWMP':
        from scipy.optimize import minimize
        #p.m['k12'], p.m['k21'], p.m['r'], p.m['s'] = \
        import numpy
        guess = numpy.array(guess)
        #test = fsolve(P_error,guess,args=(s,p)) 
        
        
#%% Plot functions
#%% 
def plot_isotherm(s,p, T_plot=281.15, added_res= 5, SingleFig=False, 
                  model_plot_points=False):
    """
    Plots isotherm and data on a single graph.
    T_plot is the specified isotherm which must be at a known data point.
    """
    # Find sets to be plotted
    Model, Data = {}, {}
    Model['P'], Model['x1'], Model['y1'], Data['P'], Data['x1'], Data['y1'] =\
    [], [], [] , [], [], []
    
    s.s['T'] = T_plot
    s.m['T'], s.c[1]['T'], s.c[2]['T'] = (s.s['T'],)*3
#%%
#    if p.m['Model'] == 'VdW standard': 
    #%% NOTE INDENT THIS CELL
    ind = 0
    for i in range(len(p.m['x1'])): # 0.02449379
    #i = 0
        if p.m['T'][i] == T_plot:
            # Append data storage
            Data['P'].append(p.m['P'][i])
            Data['x1'].append(p.m['x1'][i])
            Data['y1'].append(p.m['y1'][i])
            
            # Find difference between each data point
            nr = float(int(added_res)) # Normalize fraction to add j
                                
            #if i == 0 or i == (len(p.m['x1'])-1) :
            if p.m['x1'][i] < 1e-5 or p.m['y1'][i] < 1e-5 \
            or p.m['x1'][i] > 1.0 - 1e-5  or p.m['y1'][i] > 1.0 - 1e-5:
                x1_dif = 0.0 # Asssume at pure component => No added res
                y1_dif = 0.0 
                P_dif  = 0.0 
                printit = True
                
            else:
                x1_dif = p.m['x1'][i+1] - p.m['x1'][i]
                y1_dif = p.m['y1'][i+1] - p.m['y1'][i]
                P_dif  = p.m['P'][i+1] - p.m['P'][i]

                printit = False
                
            for j in range(int(added_res)):
                Model['x1'].append(p.m['x1'][i] + x1_dif*j/nr )
                Model['y1'].append(p.m['y1'][i] + y1_dif*j/nr )

            # Find model point
                s.s['P'] = p.m['P'][i] + P_dif*j/nr
                s.m['P'], s.c[1]['P'], s.c[2]['P'] = (s.s['P'],)*3
            
                s.c[1]['x'], s.c[1]['y'] = p.m['x1'][i] + x1_dif*j/nr \
                                         , p.m['y1'][i] + y1_dif*j/nr
                
                if printit:
                    print s.m['P']
                    print s.c[1]['x']
                    print s.c[1]['y'] 
                             
                P = P_VdW_Multicomp(s, p, P_guess = s.m['P'], T=s.m['T'], \
                                x_1=s.c[1]['x'], y_1=s.c[1]['y'], \
                                update_pure=None)
                Model['P'].append(P)
                
                    

#                    P = P_VdW_Multicomp(s, p, P_guess = s.m['P'], T=s.m['T'], \
#                                    x_1=s.c[1]['x'], y_1=s.c[1]['y'], \
#                                    update_pure=None)
#                    Model['P'].append(P)
                
#                    print 'Psat (at x1 = {}, y1 = {}) = {}'.format(p.m['x1'][i], \
#                                                               p.m['y1'][i], P)
#                    
#%% ERROR PLOTS 
#            if ind == 1 or ind == 2  or ind == 3 : # TEST ROUTINE DELETE
#                print 'Printing error plot for x_1 = {}' \
#                .format(s.c[1]['x'])
#                fug_err_plot(s,p, i)
#           #fug_err_plot(s,p, i)
#            ind += 1
                