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
          # Compound to simulate.
#         'Compounds'    : ['benzene','cyclohexane'], 
#         'Compounds'    : ['acetone','water'], 
         'Compounds'    : ['carbon_dioxide','ethane'], 
 #         'Compounds'    : ['ethane','carbon_dioxide'],
         #'Compounds'    : ['cyclohexane','benzene'], # TEST
         'Mixture model': 'DWMP', # Removed 'VdW standard', set r = s = 1
         'Model'       : 'Adachi-Lu',   # Model used in the simulation, 
                                     # options:
                                      # 'Soave'                   
                                      # 'Adachi-Lu'   
         
         # Optional inputs
         'T'           : 281.15,  # 281.15
         'P'           : 6278.150329,   # 
         
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
#%% Define paramter class
class MixParameters:
    """Store mixture and pure parameters in the same class."""
    def __init__(self):
        #self.name = name
        self.c = []    # creates a new empty list for components 
        self.c.append('nan') #  Define an empty set in index 0 

    def mixture_parameters(self,Data,I):
        """Mixture model parameters"""
        M = {'T'   : Data['T (K)'],
             'P'   : Data['P (Pa)'],
             'x1'  : Data['x1'],
             'y1'  : Data['y1'],
             'Model': I['Mixture model']
             }
        
        if I['Mixture model'] == 'VdW standard':
            M['k12'] = Data['k1_VdW']#[0]
            M['k21'] = Data['k2_VdW']#[0]
            
        if I['Mixture model'] == 'DWMP':
            M['k12'] = Data['k12'],
            M['k21'] = Data['k21'],
            M['r'] =  Data['r'],
            M['s'] = Data['s']

        for key, value in M.iteritems(): # Filter out '' values
            if not value.__class__ == float:
                M[key] = filter(lambda a: a != '', value)
        
        self.m = M

    def parameters(self,Data,I):
        """
        Move data container to parameter output dictionary and find the 
        critical Van der Waals contants if not defined.
        
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
             'name' : Data['name'],
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
        
        self.c.append(p)

#%% Define state variable class
class state:
    """Class defining state vars """
    def __init__(self):
        self.s = {} # System state vars 
        self.c = [] # creates a new empty list for components 
        self.c.append('nan') #  Define an empty set in index 0 
        
    def mixed(self):
        self.m = {} # Mixture states
        
    def pure(self):
        self.c.append({}) # Pure component states
        
#%% Define binary mixture models
def a_mix(s,p):
    """
    # Return a_mix at x_1, x_2, for specified a11, a22, a12, a21, r, s
    """  
    if p.m['Model'] == 'VdW standard':
        return   s.c[1]['a'] * s.c[1]['x']**2 \
               + 2*s.m['a12'] * s.c[1]['x'] * s.c[2]['x'] \
               + s.c[2]['a'] * s.c[2]['x']**2
        
    if p.m['Model'] == 'DWMP': 
        return (\
                s.c[1]['x'] * (s.c[1]['x']*s.c[1]['a']**p.m['s'] \
                + s.c[2]['x']*s.m['a12']**p.m['s'])**(p.m['r']/p.m['s']) \
                + s.c[2]['x'] * (s.c[1]['x']*s.m['a21']**p.m['s'] \
                + s.c[2]['x']*s.c[2]['a']**p.m['s'])**(p.m['r']/p.m['s']) \
                )**(1.0/p.m['r']) 

def a_mix_partial(s,p):
    """
    Return a_mix_partial of component 1
    """
    if p.m['Model'] == 'VdW standard':
        amix = a_mix(s,p)
        return -amix + 2*s.c[2]['x']*a12(s,p) + 2*s.c[1]['x']*s.c[1]['a']
        
    if p.m['Model'] == 'DWMP': 
        amix = a_mix(s,p)
        return (1 - 1/p.m['r'] - 1/p.m['s'])*amix \
               + (1/p.m['r'])*amix**(1-p.m['r']) \
               *( (p.m['r'] /p.m['s'] ) \
               *(s.c[1]['x']*(s.c[1]['x']*s.c[1]['a']**p.m['s'] \
               + s.c[2]['x']*a12(s,p)**p.m['s'])**(p.m['r'] /p.m['s'] - 1) \
               * s.c[1]['a']**p.m['s'] + \
                 s.c[2]['x']*(s.c[1]['x']*a21(s,p)**p.m['s'] \
               + s.c[2]['x']*s.c[2]['a']**p.m['s'])**(p.m['r'] /p.m['s'] - 1) \
               * a21(s,p)**p.m['s'])  \
               + (s.c[1]['x']*s.c[1]['a']**p.m['s'] \
               + s.c[2]['x']*a12(s,p)**p.m['s'])**(p.m['r']/p.m['s']))
#%%
def a_mix_partial_i(s,p, i = 1):
    """
    Return a_mix_partial of component i TO DO use for loops to exec sums from strings
    """
    if p.m['Model'] == 'VdW standard':
        amix = a_mix(s,p)
        
        if i == 1: k = 2 # For binary case
        if i == 2: k = 1 # For binary case
            
        return -amix + 2*s.c[k]['x']*a12(s,p) + 2*s.c[i]['x']*s.c[i]['a']
        
    if p.m['Model'] == 'DWMP': 
        amix = a_mix(s,p)
        return (1 - 1/p.m['r'] - 1/p.m['s'])*amix \
               + (1/p.m['r'])*amix**(1-p.m['r']) \
               *( (p.m['r'] /p.m['s'] ) \
               *(
               s.c[1]['x']*(s.c[1]['x']*s.c[1]['a']**p.m['s'] \
               + s.c[2]['x']*a12(s,p)**p.m['s'])**(p.m['r'] /p.m['s'] - 1) \
               * akj(s,p,i=1,j=i)**p.m['s'] + \
               s.c[2]['x']*(s.c[1]['x']*a21(s,p)**p.m['s'] \
               + s.c[2]['x']*s.c[2]['a']**p.m['s'])**(p.m['r'] /p.m['s'] - 1) \
               * akj(s,p,i=2,j=i)**p.m['s'])  \
               +(s.c[1]['x']*akj(s,p,i=i,j=1)**p.m['s']  \
               + s.c[2]['x']*akj(s,p,i=i,j=2)**p.m['s'])**(p.m['r']/p.m['s']))

#%%
def b_mix(s,p):
    return s.c[1]['x']*s.c[1]['b'] + s.c[2]['x']*s.c[2]['b']
    
def akj(s,p,i=1,j=1): # find a for specified indices
    from math import sqrt
    if i == j:
        return s.c[i]['a'] # Return pure paramater
    else: # find mixture aij i =/= j
        exec 'ij = p.m[\'k{}{}\']'.format(i,j) 
        return (1 - ij) * sqrt(s.c[i]['a'] * s.c[j]['a'])     
    
def a12(s,p): # TO DO TRY MODEL  p.m['k12'] ==  k1*x1 + k2*x2 
    from math import sqrt
    return (1 - p.m['k12']) * sqrt(s.c[1]['a'] * s.c[2]['a'])   
    
def a21(s,p): 
    from math import sqrt
    return (1 - p.m['k21']) * sqrt(s.c[1]['a'] * s.c[2]['a'])   

def b12(s,p): # NOTE: Not currently in use for VdW Standard or DWMP! 
    from math import sqrt
    return sqrt(s.c[1]['b']*s.c[2]['b']) 
#%% Partial Fugacities
def ln_fug_coeff_i_partial(s, p, i = 1, phase = 'l'):     ### NOTE ADD EXCEPTION HANDLING FOR log(Z - Beta)
    """ 
    Natural log of 
    Liquid partial fugacity coefficient of specified phase for i 
    ln(phi^ _i)
    phase = string 'l' for liquid or 'v' for vapour.
    
    NOTE: Volumes must be updated to current phase under consideration
    """ 
    from math import log
    # Mixture parameters at phase 
    amix = a_mix(s,p) # Note already calculated as s.m['a']
    bmix = b_mix(s,p) # Note already calculated as s.m['b']
    #apartial1 =s.m['a_mix_partial']
    #b1 = s.c[1]['b']
    if phase == 'l': # Use correct phase volume
        V = s.m['V_l'] 
    elif phase =='v':
        V = s.m['V_v'] 
    
    # Find partial properties of specified component
    s.m['a_mix_partial'] = a_mix_partial_i(s, p, i) #check
        
    Beta = bmix * s.s['P'] / (p.m['R'] * s.s['T'])
    Z = s.s['P'] * V / (p.m['R'] * s.s['T'])
    I = Beta/Z
    q = amix / (bmix * p.m['R'] * s.s['T'])
    #bar_q_i = q * (1 + s.m['a_mix_partial'] / amix - s.c[1]['b'] / bmix)
    bar_q_i = q * (1 + s.m['a_mix_partial'] / amix - s.c[i]['b'] / bmix)
    return (s.c[i]['b'] / bmix) * (Z - 1) - log(Z - Beta) - bar_q_i * I
    
#%% g - g_ref (add g of volume phase)
def g_R_v(s,p):
    '''g_R_v(T,P) of a single component'''
    from math import log
    return s['P']*s['V_v']/(p['R']*s['T']) - 1.0 - log(s['P']/(p['R']*s['T']))\
    - log(s['V_v'] - s['b']) - s['a']/(p['R']*s['T']*s['V_v'])
    
def g_R_l(s, p):
    '''g_R_l(T,P) of a single component'''
    from math import log
    return s['P']*s['V_l']/(p['R']*s['T']) - 1.0 - log(s['P']/(p['R']*s['T']))\
    - log(s['V_l'] - s['b']) - s['a']/(p['R']*s['T']*s['V_l'])
    
def del_g_IG_v(s1,s2):
    """del_g_IG_v(x), change in gibbs energy for mixing of ideal gasses"""
    from math import log
    if s1['y'] == 0.0: 
        return 0.0 # should be = 0 as s2['y']*log(s2['y']) = 1*log(1) = 0
    if s2['y'] == 0.0:
        return 0.0
    return s1['y']*log(s1['y']) + s2['y']*log(s2['y'])
    
def del_g_IG_l(s1,s2):
    """del_g_IG_v(x), change in gibbs energy for mixing of ideal solutions"""
    from math import log
    if s1['x'] == 0.0: 
        return 0.0 # should be = 0 as s2['x']*log(s2['x']) = 1*log(1) = 0
    if s2['x'] == 0.0:
        return 0.0
    return s1['x']*log(s1['x']) + s2['x']*log(s2['x'])   
    
def del_g_v_min_l0(s,p):
    """This function returns the analytical solution of the integral between 
    the Gibbs energy of the Vapour phase of a pure component and the liquid 
    phase pure component reference state."""
    from math import log
    return 2*s['a']/(p['R'] * s['T']) *((1/s['V_l']) - (1/s['V_v'])) \
    + s['b'] * ((1/(s['V_v'] - s['b'])) - (1/(s['V_l'] - s['b']))) \
    + log((s['V_l'] - s['b'])/(s['V_v'] - s['b']))
            
#%% dg_mix
def del_g_mix(s, p, x_1=None, update_pure=False, g_ref=None, derivative=None):
    """ 
    Returns the change in gibbs energy Gmix/RT at specified composition
    x1 = x[0], Pressure and Temperature
    """
    #%% Update new x_1 if specified
    if x_1 is not None:
        s.c[1]['x'],s.c[2]['x'] = x_1, (1.0-x_1)
        s.c[1]['y'],s.c[2]['y'] = s.c[1]['x'],s.c[2]['x']    
    #%% Find params at current composition and P,T conditions
    if update_pure:
        s.c[1]['a'] = VdW.a_T(s.c[1],p.c[1])['a']
        s.c[2]['a'] = VdW.a_T(s.c[2],p.c[2])['a']
    try: # Note: Highly non-linear models
        s.m['a12'] = a12(s,p)
        s.m['a21'] = a21(s,p)
        s.m['a']   = a_mix(s,p)
        #s.m['b12'] = b12(s,p)  # Not currently in use
        s.m['b']   = b_mix(s,p)
    except (ValueError, ZeroDivisionError): # DO NOT RAISE, SET PENALTY
        s.s['Math Error'] = True
        print 'WARNING: Math Domain error at %s Pa %s K' %(s.s['P'],s.s['T']) 
        #"""UNCOMMENT!!!!!!!!!
        s.m['a12'], s.m['a21'], s.m['a'], s.m['b12'], s.m['b'] = (0.0,)*5
    #%% Find Volume Roots ('V_v' and 'V_l') at P, T, a, b for both components.
    if update_pure:
        s.c[1], s.c[2] = VdW.V_root(s.c[1], p.c[1]), VdW.V_root(s.c[2], p.c[2])
    s.m = VdW.V_root(s.m, p.m) # 'V_v' and 'V_l' mixture volumes at x1, x2
    #%% Find Change in gibbs energies of mixing
    try:
        # del g_mix^R =  g_mix^R(T,P) - SIGMA x_i * g_i^R(T,P)
         # del g_mix^R for each phases.
        if g_ref == None:
            g_Res_Pure = - s.c[1]['x']*g_R_l(s.c[1],p.c[1]) \
                         - s.c[2]['x']*g_R_l(s.c[2],p.c[2])
        else:
            g_Res_Pure = g_ref 

        s.m['del g_mix^R l'] = g_R_l(s.m,p.m) \
                               + g_Res_Pure
#                               - s.c[1]['x']*g_R_l(s.c[1],p.c[1]) \
#                               - s.c[2]['x']*g_R_l(s.c[2],p.c[2])

               # NOTE!!! x = y = linspace for mapping del g_mix^R
        s.m['del g_mix^R v'] = g_R_v(s.m,p.m) \
                               + g_Res_Pure
#                               - s.c[1]['x']*g_R_l(s.c[1],p.c[1]) \
#                               - s.c[2]['x']*g_R_l(s.c[2],p.c[2])
                               #- s.c[1]['y']*g_R_v(s.c[1],p.c[1]) \
                               #- s.c[2]['y']*g_R_v(s.c[2],p.c[2])
        
        # Calculate pure vapour phase total Gibbs energy minus liquid reference
         # state
#        s.c[1]['g_1^V,0 - g_1^l,0'] = del_g_v_min_l0(s.c[1],p.c[1])
#        s.c[2]['g_2^V,0 - g_2^l,0'] = del_g_v_min_l0(s.c[2],p.c[2])
                                
        # del g_mix(T,P,x) = del g_mix^R(T,P,x) + del g^IG(T,P,x)                
        s.m['del g_mix l'] = s.m['del g_mix^R l'] + del_g_IG_l(s.c[1],s.c[2])
                                                        
        s.m['del g_mix v'] = s.m['del g_mix^R v'] + del_g_IG_v(s.c[1],s.c[2]) \
                             #+ s.c[1]['y']*(s.c[1]['g_1^V,0 - g_1^l,0']) \
                             #+ s.c[2]['y']*(s.c[2]['g_2^V,0 - g_2^l,0']) 

        #s.m['del g_mix^R'] = min(s.m['del g_mix^R l'], s.m['del g_mix^R v'])
        
        s.m['del g_mix'] = min(s.m['del g_mix l'], s.m['del g_mix v'])

        s.s['Math Error'] = False
        
        #%%
    except (ValueError, ZeroDivisionError):
        s.s['Math Error'] = True
        #print 'WARNING: Math Domain error in del_g_mix(s,p)!' """UNCOMMENT!!!
        s.m['del g_mix l'], s.m['del g_mix v'] = 0.0, 0.0
        s.m['del g_mix'] = 0.0

    return s
#%% dg Total TEST
def d_g_mix_T(s, p, x_1=None, y_1=None, update_pure=False, derivative=None):
    
    #%# FIND REFERENCE
    if x_1 is not None:
        s.c[1]['x'],s.c[2]['x'] = x_1, (1.0-x_1)
        s.c[1]['y'],s.c[2]['y'] = s.c[1]['x'],s.c[2]['x']    
    if update_pure:
        s.c[1]['a'] = VdW.a_T(s.c[1],p.c[1])['a']
        s.c[2]['a'] = VdW.a_T(s.c[2],p.c[2])['a']
        
    s.m['a12'] = a12(s,p)
    s.m['a21'] = a21(s,p)
    s.m['a']   = a_mix(s,p)
    s.m['b']   = b_mix(s,p)
    if update_pure:
        s.c[1], s.c[2] = VdW.V_root(s.c[1], p.c[1]), VdW.V_root(s.c[2], p.c[2])
    s.m = VdW.V_root(s.m, p.m) # 'V_v' and 'V_l' mixture volumes at x1, x2
    
    
    g_Res_Pure_l = - s.c[1]['x']*g_R_l(s.c[1],p.c[1]) \
             - s.c[2]['x']*g_R_l(s.c[2],p.c[2])
    g_Res_Pure_v = - y_1*g_R_l(s.c[1],p.c[1]) \
             - (1 - y_1)*g_R_l(s.c[2],p.c[2])   
             
             
    return (min(del_g_mix(s, p, x_1, update_pure, 
                      g_ref=g_Res_Pure_l).m['del g_mix l'],  # Sigma n_il g_il
                del_g_mix(s, p, y_1, update_pure, 
                      g_ref=g_Res_Pure_v).m['del g_mix v'])) # Sigma n_iv g_iv
      
#%% TEMPORARY (CENTRAL FD) DERIVATIVE ESTIMATE; replace with analyt-est. hybrid
def d_del_g_mix_dx1(s, p, x_1, dx=1e-6):
    return (del_g_mix(s, p, x_1 + dx).m['del g_mix'] \
            - del_g_mix(s, p, x_1 - dx).m['del g_mix'])/dx
#%%
def g_range(s, p, x_r): 
    """
    Formerly in phase_split func, solve g_mix(x) over x at a P, T
    TO DO OPTIMIZE WITH MAP FUNCTION/ COMPILE IN C / def append etc."""
    from numpy import linspace
    #% Initialize
    s.s['del g_mix l soln'], s.s['del g_mix v soln'] = [], []
    s.s['del g_mix soln'], s.s['Math Error soln'] = [], []
    s.c[1]['x_range'], s.c[2]['x_range'] = linspace(0,1,x_r), linspace(1,0,x_r)
    #% Find pure component parameters at current s['P'], s['T']
    s.c[1]['a'] = VdW.a_T(s.c[1],p.c[1])['a']
    s.c[2]['a'] = VdW.a_T(s.c[2],p.c[2])['a']
    s.c[1], s.c[2] = VdW.V_root(s.c[1], p.c[1]), VdW.V_root(s.c[2], p.c[2])
    #% Solve for x_range 
     # TO DO replace s.c[1]['x_range'] with generator function
    for i in range(len(s.c[1]['x_range']-1)): 
        s.c[1]['x'],s.c[2]['x'] = s.c[1]['x_range'][i],s.c[2]['x_range'][i]
    # To avoid future confusions, we will set extra var y1 = x1, y2 = x2
        s.c[1]['y'],s.c[2]['y'] = s.c[1]['x'],s.c[2]['x'] 
        s = del_g_mix(s, p, update_pure=False)   #.m['del g_mix l']
        s.s['del g_mix l soln'].append(s.m['del g_mix l'])
        s.s['del g_mix v soln'].append(s.m['del g_mix v'])
        s.s['del g_mix soln'].append(s.m['del g_mix'])
        s.s['Math Error soln'].append(s.s['Math Error'])
   
    return s
                    
#%%  
def pseudo_analytical_sys(x, s, p, dx_constraint = 1e-20):
    """
    Returns a system of the 4 Pseudo Analytical equations using the common 
    tangent.
    """
    [s.c[1]['x_alpha'], s.c[1]['x_beta'], m, c] = x
    #% Constraints
    P = 0.0 # Penalty
    if c > 1e-15: # If c not <= 0
        P += 1.1**(10.0*abs(c))
    if m + c > 1e-15:
        P += 1.1**(10.0*abs(abs(m) + abs(c)))
    if abs(s.c[1]['x_alpha'] - s.c[1]['x_beta'])  < dx_constraint:
        P += 1.1**(1e-5*abs(1/(s.c[1]['x_alpha'] - s.c[1]['x_beta']) ))
    #% Return sys
    return [del_g_mix(s, p, x_1=s.c[1]['x_alpha']).m['del g_mix'] \
            - (m*s.c[1]['x_alpha'] + c) + P,
            d_del_g_mix_dx1(s, p, x_1=s.c[1]['x_alpha']) - m  + P,
            del_g_mix(s, p, x_1=s.c[1]['x_beta']).m['del g_mix'] \
            - (m*s.c[1]['x_beta'] + c)  + P,
            d_del_g_mix_dx1(s, p, x_1=s.c[1]['x_beta']) - m  + P
            ]

#%%
def phase_split(s, p, method='Pseudo-Analytical', est = None, x_r = 100):
    """
    Note if est = None, phases are detectect and found otherwise only a 
    specific equilibruim point is found, if it exists.
    
    TO DO: - Implement d^2 g^R / d x_1 ^2 analysis to detect LLE (Kontogeoris,
        p. 456)
    
     - return s.m['x1'] s.m['y1'] at specified parameters and P, T 

    If est is not None:
    Returns the phase split mole fractions x1 and y1 calculated near the 
    specified data or estimated value pair est=[x1_d,y1_d].
    
    NOTE, pressures and pure component params must be updated in s before call    
    
    TO DO; put "Find equilibrium close to the estimated phase split point" cell
    into its own function so that multiple equilibrium points can be found, if
    detected.
    """
    #%% if est in None, detect phase split
    if est is None:
        from scipy.signal import argrelextrema
        from numpy import linspace, less, array
        s = g_range(s, p, x_r) # Get points to detect phase seperation
        s.s['g_min'] = argrelextrema(array(s.s['del g_mix soln']), less) 
        # TO DO; FINISH THIS SCRIPT
         # SET est to  d^2 g^R / d x_1 ^2 = 0 roots
         # SET MEthond to (Kontogeoris, p. 456)
        
    #%% Find equilibrium close to the estimated phase split point
    if est is not None:
        from scipy.optimize import fsolve

        #% Update pure volumes and paramters at current s['P'], s['T']
        s.c[1]['a'] = VdW.a_T(s.c[1],p.c[1])['a']
        s.c[2]['a'] = VdW.a_T(s.c[2],p.c[2])['a']
        s.c[1], s.c[2] = VdW.V_root(s.c[1], p.c[1]), VdW.V_root(s.c[2], p.c[2])
        
        if method == 'Pseudo-Analytical':
            #% Estimates for solver
            [s.c[1]['x_alpha'], s.c[1]['x_beta']] = est
            m_est = (del_g_mix(s, p, x_1=s.c[1]['x_beta']).m['del g_mix'] \
                    - del_g_mix(s, p, x_1=s.c[1]['x_alpha']).m['del g_mix']) \
                    /(s.c[1]['x_beta'] - s.c[1]['x_alpha'] + 1e-10)
            c_est = del_g_mix(s, p, x_1=s.c[1]['x_beta']).m['del g_mix'] \
                    - m_est*s.c[1]['x_beta']
            #% Solve
            dx_constraint = 0.2*(s.c[1]['x_beta'] - s.c[1]['x_alpha'])
            ans = fsolve(pseudo_analytical_sys,
                         [s.c[1]['x_alpha'], s.c[1]['x_beta'], m_est, c_est],
                         args=(s, p, dx_constraint) )
    #%% Find phase of equilibrium point TO DO; RECHECK THIS CELL
        s.c[1]['x_alpha'], s.c[1]['x_beta'] = ans[0], ans[1]
        g_a = del_g_mix(s, p, x_1=s.c[1]['x_alpha'])
        g_b = del_g_mix(s, p, x_1=s.c[1]['x_beta'])

        if g_a.m['del g_mix^R l'] < g_a.m['del g_mix^R v']:
            s.m['x1'] = s.c[1]['x_alpha']
            if g_b.m['del g_mix^R v'] < g_b.m['del g_mix^R l']:
                   s.m['y1'] = s.c[1]['x_beta']
            else: # Same phase!
                s.m['No Phase Penatly'] = True
                s.m['x1'], s.m['y1'] = (0.0,)*2
                
        elif g_a.m['del g_mix^R v'] < g_a.m['del g_mix^R l']:
            s.m['y1'] = s.c[1]['x_alpha']
            if g_b.m['del g_mix^R l'] < g_b.m['del g_mix^R v']:
                   s.m['x1'] = s.c[1]['x_beta']
            else: # Same phase!
                s.m['No Phase Penatly'] = True
                s.m['x1'], s.m['y1'] = (0.0,)*2
                
        elif g_a.m['del g_mix^R v'] == g_a.m['del g_mix^R l']:
                s.m['No Phase Penatly'] = True
                s.m['x1'], s.m['y1'] = (0.0,)*2
        #s.m['ans'] = ans # DEBUGGING; DELETE

    return s, p
#%% Optimization and Error Functions
def Sigma_K(Py, s, p, x_1):
    """ Returns K values used in y_i_error and Bubble_error"""
    import math 
    P = Py[0] # First entry in variable vector = pressure
    y_1 = Py[1] # Second entry in variable vector = y_1

    
    if P < p.m['P'][0]*1e-2 or P > p.m['P'][len(p.m['P'])-1] * 1.1: ### TEEEEST 
        P = s.m['P']#p.m['P']
    # Update P and update pure comp. + mixture 
    s.s['P'], s.m['P'], s.c[1]['P'], s.c[2]['P'] = (P,)*4       
    
    #%%  Find g_l^R (P,V_l,T,x) for component i
    #Update new x_1 if specified or save x1 and y1 vals.
    s.c[1]['x'], s.c[2]['x'] = x_1, (1.0-x_1)
    #s.c[1]['y'], s.c[2]['y'] = x_1, (1.0-x_1) # (find liquid values first)    
    #% Find params at current composition and P,T conditions. Note : P CHANGES
    s.c[1]['a'] = VdW.a_T(s.c[1],p.c[1])['a'] # Used to find a_mix
    s.c[2]['a'] = VdW.a_T(s.c[2],p.c[2])['a']
    try: # Note: Highly non-linear models
        s.m['a12'] = a12(s,p)
        s.m['a21'] = a21(s,p)
        s.m['a']   = a_mix(s,p)
        #s.m['b12'] = b12(s,p)  # Not currently in use
        s.m['b']   = b_mix(s,p)
    except (ValueError, ZeroDivisionError, 
            numpy.linalg.linalg.LinAlgError): # DO NOT RAISE, SET PENALTY
        s.s['Math Error'] = True
        s.m['a12'], s.m['a21'], s.m['a'], s.m['b12'], s.m['b'] = (1e-5,)*5     
        print 'WARNING: Math error in fugacity_error. Failed to calculated '+\
        'a12, a21, a_mix or b_mix for g_l^R, setting to 0.0'

    try:
        #% Find Volume Roots ('V_v' and 'V_l') at P, T, a, b for both components.
        s.m = VdW.V_root(s.m, p.m) # 'V_v' and 'V_l' mixture volumes at x1, x2
        #% Find Change in gibbs energies of mixing
        # Find the phi_1^l g_mix_x^R l'
        # Find Partial fugacities
        s.c[1]['g_mix_x_i l'] = ln_fug_coeff_i_partial(s,p,i=1,phase='l') 
        s.c[2]['g_mix_x_i l'] = ln_fug_coeff_i_partial(s,p,i=2,phase='l') 
        #s.m['phi_i^l'] = math.e**s.m['g_mix_x_1 l'] # Uses V_l at x_
        # Find phi_i^l
        s.c[1]['phi^l'] = math.e**s.c[1]['g_mix_x_i l'] # Uses V_l at x_1
        s.c[2]['phi^l'] = math.e**s.c[2]['g_mix_x_i l'] # Uses V_l at x_1

    except (ValueError, ZeroDivisionError, KeyError, 
            numpy.linalg.linalg.LinAlgError):
        s.s['Math Error'] = True
        print 'WARNING: Math error in fugacity_error. Failure to calculate '+\
        'g_mix_x^R l or phi_i^l for g_l^R, setting to 1e-5'
        s.c[1]['phi^l'] = 1e-5
        s.c[2]['phi^l'] = 1e-5
        #print s.m['g_mix_x^R l']
        try: # TEST DELETE
            print s.m['g_mix_x_i l']
        except (KeyError):
            pass
    #%%  Find g_v^R (P,V_v,T,x)
    #%%
    #Update new x_1 if specified or save x1 and y1 vals.
    s.c[1]['x'],s.c[2]['x'] = y_1, (1.0-y_1) # Because y = x in a_mix 
    # s.c[1]['y'],s.c[2]['y'] = y_1, (1.0-y_1) # (find vapour values)      
    #% Find params at current composition and P,T conditions
    s.c[1]['a'] = VdW.a_T(s.c[1],p.c[1])['a'] # Used to find a_mix 
    s.c[2]['a'] = VdW.a_T(s.c[2],p.c[2])['a'] # TO DO NO NEED TO UPDATE!!!?
    try: # Note: Highly non-linear models
        s.m['a12'] = a12(s,p)
        s.m['a21'] = a21(s,p)
        s.m['a']   = a_mix(s,p)
        #s.m['b12'] = b12(s,p)  # Not currently in use
        s.m['b']   = b_mix(s,p)
    except (ValueError, ZeroDivisionError, 
            numpy.linalg.linalg.LinAlgError): # DO NOT RAISE, SET PENALTY
        s.s['Math Error'] = True
        s.m['a12'], s.m['a21'], s.m['a'], s.m['b12'], s.m['b'] = (1e-5,)*5
        print 'WARNING: Math error in fugacity_error. Failed to calculated '+\
        'a12, a21, a_mix or b_mix for g_v^R, setting to 1e-5'

    try:
        #% Find Volume Roots ('V_v' and 'V_l') at P, T, a, b for both components.
        s.m = VdW.V_root(s.m, p.m) # 'V_v' and 'V_l' mixture volumes at y1, y2
        #% Find Change in gibbs energies of mixing
        # Find the phi_1^v g_mix_x^R v'
        # Find Partial fugacities
        s.c[1]['g_mix_x_i v'] = ln_fug_coeff_i_partial(s,p,i=1,phase='v') 
        s.c[2]['g_mix_x_i v'] = ln_fug_coeff_i_partial(s,p,i=2,phase='v') 
        # Find phi_i^l
        s.c[1]['phi^v'] = math.e**s.c[1]['g_mix_x_i v'] # Uses V_v at y_1
        s.c[2]['phi^v'] = math.e**s.c[2]['g_mix_x_i v'] # Uses V_v at y_1
    except (ValueError, ZeroDivisionError, KeyError, 
            numpy.linalg.linalg.LinAlgError):
        s.s['Math Error'] = True
        print 'WARNING: Math error in fugacity_error. Failure to calculate '+\
        'g_mix_x^R v or phi_i^v for g_l^R, setting to zero'
        #print s.m['g_mix_x_i v']
        s.c[1]['phi^v'] = 1e-5
        s.c[2]['phi^v'] = 1e-5
        try: # TEST DELETE
            print s.m['g_mix_x_i l']
        except (KeyError):
            pass

    #%%  Find Summ of all K_i x_i yi = Ki xi / SUM Ki xi  
    # Currently can be used to return s.c[#]['y'] 
    s.c[1]['K'] = s.c[1]['phi^l'] / s.c[1]['phi^v']
    s.c[2]['K'] = s.c[2]['phi^l'] / s.c[2]['phi^v']
      # Summ of all K_i x_i
    s.m['Sigma K'] =  x_1 * s.c[1]['K'] + (1 - x_1)  * s.c[2]['K'] 
    
    return s, p
    
    
def y_i_error(y, s, p, P, x_1):
    """
    TO DO if Bubble_error not converging
    This function returns the error for SIGMA K_i x_i to optimize for {y_i} for
    a specified P and {x_i}
    
    Significant improvement to Bubble_error can be gained if bounded internal
    optimzation on {y_i} is used.
    """
    Py = [P, y] # Note P should be constant for optimization
    Py[1] = y # Set y in Py vector to new y
    
    s, p = Sigma_K(Py, s, p, x_1)
    return y - ( x_1 * s.c[1]['K'] / s.m['Sigma K'] )


def Bubble_error(Py, s, p, x_1, returns = 'error'):
    """ 
    Used in Py_VdW_Multicomp
    
    Returns the error in the goal function:
    SIGMA K_i x_i - 1 == 0
    
    To optimize for P and {yi}

    Input 
    
    Returns = error in root finder / optimizer
    if returns = 'y', return "s" to get {y_i} points as { s.c[1]['y'] }
    """     
    s, p = Sigma_K(Py, s, p, x_1)
    
    #%%  Find yi = Ki xi / SUM Ki xi  
    
    #yi = scipy.optimize.brentq(y_i_error, 0.0, 1.0, args=(s, p, Py[0], x_1)) 
    # note, brent works for scalar funcs only, will need new root solver
   # Py[1] = yi
    
    #%% Find error from SUM Ki xi  - 1.0  == 0
    if returns == 'error':
        s, p = Sigma_K(Py, s, p, x_1) # Solve with correct yi
        return [s.m['Sigma K'] - 1.0 ,
                Py[1] - (x_1 * s.c[1]['K'] / s.m['Sigma K'] )]
                
    if returns == 'y': # not really used
        s.c[1]['y'] = (x_1 * s.c[1]['K']) / s.m['Sigma K']
        s.c[2]['y'] = (x_1 * s.c[2]['K']) / s.m['Sigma K']
        return s

#%% TEST PLOT RANGE
def fug_err_plot(s,p, i):
    from numpy import linspace
    x_1 = p.m['x1'][i]
    y_1 = p.m['y1'][i]
    Pr = linspace(3000, 10e7, 2000)
    Pr = linspace(1e6, 5e6, 2000)
    Pr = linspace(1e6, 20e6, 2000)
    Pr = linspace(1e5, 20e6, 2000)
    #Pr = linspace(-1e3, 20e6, 2000)
    Pr = linspace(9e4, 8e6, 2000)
    errstore = []
    for P in Pr:
        #errstore.append(fugacity_error(P, s, p, x_1 , y_1 ))
        errstore.append(\
         Bubble_error([P, y_1], s, p, x_1, returns = 'error')[0])
    from matplotlib import pyplot as plot
    plot.figure()
    #plot.plot(Pr, errstore, 'o--r')
    plot.plot(Pr, errstore, 'r')
    plot.plot([0,max(Pr)], [0,0], 'b--')
    plot.xlabel('P', fontsize=14)
    plot.ylabel(r"$\epsilon$    ", fontsize=14, rotation='horizontal')
    #plot.legend()
    return

#fug_err_plot(s,p, 0)
#%% Mulptiphase Pressure
def Py_VdW_Multicomp(s, p, P_guess=101.3e3, y_1_guess=0.0, T=None, x_1=None, \
                    update_pure=None):
    """
    Returns the Pressure and {y_i} of a (binary) phase muxture at s.m['T'] and
    s.c[1]['x'].
    
    y returns as s.c[]['y']
    
    Initial P_guess=101.3 unless otherwise specified.    
    Initial y_1_guess = x_1 unless otherwise specified.    
    Specify T=None, x1=None.
    """
    from scipy.optimize import fsolve
    #%% Save x1 and y1 if needed
    #if x_1 is None and 'x1' not in s.c[1]: # Save x
    if x_1 is None: # Vals giveN
        if 'x' in s.c[1] and 'x' in s.c[2]:
            x_1, x_2 = s.c[1]['x'], s.c[2]['x'] 
        else:
            raise IOError('Error in P_VdW_Multicomp: No x_1 specified')
    #if y_1 is None and 'y1' not in s.c[1]:
#    if y_1 is None:
#        if 'y' in s.c[1] and 'y' in s.c[2]:
#            y_1, y_2 = s.c[1]['y'], s.c[2]['y'] 
#        else:
#            raise IOError('Error in P_VdW_Multicomp: No y_1 specified')
    #%% Update Temperature if needed
    if T is not None:
        s.m['T'], s.s['T'], s.c[1]['T'], s.c[2]['T'] = (T,)*4
    elif 'T' not in s.m or 'T' not in s.c[1] or 'T' not in s.c[2]:
        raise IOError('Error in P_VdW_Multicomp: No Temperature specified')
    #%% Update Pressure and y1
    s.m['P'], P = (P_guess,)*2  
    y1, s.c[1]['y'], s.c[2]['y']  = y_1_guess, y_1_guess, 1 - y_1_guess
    #% Py vector used in solver:
    #    P   {y}
    Py = [P, y1] 
    # Find optimized Pressure from x_i * phi_i^l(P) - y_i * phi_i^v(P) == 0
    ans = fsolve(Bubble_error,Py,args=(s,p,x_1,'error'),xtol=1e-8)
    #%% TEST CONVERGENCE:
#    testPerr = Bubble_error(ans, s, p, x_1)
#    if abs(testPerr[0]) > 1e-5:
#        print 'WARNING: Poor convergence in Py_VdW_Multicomp detected, '\
#              'x_i * phi_i^l(P) - y_i * phi_i^v(P) == {}.'.format(testPerr) +\
#              ' Attempting lower intitial P guess at 0.9*P.m[\'P\']'
#        #ans = fsolve(fugacity_error,P-0.9*P,args=(s,p,x_1,y_1),xtol=1e-8)
#        if abs(testPerr[0]) > 1e-30:
#              print 'WARNING: Poor convergence in Py_VdW_Multicomp detected, '\
#              'x_i * phi_i^l(P) - y_i * phi_i^v(P) == {}'.format(testPerr) +\
#              ''
    return ans 
    
#%% Optimize parameters via pressure minimization.
def Py_error(p_set, s, p):
    """
    Returns the error function at specified paramter set p_set.
    
    TO DO: ADD PUNISHMENT FOR MATH ERRORS (Every for loop?)
    
    (Previously inside optim_mixture_parameter(s, p) )
    
    For 'VdW standard' model: p_set is scalar = p.m['k12'] 
    For 'DWMP' model: p_set is a vector [k12, k21, r, s]
    """
    import numpy
    if p.m['Model'] == 'VdW standard':
#        print p_set
        p.m['k12'] = p_set[0]
        p.m['k21'] = p.m['k12']
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
            P = Py_VdW_Multicomp(s, p, 
                                 P_guess = s.m['P'],
                                 y_1_guess = s.c[1]['y'], 
                                 T = s.m['T'], x_1=s.c[1]['x'], 
                                 update_pure=None)
            
            #print P[0]       
            #print err     
            #print '(p.m[P][i] - P[0])**2 * 1e-12 = {}'.format( 
            #       (p.m['P'][i] - P[0])**2 * 1e-12 )
            #print '2*abs(p.m[y1][i] - P[1])  = {}'.format( 
            #       2*abs(p.m['y1'][i] - P[1]) )
            err += (p.m['P'][i] - P[0])**2 * 1e-12  \
                   + 2*abs(p.m['y1'][i] - P[1]) # 
                  
        except (ValueError, ZeroDivisionError, 
                numpy.linalg.linalg.LinAlgError):

            err += 10.0*err#max(err)
            print 'WARNING: Math Error in P_error, raising goalfunc error.'
        

    
    
    #if p.m['Model'] == 'DWMP':
    #    return [err, 0, 0, 0]
    #print p_set    
    print err
    return err
    

#%% Error func (over all data points)
def optim_param_set(s, p, guess=None, constraints=None):
    """
    Optimizes the parameter set TO DO SET A PROPER OPTIMIZATION ROUTINE
    """
#    from scipy.optimize import differential_evolution
    from scipy.optimize import minimize, brute, differential_evolution

    if p.m['Model'] == 'VdW standard':
        if guess == None:
            guess = (0.01)
        from scipy.optimize import fsolve
#        p.m['k12'] = fsolve(Py_error,guess,args=(s,p)) #Scalar 
#        ans = minimize(Py_error, guess, method='SLSQP', bounds=[(-2,2)],
#                     args=(s, p) )
#        p.m['k12'] = ans.x[0]
        
        #rranges = slice(-1.0,1.0, 0.25)
        rranges = (slice(-0.1, 1.0, 0.001), slice(-0.1, 1.0, 0.001))
        ans = brute(Py_error, ranges= rranges , args=(s,p))
        print 'ANSWER = {}'.format(ans)
        
    if p.m['Model'] == 'DWMP':
        if guess == None:
            guess = (0.5, 0.5, 1, 1)
        bounds = [(-1.5,1.5),       # p.m['k12']
                  (-1.5,1.5),       # p.m['k21']
                  (-25,25),     # p.m['r']
                  (-25,25)      # p.m['s'] 
                  ]       
        bounds = [(-1.0,1.0),       # p.m['k12']
                  (-1.0,1.0),       # p.m['k21']
                  (-1,25),     # p.m['r']
                  (-1,25)      # p.m['s'] 
                  ]           


        #p.m['k12'], p.m['k21'], p.m['r'], p.m['s'] = \
        
        ans = differential_evolution(Py_error, bounds, args=(s, p))
        
        ans
        print ans
        import numpy

        
#        guess = numpy.array(guess)   

#        ans = minimize(Py_error, guess, method='SLSQP', bounds=bounds,
#                     args=(s, p) )#,
                    #constraints=cons
        #test = fsolve(P_error,guess,args=(s,p)) 
        p.m['k12'] = ans.x[0]
        p.m['k21'] = ans.x[1]
        p.m['r'] = ans.x[2]
        p.m['s'] = ans.x[3]
    return p
    
#%% Trim Azeotrope and single phase points    
def trim_azeo_and_pure(s, p):
    """Previously insdie optim_mixture_parameter(s, p) 
        (Not currently in use)"""
    p.m['x1 Pure vapour'], p.m['y1 Pure vapour'], p.m['P Pure vapour'], \
    p.m['T Pure vapour'], p.m['x1 Pure liquid'], p.m['y1 Pure liquid'], \
    p.m['P Pure liquid'], p.m['T Pure liquid'], p.m['x1 Azeo'], p.m['y1 Azeo']\
    , p.m['P Azeo'], p.m['T Azeo'] = [],[],[],[],[],[],[],[],[],[],[],[]

    def trim_conditions(p, i): 
        """Defined here to simplify proceeding loop readablity"""
        if p.m['x1'][i] > 0.9999 and p.m['y1'][i] > 0.9999:
            p.m['x1 Pure vapour'].append(p.m['x1'].pop(i))
            p.m['y1 Pure vapour'].append(p.m['y1'].pop(i))
            p.m['P Pure vapour'].append(p.m['P'].pop(i))
            p.m['T Pure vapour'].append(p.m['T'].pop(i))
        if p.m['x1'][i] < 1e-5 and p.m['y1'][i] < 1e-5:
            p.m['x1 Pure liquid'].append(p.m['x1'].pop(i))
            p.m['y1 Pure liquid'].append(p.m['y1'].pop(i))
            p.m['P Pure liquid'].append(p.m['P'].pop(i))
            p.m['T Pure liquid'].append(p.m['T'].pop(i))
        else:
            p.m['x1 Azeo'].append(p.m['x1'].pop(i))
            p.m['y1 Azeo'].append(p.m['y1'].pop(i))
            p.m['P Azeo'].append(p.m['P'].pop(i))
            p.m['T Azeo'].append(p.m['T'].pop(i))
        return p

    for i in range(len(p.m['P'])):
        # Stop loop when index > len(p.m['P']) (where p.m['P'] is trimed)
        if i > len(p.m['P'])-2:
            if abs(p.m['x1'][i] - p.m['y1'][i]) < 1e-5:
                p = trim_conditions(p, i-1)
                
            break
        # WHEN CONDITION IS MET, Should redo (same i index now shifted back)
        while abs(p.m['x1'][i] - p.m['y1'][i]) < 1e-5:
            p = trim_conditions(p, i)  

    return s, p
    
    
    
    
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
                try:
                    x1_dif = p.m['x1'][i+1] - p.m['x1'][i]
                    y1_dif = p.m['y1'][i+1] - p.m['y1'][i] # used for guess val
                    P_dif  = p.m['P'][i+1] - p.m['P'][i]
                except(IndexError):
                    pass # Use preveious difference values
                
                printit = False
                
            for j in range(int(added_res)):
                Model['x1'].append(p.m['x1'][i] + x1_dif*j/nr )
                #Model['y1'].append(p.m['y1'][i] + y1_dif*j/nr )

                # Find model point
                s.s['P'] = p.m['P'][i] + P_dif*j/nr
                s.m['P'], s.c[1]['P'], s.c[2]['P'] = (s.s['P'],)*3
            
                s.c[1]['x'], s.c[1]['y'] = p.m['x1'][i] + x1_dif*j/nr \
                                         , p.m['y1'][i] + y1_dif*j/nr
                
                if True:#printit:
                    print 'P_data = {}'.format(s.m['P'])
                    print 'x_data = {}'.format(s.c[1]['x'])
                    print 'y_data = {}'.format(s.c[1]['y'])
                    
                Py = Py_VdW_Multicomp(s, p, 
                                      P_guess = s.m['P'],
                                      y_1_guess = s.c[1]['y'], 
                                      T = s.m['T'], 
                                      x_1 = s.c[1]['x'], 
                                      update_pure=None)   
                             
                #P = P_VdW_Multicomp(s, p, P_guess = s.m['P'], T=s.m['T'], \
                #                x_1=s.c[1]['x'], y_1=s.c[1]['y'], \
                #                update_pure=None)
                Model['P'].append(Py[0])
                Model['y1'].append(Py[1])
                print 'P_model = {}'.format(Py[0])
                print 'y_model = {}'.format(Py[1])
                
                    

                    
#                    P = P_VdW_Multicomp(s, p, P_guess = s.m['P'], T=s.m['T'], \
#                                    x_1=s.c[1]['x'], y_1=s.c[1]['y'], \
#                                    update_pure=None)
#                    Model['P'].append(P)
                
#                    print 'Psat (at x1 = {}, y1 = {}) = {}'.format(p.m['x1'][i], \
#                                                               p.m['y1'][i], P)
#                    
#%% ERROR PLOTS 
            #if ind == 0 or ind == 1 or ind == 2  or ind == 3 or\
            if ind == 8 or ind == 9 or ind == 10 or ind == 11: # TEST ROUTINE DELETE
                pass
                print 'Printing error plot for x_1 = {}' \
                .format(s.c[1]['x'])
                #fug_err_plot(s,p, i)
           #fug_err_plot(s,p, i)
            ind += 1
                
#%% Add more model points
# TO DO 
#%%          
    from matplotlib import rc
    from matplotlib import pyplot as plot
    from numpy import array
    if SingleFig:
        plot.figure(9000)
    else:
        plot.figure()
    plot.plot(Data['x1'], Data['P'], 'x--b', label='Liquid data points')
    plot.plot(Data['y1'], Data['P'], 's--r', label='Vapour data points')

    x1pf, y1pf= '-b', '-r'
    if model_plot_points:
        x1pf = 'o-b'
        y1pf ='^-r'  
        
    plot.plot(Model['x1'], Model['P'], x1pf, label='Liquid model')
    plot.plot(Model['y1'], Model['P'], y1pf, label='Vapour model')
    plot.xlabel(r"$z_1$", fontsize=14)
    plot.ylabel("P", fontsize=14, rotation='horizontal')
    plot.title("{}-{} isotherm at {}".format(p.c[1]['name'][0],
                                             p.c[2]['name'][0],
                                             T_plot))
    #plot.legend()
    
    ##%% TEST PURE PRESSURE
    print 'Psat1 (at x1 = 1.0) = {}'.format( \
    VdW.Psat_V_roots(s.c[1],p.c[1],tol=1e-1)['P_sat'])
    print 'Psat2 (at x1 = 0) = {}'.format( \
    VdW.Psat_V_roots(s.c[2],p.c[2],tol=1e-1)['P_sat'])
    return 
    
#%% 
   
def plot_dg_mix(s,p):
    '''TO DO: Update'''
    from matplotlib import rc
    from matplotlib import pyplot as plot
    from numpy import array
    
    def plotprop(x, name, overall, liquid=None, vapour=None):
    	if overall is not None:
    		plot.plot(x, overall, 'g', label='Overall')
    	if liquid is not None:
    		plot.plot(x, liquid, 'b', label='Liquid')
    	if vapour is not None:
    		plot.plot(x, vapour, 'r', label='Vapour')
    	plot.xlabel(r"$x_1$", fontsize=14)
    	plot.ylabel(name, fontsize=14)
    
    valm, idx = max((valm, idx) for (idx, valm) in \
    enumerate(array(s.s['del g_mix l soln']) - array(s.s['del g_mix v soln'])))
    
    plot.rcParams.update(I['Plot options'])
    plot.figure(5)
    plotprop(s.c[1]['x_range'], r"$\Delta$g", 
             s.s['del g_mix soln'], s.s['del g_mix l soln'], 
             s.s['del g_mix v soln'])
             
#    plot.plot([s.s['Data x1'],s.s['Data y1']], [s.s['Data gx1'],
#               s.s['Data gy1']],'k')
    #plot.figure(2)
    #plotprop(s.c[1]['x_range'], 'Volumes', None, stores['V_l_m'],
    #           stores['V_l_m'])
#    plot.figure(50)
#    plotprop(s.c[1]['x_range'], '$\Delta g_{mix}^l - \Delta g_{mix}^v $',  
#             array(s.s['del g_mix l soln']) - array(s.s['del g_mix v soln']) )
#    plot.text(s.c[1]['x_range'][idx], valm+0.000006, "%f" %round(p.m['k12'],6))

def error_func(s,p):
    """Error function meshgrid plot for DWMP"""
    import numpy as np
    from scipy import optimize
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
   
    # b11 and b22 taken as the critical parameters
    s.c[1]['b'], s.c[2]['b'] = p.c[1]['b_c'], p.c[2]['b_c']  
    # Trim raw data points
    #s, p = trim_azeo_and_pure(s, p)

    # Get meshgrid for error func range
    #rrange = np.linspace(-2, 1.5)
    #srange = np.linspace(-2, 1.5)
    rrange = np.linspace(0.1e-12, 0.1e-10)
    srange = np.linspace(0.1e-12, 0.1e-10)
    rg, sg = numpy.meshgrid(rrange, srange)
    zerr = numpy.zeros((50,50))

    # Find error values
    for i in range(rg.shape[0]):
        for j in range(sg.shape[0]):
            # Do stuff to find "z"
            print 'i = {}'.format(i)
            print 'j = {}'.format(j)
            p.m['r'] = rg[i, j]
            p.m['s'] = sg[i, j]
            print 'Progress = {} %'.format(((i+1)/50.0)*100.0)
            #try:
            zer = Py_error([p.m['k12'], p.m['k21'], p.m['r'], p.m['s']],s,p)
            #maximum = max([max(zerr[i]) for i in range(len(zerr))]) 
            #if zer > 1e5*maximum:
            #    zer = 1.5*maximum
            #except:
            #    print 'WARNING: Error excepted in error_func'
            #    zer = 1.6e15 # TO DO FIND HIGHER ESTIMATE THAN MAX zerr
            zerr[i,j] = zer
            
    # Plot Meshgrid        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(rg, sg, zerr, rstride=1, cstride=1,
                           cmap=plt.cm.jet, linewidth=0, antialiased=False)
    
    ax.set_xlabel('r')
    ax.set_ylabel('s')
    ax.set_zlabel('$\epsilon$',rotation='vertical')
    #ax.set_title('Six-hump Camelback function')
    return rg, sg, zerr, surf 

def error_plot_VdW_standard(s,p): #% TEST PLOT DELETE
    from numpy import linspace
    k12r = linspace(0.0,0.2, num=1e3)
    s.m['errstore'] = []
    for k12 in k12r:
        s.m['errstore'].append(P_error(k12,s,p))
        print '{} % complete'.format((float(len(s.m['errstore'])) \
                                     /float(len(k12r)))*100.0)

    #% TEST PLOT DELETE
    from matplotlib import rc
    from matplotlib import pyplot as plot
    from numpy import array
    plot.figure()
    plot.plot(k12r, s.m['errstore'], 'k', label='Error function')
    plot.xlabel(r"$k_{12}$", fontsize=14)
    plot.ylabel(r"$\epsilon$    ", fontsize=14, rotation='horizontal')

    return
    
def g_T_surf(s,p):
    """Surface of the total Gibbs energy relative to comp. 1 pure liquid"""
    import numpy as np
    from scipy import optimize
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
   

    x_range = np.linspace(0.0, 1.0)
    y_range = np.linspace(0.0, 1.0)
    xg, yg = numpy.meshgrid(x_range, y_range)
    gerr = numpy.zeros((50,50))

    # Find error values
    for i in range(xg.shape[0]):
        for j in range(yg.shape[0]):
            # Do stuff to find "z"
            print 'i = {}'.format(i)
            print 'j = {}'.format(j)
#            p.m['r'] = xg[i, j]
#            p.m['s'] = yg[i, j]
            x = xg[i, j]
            y = yg[i, j]
            print 'Progress = {} %'.format(((i+1)/50.0)*100.0)
            #try:
            #zer = Py_error([p.m['k12'], p.m['k21'], p.m['r'], p.m['s']],s,p)
            ger = d_g_mix_T(s, p, x_1=x, y_1=y, update_pure=True, 
                            derivative=None)

            gerr[i,j] = ger
            
    # Plot Meshgrid        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(xg, yg, gerr, rstride=1, cstride=1,
                           cmap=plt.cm.jet, linewidth=0, antialiased=False)
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('$\Delta g$',rotation='vertical')
    #ax.set_title('Six-hump Camelback function')
    return xg, yg, gerr, surf 
    
#%%
if __name__ == '__main__':
    #%% Load Data
    try: # Dectect input, use local script if not defined then import data
        Compounds = I['Compounds'] # Variable to draw data in data_handling.py
        data = data_handling.ImportData()  
        data.load_pure_data(I['Compounds'])
        data.load_VLE(Compounds)
    except(KeyError,NameError): # Define local inputs if I is not found.
        I = inputs() # 
        Compounds = I['Compounds'] # Variable to draw data in data_handling.py
        data = data_handling.ImportData()  
        data.load_pure_data(I['Compounds'])
        data.load_VLE(Compounds)
        #execfile('data_handling.py') # Load data

    #%% Find pure component model parameters if not defined    
    for compound in I['Compounds']:
        if False:# TO DO TRY EXCEPTION HANDLING TO DETECT NEEDED PARAMS
            I['Compound'] = [compound]
            execfile('pure.py')
   
    #%% Initialize binary and single component paramters
    p = MixParameters() 
    p.parameters(data.c[0],I) # Defines p.c[1]
    p.parameters(data.c[1],I) # Defines p.c[2]
    p.mixture_parameters(data.VLE,I)
    p.m['R'] = p.c[1]['R'] # Use a component Universal gas constant

    #%% Initialize state variables
    s = state()
    s.mixed() # Define mix state variable, call using s.m['key']
    # Define three component state variables (use index 1 and 2 for clarity)
    s.pure(), s.pure()  # Call using s.c[1]['key'] and s.c[2]['key']
    
    #%% Plot error function
    if False:
        p.m['k12'] = 0.0#0.124#0.0 # TO DO
        p.m['k21'] = 0.0#0.124#0.0 # TO DO
        rg, sg, zerr, surf  = error_func(s,p)
        
    
    #%% Find mixture model parameters if not defined
    
    #% Test P_VdW_Multicomp
    if True:
        # b11 and b22 taken as the critical parameters
        s.c[1]['b'], s.c[2]['b'] = p.c[1]['b_c'], p.c[2]['b_c']  
    
        p.m['k12'], p.m['k21'], p.m['r'], p.m['s']= (1e-5,)*4 # First estimate
        #ans = P_VdW_Multicomp(s, p, T=283.15, x_1=0.7256, y_1=0.6626)
        #print ans
        
        if False:
            s, p = trim_azeo_and_pure(s, p) # Trim azeos and pure
        
        
        if False: # Optimize for model:
            if p.m['Model'] == 'VdW standard':
                p = optim_param_set(s, p, guess=p.m['k12']) 
                print 'Optimized paramters:'
                print 'k12 = {}'.format(  p.m['k12']    )
            if p.m['Model'] == 'DWMP':
                print 'Console Test'
                p.m['k12'], p.m['k21'], p.m['r'], p.m['s']= 0.1, 0.1, 1, 1
                p = optim_param_set(s, p, guess=[p.m['k12'], p.m['k21'], \
                                                 p.m['r'], p.m['s']]) 
                print 'Optimized paramters:'
                print 'k12 = {}'.format(  p.m['k12']    )
                print 'k21 = {}'.format(  p.m['k21']    )
                print 'r = {}'.format(    p.m['r']      )
                print 's = {}'.format(    p.m['s']      )
                
        #p.m['k12'] = 0.0247         # From Patel et. al.
        
        if False: # Plot a k12 range
            from numpy import linspace
            k12r  = linspace(0,1.0e-3,11)
            for p.m['k12'] in k12r:
                plot_isotherm(s, p, T_plot = 281.15)
                #plot_isotherm(s, p, T_plot = 293.1)
                #plot_isotherm(s, p, T_plot = 223.11)
                
                #plot_isotherm(s, p, T_plot = 263.1)
        
        if False: # Plot all isotherms
#            T_isos = [281.15, 283.15, 287.15, 293.15, 298.06, 298.15, 313.15, \
#                      318.15, 328.2, 333.15, 343.15]
            
            
            # Acetone-Water
            p.m['k12'] = 0.9#1.36758234609
            p.m['k21'] = 0.522182150663
            p.m['r'] = -3.01592078893
            p.m['s']  = -3.55163821423
            
            #p.m['k12'] = 1.36758234609
            #p.m['k21'] = 0.522182150663
            #p.m['r'] = 3.01592078893
            #p.m['s']  = 3.55163821423
            
            T_isos =  p.m['T']   
            from  more_itertools import unique_everseen
            T_isos = list(unique_everseen(T_isos))
            for T_i in T_isos:
                plot_isotherm(s, p, T_plot = T_i, added_res= 50) 

        if False: # Acetone-Water
            
            #p.m['k12'] = 0.0
            #p.m['k21'] = 0.0
            p.m['k12'] = -9.56080361137
            p.m['k21'] = 8.57589067108
            p.m['r'] = -41.8958376865
            p.m['s']  = 102.158239374

            p.m['k12'] = 9.27950399002
            p.m['k21'] = -0.437213205543
            p.m['r'] = -33.5607442337
            p.m['s']  = -177.105587853
            

            p.m['k12'] = 1.36758234609
            p.m['k21'] = 0.522182150663
            p.m['r'] = -3.01592078893
            p.m['s']  = -3.55163821423
            
            p.m['k12'] = 1.0
            p.m['k21'] = -1.0
            p.m['r'] = 6.43313519554
            p.m['s'] = -0.526444090473    
          
            p.m['k12'] = 0.0528180441074
            p.m['k21'] = 0.623384831942
            p.m['r'] = 7.53330786789
            p.m['s']  = 0.107160035705
            
            plot_isotherm(s, p, T_plot = 308.15, added_res= 50)
            plot_isotherm(s, p, T_plot = 333.15, added_res= 50)
            plot_isotherm(s, p, T_plot = 373.15, added_res= 50)
            plot_isotherm(s, p, T_plot = 423.15, added_res= 50)
            plot_isotherm(s, p, T_plot = 523.15, added_res= 50)
            
            p.m['k12'] = 0.0#0.1
            p.m['k21'] = 0.0#0.3
            p.m['r'] = 1
            p.m['s']  = 1
            
        if False: # Plot data from Kont. Graph for # CO2-Ethan
            
            #p.m['k12'] = 0.0
            #plot_isotherm(s, p, T_plot = 263.1, SingleFig=True)
            p.m['r'], p.m['s'] = 1.0, 1.0
            p.m['k12'] = 0.124
            p.m['k21'] = p.m['k12']
            plot_isotherm(s, p, T_plot = 263.1, SingleFig=True)
            p.m['k12'] =  0.0
            p.m['k21'] = p.m['k12']
            plot_isotherm(s, p, T_plot = 263.1, added_res= 50, SingleFig=True)
           

            #Benze-cyclo
            #p.m['k12'] = 0.03059005  # p.m['k12'] =5.28317633e-06 
            #plot_isotherm(s, p, T_plot = 281.15, SingleFig=True) 
            
        if False:                
            error_plot_VdW_standard(s,p)
            p.m['k12'] = 0.03922446 # benzene-cyclohexane value
        
            
        
        if False: # DWMP ERROR FUNC
            error_func(s,p)
            
        if True: # Test Gibbs curves
            p.m['r'], p.m['s'] = 1.0, 1.0
            p.m['k12'] = 0.124
            p.m['k21'] = p.m['k12']

            s.s['T'] = 263.1
            s.m['T'], s.c[1]['T'], s.c[2]['T'] = (s.s['T'],)*3
            s.s['P'] = 24e5
            s.m['P'], s.c[1]['P'], s.c[2]['P'] = (s.s['P'],)*3
            x_r = 500
            g_range(s, p, x_r)
            plot_dg_mix(s,p)
            
            #
            plot_isotherm(s, p, T_plot = 263.1, SingleFig=True)
#            
            g_T_surf(s,p)
            
            
            
    if False:
    #if I['Mixture model'] == 'VdW standard':
        #if p.m['k12'] == '':
        if False: # TESTING p.m['k12']; DELETE
            p.m['k12'] = 0.0135
            
            
            #optim_mixture_parameter(s, p)
            #s, p = optim_mixture_parameter(s, p) # Find k1 for VLE data set.

    #%% Find phase equilibrium at specified T, P point
    if False:
    #if I['Phase split']:
        #I['T'] = 313.15
        #I['P'] = 27245.75912
        I['T'] = p.m['T'][33] #263.1
        I['P'] = p.m['P'][33] #2660800.0
        if I['T'] and I['P']:
            #p.m['k12'] = 0.0135 # DEBUGGING
            #p.m['k12'] = 0.03922446 
        
            # Upate P, T state
            s.s['P'], s.s['T'] = I['P'], I['T'] 
            # b11 and b22 taken as the critical parameters
            s.c[1]['b'], s.c[2]['b'] = p.c[1]['b_c'], p.c[2]['b_c']  
            s.m['P'], s.c[1]['P'], s.c[2]['P'] = (s.s['P'],)*3
            s.m['T'], s.c[1]['T'], s.c[2]['T'] = (s.s['T'],)*3
            p.m['Plot binary'] = I['Plot binary']
            #s, p = phase_split(s, p, est=[0.556,0.55])
           # s, p = phase_split(s, p, est=[0.3551,0.3849])
            #s, p = phase_split(s, p, est=[0.265,0.3849])
            #s, p = phase_split(s, p, est=[0.71,0.85])
            s, p = phase_split(s, p, est=[ p.m['x1'][33], p.m['P'][33] ])
            #print s.m['ans']
            #print s.m['ans'][1] - s.m['ans'][0]
            print s.m['x1']
            print s.m['y1']

          
            if p.m['Plot binary']: 
                pass
                #g_range(s, p, x_r=100)
                #plot_dg_mix(s,p)
        else:
            raise IOError('No Temperature or Pressure point specified')
       
         # a_c[0] = a_c11 | a_c[1] = a_c22    
        
          #   'T'           : 281.15,  
          #   'P'           : 5663.534211
    
        # datapoint = array([0.2027,0.2567]) # @ P = 6.278150329 T = 281.15 
        # datapoint = array([0.9618,0.9297]) # @ P = 5.663534211  T = 281.15    
        # datapoint = array([0.4535,0.4605]) # @ P = 6492.799322 T = 281.15    
             
        #%% TESTS AND DEBUGGING
        #s, p = phase_split(s,p,x_r=200)
        from numpy import linspace
        if False:
            for p.m['k12'] in linspace(0.03125,0.0325,5):  
            # TESTS linspace(0.05+ 1e-30,0.07,30)
                print p.m['k12'] # TEST; DELETE
                s, p = phase_split(s,p,x_r=200)
            
        if False:
            for p.m['k12'] in linspace(0.011,0.033,5):
                print p.m['k12'] # TEST; DELETE
                s, p = phase_split(s,p,x_r=200)
                

    
    
                
    #%% Find phase equilibrium at specified T, P point

    
    #%% Plotting if True
    
    #%% !!! NB NB NB !!! Remember to revert indiced from 1,2 back to [0] [1] 
    # when saving data from s or p components

    if False: # Test plots for co2-ethane
        p.m['r'], p.m['s'] = 1.0, 1.0
        p.m['k12'] = 0.124
        p.m['k21'] = p.m['k12']
        I['T'] = p.m['T'][33] #263.1
        I['P'] = p.m['P'][33] #2660800.0
        s.s['P'], s.s['T'] = p.m['T'][33] , p.m['P'][33]
        s.c[1]['b'], s.c[2]['b'] = p.c[1]['b_c'], p.c[2]['b_c']  
        s.m['P'], s.c[1]['P'], s.c[2]['P'] = (s.s['P'],)*3
        s.m['T'], s.c[1]['T'], s.c[2]['T'] = (s.s['T'],)*3
        p.m['Plot binary'] = I['Plot binary']
        s.c[1]['x'], s.c[2]['x']  = p.m['x1'][33], (1 - p.m['x1'][33])
        s.c[1]['y'], s.c[2]['y']  = p.m['y1'][33], (1 - p.m['y1'][33])
        #s, p = phase_split(s, p, est=[ p.m['x1'][33], p.m['P'][33] ])
        #plot_dg_mix(s,p)
        g_range(s, p, x_r=100)

    
    
    

    
      