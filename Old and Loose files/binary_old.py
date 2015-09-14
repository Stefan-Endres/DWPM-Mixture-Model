#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
#%% Imports
from __future__ import division
from scipy.interpolate import interp1d
import data_handling, Van_der_Waals, numpy
VdW = Van_der_Waals.VdW()

try: #DEBUGGING; DELETE
    del s
    del p
    del Iz
except NameError:
    pass

#%% Inputs (will be called if no input container "I" is defined before exec)
def inputs():
    I = {# Model inputs
         'Compounds'    : ['benzene','cyclohexane'], # Compound to simulate.
         'Mixture model': 'VdW standard', #  'VdW standard' 'DWMP'
         'Model'       : 'Adachi-Lu',   # Model used in the simulation, 
                                     # options:
                                      # 'Soave'                   
                                      # 'Adachi-Lu'   
         
         # Optional inputs
         'T'           : 281.15,  
         'P'           : 6278.150329 ,   
         'Save results': True,
         'Plot binary' : True,
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
        M = {'k12' : Data['k12'][0],
             'k1'  : Data['k1'][0],
             'k2'  : Data['k2'][0],
             'r'   : Data['r'][0],
             's'   : Data['s'][0],
             'T'   : Data['T (K)'],
             'P'   : Data['P (Pa)'],
             'x1'  : Data['x1'],
             'y1'  : Data['y1'],
             'Model': I['Mixture model']
             }
             
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
        raise IOError('DWMP not implemented yet')
    
def b_mix(s,p):
    return s.c[1]['b']*s.c[1]['x'] + s.c[2]['b']*s.c[2]['x']
    
def a12(s,p): # TO DO TRY MODEL  p.m['k12'] ==  k1*x1 + k2*x2 
    from math import sqrt
    return (1 - p.m['k12']) * sqrt(s.c[1]['a']*s.c[2]['a'])   
    
def b12(s,p): # NOTE: Not currently in use for VdW Standard or DWMP! 
    from math import sqrt
    return sqrt(s.c[1]['b']*s.c[2]['b']) 
#%% g - g_ref (add g of volume phase)
def g_R(s,p):
    ''' NOTE: Check VOLUMES'''
    from math import log
    return s['P']*s['V_v']/(p['R']*s['T']) - 1.0 + log(s['P']/(p['R']*s['T'])) \
    - log(s['V_v'] - s['b']) - s['a']/(p['R']*s['T']*s['V_v'])
    
def g_IG(s,sc,p): # TO DO!!!!
    """ex. s.c[1]['g_R'] = g_R(s.c[1],p.c[0])"""
    #from math import log
    H_IG = 2
    S_IG = 02
    G_IG = H_IG - s.s['T'] * S_IG  # G = H - TS
    return 129655.0/(p.m['R']*s.s['T'])
    
def g(s,sc,p):
    """ex. s.c[1]['g_R'] = g_R(s.c[1],p.c[0])"""
    return 0.0#g_R(s,p) + g_IG(s,p) # g^R_i + g^IG_i
#%% dg_mix
def del_g_mix(s,p):
    """ Returns the change in gibbs energy Gmix/RT at specified composition
    x1 = x[0], Pressure and Temperature"""
    from math import log   
    #%% Find params at current composition
    s.c[1]['a'] = VdW.a_T(s.c[1],p.c[1])['a']
    s.c[2]['a'] = VdW.a_T(s.c[2],p.c[2])['a']
    try:
        s.m['a12'] = a12(s,p)
        s.m['a'] = a_mix(s,p)
        s.m['b12'] = b12(s,p) 
        s.m['b'] = b_mix(s,p)
    except (ValueError, ZeroDivisionError):
        s.s['Math Error'] = True
        print 'WARNING: Math Domain error at %s Pa %s K' %(s.s['P'],s.s['T'])
        s.m['a12'], s.m['a'], s.m['b12'], s.m['b'] = (0.0,)*4
    
    #%% Find Volume Roots ('V_v' and 'V_l') at P, T, a, b for both components.
    s.c[1], s.c[2] = VdW.V_root(s.c[1], p.c[1]), VdW.V_root(s.c[2], p.c[2])
    s.m = VdW.V_root(s.m, p.m) # 'V_v' and 'V_l' mixture volumes
    #%% Change in Volume used in main equation
    s.m['del V_l'] = s.m['V_l'] - s.c[1]['x']*s.c[1]['V_l'] \
                                - s.c[2]['x']*s.c[2]['V_l']
    s.m['del V_v'] = s.m['V_v'] - s.c[1]['x']*s.c[1]['V_v'] \
                                - s.c[2]['x']*s.c[2]['V_v']               
    #%% Find Change in gibbs energies of mixing
    try:
        s.m['del g_mix^R l'] = s.s['P']*s.m['del V_l'] / (p.m['R']*s.s['T']) \
                            \
                            + s.c[1]['x']*log( (s.c[1]['V_l'] - s.c[1]['b']) \
                                              /(s.m['V_l']    - s.m['b']   ) )\
                            + s.c[2]['x']*log( (s.c[2]['V_l'] - s.c[2]['b']) \
                                              /(s.m['V_l']    - s.m['b']   ) )\
                            \
                            + s.c[1]['x'] * s.c[1]['a'] \
                             /(p.m['R']*s.s['T'] * s.c[1]['V_l']) \
                            + s.c[2]['x'] * s.c[2]['a'] \
                             /(p.m['R']*s.s['T']*s.c[2]['V_l']) \
                            \
                            - s.m['a'] / (p.m['R']*s.s['T']*s.m['V_l'])

               
        s.m['del g_mix^R v'] = s.s['P']*s.m['del V_v'] / (p.m['R']*s.s['T']) \
                            \
                            + s.c[1]['x']*log( (s.c[1]['V_v'] - s.c[1]['b']) \
                                              /(s.m['V_v']    - s.m['b']   ) )\
                            + s.c[2]['x']*log( (s.c[2]['V_v'] - s.c[2]['b']) \
                                              /(s.m['V_v']    - s.m['b']   ) )\
                            \
                            + s.c[1]['x'] * s.c[1]['a'] \
                             /(p.m['R']*s.s['T'] * s.c[1]['V_v']) \
                            + s.c[2]['x'] * s.c[2]['a'] \
                             /(p.m['R']*s.s['T']*s.c[2]['V_v']) \
                            \
                            - s.m['a'] / (p.m['R']*s.s['T']*s.m['V_v'])
        
        s.m['del g v'] =  s.c[1]['x'] * g(s,s.c[1],p.c[1])\
                        + s.c[2]['x'] * g(s,s.c[2],p.c[2]) 
        
        
        s.m['del g_mix l'] = s.m['del g_mix^R l'] # Reference change in liquid
        s.m['del g_mix v'] = s.m['del g_mix^R v'] + s.m['del g v']
        s.m['del g_mix']   = min(s.m['del g_mix l'],s.m['del g_mix v'])
        s.s['Math Error'] = False
    except (ValueError, ZeroDivisionError):
        s.s['Math Error'] = True
        print 'WARNING: Math Domain error.'
        s.m['del g_mix l'], s.m['del g_mix v'] = 0.0, 0.0

    return s

def phase_split(s,p,x_r=100):
    """return s.m['x1'] s.m['y1']"""
    #%% Imports
    import itertools
    from numpy import linspace, less, array
    from scipy.signal import argrelextrema
    #%%
    s.m['P'], s.c[1]['P'], s.c[2]['P'] = (s.s['P'],)*3
    s.m['T'], s.c[1]['T'], s.c[2]['T'] = (s.s['T'],)*3
    
    s.c[1]['x_range'], s.c[2]['x_range'] = linspace(0,1,x_r), linspace(1,0,x_r)
    def g_range(s,p): 
        """TO DO OPTIMIZE WITH MAP FUNCTION/ COMPILE IN C / def append etc."""
        #% Initialize
        s.s['del g_mix l soln'],s.s['del g_mix v soln'] = [],[]
        s.s['del g_mix soln'], s.s['Math Error soln'] = [],[]
        #% Solve for x_range
        for i in range(len(s.c[1]['x_range']-1)):
            s.c[1]['x'],s.c[2]['x'] = s.c[1]['x_range'][i],s.c[2]['x_range'][i]
            s = del_g_mix(s,p)   #.m['del g_mix l']
            s.s['del g_mix l soln'].append(s.m['del g_mix l'])
            s.s['del g_mix v soln'].append(s.m['del g_mix v'])
            s.s['del g_mix soln'].append(s.m['del g_mix'])
            s.s['Math Error soln'].append(s.s['Math Error'])
                
        return s
        
    s = g_range(s,p) # Get points to detect phase seperation
    
    s.s['g_min'] = argrelextrema(array(s.s['del g_mix soln']), less) 
                    # (array([52]),)
    
    
    #%% Data point tests
    s.s['Data x1'] = 0.2027
    s.s['Data y1'] = 0.2567
    s.c[1]['x'],s.c[2]['x'] = s.s['Data x1'],(1.0 -s.s['Data x1'])
    s.s['Data gx1']  = del_g_mix(s,p).m['del g_mix']
    
    s.c[1]['x'],s.c[2]['x'] = s.s['Data y1'],(1.0 -s.s['Data y1'])
    s.s['Data gy1']  = del_g_mix(s,p).m['del g_mix']
    #%%
    
    if p.m['Plot binary']: 
        plot_dg_mix(s,p)
    
    


    return s, p

def optim_mixture_parameter():
    # Error func (over all data points)
    pass

#%% Plot functions
def plot_dg_mix(s,p):
    ''''''
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
    plot.figure(1)
    plotprop(s.c[1]['x_range'], r"$\Delta$g$_{mix}$", 
             s.s['del g_mix soln'], s.s['del g_mix l soln'], 
             s.s['del g_mix v soln'])
             
    plot.plot([s.s['Data x1'],s.s['Data y1']], [s.s['Data gx1'],
               s.s['Data gy1']],'k')
    #plot.figure(2)
    #plotprop(s.c[1]['x_range'], 'Volumes', None, stores['V_l_m'],
    #           stores['V_l_m'])
    plot.figure(50)
    plotprop(s.c[1]['x_range'], '$\Delta g_{mix}^l - \Delta g_{mix}^v $',  
             array(s.s['del g_mix l soln']) - array(s.s['del g_mix v soln']) )
    plot.text(s.c[1]['x_range'][idx], valm+0.000006, "%f" %round(p.m['k12'],6))

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
    p.parameters(data.c[0],I)
    p.parameters(data.c[1],I)
    p.mixture_parameters(data.VLE,I)
    p.m['R'] = p.c[1]['R'] # Use a component Universal gas constant
    
    #p.m['k12'] = -0.002403649 #DEBUGGING; DELETE
    
    #%% Initialize state variables
    s = state()
    s.mixed() # Define mix state variable, call using s.m['key']
    # Define three component state variables (use index 1 and 2 for clarity)
    s.pure(), s.pure(), s.pure()  # Call using s.c[1]['key'] and s.c[2]['key']

    #%% TEST PHASE SPLIT
    s.s['P'], s.s['T'] = I['P'], I['T'] # Upate P, T state
    # b11 and b22 taken as the critical parameters
    s.c[1]['b'], s.c[2]['b'] = p.c[1]['b_c'], p.c[2]['b_c']  
     # a_c[0] = a_c11 | a_c[1] = a_c22

    p.m['Plot binary'] = I['Plot binary']
    
    
      #       'T'           : 281.15,  
      #   'P'           : 5663.534211 , 5663.534211 

    # datapoint = array([0.2027,0.2567]) # @ P = 6.278150329 T = 281.15 
    # [-0.0024400000201203567, -0.0024000000201203566]
    from numpy import linspace
    for p.m['k12'] in linspace(-1,1.0,30): # k12 = [-0.5, 0.5] #  (0.5,0.51
        #p.m['k12'] = p.m['k12']/1200.0
        print p.m['k12']
        s, p = phase_split(s,p,x_r=200)
    #%%
    
    
    
    
    
    
    
    
    
    
    # datapoint = array([0.9618,0.9297]) # @ P =5.663534211  T = 281.15    
    
    #%% Find mixture model parameters if not defined
    
       
    
    #%% Find phase equilibrium at specified T, P point
    if I['T'] and I['P']:
        pass
    
    #%% Plotting if True
    
    #%% !!! NB NB NB !!! Remember to revert indiced from 1,2 back to [0] [1] 
    # when saving data from s or p components

    

    #%% TESTS
    '''
    p.m['k1'] = 0.001
    
    
    s = state()
    s.mixed() # Define mix state variable
    
    #s.pure() 
    s.pure(), s.pure(), s.pure() 
    
    
    s.m['test'] = 3
    s.m['test'] 
    
    
    
    s.s['P'] = I['P']#p.m['P'][6]
    s.s['T'] = I['T']#p.m['T'][6]
    s.s['V'] = 5e-3
    s.s['b'] = 3e-5
    s.s['a'] = 13e-1
    #s.s['g_R'] = g_R(s.s,p.c[1])
    
    s.m['P'] = I['P']#p.m['P'][6]
    s.m['T'] = I['T']#p.m['T'][6]
    s.m['V'] = 5e-3
    s.m['b'] = 3e-5
    s.m['a'] = 13e-1
    #s.m['g_R'] = g_R(s.s,p.c[1])
    
    s.c[1]['x'],s.c[2]['x'] = 0.5,0.5
    s.c[1]['P'] = I['P']#p.m['P'][6]
    s.c[1]['T'] = I['T']#p.m['T'][6]
    s.c[1]['V'] = 5e-3
    s.c[1]['b'] = 3e-5
    s.c[1]['a'] = 13e-1
    #s.c[1]['g_R'] = g_R(s.c[1],p.c[1])
    
    
    s.c[2]['P'] = I['P']#p.m['P'][6]
    s.c[2]['T'] = I['T']#p.m['T'][6]
    s.c[2]['V'] = 5e-3
    s.c[2]['b'] = 3e-5
    s.c[2]['a'] = 13e-1
    #s.c[2]['g_R'] = g_R(s.c[2],p.c[2])
    
    
    
    s.m['a12'] = a12(s,p)

    s, p = del_g_mix(s,p)
    '''
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    