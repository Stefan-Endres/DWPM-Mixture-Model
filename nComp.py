#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to simulate phase equilibria of multicomponent systems.
"""

# %% Imports
from __future__ import division
from scipy.interpolate import interp1d
import data_handling, Van_der_Waals, numpy
VdW = Van_der_Waals.VdW()

try: #DEBUGGING; DELETE
    del s
    del p
    del I
except NameError:
    pass

# %% Inputs (will be called if no input container "I" is defined before exec)
def inputs():
    I = {# Model inputs
          # Compounds to simulate.
         'Compounds'    : ['acetone', 'water', 'phenol'], 
         'Mixture model': 'DWMP', #  'VdW standard' 'DWMP'
         'Model'        : 'Adachi-Lu',   # Model used in the simulation, 
                                     # options:
                                      # 'Soave'                   
                                      # 'Adachi-Lu' 
         
         'Valid phases' : ['x', 'y'], # List of valid phases in equilibrium
                                       # ex. for VLE use ['x', 'y']
         
         # Optional inputs
         'T'           : 281.15,  # 281.15
         'P'           : 6278.150329,   # 
         'Phase split' : True, # Find phase split at specified T, P.
         'Save results': True,
         'Fig. number' : 2,
         'Plot pure'   : False,       
         'Plot options': {'text.usetex' : True, # Options for all plots
                          'font.size' : 11,            
                          'font.family' : 'lmodern',
                          'text.latex.unicode': True

                          },
         }

    return I
    
# %% Define paramter class
class MixParameters:
    """Store mixture and pure parameters in the same class."""
    def __init__(self):
        #self.name = name
        self.c = []    # creates a new empty list for components 
        self.c.append('nan') #  Define an empty set in index 0 

    def mixture_parameters(self,Data,I):
        """Mixture model parameters"""
        M = {'T'      : Data['T (K)'], # Temperature Pressure data
             'P'      : Data['P (Pa)'],
             'n'      : len(I['Compounds']),
             'phases' : len(I['Valid phases']),
             'Model'  : I['Mixture model']
             }
        
        # Define phases
        for i in range(len(I['Valid phases'])):
            M[I['Valid phases'][i]] = ['nan'] 
            
        # Define equilibria for each component in phase    
        for j in range(1, M['n'] + 1): 
            for i in range(len(I['Valid phases'])):
                exec 'M[ I[\'Valid phases\'][i] ]' + \
                '.append(Data[ I[\'Valid phases\'][i] + \'{}\'])'.format(j)

                # NOTE: This routine will change component strings the for the 
                # equilibrium of each phase into a list simple list for each 
                # phase ex. Data['x1'], Data['x2'] becomes M['x'][1], M['x'][2] 

               
        # Define model paramters
        # Empty lists for model interaction paramters
        M['k'] = []
        [M['k'].append(['nan']) for j in range(M['n'] + 1)]
        
#        if I['Mixture model'] == 'VdW standard':
#            for i in len(I):
#                M['k12'] = Data['k1_VdW']
#                M['k21'] = Data['k2_VdW']
            
        if I['Mixture model'] == 'DWMP':
            # Find the interaction paramters between and put them into 
            # component lists (ex. Data['k12']  --> M['k'][1][2])
        
            for j in range(1, M['n'] + 1): 
                for i in range(1, M['n'] + 1): 
                    # Define empty list
                    exec 'M[\'k\'][{J}].append(\'nan\')'.format(J = j, I = i)
                    if i != j: # Define interaction paramter
                        exec 'M[\'k\'][{J}][{I}] = '.format(J = j, I = i) + \
                        'Data[\'k{J}{I}\']'.format(J = j, I = i)
            
            M['r'] =  Data['r'],
            M['s'] = Data['s']

        for key, value in M.iteritems(): # Filter out '' values
            if not value.__class__ == float:
              #if key != 'x' and key != 'y' and key != 'k' and key != 'phases':
                try:
                    M[key] = filter(lambda a: a != '', value)
                except(TypeError):
                    pass
        
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

# %% Define state variable class
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
        
        
# %%
if __name__ == '__main__':

    # %% Load Data
    try: # Dectect input, use local script if not defined then import data
        data = data_handling.ImportData()  
        data.load_pure_data(I['Compounds'])
        data.load_E(I['Compounds'])
    except(KeyError,NameError): # Define local inputs if I is not found.
        I = inputs() # 
        data = data_handling.ImportData()  
        data.load_pure_data(I['Compounds'])
        data.load_E(I['Compounds'])
        #execfile('data_handling.py') # Load data
        

    # %% Find pure component model parameters if not defined    
    for compound in I['Compounds']:
        if False:# TO DO TRY EXCEPTION HANDLING TO DETECT NEEDED PARAMS
            I['Compound'] = [compound]
            execfile('pure.py')
 
   # %% Initialize binary and single component paramters
    p = MixParameters() 
    p.mixture_parameters(data.VLE,I)
    #p.m['n'] = len(I['Compounds']) # Define system size

    for i in range(p.m['n']): # Set params for all compounds
        p.parameters(data.c[i],I) # Defines p.c[i]

    p.m['R'] = p.c[1]['R'] # Use a component Universal gas constant

    
    # %% Initialize state variables
#    s = state()
#    s.mixed() # Define mix state variable, call using s.m['key']
#    # Define three component state variables (use index 1 and 2 for clarity)
#    s.pure(), s.pure()  # Call using s.c[1]['key'] and s.c[2]['key']