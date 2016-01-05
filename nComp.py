#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to simulate phase equilibria of multicomponent systems.
"""

# %% Imports
from __future__ import division
import data_handling
import Van_der_Waals
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
         'Mixture model': 'DWPM',  # Removed 'VdW standard', set r = s = 1
         'Model'        : 'Adachi-Lu',  # Activity coefficient Model used in 
                                        # the simulation, 
                                         # options:
                                          # 'Adachi-Lu' 
                                          # 'Soave'

         'Valid phases' : ['x'],#, 'y'], # List of valid phases in equilibrium
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


# %% Define paramter class
class MixParameters:
    """
    Store mixture and pure parameters in the same class.

    Parameters
    ----------
    Data : Dictionary containing data loaded from the stored .csv file.

    I : Input dictionary, must contain key 'Model' defining the a parameter
        dependancy model.
    """
    def __init__(self):
        self.c = []  # creates a new empty list for components
        self.c.append('nan')  # Define an empty set in index 0

    def mixture_parameters(self, Data, I):
        """Mixture model parameters"""
        M = {'T'      : Data['T (K)'], # Temperature Pressure data
             'P'      : Data['P (Pa)'],
             'n'      : len(I['Compounds']),
             'phases' : len(I['Valid phases']),
             'Model'  : I['Mixture model'],
             'Valid phases' : I['Valid phases']
             }

        # Define phases
        for i in range(len(I['Valid phases'])):
            M[I['Valid phases'][i]] = ['nan']

        # Define equilibria for each component in phase
        for j in range(1, M['n'] + 1):
            for i in range(len(I['Valid phases'])):
                M[I['Valid phases'][i]].append(Data[I['Valid phases'][i] +
                                               '{}'.format(j)])

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

        if I['Mixture model'] == 'DWPM':
            # Find the interaction paramters between and put them into
            # component lists (ex. Data['k12']  --> M['k'][1][2])

            for j in range(1, M['n'] + 1):
                for i in range(1, M['n'] + 1):
                    # Define empty list
                    M['k'][j].append('nan')
                    if i != j:  # Define interaction paramter
                        M['k'][j][i] = Data['k{J}{I}'.format(J=j, I=i)][0]

            M['r'] = Data['r']
            M['s'] = Data['s']

        for key, value in M.iteritems():  # Filter out '' values
            if not value.__class__ == float:
            #  if key != 'x' and key != 'y' and key != 'k' and key != 'phases':
                try:
                    M[key] = filter(lambda a: a != '', value)
                except(TypeError):
                    pass

        self.m = M

    def parameters(self, Data, I):
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

        if I['Model'] == 'Adachi-Lu':  # Find model params if not defined
            p['m'] = Data['m (Adachi-Lu)'][0]
        elif I['Model'] == 'Soave':
            p['m'] = Data['m (Soave)'][0]

        if p['a_c'] == '' or p['b_c'] == '':
            p['b_c'] = p['R']*p['T_c']/(8*p['P_c'])
            p['a_c'] = 27*(p['R']**2)*(p['T_c']**2)/(64.0*p['P_c'])
        else:
            pass

        for key, value in p.iteritems():  # Filter out '' values
            if not value.__class__ == float:
                p[key] = filter(lambda a: a != '', value)

        self.c.append(p)


# %% Define state variable class
class state:
    """Class defining state variables """
    def __init__(self):
        self.s = {}  # System state vars
        self.c = []  # creates a new empty list for components
        self.c.append('nan')  # Define an empty set in index 0
       

    def mixed(self):
        self.m = {}  # Mixture states
        self.m['a_mix'] = {}  # Mixture activity coefficient states
        self.m['b_mix'] = {}  # Mixture co-volume coefficient states

    def pure(self, p, i):
        self.c.append({})  # Pure component states
        self.c[i]['b'] = p.c[i]['b_c']  # Invariant co volume parameter.
        
    # % Calculate state at P, T, {composition} for all phases
        
    def update_state(self, s, p, P=None, T=None, phase=['All'], X=None,
                     Force_Update=False):
        """
        This function calculates the new state variables of the system at the 
        specified pressure, temperature and composition vector.
        
        Parameters
        ----------
        s : class
            Contains the dictionaries with the state of each component.
    
        p : class
            Contains the dictionary describing the mixture parameters.
            
        P : scalar, optional
            Pressure (Pa), if None the current state pressure will be used.
    
        T : scalar, optional
            Temperature (K), if None the current state temperature will be 
            used.
            
        phase : string inside a list, optional
                Phase to be updated, if the default value 'All' is used then
                all viable phases will be updated.
            
        X : vectors embedded inside a list, optional
            Entries of list conaint the compositon vector with n - 1 components
            List must be in correct phase order specified in 
            p.m['Valid phases'] and components in correct component number as 
            specified in the data. If None the current state composition will 
            be used.
            
            Example specification with 4 independent components:
                    #     x_1  x_2  x_3    y_1  y_2  y_3
                    X = [[0.1, 0.3, 0.6], [0.3, 0.4, 0.3]]
                    
            If only less than 'All' phases are being updated, specify only the
            needed composition vector.
            
            Example for phase='x' with 4 independent components:
                    #     x_1  x_2  x_3
                    X = [[0.1, 0.3, 0.6]]
                    
        Force_Update : booleean, optional
                       If True this will force an update of a single vector
                       input for X for all phases regardless of the number
                       of phases specified with "phase". X must be a single
                       vector/list to avoid failure.

        Dependencies
        ------------
        numpy.array
    
        Returns
        -------
        s : class output.
            Contains the new updated state. Values changed:
                s.c[i]['a'] for all components p.m['n']
                s.m['a']
                s.b['a']
                if P: 
                    s.c[n]['P'] for all components p.m['n']
                if T: 
                    s.c[n]['T'] for all components p.m['n']  
                if X: 
                    s.c[i][ph] for all components p.m['n'] for all phases.
                    
        Examples
        --------
        s = s.update_state(s, p, 
                           P=I['P'], 
                           T=I['T'], 
                           phase=['x','y'], 
                           X=[[0.5,0.2],[0.2,0.2]])
            
        """
        from numpy import array, size
        # Independent updates
        if P is not None:  # Update pressure
            s.m['P'] = P
            for i in range(1, p.m['n'] + 1): 
                s.c[i]['P'] = P

        if T is not None:  # Update temperature
            s.m['T'] = T
            for i in range(1, p.m['n'] + 1): 
                s.c[i]['T'] = T

        if X is not None: # Update new compositions
            if Force_Update:
                for ph in p.m['Valid phases']:
                    Sigma_x_dep = 0.0 # Sum of dependent components
                    for i in range(1, p.m['n']):
                        if size(X) == 1:   # Ugly fix because of the zip 
                            X = array([X])  # conversion to float or list
                            
                        s.c[i][ph] = X[i-1] 
                        Sigma_x_dep += X[i-1]
                    #  Dependent component
                    s.c[p.m['n']][ph] = 1.0 - Sigma_x_dep
                        
            elif phase[0] is 'All':
                if len(X) != len(p.m['Valid phases']): # Check for dimensions
                    raise IOError('The array of X specified does not match'
                                   +' the number of phases expected. len(X) = '
                                   '{} .'.format(len(X))
                                   +'len(p.m[\'Valid phases\']) = {}'
                                    .format(len(p.m['Valid phases'])))
                                   
                for xp, ph in zip(X, p.m['Valid phases']):
                    Sigma_x_dep = 0.0 # Sum of dependent components
                    for i in range(1, p.m['n']): #  Independent components
                        if size(xp) == 1:   # Ugly fix because of the zip 
                            xp = array([xp])  # conversion to float or list

                        s.c[i][ph] = xp[i-1] 
                        Sigma_x_dep += xp[i-1] 
                        
                    #  Dependent component
                    s.c[p.m['n']][ph] = 1.0 - Sigma_x_dep
                    
            elif phase is not None:
                if len(X) != len(phase): # Check for dimensions
                    raise IOError('The array of X specified does not match'
                                   +' the number of phases expected. len(X) = '
                                   '{} .'.format(len(X))
                                   +'len(phase) = {}'
                                    .format(len(phase)))
                                    
                for xp, ph in zip(X, phase):
                    Sigma_x_dep = 0.0 # Sum of dependent components
                    for i in range(1, p.m['n']): #  Independent components
                        if size(xp) == 1:   # Ugly fix because of the zip 
                            xp = array([xp])  # conversion to float
                        
                        s.c[i][ph] = xp[i-1] 
                        Sigma_x_dep += xp[i-1] 
                        
                    #  Dependent component
                    s.c[p.m['n']][ph] = 1.0 - Sigma_x_dep
                    
        # Dependent updates
        for i in range(1, p.m['n']+1): # Update pure activity coefficents
            s.c[i]['a'] = VdW.a_T(s.c[i],p.c[i])['a'] # a(T)
        
        try:  # Note: Highly non-linear models
            if phase[0] is 'All':
                for ph in p.m['Valid phases']:
                    s.m['a_mix'][ph] = a_mix(s, p, phase=ph)
                    s.m['b_mix'][ph] = b_mix(s, p, phase=ph)
                
            elif phase is not None:
                for ph in phase:
                    s.m['a_mix'][ph] = a_mix(s, p, phase=ph)
                    s.m['b_mix'][ph] = b_mix(s, p, phase=ph)
                    
        except(ValueError, ZeroDivisionError): # DO NOT RAISE, SET PENALTY
            s.s['Math Error'] = True
            print('WARNING: Math Domain error in s.update_state'
                  +'at %s Pa %s K' %(s.s['P'],s.s['T']) )
        
        try: # Default optimization state z
            s.m['a'] = a_mix(s, p, phase='x') 
            s.m['b'] = b_mix(s, p, phase='x')   
            
        except(KeyError):
            raise IOError('Specify at least one viable phase as \'x\' for' 
                           'optimization routines')
            
        # Find Volume Roots ('V_v' and 'V_l') at P, T, a, b for all components
        # and phases
        for i in range(1, p.m['n']+1): 
            s.c[i] = VdW.V_root(s.c[i], p.c[i])
            
        s.m = VdW.V_root(s.m, p.m) # 'V_v' and 'V_l' mixture volumes at x1, x2
  
        return s
            

# %% Define mixture models
def a_ij(s, p, i=1, j=1):  # (Validated)
    """
    Returns the temperature dependent a_ij parameter at specified indices.

    Parameters
    ----------
    s : class
        Contains the dictionaries with the state of each component.
        Note: s.c[i]['a'] and s.c[j]['a'] MUST be updated to the current state
              temperature value (s.s['T']) before calling.

    p : class
        Contains the dictionary describing the mixture parameters.
        Holds p.m['k'][i][j]

    i, j : int, optional
           The first component indices.

    phase : string, optional
            Phase to be calculated, ex. liquid phase 'x'.

    Dependencies
    ------------
    math.sqrt

    Returns
    -------
    a_ij : scalar output.
    """
    from math import sqrt
    if i == j:
        return s.c[i]['a']  # Return pure paramater
    else:  # find mixture aij i =/= j
        return (1 - p.m['k'][i][j]) * sqrt(s.c[i]['a'] * s.c[j]['a'])


def a_mix(s, p, phase='x'):  # (Validated)
    """
    Returns the calculated DWPM mixture parameter at current system state for
    the specified phase.

    Parameters
    ----------
    s : class
        Contains the dictionaries with the state of each component.

    p : class
        Contains the dictionary describing the mixture parameters.
        Holds p.m['n'], p.m['r'] and p.m['s']

    phase : string, optional
            Phase to be calculated, ex. liquid phase 'x'.

    Dependencies
    ------------
    nComp.a_ij

    Returns
    -------
    a_mix : scalar output.
    """
    # SUM^n_i x_i * [Sigma^n_j=1 (x_k a_ij^s)^(r/s)]
    Sigma2 = 0
    for i in range(1, p.m['n']+1):
        # SUM^n_j=1 (x_k a_ij^s)
        Sigma1 = 0
        for j in range(1, p.m['n']+1):
            Sigma1 += s.c[j][phase] * a_ij(s, p, i, j)**p.m['s']

        # SUM^n_j=1 (x_k a_ij^s)^(r/s)
        Sigma1rs = Sigma1**(p.m['r']/p.m['s'])

        Sigma2 += s.c[i][phase] * Sigma1rs

    # a_mix = (SUM^n_i x_i * [Sigma^n_j=1 (x_k a_ij^s)^(r/s)])^(1/r)
    return Sigma2**(1/p.m['r'])


def a_mix_partial_k(s, p, k=1, phase='x'):  # (Validated)
    """
    Returns a_mix_partial of component k in the specified phase.


    Parameters
    ----------
    s : class
        Contains the dictionaries with the state of each component.

    p : class
        Contains the dictionary describing the mixture parameters.
        Holds p.m['n'], p.m['r'] and p.m['s']

    phase : string, optional
            Phase to be calculated, ex. liquid phase 'x'.

    Dependencies
    ------------
    nComp.a_ij

    Returns
    -------
    a_mix_partial_k : scalar output.
    """
    amix = a_mix(s, p, phase)

    Term1 = (1 - 1/p.m['r'] - 1/p.m['s']) * amix
    Term2 = 1/p.m['r'] * amix**(1 - 1/p.m['r'])

    # CT1 = SUM^n_i=1 x_i * a_ik^s * [SUM^n_j=1 (x_j * a_ij^s)^(r/s - 1)]
    Sigma2 = 0
    for i in range(1, p.m['n']+1):
        # SUM^n_j=1 x_j * a_ij ^ s
        Sigma1 = 0
        for j in range(1, p.m['n']+1):
            Sigma1 += s.c[j][phase] * a_ij(s, p, i, j)**p.m['s']

        # SUM^n_j=1 (x_j * a_ij^s)^(r/s - 1)
        Sigma1rs = Sigma1**(p.m['r']/p.m['s'] - 1)

        # x_i * a_ik^s * [SUM^n_j=1 (x_j * a_ij^s)^(r/s - 1)]
        Sigma2 += s.c[i][phase] * (a_ij(s, p, i, k)**p.m['s']) * Sigma1rs

    CT1 = (p.m['r']/p.m['s']) * Sigma2

    # CT2 = (SUM^n_j=1 (x_j * a_ij^s)^(r/s)
    Sigma3 = 0
    for j in range(1, p.m['n']+1):
        Sigma3 += s.c[j][phase] * a_ij(s, p, k, j)**p.m['s']

    CT2 = Sigma3**(p.m['r']/p.m['s'])

    return Term1 + Term2 * (CT1 + CT2)


def b_mix(s, p, phase='x'):  # (Validated)
    """
    Returns the calculated volume mixture parameter at current system state for
    the specified phase.

    Parameters
    ----------
    s : class
        Contains the dictionaries with the state of each component.

    p : class
        Contains the dictionary describing the parameters.

    phase : string, optional
            Phase to be calculated, ex. liquid phase 'x'.

    Returns
    -------
    b_mix : scalar output.
    """
    b_mix = 0.0
    for i in range(1, p.m['n']+1):
        b_mix += s.c[i][phase]*s.c[i]['b']
    return b_mix

def d_b_mix_d_x(s, p, d=1, z=1, m=1, phase='x'):  # (Validated)
    """
    Returns the calculated differential of the volume mixture parameter at 
    current system state for the specified component phase.

    Parameters
    ----------
    s : class
        Contains the dictionaries with the state of each component.

    p : class
        Contains the dictionary describing the parameters.

    d : int, optional
        The differential order. Should be 1 for Jacobian entires and 2 for 
        Hessian entires.

    z : int, optional
        The n - 1 independent component number of the first differential.

    m : int, optional
        The n - 1 independent component number of the second differential.

    phase : string, optional
            Phase to be calculated, ex. liquid phase 'x'.
            
    Returns
    -------
    d_b_mix_dx : scalar output.
    """
    if d == 1: 
        return s.c[z]['b'] -  s.c[p.m['n']]['b']
        #  Note s.c[p.m['n']]['b'] refers to the volume parameter of the 
        #       independent component
    if d == 2:
        return 0.0  # For all x_z
    
# %% Define Gibbs energy functions for VdW EoS (cubic with 2 volume roots)
def g_R_k_i(s, p, k='x', i=1):  # (Validated)
    """
    Residual Gibbs energy g_R_k_i((T,P) of a single component where k is the 
    specified phase and i is the component in phase k. 
    
    Parameters
    ----------
    s : class
        Contains the dictionaries with the state of the specified component.

    p : class
        Contains the dictionary describing the component parameters.

    k : string, optional
        Phase to be calculated, ex. liquid phase 'x'.
        
    i : integer, optional
        The component number in phase k to calculate. ex. 1

    Dependencies
    ------------
    math
    
    Returns
    -------
    g_R_k_i : scalar output.
    """
    from math import log
    if k == 'y':  # 'y' = Vapour phase standard
        V = s.c[i]['V_v'] 
    else:  # Assume all phases other than vapour are liquid, ex. 'x'
        V = s.c[i]['V_l'] 
        
    return (s.c[i]['P'] * V / (p.c[i]['R'] * s.c[i]['T']) - 1.0
            - log(s.c[i]['P'] / (p.c[i]['R'] * s.c[i]['T']))
            - log(V - s.c[i]['b']) 
            - s.c[i]['a'] / (p.c[i]['R'] * s.c[i]['T'] * V))
            
def g_R_mix_i(s, p, k='x'):  # (Validated)
    """
    Residual Gibbs energy g_R_mix_i((T,P) of a mixture where k is the specified 
    phase.

    
    Parameters
    ----------
    s : class
        Contains the dictionaries with the state of the specified component.

    p : class
        Contains the dictionary describing the component parameters.

    k : string, optional
        Phase to be calculated, ex. liquid phase 'x'.
        
    Dependencies
    ------------
    math
    
    Returns
    -------
    g_R_k_i : scalar output.
    """
    from math import log
    if k == 'y':  # 'y' = Vapour phase standard
        V = s.m['V_v'] 
    else:  # Assume all phases other than vapour are liquid, ex. 'x'
        V = s.m['V_l'] 
        
    return (s.m['P'] * V / (p.m['R'] * s.m['T']) - 1.0
            - log(s.m['P'] / (p.m['R'] * s.m['T']))
            - log(V - s.m['b']) 
            - s.m['a'] / (p.m['R'] * s.m['T'] * V))
    
def g_IG_k(s, k='x'):  # (Validated)
    """
    Change in gibbs energy for mixing of ideal gasses.
    
    Parameters
    ----------
    s : class
        Contains the dictionaries with the system state information.

    k : string, optional
        Phase to be calculated, ex. liquid phase 'x'. 
            
    Dependencies
    ------------
    math
    
    Returns
    -------
    g_IG_k : scalar output.
    """
    from math import log
          
    Sigma_g_IG_k = 0.0 # Sum of ideal gas terms
    for i in range(1, p.m['n']+1): 
        if s.c[i][k] == 0.0:  # Prevent math errors from zero log call.
            pass  # should be = 0 as s2['y']*log(s2['y']) = 1*log(1) = 0
        else:
            Sigma_g_IG_k += s.c[i][k] * log(s.c[i][k])
    return Sigma_g_IG_k

    
def g_mix(s, p, k=None, ref='x', update_system=False):  # (Validated)
    """
    TO DO: -Add float conversion for 1x1 arrays.
           -Force 'k' lowest phase
           -
    
    Returns the gibbs energy at specified composition relative to a reference
    phase pure component gibbs energy.
    
    Parameters
    ----------
    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...
    
    p : class
        Contains the dictionary describing the parameters.
        
    k : string, optional # TO DO UNFINISHED
        Force phase to be calculated, ex. liquid phase 'x'. 
        
    ref : string, optional
          Selected reference phase. Note that is common practice to choose a
          liquid phase 'x' as the reference phase.
          
    update_system : boolean, optional # TO DO UNFINISHED
                    This updates the system state to the current P, T and 
                    composition conditions. Only use if the system dictionary 
                    has not been updated to the current independent variables.

    Dependencies
    ------------
    math
    
    Returns
    -------
    s : class output.
        Contains the following values (or more if more phases are chosen):
          s.m['g_mix']['t'] : scalar, Total Gibbs energy of mixing.
          s.m['g_mix']['x'] : scalar, Gibbs energy of mixing for liquid 
                                          phase.
          s.m['g_mix']['y'] : scalar, Gibbs energy of mixing for vapour 
                                          phase.
          s.s['Math Error'] : boolean, if True a math error occured during
                                       calculations. All other values set to 0.
    """
    if update_system: # TO DO TEST
        Xvec = [[]] # Construct update vector
        for i in range(1, p.m['n']): # for n-1 independent components
            Xvec[0].append(s.c[i][k])
        
        s = s.update_state(s, p, P = s .m['P'], T = s.m['T'], phase = k,
                               X = Xvec)
    
    s = s.update_state(s, p) # Update Volumes and activity coeff.
    
    try:
        Sigma_g_ref = 0.0
        for i in range(1, p.m['n'] + 1): 
            Sigma_g_ref -= s.c[i][ref] * g_R_k_i(s, p, k = ref, i=i)

        s.m['g_mix^R'] = {}
        for ph in p.m['Valid phases']:
                s.m['g_mix^R'][ph] = g_R_mix_i(s, p, k = ph) + Sigma_g_ref

        s.m['g_mix'] = {}
        g_min = []
        for ph in p.m['Valid phases']:
                s.m['g_mix'][ph] = s.m['g_mix^R'][ph] + g_IG_k(s, k=ph)
                g_min.append(s.m['g_mix'][ph])

        s.m['g_mix']['t'] = min(g_min)

        s.s['Math Error'] = False

    except(ValueError, ZeroDivisionError):
        s.m['g_mix'] = {}
        g_min = []
        s.s['Math Error'] = True
        print 'WARNING: Math Domain error in g_mix(s,p)!'
        for ph in p.m['Valid phases']:
                s.m['g_mix'][ph] = 0.0

        s.m['g_mix']['t'] = 0.0
            
    return s

# %% Duality formulation
def ubd(Lambda, g_x_func, X_d, Z_0, s, p, X_bounds, k=['All']): # 
    """
    TO DO: Improve bounds
    
    Returns the upper bounding problem of the dual extremum. Return is negative
    to change to minimization problem.
    
    Parameters
    ----------
    Lambda : vector (1xn array)
             Contains the langrange (or diality) multipliers Lambda \in R^m
             to be optimised to the maximum value of the ubd.

    g_x_func : function
               Returns the gibbs energy at a the current composition 
               point. Should accept s, p as first two arguments.
               Returns a class containing scalar value .m['g_mix']['t']
               
    X_d : vector (1xn array)
          Contains the current composition point in the overall dual 
          optimisation. Constant for the upper bounding problem.
                
    Z_0 : vector (1xn array)
          Feed composition. Constant.
    
    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...
    
    p : class
        Contains the dictionary describing the parameters.
        
    k : list, optional
        List contain valid phases for the current equilibrium calculation.
        ex. k = ['x', 'y']
        If default value None is the value in p.m['Valid phases'] is retained.  
        
    Dependencies
    ------------
    numpy.array
    math.e

    Returns
    -------
    ubd : scalar
          Value of the upper bounding problem at Lamda.
    """
    from math import e
    from numpy import array
    X_d = array(X_d)  # Prevent float converstion of 1x1 arrays
    # Reset system from changed composition in bound calculation each call
    # Comp. invariant in UBD
    s = s.update_state(s, p,  X = X_d, phase = k, Force_Update=True)
    UBD = g_x_func(s, p).m['g_mix']['t'] + sum(Lambda * (Z_0 - X_d)) 

    # NOTE: G_p = g_x_func(s, p) with s.update_state(s, p,  X = Z_0)
    #       must strictly be added as a constraint to this problem.

    s = s.update_state(s, p, phase = k,  X = Z_0, Force_Update=True)  
    # G_p (Primal problem Z_0_i - x_i = 0 for all i)
    P = 0
    G_P = g_x_func(s,p).m['g_mix']['t']  
    if UBD > G_P:
        P +=  e**(abs(G_P - UBD)*1e2) # Penalty
        #print "WARNING: Primal problem > UBD"
        #print ('WARNING: Primal problem > UBD: abs(G_P - UBD)*1e3 = {}'
        #       .format(P))
               
    # Lamda bounds
    # Upper
    s = s.update_state(s, p,  X = X_bounds[0], phase = k, Force_Update=True)  
    G_upper = g_x_func(s, p).m['g_mix']['t']
    # Lower
    s = s.update_state(s, p,  X = X_bounds[1], phase = k, Force_Update=True) 
    G_lower = g_x_func(s, p).m['g_mix']['t']
    for i in range(p.m['n']-1):
        UB = ((UBD - G_upper)/(Z_0[i] - X_bounds[0][i]))
        LB = ((UBD - G_lower)/(Z_0[i] - X_bounds[1][i]))
        if Lambda[i] > UB:
            P +=  e**abs(Lambda[i] - UB)*1e2
        
        if Lambda[i] < LB:
            P +=  e**abs(Lambda[i] - LB)*1e2

    s = s.update_state(s, p,  X = Z_0, Force_Update=True) 
    return -UBD + P # -UBD to minimize max problem



def lbd(X, g_x_func, Lambda_d, Z_0, s, p, k=['All']):
    """
    Returns the lower bounding problem of the dual extremum.
    
    Parameters
    ----------
    X_d : vector (1xn array)
          Contains the current composition point in the overall dual 
          optimisation to be optimised to the minimum value of the lbd.

    g_x_func : function
               Returns the gibbs energy at a the current composition 
               point. Should accept s, p as first two arguments.
               Returns a class containing scalar value .m['g_mix']['t']
                        
    Lambda_d : vector (1xn array)
               Contains the diality multipliers Lambda \in R^m.
               Constant for the lower bounding problem.
             
    Z_0 : vector (1xn array)
          Feed composition. Constant.
    
    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...
    
    p : class
        Contains the dictionary describing the parameters.
        
    k : list, optional
        List contain valid phases for the current equilibrium calculation.
        ex. k = ['x', 'y']
        If default value None is the value in p.m['Valid phases'] is retained.  
        
    Dependencies
    ------------
    numpy.array
    math.e

    Returns
    -------
    lbd : scalar
          Value of the lower bounding problem at X.
    """
    Xn = X
    # Update system to new composition.
    s.update_state(s, p,  X = Xn, phase=k, Force_Update=True)  

    return g_x_func(s, p).m['g_mix']['t'] + sum(Lambda_d * (Z_0 - X))

     
def dual_equal(s, p, g_x_func, Z_0, k=None, P=None, T=None, 
               tol=1e-2):
    """
    Dev notes and TODO list
    -----------------------
    TO DO: -Add valid phases option.
           -Add constraints to x_i =/= 0 or 1 for all i to prevent vertical
            hyperplanes.
           -Construct bounds which is a list of tuples in [0,1] \in R^n
           -UBD X_d[0] --> X_d
           -Lambda_Sol: Allow for multiple solutions and impliment bounds to
                        find any other solutions if they exist.
           -Lambda_Sol bounds method: The current implementation is only viable
                                      for binary systems. Improve.
                                      
    NOTES: -Strictly the composition in all phases in should be specified in
            X_d, refer to older versions of this script when different comp.
            spaces need to be used.
    -----------------------
    Find the phase equilibrium solution using the daul optimization algorithm.
    Ref. Mitsos and Barton (2007)
    
    Parameters
    ----------
    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...
    
    p : class
        Contains the dictionary describing the parameters.

    g_x_func : function
               Returns the gibbs energy at a the current composition 
               point. Should accept s, p as first two arguments.
               Returns a class containing scalar value .m['g_mix']['t']

    k : list, optional
        List contain valid phases for the current equilibrium calculation.
        ex. k = ['x', 'y']
        If default value None is the value in p.m['Valid phases'] is retained.  

    P : scalar, optional
        Pressure (Pa), if unspecified the current state pressure will be used.

    T : scalar, optional
        Temperature (K), if unspecified  the current state temperature will be 
        used.

    Z_0 : vector
          Contains the feed composition point (must be and unstable point to 
          find multiphase equilibria).
          
    tol : scalar, optional
          Tolerance, if epsilon >= UBD - LBD that will terminate the routine.
          
    Dependencies
    ------------
    numpy

    Returns
    -------
    s.m['Z_eq'] : vector
                  Contains the optimised equilibrium point. 
    s.m['Lambda_d'] : vector
                      Contains the langrange (or diality) multipliers Lambda at
                      the equilibrium solution.
                      
    """
    from numpy import array, zeros_like, zeros
    from scipy.optimize import minimize, differential_evolution
    from tgo import tgo
    
    def x_lim(X): # limiting function used in TGO
        import numpy
        return numpy.sum(X, axis=1) - 1.0 <= 0.0
    
    if k == None:
        k = p.m['Valid phases']
        
    # Initialize
    Z_0 = array(Z_0)
    LBD = - 1e300  # -inf
    s.update_state(s, p,  X = Z_0, phase = k, Force_Update=True) 
        
    # G_p (Primal problem Z_0_i - x_i = 0 for all i):
    UBD = g_x_func(s, p).m['g_mix']['t']  
    Lambda_d = zeros_like(Z_0)

    # X bounds used in UBD optimization
    X_bounds = [zeros(shape=(p.m['n']-1)), # Upper bound
                zeros(shape=(p.m['n']-1))  # Lower bound
                ]  
    for i in range(p.m['n']-1): 
        Sigma_ind = 0.0  # Sum of independent components excluding i
                         # (lever rule) 
        for j in range(p.m['n']-1):
            if j != i:
                Sigma_ind += Z_0[j]

        X_bounds[0][i] = 1.0 - Sigma_ind 
        
    # Construct physical bounds x \in [0, 1] for all independent components
    # (used in differential_evolution)
    Bounds = []
    for i in range(p.m['n'] - 1):
        #Bounds.append((0, 1))
        Bounds.append((1e-5, 0.99999))
    
    # Construct composition container from X_d for all phases. (REMOVED)
    X_d = Z_0
#    X_d = []
#    for ph in k:
#        X_d.append(Z_0)

    #%% Normal calculation of daul problem if Z_0 is unstable.
    while abs(UBD - LBD) >= tol:
        # Solve UBD
         # Update system to new composition.
         # (Comp. invariant in UBD)
        s = s.update_state(s, p,  X = X_d , phase = k, Force_Update=True) 
        Lambda_sol = minimize(ubd, Lambda_d, method='L-BFGS-B', 
                              args=(g_x_func, X_d[0], Z_0, s, p, X_bounds,
                                    k))['x']
                                    
        Lambda_d = array(Lambda_sol)  # If float convert back to 1x1 array
    
        # NOTE: NEGATIVE THE MAX DEFINED PROBLEM:
        UBD = -ubd(Lambda_d, g_x_func, X_d[0], Z_0, s, p, X_bounds, k)

        X_sol = tgo(lbd, Bounds, args=(g_x_func, Lambda_d, Z_0, s, p, k),
                                 g_func=x_lim,
                                 n = 100,
                                 skip=2)
        # Solve LBD                  
        LBD = lbd(X_sol, g_x_func, Lambda_d, Z_0, s, p, k) 
        X_d = X_sol
        # End
       

    if False:  # Print results optional
        print 'Final UBD = {}'.format(UBD)
        print 'Final LBD = {}'.format(LBD)
        print 'Final UBD - LBD = {}'.format(UBD - LBD)
        print 'Final Z_eq = {}'.format(X_d)
        print 'Final Lambda_d = {}'.format(Lambda_d)
    # Returns
    s.m['Z_eq'] = X_d
    s.m['Lambda_d'] = Lambda_d
            
    return s



def Eq_sol(X, g_x_func, Lambda_d, X_I, s, p, k=['All']):
    """
    
    TO DO WRITE THIS
    
TO DO, Why does returning the natiove present an objective function with the
    second minimizer of the equilibrium problem as the global minima?
    Is this true for all systems?

    Returns the lower bounding problem of the dual extremum.
    
    Parameters
    ----------
    X : vector (1xn array)
        Contains the current composition point the solution of which is the
        point in equilibrium with X_I. 

    g_x_func : function
               Returns the gibbs energy at a the current composition 
               point. Should accept s, p as first two arguments.
               Returns a class containing scalar value .m['g_mix']['t']
                        
    Lambda_d : vector (1xn array)
               Contains the solution diality multipliers Lambda \in R^m.
               Constant for the lower bounding problem.
             
    X_I : vector (1xn array)
          First equilibrium solution point. Constant.
    
    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...
    
    p : class
        Contains the dictionary describing the parameters.
        
    k : list, optional
        List contain valid phases for the current equilibrium calculation.
        ex. k = ['x', 'y']
        If default value None is the value in p.m['Valid phases'] is retained.  
        
    Dependencies
    ------------
    numpy.array
    math.e

    Returns
    -------
    lbd : scalar
          Value of the lower bounding problem at X.
    """
    import numpy
    Xn = X_I
    # Update system to new composition.
    s.update_state(s, p,  X = Xn, phase=k, Force_Update=True)  
    G_I = g_x_func(s, p).m['g_mix']['t'] 

    Xn = numpy.array(X)
    s.update_state(s, p,  X = Xn, phase=k, Force_Update=True)  

#    return abs(G_I + sum(Lambda_d * (X - X_I))  # Root prblem
#               - g_x_func(s, p).m['g_mix']['t'] ) 

    return -(G_I + Lambda_d * (X - X_I)
               - g_x_func(s, p).m['g_mix']['t'] )
          
# %% Parameter goal Functions

# %% Plots
    
def plot_g_mix(s, p, g_x_func,  Tie=None, x_r=1000):
    """
    Plots the surface of the Gibbs energy of mixing function. Binary and 
    Trenary plots only.
    
    Parameters
    ----------
    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...
    
    p : class
        Contains the dictionary describing the parameters.

    g_x_func : function
               Returns the gibbs energy at a the current composition 
               point. Should accept s, p as first two arguments.
               Returns a class containing scalar value .m['g_mix']['m']
               
    x_r : int, optional
          Number of composition points to plot.
    
    Tie : list of vectors, optional
          Equilibrium tie lines (of n - 1 independent components) to be added 
          to the plots.

    Dependencies
    ------------
    numpy
    matplotlib

    """
    from matplotlib import rc
    from matplotlib import pyplot as plot
    from numpy import array, linspace, float64 
    #% Binary
    if p.m['n'] == 2:
        #% Initialize
        s.m['g_mix range'] = {} # Contains solutions for all phases
        s.m['g_mix range']['t'] = []
        for ph in p.m['Valid phases']:
            s.m['g_mix range'][ph] = []
            
        s.c[1]['x_range'] = linspace(0, 1, x_r)
        for i in range(len(s.c[1]['x_range']-1)):  # Calculate range of G_mix
            X = [s.c[1]['x_range'][i]]
            s = s.update_state(s, p, X = X, Force_Update=True)
            s = g_x_func(s, p)
            s.m['g_mix range']['t'].append(float64(s.m['g_mix']['t']))
            for ph in p.m['Valid phases']:
                s.m['g_mix range'][ph].append(float64(s.m['g_mix'][ph]))

        plot.figure()

        for ph in p.m['Valid phases']:
            plot.plot(s.c[1]['x_range'], s.m['g_mix range'][ph], label=ph)

        plot.plot(s.c[1]['x_range'], s.m['g_mix range']['t'])
        if Tie is not None: # Add tie lines
            for point in Tie:
                X = [point[0]]
                s = s.update_state(s, p, X = X, Force_Update=True)
                s = g_x_func(s, p)
                G1 = float64(s.m['g_mix']['t'])
                X = [point[1]]
                s = s.update_state(s, p, X = X, Force_Update=True)
                s = g_x_func(s, p)
                G2 = float64(s.m['g_mix']['t'])
                plot.plot(point,[G1, G2])
                Slope = (G1 - G2)/ (point[0] - point[1])
                plot.annotate('slope = {}'.format(Slope), 
                              xy=(point[1], G2 + G2*0.5) )
                
        plot.xlabel(r"$x_1$", fontsize=14)
        plot.ylabel(r"$\Delta$g", fontsize=14)
        plot.title("{}-{} Gibbs energy of mixing at T = {}, P = {}".format(
                    p.c[1]['name'][0],
                    p.c[2]['name'][0],
                    s.m['T'],
                    s.m['P']))
        plot.legend()
        return
    
    #% Trenary
    if p.m['n'] == 3:
        pass
    else:
        raise IOError('Too many independent components to plot hypersurface'
                      +'in R^3.')

def plot_ep(func, X_r, s, p, args=()):
    """
    Plot the speficied error function over a range X_r
    """
    from matplotlib import rc
    from matplotlib import pyplot as plot
   # from numpy import float32
    #% Binary
    if p.m['n'] == 2:
        ep = []
        for X in X_r:
            s.update_state(s, p,  X = X, phase = ['All'], Force_Update=True) 
            ep.append(func(X, *args))
            
        plot.figure()
        plot.plot(X_r, ep)
        plot.xlabel(r"$x_1$", fontsize=14)
        plot.ylabel(r"$\epsilon$", fontsize=16)
        plot.title("{}-{} at T = {}, P = {}".format(
                    p.c[1]['name'][0],
                    p.c[2]['name'][0],
                    s.m['T'],
                    s.m['P']))
    return
    
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
                                             
#%
def g_range_test(s, p, x_r):  # NOTE THESE ARE FOR BINARY FUNCS
    from numpy import linspace
    #% Initialize
    s.s['del g_mix l soln'] = []
    s.c[1]['x_range'], s.c[2]['x_range'] = linspace(0,1,x_r), linspace(1,0,x_r)
    for i in range(len(s.c[1]['x_range']-1)): 
        s.c[1]['x'],s.c[2]['x'] = s.c[1]['x_range'][i],s.c[2]['x_range'][i]

        s.s['del g_mix l soln'].append(g_x_test_func(s, p).m['g_mix']['t'])

        
# %% Define Gibbs energy plotting functions.
def g_range(s, p, x_r):  # NOTE THESE ARE FOR BINARY FUNCS
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
        s = g_mix(s, p)   #.m['del g_mix l']
        s.s['del g_mix l soln'].append(s.m['g_mix']['x'])
        s.s['del g_mix v soln'].append(s.m['g_mix']['y'])
        s.s['del g_mix soln'].append(s.m['g_mix']['t'])
        s.s['Math Error soln'].append(s.s['Math Error'])
   
    return s

#  Numerical FD estimates for validation
def CDM(f, s, p, d=1, z=1, m=1, dx=1e-6, gmix=False):  
    """"
    Central difference estimate of a function f(s, p) used to validate the
    symbolic expressions.

    Parameters
    ----------
    f : function
        Function to differentiate
    
    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...
    
    p : class
        Contains the dictionary describing the parameters.
        
    d : int. optional
        Differential order (1-2)
        
    z : int, optional
        First component to differentiate to.
        
    m : int, optional
        Second component to differentiate to.
        
    dx : float, optional
         Spatial discretization size.
         
    fmix : booleaan, optional
           Specify True when using a Gibbs surface function with a class return
           and _.m['g_mix']['t'] float return.
    """
    if d == 1:
        s.c[z]['x'] += 0.5*dx
        X_d = []
        for i in range(1, p.m['n']):
            X_d.append(s.c[i]['x'])
            
        s = s.update_state(s, p, X = X_d, Force_Update=True) 
        if gmix:
            f1 = f(x,p).m['g_mix']['t']
        else:
            f1 = f(s, p)
            
        s.c[z]['x'] -= 1.0*dx
        X_d = []
        for i in range(1, p.m['n']):
            X_d.append(s.c[i]['x'])
            
        s = s.update_state(s, p, X = X_d, Force_Update=True) 
        if gmix:
            f2 = f(x,p).m['g_mix']['t']
        else:
            f2 = f(s, p)
        
        return (f1 - f2)/dx
        
    if d == 2:
        if gmix:
            f1 = f(x,p).m['g_mix']['t']
        else:
            f1 = f(s, p)
            
        s.c[z]['x'] += 1.0*dx
        X_d = []
        for i in range(1, p.m['n']):
            X_d.append(s.c[i]['x'])
            
        s = s.update_state(s, p, X = X_d, Force_Update=True) 
        if gmix:
            f2 = f(x,p).m['g_mix']['t']
        else:
            f2 = f(s, p)
            
        s.c[z]['x'] -= 2.0*dx
        X_d = []
        for i in range(1, p.m['n']):
            X_d.append(s.c[i]['x'])
            
        s = s.update_state(s, p, X = X_d, Force_Update=True) 
        if gmix:
            f3 = f(x,p).m['g_mix']['t']
        else:
            f3 = f(s, p)
        
        return (f2 - 2 * f1 + f3)/(dx*dx)
                
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
        if True: #%% TEST CURVE Mitsos et al. (2007)  ##  True: Validated 
            # Single phase test
            X_d = array([0.4]) 
            Z_0 = array([0.5]) 
            # Approx. Solution  Lambda = array([-0.03244]) X_d = array([0.004557])      
            s = s.update_state(s, p, P=101e3, T=300.0,  X = X_d) 
            s = dual_equal(s, p, g_x_test_func, Z_0 , tol=1e-9)
            X_I = s.m['Z_eq']
            
            Args= (g_x_test_func, s.m['Lambda_d'], X_I, s, p, ['All'])
            from tgo import tgo
            X_II = tgo(Eq_sol, [(1e-5, 0.99999)], args=Args, n=10)
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
        if False: 
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
                        
    #%% Jacobian and Hessian tests.
    if False: # Binary
        # b_mix

        X_d = [[0.6], [0.6]]#array([0.6, 0.6]) 
        s = s.update_state(s, p, P=24e5, T=263.1,  X = X_d) 
        

                
        print 'd_b_mix_d_x = {}'.format(d_b_mix_d_x(s, p))
        print 'FDM est. = {}'.format( d1(b_mix, s, p, dx=1e-1))
    
        print 'd2_b_mix_d_x2 = {}'.format(d_b_mix_d_x(s, p, d=2))
        print 'FDM est. = {}'.format( d1(b_mix, s, p, d=2, dx=1e-1))
    
    
    