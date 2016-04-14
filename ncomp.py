#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Script to simulate phase equilibria of multicomponent systems.
"""

# %% Imports
from __future__ import division
import data_handling
from models import van_der_waals
VdW = van_der_waals.VdW()
import data_handling
import plot


# %% Define state variable class
class state:
    """Class defining state variables """
    def __init__(self):
        self.s = {}  # System state vars
        self.c = []  # creates a new empty list for components
        self.c.append('nan')  # Define an empty set in index 0

        # Define the dynamically
        self.update_dp = None


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
                       vector/list to avoid failure. Use conservatively.

        Dependencies
        ------------
        numpy
        logging
    
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
        import logging
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
                            X = array([X])  # c
                            # onversion to float or list
                            
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
            logging.warn('Math Domain error in s.update_state'
                          +'at {} Pa {} K'.format(s.s['P'], s.s['T'])
                             )
                             
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

# Define multi component initialisation
def n_comp_init(data):
    """
    Initiate state and parameter class objects based on specified data
    input.
    """
    # Define parameter class
    p = data_handling.MixParameters()
    p.mixture_parameters(data.VLE, data)
    p.m['n'] = len(data.comps)  # Define system size
    for i in range(p.m['n']):  # Set params for all compounds
        p.parameters(data.c[i])  # Defines p.c[i]
        #p.parameter_build(data.c[i])
    p.m['R'] = p.c[1]['R']  # Use a component Universal gas constant

    # %% Initialize state variables
    s = state()
    s.mixed()  # Define mix state variable, call using s.m['key']
    # Define three component state variables (use index 1 and 2 for clarity)
    for i in range(1, p.m['n']+1):
        s.pure(p, i)  # Call using ex. s.c[1]['key']

    p.m['R'] = p.c[1]['R']  # Use a component Universal gas constant

    # %% Initialize state variables
    s = state()
    s.mixed()  # Define mix state variable, call using s.m['key']
    # Define three component state variables (use index 1 and 2 for clarity)
    for i in range(1, p.m['n']+1):
        s.pure(p, i)  # Call using ex. s.c[1]['key']

    return s, p


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
    from numpy import log
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
    from numpy import log
    if k == 'y':  # 'y' = Vapour phase standard
        V = s.m['V_v'] 
    else:  # Assume all phases other than vapour are liquid, ex. 'x'
        V = s.m['V_l'] 
        
    return (s.m['P'] * V / (p.m['R'] * s.m['T']) - 1.0
            - log(s.m['P'] / (p.m['R'] * s.m['T']))
            - log(V - s.m['b']) 
            - s.m['a'] / (p.m['R'] * s.m['T'] * V))
    
def g_IG_k(s, p, k='x'):  # (Validated)
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
    from numpy import log
          
    Sigma_g_IG_k = 0.0 # Sum of ideal gas terms
    for i in range(1, p.m['n']+1): 
        if s.c[i][k] == 0.0:  # Prevent math errors from zero log call.
            pass  # should be = 0 as s2['y']*log(s2['y']) = 1*log(1) = 0
        else:
            Sigma_g_IG_k += s.c[i][k] * log(s.c[i][k])
    return Sigma_g_IG_k

    
def g_mix(s, p, k=None, ref='x', update_system=False):  # (Validated)
    """
    Returns the gibbs energy at specified composition relative to a reference
    phase pure component gibbs energy.
    
    Parameters
    ----------
    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...
    
    p : class
        Contains the dictionary describing the parameters.
        
    k : string, optional # TODO UNFINISHED
        Force phase to be calculated, ex. liquid phase 'x'. 
        
    ref : string, optional
          Selected reference phase. Note that is common practice to choose a
          liquid phase 'x' as the reference phase.
          
    update_system : boolean, optional # TODO UNFINISHED
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
          s.m['g_mix']['x'] : scalar, Gibbs energy of mixing for liquid phase.
          s.m['g_mix']['y'] : scalar, Gibbs energy of mixing for vapour phase.
          s.m['g_mix']['ph min'] : string, contains the phase/volume root of 
                                           with lowest Gibbs energy.
          s.s['Math Error'] : boolean, if True a math error occured during
                                       calculations. All other values set to 0.
    """
    import logging
    if update_system: # TODO TEST
        Xvec = [[]] # Construct update vector
        for i in range(1, p.m['n']): # for n-1 independent components
            Xvec[0].append(s.c[i][k])
        
        s.update_state(s, p, P = s .m['P'], T = s.m['T'], phase = k,
                               X = Xvec)
    
    s.update_state(s, p) # Update Volumes and activity coeff.
    
    try:
        Sigma_g_ref = 0.0
        for i in range(1, p.m['n'] + 1): 
            Sigma_g_ref -= s.c[i][ref] * g_R_k_i(s, p, k = ref, i=i)

        s.m['g_mix^R'] = {}
        for ph in p.m['Valid phases']:
                s.m['g_mix^R'][ph] = g_R_mix_i(s, p, k = ph) + Sigma_g_ref

        s.m['g_mix'] = {}
        g_min = []
        g_abs_min = long  # "inf" large int
        for ph in p.m['Valid phases']:
                s.m['g_mix'][ph] = s.m['g_mix^R'][ph] + g_IG_k(s, p, k=ph)
                g_min.append(s.m['g_mix'][ph])
                if s.m['g_mix'][ph] < g_abs_min: # Find lowest phase string
                    s.m['g_mix']['ph min'] = ph
                    g_abs_min = s.m['g_mix'][ph]


        s.m['g_mix']['t'] = min(g_min)

        s.s['Math Error'] = False

    except(ValueError, ZeroDivisionError):
        s.m['g_mix'] = {}
        s.s['Math Error'] = True
        logging.error('Math Domain error in g_mix(s,p)')
        for ph in p.m['Valid phases']:
                s.m['g_mix'][ph] = 0.0

        s.m['g_mix']['t'] = 0.0
            
    return s

# %% Duality formulation
def ubd(X_D, Z_0, g_x_func, s, p, k=None):
    """
    Returns the arguments to be used in the optimisation of the upper bounding
    problem with scipy.optimize.linprog.

    used
    Parameters
    ----------
    X_D : vector (1xn array)
          Contains the current composition point in the overall dual
          optimisation. Constant for the upper bounding problem.

    Z_0 : vector (1xn array)
          Feed composition. Constant.

    g_x_func : function
               Returns the gibbs energy at a the current composition
               point. Should accept s, p as first two arguments.
               Returns a class containing scalar value .m['g_mix']['t']

    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...

    p : class
        Contains the dictionary describing the parameters.

    k : list, optional (TODO)
        List contain valid phases for the current equilibrium calculation.
        ex. k = ['x', 'y']
        If default value None is the value in p.m['Valid phases'] is retained.

    Returns
    -------
    c : array_like
        Coefficients of the linear objective function to be minimized.

    A : A_eq : array_like, optional
        2-D array which, when matrix-multiplied by x, gives the values of the
        upper-bound inequality constraints at x.

    b : array_like, optional
        1-D array of values representing the upper-bound of each inequality
        constraint (row) in.
    """
    import numpy
    # Coefficients of UBD linear objective function
    c = numpy.zeros([p.m['n']]) # linrpog maximize/minimize? D
                 # Documentation is contradictory across version; check
    c[p.m['n'] - 1] = -1.0  # -1.0 change max --> min problem

    # Coefficients of Lambda inequality constraints

    A = numpy.zeros([len(X_D) + 1, # rows = for all X_D + Global ineq
                     p.m['n']]     # cols = n comps + eta
                    )
    b = numpy.zeros(len(X_D) + 1)


    # Global problem bound (Fill last row of A and last element in b
    # G_p (Primal problem Z_0_i - x_i = 0 for all i)
    s = s.update_state(s, p, X = Z_0, Force_Update=True)
    G_P = g_x_func(s, p).m['g_mix']['t']
    A[len(X_D), p.m['n'] - 1] = 1  # set eta to 1
    b[len(X_D)] = G_P

    # Bounds for all X_d in X_D
    A[:, p.m['n'] - 1] = 1  # set all eta coefficients = 1
    for X, k in zip(X_D, range(len(X_D))):
        # Find G(X_d)
        s = s.update_state(s, p, X = X, Force_Update=True)
        G_d = g_x_func(s, p).m['g_mix']['t']
        b[k] = G_d
        for i in range(p.m['n'] - 1):
            A[k, i] = -(Z_0[i] - X_D[k][i])

    if False:
        print 'c shape = {}'.format(numpy.shape(c))
        print 'A shape = {}'.format(numpy.shape(A))
        print 'b shape = {}'.format(numpy.shape(b))
        print 'c = {}'.format(c)
        print 'A = {}'.format(A)
        print 'b = {}'.format(b)

    return c, A, b


def lbd(X, g_x_func, Lambda_d, Z_0, s, p, k=['All']):
    """
    Returns the lower bounding problem of the dual extremum.
    
    Parameters
    ----------
    X : vector (1xn array)
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
    # Update system to new composition.
    s.update_state(s, p,  X = X, phase=k, Force_Update=True)

    return g_x_func(s, p).m['g_mix']['t'] + sum(Lambda_d * (Z_0 - X))

     
def dual_equal(s, p, g_x_func, Z_0, k=None, P=None, T=None, tol=1e-9, n=100):
    """
    Dev notes and TODO list
    -----------------------
    TODO: -The material bounds is too high since (Mitsos')  \hat{X} is the true
            upper limit for a given feedpoint
          -Look into using X_bounds scheme of Pereira instead for low mole Z_0
          -Add valid phases option.
          -Add constraints to x_i =/= 0 or 1 for all i to prevent vertical
            hyperplanes.
          -Construct bounds which is a list of tuples in [0,1] \in R^n

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

    n : scalar, optional
        Number of sampling points used in the tgo routine in solving LBD of the
        dual problem.
        Note: It is recommended to use at least ``100 + p.m['n'] * 100``
          
    Dependencies
    ------------
    numpy

    Returns
    -------
    X_sol : vector
            Contains the first optimised equilibrium point of the dual problem

    Lambda_sol : vector
                 Contains the optimised lagrange multipliers (partial chemical
                 potential) sol. of the dual solution hyperplane

    d_res : optimisation object return
            Contains the final solution to the dual problem with the
            following values:
                d_res.fun : lbd plane solution at equil point
                d_res.xl : Other local composition solutions from final tgo
                d_res.funl : lbd plane at local composition solutions

    """
    import numpy
    from scipy.optimize import linprog
    from tgo import tgo

    def x_lim(X): # limiting function used in TGO defining material constraints
        import numpy
        return -numpy.sum(X, axis=-1) + 1.0
    
    if k == None:
        k = p.m['Valid phases']
        
    # Initialize
    Z_0 = numpy.array(Z_0)
    LBD = - 1e300  # -inf
    s.update_state(s, p,  X = Z_0, phase = k, Force_Update=True) 
        
    # G_p (Primal problem Z_0_i - x_i = 0 for all i):
    UBD = g_x_func(s, p).m['g_mix']['t']  
    Lambda_d = numpy.zeros_like(Z_0)

    # X bounds used in UBD optimization
    X_bounds = [[],  # Upper bound (bar x)
                []   # Lower bound (hat x)
                ]

    for i in range(p.m['n']-1):
        # Append an independent coordinate point for each bound
        X_bounds[0].append(numpy.zeros(shape=(p.m['n']-1)))
        X_bounds[1].append(numpy.zeros(shape=(p.m['n']-1)))

        # Set upper bound coordinate point i
        Sigma_ind = 0.0  # Sum of independent components excluding i
                         # (lever rule)
        for k_ind in range(p.m['n']-1): # Note: k var name is used as phase
            # Set k != i  (k==i changed at end of loop)
            X_bounds[0][i][k_ind] = Z_0[k_ind]
            if k_ind != i:
                Sigma_ind += Z_0[k_ind] # Test, use numpy.sum if working

                # Set Lower bound coordinate point i
                X_bounds[1][i][k_ind] = Z_0[k_ind]
                # (Remaining bound coordinate points kept at zero)

        X_bounds[0][i][i] = 1.0 - Sigma_ind # change from Z_0[k]

    # Construct physical bounds x \in [0, 1] for all independent components
    Bounds = []
    L_bounds = [] # Lambda inf bounds used in linprog.
    for i in range(p.m['n'] - 1):
        Bounds.append((1e-10, 0.99999999))
        L_bounds.append((-numpy.inf, numpy.inf))

    L_bounds.append((-numpy.inf, numpy.inf)) # Append extra bound set for eta

    # Update state to random X to initialise state class.
    X_sol = numpy.array(X_bounds[1][0])  # set X_sol to lower bounds
    s.update_state(s, p,  X = X_sol , phase = k, Force_Update=True)

    X_D = []  # set empty list
    # Add every X_bounds to X_D list to use in linprog
    for i in range(p.m['n']-1):
        X_D.append(X_bounds[0][i])
        X_D.append(X_bounds[1][i])

    #%% Normal calculation of daul problem if Z_0 is unstable.
    while abs(UBD - LBD) >= tol:
        # Solve UBD
        # Find new bounds for linprog
        c, A, b = ubd(X_D, Z_0, g_x_func, s, p)
        # Find mulitpliers with max problem.
        lp_sol = linprog(c, A_ub=A, b_ub=b, bounds=L_bounds)
        Lambda_sol = numpy.delete(lp_sol.x, numpy.shape(lp_sol.x)[0] - 1)
        # If float convert back to 1x1 array
        Lambda_sol = numpy.array(Lambda_sol)

        UBD = -lp_sol.fun  # Final func value is neg. of minimised max. problem

        # Solve LBD
        d_res = tgo(lbd, Bounds, args=(g_x_func, Lambda_sol, Z_0, s, p, k),
                                 g_func=x_lim,
                                 n = n,
                                 #n = 100 + 100*(p.m['n'] - 1),
                                 skip=2)

        X_sol = d_res.x
        # Calculate LBD
        LBD = lbd(X_sol, g_x_func, Lambda_sol, Z_0, s, p, k)
        X_D.append(X_sol)
        # End
       
    if False:  # Print results optional
        print 'Final UBD = {}'.format(UBD)
        print 'Final LBD = {}'.format(LBD)
        print 'Final UBD - LBD = {}'.format(UBD - LBD)
        print 'Final Z_eq = {}'.format(X_sol)
        print 'Final Lambda_d = {}'.format(Lambda_d)
        
    # Returns
    return X_sol, Lambda_sol, d_res


# %% Dual plane eta for optim visualization
def dual_plane(X, Lambda_sol, Z_0, g_x_func, s, p, k=['All']):
    """          Lambda_d, G_sol,
    Returns the scalar output of the dual solution hyperplane at X

    Parameters
    ----------
    X : vector
        Contains an input composition point to calculate the plane scalar
        output at X.

    Z_0 : vector
    Contains the feed composition point (must be and unstable point to
    find multiphase equilibria).

#TODO: Finish doc strings

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
    numpy

    Returns
    -------
    """
    s.update_state(s, p,  X = X, phase=k, Force_Update=True)
    return g_x_func(s, p).m['g_mix']['t'] + sum(Lambda_sol * (Z_0 - X))

# %% Multi-minimiser approach to phase equil. calculation
def phase_equilibrium_calculation(s, p, g_x_func, Z_0, k=None, P=None, T=None,
                                  tol=1e-9, gtol=1e-2, n=100, phase_tol=1e-3,
                                  Print_Results=False,
                                  Plot_Results=False,
                                  Sampling_Stepping=True):
    # TODO: Remove the roman numeral comp. returns and change the way
    # phase_seperation_detection works to iterate along all points found
    # the plane instead.

    """
    This function uses the duality formulation to calculate all equilibrium
    points from an unstable feed point Z_0.

    Parameters
    ----------
    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state before call at P, T, {x}, {y}
              if "P" and "T "are not specified.

    p : class
        Contains the dictionary describing the parameters.

    g_x_func : function
               Returns the gibbs energy at a the current composition
               point. Should accept s, p as first two arguments.
               Returns a class containing scalar value .m['g_mix']['t']

    Z_0 : 1xn vector
          Contains the feed composition point (must be and unstable point to
          find multiphase equilibria).

    k : list, optional
        List contain valid phases for the current equilibrium calculation.
        ex. k = ['x', 'y']
        If default value None is the value in p.m['Valid phases'] is retained.

    P : scalar, optional
        Pressure (Pa), if unspecified the current state pressure will be used.

    T : scalar, optional
        Temperature (K), if unspecified  the current state temperature will be
        used.

    tol : scalar, optional
          Tolerance used in ``dual_equal``, if epsilon >= UBD - LBD that will
          terminate the routine.

    gtol : scalar, optional
          Minimum tolerance between hyperplane solution
          Note: The Dual solution is not perfect so a low tolerance is
          required, but a too low tol could potentially include points that do
          not truly lie on the equilibrium plane within the considered
          instability region.

    n : scalar, optional
        Number of sampling points used in the tgo routine in solving LBD of the
        dual problem.
        Note: It is recommended to use at least ``100 + p.m['n'] * 100``

    phase_tol : scalar, optional
                The minimum seperation between equilibrium planes required to
                be considered a phase. Defaults to 0.001

    Print_Results : boolean, optional
                    If True the results of the calculation will be printed in
                    the console.

    Plot_Results : boolean, optional
                   If True the g_mix curve with tie lines will be plotted for
                   binary and ternary systems.

    Sampling_Stepping : boolean, optional
                        If True the number of sampling points will be increased
                        by a factor 10 when an equilibrium point is not detect-
                        ed.


    Returns
    -------
    X_eq : list of vectors
           Contains all points in equilibrium in the region of the dual plane.

    g_eq : dict containing scalars and vectors
           Contains the Gibbs evaluation at equil. each point including keys:
           'ph min' : phase key returning the minimum phase string information.
           't' : Gibbs energy key returning the minimum Gibbs scalar value.

    phase_eq : list of strings
               Contains the information on the phase of the equilibrium point
    '"""
    import numpy
    import logging
    # init
    X_eq  = []  # Empty list storing equilibrium points
    g_eq = []
    phase_eq = []  # Empty list storing the phase strings

    # Calculate hyperplane at Z_0
    s.update_state(s, p, P=P, T=T,  X = Z_0, Force_Update=True)
    X_sol, Lambda_sol, d_res = dual_equal(s, p, g_x_func, Z_0, tol=tol, n=n)
    # X_sol : the dual solution used as a reference point
    # Usable returns:
     # X_sol  # Equilibrium solution for calculated hyperplane
     # Lambda_sol  # Lambda (chemical potential) sol.
     # d_res.fun  # lbd plane solution at equil point
     # d_res.xl  # local composition solutions from final tgo
     # d_res.funl  # lbd plane at local composition solutions

    X_eq.append(X_sol)  # Add first point to solution set

    #TODO: Of len(d_res.xl) = 1 then we only need to optimize the solution
    # plane again
    if (len(d_res.xl) < 2) and Sampling_Stepping:
        logging.basicConfig(level=logging.DEBUG)
        logging.warn('Less than 2 equilibrium points found in dual, increasing'
                     ' duality sampling points to n=n * p.m[\'n\'] * 10')
        # Change n and tol
        X_sol, Lambda_sol, d_res = dual_equal(s, p, g_x_func, Z_0,
                                              tol= tol,
                                              n=n * p.m['n'] * 10)
        if len(d_res.xl) < 2:
            logging.warn('Less than 2 equilibrium points found in dual')
            return  X_eq, g_eq, phase_eq
        else:
            logging.info('Succesfully converged to new equilibrium point')
            #print d_res.xl


    # Exclude any Sigma X_i > 1 (happens with unbounded local solvers etc.)
    Flag = []
    for i in range(len(d_res.xl)):
        if sum(d_res.xl[i]) > 1.01:  # 1.1 is an approx. tol, should be 1.0
            Flag.append(i)

    d_res.xl = numpy.delete(d_res.xl, Flag, axis=0)
    d_res.funl  = numpy.delete(d_res.funl , Flag, axis=0)


    # Exclude same composition points within tolerance.
    # TODO: Find more efficient way to exclude minima points within a tolerance
    def x_plane_tol(xl, phase_tol):
        Flag = []  # Flag index for deletion from solution array
        for i in range(len(d_res.xl)):
            for j in range(i + 1, len(xl)):
                if numpy.allclose(xl[i], xl[j],
                                  rtol=phase_tol, atol=phase_tol):
                    Flag.append(j)  # Flag index j for deletion if points are
                                    # within a the specified phase tolerance from
                                    # each other

        return Flag

    Flag = x_plane_tol(d_res.xl, phase_tol)

    # Check of all points are flagged for deletion with no equilibrium point.
    if len(Flag) == len(d_res.xl):
        logging.warn('All equilibrium points excluded with current phase'
                     'tolerance, temporarily lowering phase_tol')
        Flag = x_plane_tol(d_res.xl, phase_tol * 1e-1)
        if len(Flag) == len(d_res.xl):
            Flag = x_plane_tol(d_res.xl, phase_tol * 1e-2)
            if len(Flag) == len(d_res.xl):
                logging.warn('Failed to find equilibrium point within '
                             'phase_tol')

                return X_eq, g_eq, phase_eq
            else:
                logging.info('Succesfully converged to new equilibrium point '
                             'within phase_tol * 1e-2 tolerance')
        else:
            logging.info('Succesfully converged to new equilibrium point'
                         'within phase_tol * 1e-1 tolerance')

    # Delete composition rows (lists are converted to arrays in func)
    d_res.xl = numpy.delete(d_res.xl, Flag, axis=0)
    # Delete the corresponding dual plane function values
    d_res.funl = numpy.delete(d_res.funl, Flag, axis=0)


    # Find the differences in plane solutions at each minima and add the
    # solutions withing tolerance to the solution set.
    def x_plane(X_eq, fun, funl, xl, gtol):
        for i in range(1, len(funl)):
            if abs(fun - funl[i]) < gtol:
                X_eq.append(xl[i])

        return X_eq

    X_eq = x_plane(X_eq, d_res.fun, d_res.funl, d_res.xl, gtol)

    if len(X_eq) < 2:
        logging.warn('Less than 2 equilibrium points found within dual plane '
                     'within surface tolerance, temporarily increasing gtol')

        X_eq = x_plane(X_eq, d_res.fun, d_res.funl, d_res.xl, gtol*1e2)

        if len(X_eq) < 2:
            #print d_res.fun
            #print d_res.funl
            return X_eq, g_eq, phase_eq
        else:
            logging.info('Succesfully converged to new equilibrium point'
                         'within gtol*1e2')

    # Find Gibbs and phase info for each point
    # Note: We need the phase information at each point which requires a
    #       calculation of the Gibbs energy to find the minimum phase.
    for i in range(len(X_eq)):
        s.update_state(s, p, P=P, T=T,  X = X_eq[i], Force_Update=True)
        g_eq.append(g_x_func(s, p).m['g_mix'])
        phase_eq.append(g_eq[i]['ph min'])


    if Plot_Results:
        # Gibbs mix func with Tie lines
        from scipy import linspace
        print 'Feeding Lamda_d = {} to ep. func.'.format(Lambda_sol)
        print [[X_eq[1], X_eq[0]]]  # TODO: Allow for more points?
        if p.m['n'] == 2: # Plot binary tie lines
            plot.plot_g_mix(s, p, g_x_func, Tie =[[X_eq[1], X_eq[0]]]
                            , x_r=1000)

            from scipy import linspace
            X_r = linspace(1e-5, 0.9999, 4000)
            plane_args = (Z_0, Lambda_sol, g_x_func, s, p, ['All'])
            plot.plot_ep(dual_plane, X_r, s, p, args=plane_args)

        if p.m['n'] == 3: # Plot ternary tie lines
            s.update_state(s, p, P=P, T=T,  X = X_sol, Force_Update=True)
            G_P = g_x_func(s, p).m['g_mix']['t']
            Tie = [[G_P,                           # G_P
                    numpy.array(X_eq[0][0]),       # x_1
                    Lambda_sol[0],                 # lambda_1
                    numpy.array(X_eq[0][1]),       # x_2
                    Lambda_sol[1]]                 # lambda_2
                    ]

            plot.plot_g_mix(s, p, g_x_func, Tie = Tie, x_r=100)

    return X_eq, g_eq, phase_eq

# Phase seperation detection
def phase_seperation_detection(g_x_func, s, p, P, T, n=100, LLE_only=False,
                               VLE_only=False, tol=1e-9, gtol=1e-3, n_dual=100,
                               phase_tol=1e-3, Print_Results=False,
                               Plot_Results=False):
    """
    Detect and calculate phase separations in hte composition space at the
    current thermodynamic state.

    Parameters
    ----------
    g_x_func : function
               Returns the gibbs energy at a the current composition
               point. Should accept s, p as first two arguments.
               Returns a class containing scalar value .m['g_mix']['t']

    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...

    p : class
        Contains the dictionary describing the parameters.

    P : scalar  # TODO: Add optional specification or update
        Pressure (Pa), if unspecified the current state pressure will be used.

    T : scalar
        Temperature (K), if unspecified  the current state temperature will be
        used.

    n : int, optional
        Number of sampling points to be tested for in the R^(p.m['n'] - 1)
        Note. For higher component systems higher numbers of n are recommended.

    LLE_only : boolean, optional
               If True then only phase seperation of same volume root
               instability will be calculated.

    VLE_only : boolean, optional
               If True then phase seperation of same volume root instability
               will be ignored.

    tol : scalar, optional
          Tolerance used in ``dual_equal``, if epsilon >= UBD - LBD that will
          terminate the routine.

    gtol : scalar, optional
          Minimum tolerance between hyperplane solution
          Note: The Dual solution is not perfect so a low tolerance is
          required, but a too low tol could potentially include points that do
          not truly lie on the equilibrium plane within the considered
          instability region.

    n_dual : scalar, optional
            Number of sampling points used in the tgo routine in solving LBD
            of the dual problem.
            Note: It is recommended to use at least ``100 + p.m['n'] * 100``

    phase_tol : scalar, optional
                The minimum seperation between equilibrium planes required to
                be considered a phase. Defaults to 0.001

    Print_Results : boolean, optional
                    If True the results of the calculation will be printed in
                    the console.

    Plot_Results : boolean, optional
                   If True the g_mix curve with tie lines will be plotted for
                   binary and ternary systems.

    Dependencies
    ------------
    numpy
    tgo

    Returns
    -------
    ph_eq : dict containing keys for each phase in p.m['Valid phases'], ex:
        ph_eq[ph] : list containing composition vectors
                    Contains a list of equilibrium points of phase (ph)
                    seperations in the same volume root of the EOS
                    (ex. LLE type)

    mph_eq : list containing composition vectors
             contains a list of equilibrium points of phase
             seperations in different volume roots of the EOS (mph)
             (ex. VLE type)

    mph_ph : list containing strings
             containts the phase string of the corresponding ``mph_eq``
             equilibrium point
    """
    # TODO: Update this documentation
    import numpy
    from UQToolbox.sobol_lib import i4_sobol_generate
    from tgo import tgo
    # init returns
    ph_eq = {}
    mph_eq = []
    mph_ph = []

    # Generate sampling points.
    m = p.m['n'] - 1
    skip = 4
    Points = i4_sobol_generate(m, n, skip)
    Points = numpy.column_stack([Points[i] for i in range(m)])
    Points = Points[numpy.sum(Points, axis=-1) <= 1.0]
    S = numpy.empty(n, dtype=bool)

    # Update P, T to specified value
    s = s.update_state(s, p,  P=P, T=T, X = Points[0], Force_Update=True)

    def subset_eqp(Points, EQ):
        # Returns a subset of "Points" outside EQP
        import numpy
        for i in range(p.m['n']-1):
            P_new_low = Points[Points[:,i] <
                               #min(X_I[i], X_II[i])]
                             min(EQ[j][i] - phase_tol for j in range(len(EQ)))]

            P_new_high = Points[Points[:,i] >
                             max(EQ[j][i] + phase_tol for j in range(len(EQ)))]

            return numpy.append(P_new_low, P_new_high, axis=0)


    # Detect instability in a same volume root phase:
    if not VLE_only:
        # define LLE instability func
        def instability_point_calc(Points, g_x_func, s, p, n, k, P=P, T=T):
            #  Find an instability point, calculated equilibrium and return
            #  new feasible subset.
            Stop = False # Boolean to run main while loop
            for i, X in zip(range(n), Points):
                # Test for instability at current equilibrium point.
                S[i] = stability(X, g_x_func, s, p, k=ph)
                if not S[i]: # If point is unstable find equilibrium point.
                    # noinspection PyTupleAssignmentBalance
                    X_eq, g_eq, phase_eq = phase_equilibrium_calculation(s, p,
                                                   g_x_func, X, k=k, P=P, T=T,
                                                   tol=tol, gtol=gtol,
                                                   phase_tol=phase_tol,
                                                   Print_Results=Print_Results,
                                                   Plot_Results=Plot_Results)

                    #s.m['ph equil P'] = [s.m['X_I'], s.m['X_II']]
                    ph_eq_P = X_eq
                    # TODO: Improve finding feasible subspace of points.
                    P_new = Points[(i+1):]
                    P_new = subset_eqp(P_new, X_eq)

                    # Stop if no values in subset
                    if numpy.shape(P_new)[0] == 0:
                        Stop = True

                    return P_new, ph_eq_P, Stop

            # If no instability was found, stop the main for loop and set eq.
            #  point to None.
            ph_eq_P = None
            Stop = True
            return Points, ph_eq_P, Stop

        # Main looping
        ph_eq = {} # Range of equilibrium points.
        for ph in p.m['Valid phases']:
            Stop = False
            ph_eq[ph] = []

            while not Stop:
                Points, ph_eq_P, Stop = instability_point_calc(Points,g_x_func,
                                                               s, p, n, ph)

                # Save an equilibrium point to the range of points in the
                # current phase if found.
                if ph_eq_P is not None:
                    ph_eq[ph].append(ph_eq_P)


    # Detect phase seperation accross volume root phases:
    if not LLE_only:
        # Define difference function
        def g_diff(X, g_x_func, s, p, ph1, ph2, ref):
            # Returns difference between Gibbs energy of phases 'ph1' & 'ph2'
            # Note, all phases must be at same composition for meaningful
            # comparison
            s = s.update_state(s, p,  X = X, Force_Update=True)
            return (g_x_func(s, p, k=ph1, ref=ref).m['g_mix'][ph1]
                    - g_x_func(s, p, k=ph2, ref=ref).m['g_mix'][ph2])


        # Define objective function for feed search
        def g_diff_obj(X, g_x_func, s, p, ph1, ph2, ref):
            # Returns difference between Gibbs energy of phases 'ph1' & 'ph2'
            # Note, all phases must be at same composition for meaningful
            # comparison
            s = s.update_state(s, p,  X = X, Force_Update=True)
            return abs(g_x_func(s, p, k=ph1, ref=ref).m['g_mix'][ph1]
                       - g_x_func(s, p, k=ph2, ref=ref).m['g_mix'][ph2])

        # Calculated difference of Gibbs energies between all phases at all
        # sampling points.
        mph_eq = []
        mph_ph = []
        for i in range(len(p.m['Valid phases'])):
            for j in range(i + 1, len(p.m['Valid phases'])):
                ph1 = p.m['Valid phases'][i]
                ph2 = p.m['Valid phases'][j]
                Fd = numpy.empty(n)
                for l, X in zip(range(n), Points):
                    Fd[l] = g_diff(X, g_x_func, s, p, ph1, ph2, ph1)

                # Look for sign cross phase seperation
                if not numpy.all(Fd > 0) or numpy.all(Fd < 0):
                    # (if all values are not greater than or less than zero)
                    Bounds = [(1e-6, 0.99999)]
                    Args=(g_x_func, s, p, ph1, ph2, ph1)
                    #TODO: Reterieve local minima instead and loop over Z_0
                    # while eliminating subspace
                    Z_0 = tgo(g_diff_obj, Bounds, args=Args, n=1000, k_t = 5).x

                    X_eq, g_eq, phase_eq = phase_equilibrium_calculation(s, p,
                                                   g_x_func, Z_0, P=P, T=T,
                                                   tol=tol, gtol=gtol, n=n,
                                                   phase_tol=phase_tol,
                                                   Print_Results=Print_Results,
                                                   Plot_Results=Plot_Results)

                    mph_eq.append(X_eq)
                    mph_ph.append(phase_eq)


    return ph_eq, mph_eq, mph_ph


# %% Stability and phase seperation
def stability(X, g_x_func, s, p, k):
    """
    Tests a specified point "X" in phase "k", returns True if stable.

    Parameters
    ----------
    X : vector (1xn array)
        Contains the current composition point to test for stability.

    g_x_func : function
               Returns the gibbs energy at a the current composition
               point. Should accept s, p as first two arguments.
               Returns a class containing scalar value .m['g_mix']['t']

    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...

    p : class
        Contains the dictionary describing the parameters.

    Dependencies
    ------------
    numpy

    Returns
    -------
    SB : boolean
         A True boolean is return if the point is stable, else False.
    """
    import numpy
    s.update_state(s, p, X = X, phase = k, Force_Update=True)
    H = hessian(g_x_func, s, p, dx=1e-6, gmix=True, k=k)
    Heig = numpy.linalg.eig(H)[0]
    HBeig = (Heig > 0.0)
    return numpy.all(HBeig)


# %%  Numerical FD estimates for validation
def FD(f, s, p, d=1, z=1, m=1, dx=1e-6, gmix=False, k=['All']):
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
        git
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

    k : list, optional
        List contain valid phases for the current equilibrium calculation.
        ex. k = ['x']
        If default value ['All'] is used, the minimum Gibbs phse ['t'] is used.
        And the composition values are assumed to be in ['x']

    Dependencies
    ------------
    numpy

    Returns
    -------
    FD : float
         Output of CFD estimate of f(s, p).

    """
    if k == ['All']:
        ph = 't'
        cph = 'x'
    else:
        ph = k[0]
        cph = k[0]

    if d == 1:
        s.c[z][cph] += 0.5*dx
        X_d = []
        for i in range(1, p.m['n']):
            X_d.append(s.c[i][cph])

        s = s.update_state(s, p, X = X_d, Force_Update=True)
        if gmix:
            f1 = f(s, p).m['g_mix'][ph]
        else:
            f1 = f(s, p)

        s.c[z][cph] -= 1.0*dx
        X_d = []
        for i in range(1, p.m['n']):
            X_d.append(s.c[i][cph])

        s = s.update_state(s, p, X = X_d, Force_Update=True)
        if gmix:
            f2 = f(s, p).m['g_mix'][ph]
        else:
            f2 = f(s, p)

        return (f1 - f2)/dx

    if d == 2:
        s.c[z][cph] += 1.0*dx
        s.c[m][cph] += 1.0*dx
        X_d = []
        for i in range(1, p.m['n']):
            X_d.append(s.c[i][cph])

        s = s.update_state(s, p, X = X_d, Force_Update=True)
        if gmix:
            f1 = f(s, p).m['g_mix'][ph]
        else:
            f1 = f(s, p)

        s.c[m][cph] -= 2.0*dx
        X_d = []
        for i in range(1, p.m['n']):
            X_d.append(s.c[i][cph])

        s = s.update_state(s, p, X = X_d, Force_Update=True)
        if gmix:
            f2 = f(s, p).m['g_mix'][ph]
        else:
            f2 = f(s, p)

        s.c[z][cph] -= 2.0*dx
        s.c[m][cph] += 2.0*dx
        X_d = []
        for i in range(1, p.m['n']):
            X_d.append(s.c[i][cph])

        s = s.update_state(s, p, X = X_d, Force_Update=True)
        if gmix:
            f3 = f(s, p).m['g_mix'][ph]
        else:
            f3 = f(s, p)

        s.c[m][cph] -= 2.0*dx
        X_d = []
        for i in range(1, p.m['n']):
            X_d.append(s.c[i][cph])

        s = s.update_state(s, p, X = X_d, Force_Update=True)
        if gmix:
            f4 = f(s, p).m['g_mix'][ph]
        else:
            f4 = f(s, p)

        return (f1 - f2 - f3 + f4)/(4.0*dx*dx)


def hessian(f, s, p, dx=1e-6, gmix=False, k =['All']):
    """
    Returns the Hessian of the function as a numpy array.

    Parameters
    ----------
    f : function
        Input function used to calculate Hessian.

    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...

    p : class
        Contains the dictionary describing the parameters.

    d : int. optional
        Differential order (1-2)

    dx : float, optional
         Spatial discretization size.

    fmix : booleaan, optional
           Specify True when using a Gibbs surface function with a class return
           and _.m['g_mix']['t'] float return.

    Dependencies
    ------------
    numpy

    Returns
    -------
    H : array
        Ouput hessian matrix.

    """
    import numpy
    N = (p.m['n'] - 1)
    H = numpy.zeros(shape=(N,N))
    for m in range(1, N + 1):
        for z in range(1, N + 1):
            H[m - 1, z - 1] = FD(f, s, p, 2, z, m, dx, gmix, k)

    return H

# Equilibrium range function
def equilibrium_range(g_x_func, s, p, Data_Range=False, PT_Range=None, n=100,
                      LLE_only=False, VLE_only=False, res=100, tol=1e-9,
                      gtol=1e-2, n_dual=100, phase_tol=1e-3,
                      Print_Results=False, Plot_Results=False):
    #TODO: Define a new function with a full phase envelope calculation method
    #      which detects the limits ex. Henderson 2004
    """
    Calculates the equilibrium composition points over a range of P, T values,
    if ``Data_Range`` is specified as True then the equilibrium only at each
    P, T datapoint wil be calculated, otherwise a range of P, T values will be
    calculated over the min-max pairs of the data or ``PT_Range`` can be used
    to specify a specific range of points.

    Parameters
    ----------
    g_x_func : function
               Returns the gibbs energy at a the current composition 
               point. Should accept s, p as first two arguments.
               Returns a class containing scalar value .m['g_mix']['t']
    
    s : class
        Contains the dictionaries with the system state information.
        NOTE: Must be updated to system state at P, T, {x}, {y}...
    
    p : class
        Contains the dictionary describing the parameters.

    Data_Range : boolean, optional
                 Specify true


    PT_Range : list of tuples
               Specify ranges as [(P_min, P_max), (T_min, T_max)]
               ex. [(1, 101e3), (100, 273.15)]

    n : int, optional
        Number of sampling points to be tested for in the R^(p.m['n'] - 1)
        Note. For higher component systems higher numbers of n are recommended.
        
    LLE_only : boolean, optional
               If True then only phase seperation of same volume root 
               instability will be calculated.
               
    VLE_only : boolean, optional
               If True then phase seperation of same volume root instability 
               will be ignored.

    res : integer
          Specifies the number of data points to be simulated within the
          specified range.

    tol : scalar, optional
          Tolerance used in ``dual_equal``, if epsilon >= UBD - LBD that will
          terminate the routine.

    gtol : scalar, optional
          Minimum tolerance between hyperplane solution
          Note: The Dual solution is not perfect so a low tolerance is
          required, but a too low tol could potentially include points that do
          not truly lie on the equilibrium plane within the considered
          instability region.

    n_dual : scalar, optional
            Number of sampling points used in the tgo routine in solving LBD
            of the dual problem.
            Note: It is recommended to use at least ``100 + p.m['n'] * 100``

    phase_tol : scalar, optional
                The minimum seperation between equilibrium planes required to
                be considered a phase. Defaults to 0.001

    Print_Results : boolean, optional
                    If True the results of the calculation will be printed in
                    the console.

    Plot_Results : boolean, optional
                   If True the g_mix curve with tie lines will be plotted for
                   binary and ternary systems.

    Dependencies
    ------------
    numpy
    scipy

    Returns
    -------

    P_range:

    T_range:

    r_ph_eq : list containing ph_eq returns:
        ph_eq : dict containing keys for each phase in p.m['Valid phases'], ex:
            ph_eq[ph] : list containing composition vectors
                        Contains a list of equilibrium points of phase (ph)
                        seperations in the same volume root of the EOS
                        (ex. LLE type)

    r_mph_eq : list containing mph_eq returns:
        mph_eq : list containing composition vectors
                 contains a list of equilibrium points of phase
                 seperations in different volume roots of the EOS (mph)
                 (ex. VLE type)

    r_mph_ph  : list containing mph_ph returns:
        mph_ph : list containing strings
                 containts the phase string of the corresponding ``mph_eq``
                 equilibrium point
    """

    import numpy
    import scipy
    Store = numpy.array([numpy.shape(p.m['P'])[0], (p.m['n'] -1 )])

    # Find limits
    if Data_Range:
        P_range, T_range = p.m['P'], p.m['T']
    elif PT_Range is not None:
        P_range = scipy.linspace(PT_Range[0][0], PT_Range[0][1], res)
        T_range = scipy.linspace(PT_Range[1][0], PT_Range[1][1], res)
    else:
        P_range = scipy.linspace(min(p.m['P']), max(p.m['P']), res)
        T_range = scipy.linspace(min(p.m['T']), max(p.m['T']), res)

    # init stores
    r_ph_eq = []
    r_mph_eq = []
    r_mph_ph = []

    for P, T in zip(P_range, T_range):
        #s.update_state(s, p, P=P, T=T)  # Updated in psd; not need
        try:
            ph_eq, mph_eq, mph_ph = phase_seperation_detection(g_x_func, s, p,
                                                               P=P, T=T, n=n,
                                           LLE_only=LLE_only, VLE_only=VLE_only,
                                           tol=tol, gtol=gtol, n_dual=n_dual,
                                           phase_tol=phase_tol,
                                           Print_Results=Print_Results,
                                          # Plot_Results=True)
                                           Plot_Results=Plot_Results)

        except(numpy.linalg.linalg.LinAlgError):
            pass
            ph_eq = []
            mph_eq = []
            mph_ph = []

        #try:
        r_ph_eq.append(ph_eq)
        r_mph_eq.append(mph_eq)
        r_mph_ph.append(mph_ph)
        # except(UnboundLocalError):
        #     r_ph_eq.append(None)
        #     r_mph_eq.append(None)
        #     r_mph_ph.append(None)

    return P_range, T_range, r_ph_eq, r_mph_eq, r_mph_ph

# Sort points into ordered dicts

# %% Parameter goal Functions
def parameter_goal_func(Params, g_x_func, s, p, n=100, LLE_only=False,
                        VLE_only=False, res=100, tol=1e-5, gtol=1e-2,
                        n_dual=100, phase_tol=1e-3, Print_Results=False,
                        Plot_Results=False):
    """
    Calculates the error function over a range of data with at the current
    parameter set ``Params``.

    Parameters
    ----------
    Params : list or 1xn numpy array
             DWPM-EOS: [p.m['r'], p.m['s'],
                        p.m['k'][1][i], ... for all i \in N <= p.m['n']
                        p.m['k'][2][i], ... for all i \in N <= p.m['n']
                        ...
                        p.m['k'][n][i], ... for all i \in N <= p.m['n']
                        ]
    g_x_func
    s
    p
    n
    LLE_only
    VLE_only
    res
    tol
    gtol
    n_dual
    phase_tol
    Print_Results
    Plot_Results

    Returns
    -------

    """
    # Set params to new state
    p.m['r'], p.m['s'] = Params[0], Params[1]
    pint = 2
   # for i in range(1, p.m['n'] + 1):
   #     for j in range(1, p.m['n'] + 1):
   #         if i == j:
   #             pass
   #         else:
   #             p.m['k'][i][j] = Params[pint]
   #             pint += 1


    if abs(p.m['r']) <= 1e-10:  # Avoid singularities
        p.m['r'] = 1e-3
    if abs(p.m['s']) <= 1e-10:
        p.m['s'] = 1e-3

    print 'Params = '
    print Params
    print p.m['r']
    print p.m['s']
    print p.m['k'][1][2]
    print p.m['k'][2][1]

    # Find points in data range
    try:
        P_range, T_range, r_ph_eq, r_mph_eq, r_mph_ph = \
            equilibrium_range(g_mix, s, p,
                              Data_Range=True,
                              PT_Range=None,
                              n=n,
                              LLE_only=LLE_only,
                              VLE_only=VLE_only,
                              res=res,
                              tol=tol,
                              gtol=gtol,
                              n_dual=n_dual,
                              phase_tol=phase_tol,
                              Print_Results=Print_Results,
                              Plot_Results=Plot_Results)
    except(IndexError):
        return len(p.m['x'][1])
    # Find data error for VLE (r_mph_eq)
    epsilon = 0
    for ph in p.m['Valid phases']:
        for c in range(1, p.m['n']): # Cycle through
            for ph_c_data_i, i in zip(p.m[ph][c], range(len(p.m[ph][c]))):
                #TODO: Cycle through each point found instead of r_mph_ph[i][0]
                # and then use the minimum error
                if len(r_mph_ph[i]) > 0:
                    for nph, j in zip(r_mph_ph[i][0],
                                      range(len(r_mph_ph[i][0]))):
                        if nph == ph:
                            epsilon += abs(ph_c_data_i - r_mph_eq[i][0][j]) ** 2.0

                else:
                    epsilon += 1.0 # Maximum error when no equilibrium point
                                   # was found

            #if x_data == 0.0 or x_data == 1.0: # Ignore purs
            #    break

           # for i in
    print 'Params = '
    print Params
    print 'epsilon = '
    print epsilon
    return epsilon

if __name__ == '__main__':
    pass