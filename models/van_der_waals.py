#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

def a_m_sol(state, parameters):
    """
    Return explicit solution for 'm', at 'a', 'T'

    Parameters
    ----------
    state : dictionary
        Contains the current temperature 'T' and attraction paramter 'a'.

    parameters : dictionary
        Contains the critical parameters 'T_c', 'a_c'.

    Dependencies
    ------------
    math
    """
    from math import sqrt, log

    T_c = parameters['T_c']
    a_c = parameters['a_c']
    T = state['T']
    a = state['a']

    if parameters['Model'][0] == "Soave":
        m = (sqrt(a/a_c) - 1) / (1 - sqrt(T/T_c))
    elif parameters['Model'][0] == 'Adachi-Lu':
        m = T_c*log(a_c/a)/(T - T_c)
    else:
        raise ValueError("Unknown model")
    return m

def a_T(state, parameters):
    """
    Return explicit solution for 'a' at 'm', 'T'

    Parameters
    ----------
    state : dictionary
        Contains the current state variable temperature 'T'

    parameters : dictionary
        Contains the critical paramters 'T_c
        ', 'a_c' and the a dependancy
        model 'Model' with parameter 'm'

    Dependencies
    ------------
    math
    """
    from math import e

    model = parameters['Model'][0] if (len(parameters['Model']) == 1) else parameters['Model']

    T_c = parameters['T_c']
    m = parameters['m']
    a_c = parameters['a_c']
    T = state['T']

    if model == "Soave":
        a = a_c * (1.0 + m*(1 - (T/T_c)**0.5))**2
    elif model == 'Adachi-Lu':
        a = a_c*e**(m*(1 - T/T_c))
    else:
        raise ValueError("Unknown model")

    state['a'] = a

    return state # = s['a']

def a_maxwell(state, parameters):
    """
    Explicit solution of 'a' parameter of the VdW EoS Maxwell integral at
    Psat.

    Parameters
    ----------
    state : dictionary
        Contains the current temperature 'T' and phase volumes 'V_V','V_l'.

    parameters : dictionary
        Parameter dictionary container containing the universal gas
        constant paramter 'R'.

    Dependencies
    ------------
    math
    """
    from math import log
    V_v = state['V_v']
    V_l = state['V_l']
    P = state['P']
    b = state['b']
    T = state['T']
    R = parameters['R']

    a = -(R*T*V_l*V_v*log(V_v/(V_l - b) - b/(V_l - b)) +
          (V_l**2*V_v - V_l*V_v**2)*P) / (V_l - V_v)

    state['a'] = a

    return state

def V_root(state, parameters):
    """
    Calculates the volume roots of the van der Waals equation using the
    analytic solution at specified values of P (or Psat), T, a(T) and b. If
    an analytical solution does not exist a numerical estimate is used.

    Parameters
    ----------
    state : dictionary
        Contains the current temperature state variables 'T', pressure 'P'
        and the VdW coefficients 'a' and 'b'.

    parameters : dictionary
        Contains the critical paramters 'T_c', 'a_c', 'R'.

    Dependencies
    ------------
    numpy, math
    """
    import logging
    from math import sqrt, acos, cos, pi
    if state['P'] == 0:
        state['P'] = 3.0
    try:
        # Coefficients of V^3 + (C_1)V^2 + (C_2)V + C_3 = 0
        C = [-(parameters['R'] * state['T'] / state['P'] + state['b']),  # Coefficient C_1
             state['a'] / state['P'],  # Coefficient C_2
             - state['a'] * state['b'] / state['P']  # Coefficient C_3
             ]
        # Substitutions (see solution of Cubic equations:
        #                         mathworld.wolfram.com/CubicFormula.html )
        w = (3*C[1] - C[0]**2)/3.0
        q = (27*C[2] - 9*C[0]*C[1] + 2*C[0]**3)/27.0
        R_t = (w/3.0)**3 + (q/3.0)**2.0
        q_s = q/abs(q) #math.copysign(1, q)

        if  R_t < 0:
            Ratio = sqrt(((q/2.0)**2)/(-(w/3.0)**3))
            if abs(Ratio) < 1.0:
                phi = acos(Ratio)
            else:
                raise ValueError # Raise Math error if no solution
                #phi = math.degrees(math.acos(Ratio-2)+math.pi)
                     # math.degrees(math.acos(Ratio-2))
        else:
            raise ValueError
        # Analytical expressions for Volume roots:
        V_roots = [-2*q_s*sqrt(-w/3.0)*cos(phi/3.0           ) - C[0]/3.0,
                   -2*q_s*sqrt(-w/3.0)*cos(phi/3.0 + 2*pi/3.0) - C[0]/3.0,
                    2*q_s*sqrt(-w/3.0)*cos(phi/3.0 + 4*pi/3.0) - C[0]/3.0
                  ]
        # Find physical volume roots
        state['V_v'], state['V_l'] = max(V_roots), min(V_roots)

    except ValueError:
        import numpy
        # Coefficients of (C_0)V^3 + (C_1)V^2 + (C_2)V + C_3 = 0
        C = [1.0,  # Coefficient C_0
             - (parameters['R'] * state['T'] / state['P'] + state['b']),  # Coefficient C_1
             state['a'] / state['P'],  # Coefficient C_2
             - state['a'] * state['b'] / state['P']  # Coefficient C_3
             ]
        #try:
        V_roots = numpy.roots(C)
        state['V_v'], state['V_l']  = max(V_roots.real), min(V_roots.real)

        #except(numpy.linalg.linalg.LinAlgError):  # Nan's in roots
        #    V_roots = numpy.array([0.0, 0.0, 0.0])
        #    s['V_v'], s['V_l'] = numpy.nan, numpy.nan

    if abs(V_roots[0].imag) > 0.1 or abs(V_roots[1].imag) > 0.1 \
                                  or abs(V_roots[2].imag) > 0.1:

        logging.warn('large imaginary roots in VdW.VRoot = '
                      + '{}, {}, {}'.format(V_roots[0].imag,
                                            V_roots[1].imag,
                                            V_roots[2].imag)
                     )
    return state

#%%
def Psat_V_roots(state, parameters, tol=1e-20, estfactor=1e-3):
    """
    Calculates the saturation pressure and volume roots of the Van der
    Waals equation at a specified temperature and pressure.

    Parameters
    ----------
    state : dictionary
        Contains the current state variable temperature 'T', pressure 'P'
        and the VdW coefficients 'a' and 'b'

    parameters : dictionrary
        Contains the critical paramters 'T_c', 'a_c', 'R' and the a
        dependancy model 'Model' with parameter 'm'

    tol : float, optional.
          Tolerance used in the Maxwell integral solution

    estfactor: float, optional
               If the first iteration attempt fails, a second saturation
               pressure is estimated from the an exponential scaling
               function based on the reduced temperature and critical
               pressure. The 'estfactor' will be multiplied with this
               factor, if convergence still fails even lower values should
               be attempted.

    Dependencies
    ------------
    numpy, math

    """
    import logging
    from math import log # NOTE: math.log is the natural logarithm, not b10
    from scipy.optimize import fsolve
    # Update s['a'] at specified T for given p['m']
    state['a'] = a_T(state, parameters)['a']
    state = V_root(state, parameters) # Update s['V_v'] and s['V_l']]
    # The phase volumes at the specified pressure:
    state['V_v_P'], state['V_l_P'] = state['V_v'], state['V_l']

    def P_maxwell(P, istate, iparameters):
        istate['P'] = P
        istate = a_T(istate, iparameters)
        istate = V_root(istate, iparameters) # Update s['V_v'] and s['V_l']
        try:
            b = istate['b']
            V_l = istate['V_l']
            V_v = istate['V_v']
            T = istate['T']
            a = istate['a']
            P = istate['P']
            R = iparameters['R']

            return P*(V_l - V_v) + a/V_v - a/V_l + R*T*log((V_v - b)/(V_l - b))

        except ValueError:
            raise IOError('Math error in P_maxwell in Psat_V_roots, try'+\
            ' to use a lower starting value s[\'P\'] before executing the'\
            +' function')
            logging.warn('Value error in P_maxwell, P = '
                         + '{}'.format(istate['P'])
                         )

    try:
        state['P_sat'] = fsolve(P_maxwell, state['P'], args=(state, parameters), xtol=tol)
    except IOError:
        try: # Scale to approx. P near P_sat
            state['P'] = parameters['P_c'] ** (state['T'] / parameters['T_c']) * estfactor
            state = V_root(state, parameters) # Update s['V_v'] and s['V_l']
            state['P_sat'] = fsolve(P_maxwell, state['P'], args=(state, parameters), xtol=tol)
        except IOError:
            raise IOError('Math error in P_maxwell in Psat_V_roots, try'+\
            ' to use a lower starting value s[\'P\'] or a lower'+\
            '\'downscale\' argument before executing the function.')

    from numpy import float64 # Convert 1x1 np arrays to floats
    for property in ['V_v', 'V_l', 'P_sat', 'P']:
        state[property] = float64(state[property])

    return state

    # Outpus are:  'V_v_P'   = Vapour phase volume at specified Pressure
    #              'V_l_P'   = Liquid phase volume at specified Pressure
    #              'V_v'     = Vapour phase volume at saturation Pressure
    #              'V_l'     = Liquid phase volume at saturation Pressure
    #              'P_sat'   = Saturation Pressure at specified Temperature

