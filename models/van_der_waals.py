#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""


class VdW:
    """
    Class containing all Van der Waals equation related functions.
    """
    def __init__(self):
        pass
    
    def a_m_sol(self, s, p):
        """
        Return explicit solution for 'm', at 'a', 'T'
        
        Parameters
        ----------
        s : dictionary
            Contains the current temperature 'T' and attraction paramter 'a'.
        
        p : dictionary
            Contains the critical parameters 'T_c', 'a_c'.
    
        Dependencies
        ------------
        math
        """
        from math import sqrt, log
        if p['Model'][0] == "Soave":
            m = (sqrt(s['a']/p['a_c']) - 1) / (1 - sqrt(s['T']/p['T_c']))
        elif p['Model'][0] == 'Adachi-Lu':
            m = p['T_c'] * log(p['a_c']/s['a']) / (s['T'] - p['T_c'])
        else:
            raise ValueError("Unknown model")
        return m
        
    def a_T(self, s, p):
        """
        Return explicit solution for 'a' at 'm', 'T'
        
        Parameters
        ----------
        s : dictionary
            Contains the current state variable temperature 'T' 
        
        p : dictionary
            Contains the critical paramters 'T_c
            ', 'a_c' and the a dependancy
            model 'Model' with parameter 'm'
            
        Dependencies
        ------------
        math
        """
        from math import e

        model = p['Model'][0] if (len(p['Model']) == 1) else p['Model']

        if model == "Soave":
            s['a'] = p['a_c']*(1.0 + p['m']*(1 - (s['T']/p['T_c'])**0.5))**2
        elif model == 'Adachi-Lu':
            s['a'] = p['a_c']*e**(p['m']*(1 - s['T']/p['T_c']))
        else:
            raise ValueError("Unknown model")

        return s # = s['a']

    def a_maxwell(self, s, p):
        """
        Explicit solution of 'a' parameter of the VdW EoS Maxwell integral at 
        Psat.
        
        Parameters
        ----------
        s : dictionary
            Contains the current temperature 'T' and phase volumes 'V_V','V_l'.
        
        p : dictionary
            Parameter dictionary container containing the universal gas 
            constant paramter 'R'.
    
        Dependencies
        ------------
        math
        """
        from math import log
        s['a'] = -(p['R']*s['T']*s['V_l']*s['V_v']*log(s['V_v'] \
                    /(s['V_l'] - s['b']) \
                 - s['b']/(s['V_l'] - s['b'])) + ((s['V_l']**2) * s['V_v'] \
                 - s['V_l'] * (s['V_v']**2)) * s['P']) / (s['V_l'] - s['V_v']) 
        return s
        
    def V_root(self, s, p): 
        """
        Calculates the volume roots of the van der Waals equation using the 
        analytic solution at specified values of P (or Psat), T, a(T) and b. If
        an analytical solution does not exist a numerical estimate is used.
        
        Parameters
        ----------
        s : dictionary
            Contains the current temperature state variables 'T', pressure 'P' 
            and the VdW coefficients 'a' and 'b'.
        
        p : dictionary
            Contains the critical paramters 'T_c', 'a_c', 'R'.
    
        Dependencies
        ------------
        numpy, math
        """   
        import logging
        from math import sqrt, acos, cos, pi
        try:
            # Coefficients of V^3 + (C_1)V^2 + (C_2)V + C_3 = 0
            C = [- (p['R']*s['T']/s['P'] + s['b']),   # Coefficient C_1
                 s['a']/s['P'],                       # Coefficient C_2
                 - s['a']*s['b']/s['P']               # Coefficient C_3
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
            s['V_v'], s['V_l'] = max(V_roots), min(V_roots) 
    
        except ValueError:
            import numpy
            # Coefficients of (C_0)V^3 + (C_1)V^2 + (C_2)V + C_3 = 0
            C = [ 1.0,                                # Coefficient C_0
                 - (p['R']*s['T']/s['P'] + s['b']),   # Coefficient C_1
                 s['a']/s['P'],                       # Coefficient C_2
                 - s['a']*s['b']/s['P']               # Coefficient C_3
                ]
            #try:
            V_roots = numpy.roots(C)
            s['V_v'], s['V_l']  = max(V_roots.real), min(V_roots.real)

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
        return s
        
    #%%  
    def Psat_V_roots(self, s, p, tol = 1e-20, estfactor=1e-3):
        """
        Calculates the saturation pressure and volume roots of the Van der 
        Waals equation at a specified temperature and pressure.
        
        Parameters
        ----------
        s : dictionary
            Contains the current state variable temperature 'T', pressure 'P' 
            and the VdW coefficients 'a' and 'b' 
        
        p : dictionrary
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
        s['a'] = self.a_T(s, p)['a']
        s = self.V_root(s,p) # Update s['V_v'] and s['V_l']]
        # The phase volumes at the specified pressure:
        s['V_v_P'], s['V_l_P'] = s['V_v'], s['V_l']
                
        def P_maxwell(P,s,p):
            s['P'] = P
            s = self.a_T(s,p)
            s = self.V_root(s,p) # Update s['V_v'] and s['V_l']
            try:
                return s['P']*(s['V_l'] - s['V_v']) \
                   + s['a']/s['V_v'] - (s['a']/s['V_l']) \
                   + p['R'] * s['T'] * log((s['V_v'] - s['b']) \
                                          /(s['V_l'] - s['b']))

            except ValueError:
                raise IOError('Math error in P_maxwell in Psat_V_roots, try'+\
                ' to use a lower starting value s[\'P\'] before executing the'\
                +' function')
                logging.warn('Value error in P_maxwell, P = '
                             + '{}'.format(s['P'])
                             )

        try:
            s['P_sat'] = fsolve(P_maxwell, s['P'], args=(s,p), xtol=tol)
        except IOError:
            try: # Scale to approx. P near P_sat
                s['P'] = p['P_c']**(s['T']/p['T_c'])*estfactor
                s = self.V_root(s,p) # Update s['V_v'] and s['V_l']
                s['P_sat'] = fsolve(P_maxwell, s['P'], args=(s,p), xtol=tol)
            except IOError:
                raise IOError('Math error in P_maxwell in Psat_V_roots, try'+\
                ' to use a lower starting value s[\'P\'] or a lower'+\
                '\'downscale\' argument before executing the function.')
                
        from numpy import float64 # Convert 1x1 np arrays to floats
        s['V_v'], s['V_l'], s['P_sat'], s['P'] = float64(s['V_v']), \
        float64(s['V_l']), float64(s['P_sat']), float64(s['P'])
        return s    
        
        # Outpus are:  'V_v_P'   = Vapour phase volume at specified Pressure
        #              'V_l_P'   = Liquid phase volume at specified Pressure
        #              'V_v'     = Vapour phase volume at saturation Pressure      
        #              'V_l'     = Liquid phase volume at saturation Pressure   
        #              'P_sat'   = Saturation Pressure at specified Temperature

