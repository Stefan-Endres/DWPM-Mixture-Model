#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
#from __future__ import division, print_function, absolute_import
from sobol_lib import i4_sobol_generate
import numpy
import scipy
#import numpy as np
#from scipy.optimize import minimize

def tgo(func, bounds, args=(), g_func=None, n=100, k=None, callback=None, 
         minimizer_kwargs=None, disp=False):
    """
    Global optimization using the Topographical Global Optimization (TGO)
    method first proposed by Törn (1990) with [ TODO: EXPLAIN k-t matrix use]
    Hendorson et. al. (2015).
    
    The TGO is a clustering method that uses graph theory to generate good 
    starting points for local search methods from points distributed uniformly 
    in the interior of the feasible set. These points are generated using the 
    Sobol (1967) sequence.
    
    Parameters TO DO: REVIEW DOC
    ---------- 
    func : callable
        The objective function to be minimized.  Must be in the form
        ``f(x, *args)``, where ``x`` is the argument in the form of a 1-D array
        and ``args`` is a  tuple of any additional fixed parameters needed to
        completely specify the function.

    bounds : sequence
        Bounds for variables.  ``(min, max)`` pairs for each element in ``x``,
        defining the lower and upper bounds for the optimizing argument of
        `func`. It is required to have ``len(bounds) == len(x)``.
        ``len(bounds)`` is used to determine the number of parameters in ``x``.
        
    args : tuple, optional
        Any additional fixed parameters needed to
        completely specify the objective function.
        
    g_func : callable, optional
        Function used to define a limited subset to defining the feasible set 
        of solutions in R^n in the form g(x) <= 0 applied as g : R^n -> R^m
        ex. x[0] + x[1] + x[2] - 1 <= 0
        
    n : int, optional
        Number of sampling points used in the construction of the topography
        matrix.
        
    k : int, optional 
        Defines the number of columns constructed in the k-t matrix. The higher
        k is the lower the amount of minimizers will be used for local search
        routines. If None the empircal model of Henderson et. al. (2015) will 
        be used. (Note: Lower k values decrease performance, but could
        potentially be more robust due to testing more local minimizers in the
        function hypersuface)

TODO:    disp : bool, optional
        Display status messages

TODO:    callback : callable, `callback(xk, convergence=val)`, optional:
        A function to follow the progress of the minimization. ``xk`` is
        the current value of ``x0``. ``val`` represents the fractional
        value of the population convergence.  When ``val`` is greater than one
        the function halts. If callback returns `True`, then the minimization
        is halted (any polishing is still carried out).
        
TODO:    minimizer_kwargs : dict, optional
        Extra keyword arguments to be passed to the minimizer
        ``scipy.optimize.minimize()`` Some important options could be:

            method : str
                The minimization method (e.g. ``"L-BFGS-B"``)
            args : tuple
                Extra arguments passed to the objective function (``func``) and
                its derivatives (Jacobian, Hessian).

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a `OptimizeResult` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes. If `polish`
        was employed, then OptimizeResult also contains the `jac` attribute.
        
    References
    ----------
    .. Henderson, N, de Sá Rêgo, M, Sacco, WF, Rodrigues, RA Jr. (2015) "A new
       look at the topographical global optimization method and its application
       to the phase stability analysis of mixtures", Chemical Engineering 
       Science, 127, 151-174
    .. Sobol, IM (1967) "The distribution of points in a cube and the 
       approximate evaluation of integrals. USSR Comput. Math. Math. Phys. 7, 
       86-112.
    .. Törn, A (1990) "Topographical global optimization", Reports on Computer
       Science and Mathematics Ser. A, No 199, 8p. Abo Akademi University, 
       Sweden

    """
    m = len(bounds)  # Dimensions # TO DO Test if func output matches dims.
    B = i4_sobol_generate(m, n)  # Generate uniform sample points in R^m
    C = numpy.column_stack([B[i] for i in range(m)])
    pass

if __name__ == '__main__':
    # Run tests
    from tgo_tests import Bounds1
    
    print len(Bounds1)
    
    
 