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

def tgo(func, bounds, args=(), g_func=None, g_args=(), n=100, skip=1, k=None, 
        callback=None, minimizer_kwargs=None, disp=False):
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
        Any additional fixed parameters needed to completely specify the 
        objective function.
        
    g_func : callable, optional
        Function used to define a limited subset to defining the feasible set 
        of solutions in R^n in the form g(x) <= 0 applied as g : R^n -> R^m
        ex. To impose the constraint  x[0] + x[1] + x[2] - 1 <= 0
            Use the function definition
            
            def g_func(A):
                return numpy.sum(A, axis=1) - 1.0# 

TODO: Improve the g_func usage to allow for numpy manipulations of scalar
      functional definitiions.

    g_args : tuple, optional
        Any additional fixed parameters needed to completely specify the 
        feasible set function.
        
    n : int, optional
        Number of sampling points used in the construction of the topography
        matrix.
        
    skip : int, optional
        Tthe number of initial points to skip in the Sobol generation sequence.
        
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
    # %% Define funcs # TO DO: MOve
    def t_matrix(H, F):
        """
        Returns the topographical matrix with True boolean values indicating
        positive entries and False ref. values indicating negative values.
        """ 
        H2 = numpy.empty_like(H, dtype=bool)
        for i in range(numpy.shape(H)[0]):
            H2[i,:] = (H[i,:] > F[i])
        
        return H2
        
    def k_t_matrix(T, k):
        """Returns the k-t topograph matrix""" 
        return numpy.delete(T, numpy.s_[k:numpy.shape(T)[1]], axis=-1)
    
    def Minimizers(K): 
        """Returns the minimizer indexes of a k-t matrix""" 
        Minimizers = numpy.all(K, axis=-1)         
        # Find data point indexes of minimizers:
        #Minimizers_indices = numpy.where(Minimizers)[0]
        return numpy.where(Minimizers)[0]
        
    def K_optimal(T):
        """
        Returns the optimimal k-t topograph with the semi-empircal correlation
        proposed by Henderson et. al. (2015)
        """
        K_1 = k_t_matrix(T, 1)  # 1-t topograph
        k_1 = len(Minimizers(K_1))
        k_i = k_1
        i = 2
        while k_1 == k_i:
            K_i = k_t_matrix(T, i)
            k_i = len(Minimizers(K_i))
            i += 1
            
        ep = i * k_i / (k_1 - k_i)
        k_c = numpy.floor((-(ep - 1) + numpy.sqrt((ep - 1.0)--2 + 80.0 * ep))
                          / 2.0)
        
        k_opt = int(k_c + 1)
        if k_opt > numpy.shape(T)[1]:  # If size of k_opt exceeds t-graph size.
            k_opt = int(numpy.shape(T)[1])
        
        K_opt = k_t_matrix(T, k_opt) 
        return K_opt
        
        
    # %% Generate sampling points.
    m = len(bounds)  # Dimensions # TO DO Assert if func output matches dims.
    #print m
    B = i4_sobol_generate(m, n, skip)  # Generate uniform sample points in R^m
    C = numpy.column_stack([B[i] for i in range(m)])
    
    # Distribute over bounds
    # TO DO: Find a better way to do this
    for i in range(len(bounds)):
        C[:,i] = C[:,i] * (bounds[i][1] - bounds[i][0]) + bounds[i][0]
    
    if g_func is not None: # TO DO: Improve
        C =  C[g_func(C)<= 0.0]  # Subspace of usable points.
    
    Y = scipy.spatial.distance.cdist(C, C, 'euclidean')
    Z = numpy.argsort(Y, axis=-1)
    A = numpy.delete(Z, 0, axis=-1)  # Topographical matrix without signs
    
    # %% Obj. function returns to be used as reference table.:
    F = numpy.zeros(numpy.shape(C)[0]) 
    for i in range(numpy.shape(C)[0]):
        F[i] = func(C[i,:], *args)
    # To Do see scipy.spatial.KDTree for F lookup

    # %% Create float value and bool topograph:
    H = F[A] # This replaces all index values in A with the function result
    #
    T = t_matrix(H, F)  # Topograph with Boolean entries
    # %% Find the optimial k+ topograph
    # Find epsilon_i parameter for current system
    K1 = k_t_matrix(T, 1)  # 1-t topograph
    K_opt = K_optimal(T)
    
    # %% Local Search: Find the minimzer float values and 
    """
    TO DO IMPROVE
    """
    Min_ind = Minimizers(K_opt)
    x_vals = []
    Func_min = numpy.zeros_like(Min_ind)
    for i, ind in zip(range(len(Min_ind)), Min_ind):
        # Find minimum x vals
        x_min = scipy.optimize.minimize(func, C[ind,:], method='L-BFGS-B', 
                              args=args)['x']
        x_vals.append(x_min)
        # Find func float vals
        Func_min[i] = func(x_min, *args)
    
    # Find global of all minimizers
    i_glob = numpy.argsort(Func_min)[0]
    x_global_min = x_vals[i_glob]
    return x_global_min

    
if __name__ == '__main__':
    # Run tests
    import os 
    __location__ = os.path.realpath(
        os.path.join(os.getcwd(), os.path.dirname(__file__)))
    #f = open(os.path.join(__location__, 'tgo_tests.py'));
#    runfile(os.path.join(__location__, 'tgo_tests.py'))
    #from tgo_tests import Bounds1
    
    
    
 