#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""

from sobol_lib import i4_sobol_generate
import numpy
import scipy
#import numpy as np
#from scipy.optimize import minimize

def tgo(func, bounds, args=(), g_func=None, N=100, k=None, callback=None, 
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
        
    N : int, optional
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
    M = len(bounds)
             
    pass

if __name__ == '__main__':
    #from nComp import
    
    
    #tgo(func, bounds, args=(), g_func=None, N=100, k=None, callback=None, 
    #     minimizer_kwargs=None, disp=False)
    
    m = 3 # dimensions
    n = 1000 # points
    skip = 1 # skip to this point in the sequence

    B = i4_sobol_generate(m, n, skip)
    
    from matplotlib import rc
    from matplotlib import pyplot as plot
    
    if False:
        plot.figure(1)
        plot.plot(B[0], B[1], 'x')
        plot.plot([0,1], [1,0], 'r--')
    
    #C = numpy.column_stack((B[0],B[1],B[2]))
    #C = numpy.column_stack((B[0],B[1]))
    C = numpy.column_stack([B[i] for i in range(m)])
    # = numpy.hstack((B[0],B[1],B[2]))
    
    condition = sum(C) >=1
    
#    F = numpy.sum(C, axis=1)
    G = C[numpy.sum(C, axis=1) <= 1.0]
    H = numpy.sum(G, axis=1)

    if False:
        plot.figure(2)
        plot.plot(G[:,0], G[:,1], 'x')
        plot.plot([0,1], [1,0], 'r--')
    
    #numpy.inner(G[0],G[5])
    
    #d = numpy.tensordot(G,G,axes=(-1,-1))
    # If x is MxN, d will be an MxM array of dot products where d[i,j] = dot(x[i], x[j]).
    
    #ix = numpy.argsort(d,axis=-1)
    #will return an array of integers with the same shape as d such that
    #d[i,ix[i,0]], d[i,ix[i,1]], d[i,ix[i,2]], ...
    #are in ascending order.
    
    # scipy.spatial.distance.cdist( X, Y ) g
    
    # from scipy import spatial
    # A = scipy.spatial.distance.euclidean
    #scipy.spatial.distance.euclidean(G[0],G[1])
    
    # scipy.spatial.distance.pdist
    
    
    #A = scipy.spatial.distance.pdist(G, 'euclidean')
    #A = scipy.spatial.distance.cdist(G, G, 'euclidean')
    #%%
    T = numpy.array([[2, 5],   # P1
                     [1, 2],   # P2
                     [3, 4],   # P3
                     [0, 1],   # P4
                     [5, 0],   # P5
                     [4, 2]])  # P6
   
    Y = scipy.spatial.distance.cdist(T, T, 'euclidean')
    
    Y = scipy.spatial.distance.cdist(G, G, 'euclidean')
    
    T = G
    #  scipy.spatial.cKDTree
    
    
    Z = numpy.argsort(Y,axis=-1)
    
    A =  numpy.delete(Z, 0, axis=-1) # Topololy matrix without signs
    
    #Y = scipy.spatial.distance.cdist(C, C, 'euclidean')
    #Z = numpy.argsort(Y,axis=-1)
    #A =  numpy.delete(Z, 0, axis=-1)
    
    def func(x, p, r):  # Test function, bounds: -1 =< x_i =< 6
        return x[0]**2 + x[1]**2
        
    args = ([1, 2], True)
    
    #B = map(func, T, args=(3)) # Rurn
    #vfunc = numpy.vectorize(func , otypes=[numpy.float]) # essentially a for loop
    # according to docs
    
    #p = 3
    #B = vfunc(T)
    # Obj. function returns to be used as reference table.:
    F = numpy.zeros(numpy.shape(T)[0]) 
    for i in range(numpy.shape(T)[0]):
        F[i] = func(T[i,:], *args)
    # To Do see scipy.spatial.KDTree for F lookup
        
    H = F[A] # This replaces all index values in A with the function result
    #H2 = F[Z]
    #(H2 < 3)
    #(H2 < H2[:,0])
    #(H < H2[:,0])
    H2 = numpy.empty_like(H, dtype=bool)
    for i in range(numpy.shape(T)[0]):
        H2[i,:] = (H[i,:] > F[i])
    
    Minimizers = numpy.all(H2, axis=-1)
    # Find data point indexes of minimizers
    Minimizers_indices = numpy.where(Minimizers  == True)[0]
    Minimizers_indices = numpy.where(Minimizers)[0]
    for i in Minimizers_indices:
        print i
