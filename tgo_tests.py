#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
#from __future__ import division, print_function, absolute_import
from UQToolbox.sobol_lib import i4_sobol_generate
import numpy
import scipy
from tgo import tgo
import scipy.spatial
import scipy.optimize

def f_test_1(x, r, s): # Test function, bounds: -1 =< x_i =< 6
    return x[0]**2 + x[1]**2
    
Bounds1 = [(-1, 6), (-1, 6)]

Bounds2 = [(0, 1), (0, 1)]

def g_test_1(C):
    return numpy.sum(C, axis=1) - 6.0 <= 0.0# 


def f_test_3(x): # Test function, bounds: -1 =< x_i =< 6
    """
    Hock and Schittkowski 19 problem (HS19). Hoch and Schittkowski (1991)
    
    Approx. Answer:
        f_test_3([14.095, 0.84296]) = -6961.814744487831
    """
    return (x[0] - 10.0)**3.0 + (x[1] - 20.0)**3.0
    
Bounds3 = [(13.0, 100.0), (0.0, 100.0)]

def g_test_3(C):
     return ((-(C[:,0] - 5)**2 - (C[:,1] - 5)**2  - 100.0 <= 0.0)
             & ((C[:,0] - 6)**2 - (C[:,1] - 5)**2  - 82.81 <= 0.0))
             
def Rosen(x): # Rosenbrock's function # Ans x1 = 1, x2 = 2, f = 0
        return (1.0 - x[0])**2.0 + 100.0 * (x[1] - x[0]**2.0)**2.0 
        
BoundsR = [(-3.0, 3.0), (-3.0, 3.0)]

def plot_2D_sequance(B):
    """Plot the generated sequence to visualize uniformity of distrubtion."""
    from matplotlib import pyplot as plot
    plot.figure()
    #plot.plot(B[0], B[1], 'x')
    #plot.plot([0,1], [1,0], 'r--')
    #plot.figure(2)
    plot.plot(B[:,0], B[:,1], 'x')
    plot.plot([0,1], [1,0], 'r--')

    return

if __name__ == '__main__':
#    a = -1
#    b = 6
#    A = i4_sobol_generate(m, n, skip)  * (b - a) + a 

    test_atol = 1e-5

    r = [1, 2, 3] # random args for test func tuple
    s = True
    x1 = tgo(f_test_1, Bounds1, args=(r,s), g_func=g_test_1, n=500, 
            skip=1, k_t=None, 
        callback=None, minimizer_kwargs=None, disp=False)
    numpy.testing.assert_allclose(x1, [0., 0.], atol=test_atol)

#    x3 = tgo(f_test_3, Bounds3, args=(), g_func=g_test_3, n=500, 
#            skip=1, k=None, 
#        callback=None, minimizer_kwargs=None, disp=False)

    # OverflowError: Python int too large to convert to C long
    #   Func_min[i] = func(x_min, *args)
    # Why?
    # To do implement bounds in local search function  
    # >>> f_test_3([ -1.04572783e+08,-3.42296527e+08])
    # -4.12493867624096e+25

    xR = tgo(Rosen, BoundsR, args=(), g_func=None, n=500, 
            skip=1, k_t=None, 
        callback=None, minimizer_kwargs=None, disp=False)
    numpy.testing.assert_allclose(xR, [1., 1.], atol=test_atol)

    #%% Old
    pass
    
    if False: # (Old practice funcs an alternatives)
        m = 3 # dimensions
        n = 1000 # points
        skip = 1 # skip to this point in the sequence
        B = i4_sobol_generate(m, n, skip)  # Generate uniform sample points in R^m
        C = numpy.column_stack([B[i] for i in range(m)])

        # = numpy.hstack((B[0],B[1],B[2]))
        
        condition = sum(C) >=1
        
    #    F = numpy.sum(C, axis=1)
        G = C[numpy.sum(C, axis=1) <= 1.0]
        H = numpy.sum(G, axis=1)
    
        #numpy.inner(G[0],G[5])
        
        #d = numpy.tensordot(G,G,axes=(-1,-1))
        # If x is MxN, d will be an MxM array of dot products where d[i,j] 
        # = dot(x[i], x[j]).
        
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
        #%
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
        
        #args = ([1, 2], True)
        
        #B = map(func, T, args=(3)) # Rurn
        #vfunc = numpy.vectorize(func , otypes=[numpy.float]) # essentially a
                                                              # for loop
        # according to docs
        
        func = f_test_1
        args=(1,2)
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
    #    for i in Minimizers_indices:
    #        print i
        

        
        
        
        
        
        
        
        































        