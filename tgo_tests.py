#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
#from __future__ import division, print_function, absolute_import
from sobol_lib import i4_sobol_generate
import numpy
import scipy
from tgo import tgo

Bounds1 = [(0, 1), (0, 1)]

def plot_2D_sequance(B):
    """Plot the generated sequence to visualize uniformity of distrubtion."""
    from matplotlib import pyplot as plot
    plot.figure(1)
    plot.plot(B[0], B[1], 'x')
    plot.plot([0,1], [1,0], 'r--')
    return
    

if __name__ == '__main__':
    pass

   #tgo(func, bounds, args=(), g_func=None, N=100, k=None, callback=None, 
    #     minimizer_kwargs=None, disp=False)
    
    m = 3 # dimensions
    n = 1000 # points
    skip = 1 # skip to this point in the sequence

    B = i4_sobol_generate(m, n, skip)

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
