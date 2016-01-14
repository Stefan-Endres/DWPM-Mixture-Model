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
import unittest

class TestFunction(object):
    def __init__(self, bounds, expected):
        self.bounds = bounds
        self.expected = expected

class Test1(TestFunction):
    def f(self, x, r, s):
        return x[0]**2 + x[1]**2

    def g(self, C):
        return numpy.sum(C, axis=1) - 6.0 <= 0.0

test1_1 = Test1(bounds=[(-1, 6), (-1, 6)],
                expected=[0, 0])
test1_2 = Test1(bounds=[(0, 1), (0, 1)],
                expected=[0, 0])

class Test3(TestFunction):
    """
    Hock and Schittkowski 19 problem (HS19). Hoch and Schittkowski (1991)

    Approx. Answer:
        f_test_3([14.095, 0.84296]) = -6961.814744487831
    """
    def f(self, x):
        return (x[0] - 10.0)**3.0 + (x[1] - 20.0)**3.0

    def g(self, C):
        return ((-(C[:,0] - 5)**2 - (C[:,1] - 5)**2  - 100.0 <= 0.0)
                & ((C[:,0] - 6)**2 - (C[:,1] - 5)**2  - 82.81 <= 0.0))


#FIXME: The bounds appear not to include the expected value
test3 = Test3(bounds=[(13.0, 100.0), (0.0, 100.0)],
              expected=[14.095, 0.84296])

class Rosenbrock(TestFunction):
    """ Rosenbrock's function  Ans x1 = 1, x2 = 1, f = 0 """
    g = None
    def f(self, x):
        return (1.0 - x[0])**2.0 + 100.0 * (x[1] - x[0]**2.0)**2.0

rosen = Rosenbrock(bounds=[(-3.0, 3.0), (-3.0, 3.0)],
                   expected=[1, 1])

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

test_atol = 1e-5

def run_test(test, args=()):
    x = tgo(test.f, test.bounds, args=args, g_func=test.g, n=500,
            skip=1, k_t=None,
            callback=None, minimizer_kwargs=None, disp=False)

    numpy.testing.assert_allclose(x, test.expected, atol=test_atol)

class TestFunctions(unittest.TestCase):
    def test_1(self):
        r = [1, 2, 3] # random args for test func tuple
        s = True
        run_test(test1_1, args=(r, s))

    @unittest.skip("OverflowError")
    def test_3(self):
        run_test(test3)

        # OverflowError: Python int too large to convert to C long
        #   Func_min[i] = func(x_min, *args)
        # Why?
        # TODO: implement bounds in local search function
        # >>> test3.f([ -1.04572783e+08,-3.42296527e+08])
        # -4.12493867624096e+25

    def test_rosen(self):
        run_test(rosen)


if __name__ == '__main__':
#    a = -1
#    b = 6
#    A = i4_sobol_generate(m, n, skip)  * (b - a) + a 
    unittest.main()

    #%% Old
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
