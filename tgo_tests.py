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
        # TO DO: implement bounds in local search function
        # >>> test3.f([ -1.04572783e+08,-3.42296527e+08])
        # -4.12493867624096e+25

    def test_rosen(self):
        run_test(rosen)


if __name__ == '__main__':
    unittest.main()
