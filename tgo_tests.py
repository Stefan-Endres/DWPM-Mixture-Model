#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test Run examples:

$ python2 -m unittest -v tgo_tests.TestTgoFuncs
$ python2 -m unittest -v tgo_tests.TestTgoSubFuncs
$ python2 -m unittest -v tgo_tests.TestTgoSubFuncs.test_t1
"""
import unittest
import numpy
from tgo import *

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

    def f(self, x):     # TODO: Add f bounds from original problem
        return (x[0] - 10.0)**3.0 + (x[1] - 20.0)**3.0

    def g(self, C):
        return ((-(C[:, 0] - 5)**2 - (C[:, 1] - 5)**2 - 100.0 <= 0.0)
                & ((C[:, 0] - 6)**2 - (C[:, 1] - 5)**2 - 82.81 <= 0.0))


# FIXME: The bounds appear not to include the expected value
test3 = Test3(bounds=[(13.0, 100.0), (0.0, 100.0)],
              expected=[14.095, 0.84296])


class Rosenbrock(TestFunction):
    """ Rosenbrock's function  Ans x1 = 1, x2 = 1, f = 0 """
    g = None

    def f(self, x):
        return (1.0 - x[0])**2.0 + 100.0*(x[1] - x[0]**2.0)**2.0


rosen = Rosenbrock(bounds=[(-3.0, 3.0), (-3.0, 3.0)],
                   expected=[1, 1])

test_atol = 1e-5


def run_test(test, args=()):
    x = tgo(test.f, test.bounds, args=args, g_func=test.g, n=500)
    numpy.testing.assert_allclose(x, test.expected, atol=test_atol)

# python2 -m unittest -v tgo_tests.TestTgoFuncs
class TestTgoFuncs(unittest.TestCase):
    """Global optimisation tests:"""
    def test_f1(self):
        r = [1, 2, 3]  # random args for test func tuple
        s = True
        run_test(test1_1, args=(r, s))

    @unittest.skip("OverflowError")
    def test_f3(self):
        """HS19 optimisation:"""
        run_test(test3)

        # OverflowError: Python int too large to convert to C long
        #   Func_min[i] = func(x_min, *args)
        # Why?
        # TODO: implement bounds in local search function
        # >>> test3.f([ -1.04572783e+08,-3.42296527e+08])
        # -4.12493867624096e+25

    def test_rosen(self):
        """Rosenbrock optimisation:"""
        run_test(rosen)

# python2 -m unittest -v tgo_tests.TestTgoSubFuncs
class TestTgoSubFuncs(unittest.TestCase):
    """TGO subfunction tests using known solution (test_f1)"""
    # int bool solution for known sampling points
    T_Ans = numpy.array([[0, 0, 0, 0, 0],
                         [0, 1, 1, 1, 1],
                         [1, 0, 0, 0, 0],
                         [1, 1, 1, 1, 1],
                         [0, 0, 1, 0, 1],
                         [1, 1, 0, 1, 0]])

    # Known order of sampling points
    A = numpy.array([[2, 1, 5, 3, 4],
                     [3, 2, 5, 0, 4],
                     [0, 5, 1, 3, 4],
                     [1, 5, 2, 0, 4],
                     [5, 1, 2, 3, 0],
                     [2, 4, 1, 0, 3]])

    # function values at test points
    F = numpy.array([29, 5, 25.81, 1, 25, 20])

    H = F[A]
    T = t_matrix(H, F).astype(int)
    def test_t1(self):
        """t-matrix construction:"""
        numpy.testing.assert_array_equal(self.T, self.T_Ans)

        #self.assertEqual(B, self.A)

    def test_t2(self):
        """k-1 topograph"""
        K_1 = k_t_matrix(self.T, 1).T[0] #
        numpy.testing.assert_array_equal(K_1 , self.T_Ans[:,0])

    def test_t3(self):
        """k-3 topograph"""
        K_3 = k_t_matrix(self.T, 3)
        Ans = numpy.delete(self.T_Ans, numpy.s_[3:numpy.shape(self.T_Ans)[1]]
                     , axis=-1)
        numpy.testing.assert_array_equal(K_3, Ans)

    def test_t3(self):
        """Minimizer function"""
        self.assertEqual(numpy.float32(minimizers(self.T_Ans)), 3)

    def K_optimal(T):
        """K_optimal"""
        numpy.testing.assert_array_equal(K_optimal(self.T), A)

def tgo_suite():
    """
    Gather all the TGO tests from this module in a test suite.
    """
    tgo_test_suite = unittest.TestSuite()
    tgo_suite1 = unittest.makeSuite(TestTgoFuncs)
    tgo_suite2 = unittest.makeSuite(TestTgoSubFuncs)
    tgo_test_suite.addTest(tgo_suite1)
    tgo_test_suite.addTest(tgo_suite2)
    return tgo_test_suite



if __name__ == '__main__':
    tgo_test_suite=tgo_suite()
    unittest.TextTestRunner(verbosity=2).run(tgo_test_suite)




