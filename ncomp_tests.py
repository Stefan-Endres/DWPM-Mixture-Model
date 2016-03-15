#!/usr/bin/env python
# TODO: Add list sorting to ensure right component is tests in phase detection.


import unittest
import data_handling
import Van_der_Waals
import numpy
import main
from nComp import *
VdW = Van_der_Waals.VdW()

#%% TEST FUNCTION  Binary NRTL
def g_x_test_func(s, p, k=None, ref='x'):
    """
    This is the test function of a binary NRTL Model of the water-butyl-acetate
    system from Misos et. al. (2007) using the parameters referenced in the
    paper.

    x_1^0 = 0.5 is an unstable point.
    """
    from math import log, e
    t_12 = 3.00498# tau paramters
    t_21 = 4.69071
    a_12 = 0.391965 # Alpha paramter
    a_21 = 0.391965#**(-1.0) # Checked. Should be a_12 in Mitsos, see SvA p 448

    for i in range(1, p.m['n']+1):
        if s.c[i]['x'] <= 1e-20:  # Prevent math errors from zero log call.
            s.m['g_mix'] = {}
            s.m['g_mix']['t'] = 0.0
            s.m['g_mix']['x'] = s.m['g_mix']['t']
            s.m['g_mix']['ph min'] = 'x'
            return s  # should be = 0 as s2['y']*log(s2['y']) = 1*log(1) = 0

    s.m['g_mix'] = {}
    s.m['g_mix']['t'] = ( s.c[1]['x'] * log(s.c[1]['x'])
                        + s.c[2]['x'] * log(s.c[2]['x'])
                        + s.c[1]['x'] * (s.c[2]['x'])
                        * ((t_12 * e**(-a_12 * t_12))
                        / (s.c[2]['x'] + s.c[1]['x'] * e**(- a_12 * t_12))
                        + (t_21 * e**(-a_21 * t_21))
                        / (s.c[1]['x'] + s.c[2]['x'] * e**(- a_21 * t_21))))

    s.m['g_mix']['x'] = s.m['g_mix']['t']
    s.m['g_mix']['ph min'] = 'x'
    return s

def g_x_test_func2(s, p, k=None, ref='x'):
    """
    This is the test function of a trenary NRTL Model of the toluene_water_
    aniline system from Misos et. al. (2007) using the parameters referenced
    in the paper.

    x_1^0 = [0.3, 0.2] is an unstable point.
    """
    from math import log, e
    t_12 = 4.93035  # tau paramters
    t_21 = 7.77063
    t_13 = 1.59806
    t_31 = 0.03509
    t_23 = 4.18462
    t_32 = 1.27932

    a_12 = 0.2485  # Alpha paramters
    a_21 = a_12
    a_13 = 0.3000
    a_31 = a_13
    a_23 = 0.3412
    a_32 = a_23

    for i in range(1, p.m['n']+1):
        if s.c[i]['x'] <= 1e-20:  # Prevent math errors from zero log call.
            s.m['g_mix'] = {}
            s.m['g_mix']['t'] = 0.0
            s.m['g_mix']['x'] = s.m['g_mix']['t']
            s.m['g_mix']['ph min'] = 'x'
            return s  # should be = 0 as s2['y']*log(s2['y']) = 1*log(1) = 0


    s.m['g_mix'] = {}
    s.m['g_mix']['t'] = ( s.c[1]['x'] * log(s.c[1]['x'])
                        + s.c[2]['x'] * log(s.c[2]['x'])
                        + s.c[3]['x'] * log(s.c[3]['x'])

                        + s.c[1]['x'] *
                        (t_21 * e**(-a_21 * t_21) * s.c[2]['x']
                         + t_31 * e**(-a_31 * t_31) * s.c[3]['x'])
                        / (s.c[1]['x']
                           + e**(-a_21 * t_21) * s.c[2]['x']
                           + e**(-a_31 * t_31) * s.c[3]['x']
                           )

                        + s.c[2]['x'] *
                        (t_12 * e**(-a_12 * t_12) * s.c[1]['x']
                         + t_32 * e**(-a_32 * t_32) * s.c[3]['x'])
                        / (e**(-a_12 * t_12) * s.c[1]['x']
                           + s.c[2]['x']
                           + e**(-a_32 * t_32) * s.c[3]['x']
                           )

                        + s.c[3]['x'] *
                        (t_13 * e**(-a_13 * t_13) * s.c[1]['x']
                         + t_23 * e**(-a_23 * t_23) * s.c[2]['x'])
                          / (e**(-a_13 * t_13) * s.c[1]['x']
                           + e**(-a_23 * t_23) * s.c[2]['x']
                           + s.c[3]['x']
                           )
                        )

    s.m['g_mix']['x'] = s.m['g_mix']['t']
    s.m['g_mix']['ph min'] = 'x'
    return s


class TestNcompFuncsBin(unittest.TestCase):
    """
    Test ncomp functions
    """
    data = data_handling.ImportData()
    data.comps = ['carbon_dioxide', 'ethane']
    data.phases = ['x', 'y']
    data.eos = 'DWPM'
    data.model = 'Adachi-Lu'
    data.r = None
    data.s = None
    if len(data.comps) > 1:  # multi component simulation.
        # Load all pure dictionaries data.c[i]
        data.load_pure_data()
        # Load VLE and mixture parameter data
        data.load()

        s, p = n_comp_init(data)

    def test_b1(self):
        """
        State and data class defs and funcs
        """
        self.p.m['r'], self.p.m['s'] = 1.0, 1.0
        self.p.m['k'][1][2] = 0.124
        self.p.m['k'][2][1] = self.p.m['k'][1][2]
        self.s.update_state(self.s, self.p, P=24e5, T=263.1,
                            X=[[0.25], [0.25]])

        self.assertTrue(type(self.s.c) is list)
        self.assertTrue(type(self.s.c[1]) is dict)
        self.assertTrue(type(self.p.c) is list)
        self.assertTrue(type(self.p.c[1]) is dict)
        self.assertTrue(type(self.p.m['k'][1][2]) is float )
        part = a_mix_partial_k(self.s, self.p, k=1, phase='x')
        #self.assertTrue(type(part) is float )
        #TODO: Type is float64, find test.

    def test_b2(self):
        """
        Equil. CO2-Ethane Equilibrium
        """
        self.p.m['r'], self.p.m['s'] = 1.0, 1.0
        self.p.m['k'][1][2] = 0.124
        self.p.m['k'][2][1] = self.p.m['k'][1][2]
        Z_0 = numpy.array([0.25])

        X_eq, g_eq, phase_eq = phase_equilibrium_calculation(self.s, self.p,
                                                           g_mix,
                                                           Z_0,
                                                           k=None,
                                                           P=24e5, T=263.1,
                                                           tol=1e-9,
                                                           gtol=1e-2,
                                                           phase_tol=1e-4,
                                                           Print_Results=False,
                                                           Plot_Results=False)

        numpy.testing.assert_allclose([X_eq[0], X_eq[1]],
                                      #[0.28226453, 0.25],
                                      [[0.193647],  [0.308676]],
                                      #More accurate, see plot
                                      rtol=1e-03,
                                      atol=1e-03)


    def test_b3(self):
        """
        Phase sep. CO2-Ethane Equilibrium
        """
        self.p.m['r'], self.p.m['s'] = 1.0, 1.0
        self.p.m['k'][1][2] = 0.124
        self.p.m['k'][2][1] = self.p.m['k'][1][2]
        ph_eq, mph_eq, mph_ph = phase_seperation_detection(g_mix, self.s,
                                                           self.p,
                                                           P=24e5, T=263.1,
                                                           n=100,
                                                           VLE_only=True,
                                                           Plot_Results=False)

        numpy.testing.assert_allclose(mph_eq[0],
                                      [numpy.array([ 0.19469983]),
                                      numpy.array([ 0.30628315])],
                                      rtol=9e-03)


    def test_b4(self):
        """
        Equil. Mitsos et al. (2007) test 1 bin
        """
        Z_0 = numpy.array([0.5])
        self.p.m['Valid phases'] = ['x']
        # noinspection PyTupleAssignmentBalance
        X_eq, g_eq, phase_eq = phase_equilibrium_calculation(self.s, self.p,
                                                           g_x_test_func,
                                                           Z_0,
                                                           k=None,
                                                           tol=1e-9,
                                                           gtol=1e-6,
                                                           phase_tol=1e-5,
                                                           Print_Results=False,
                                                           Plot_Results=False)

        #numpy.testing.assert_allclose([s.m['X_I'][0], s.m['X_II'][0]],

                                      
    def test_b5(self):
        """
        Phase sep. Mitsos et al. (2007) test 1 bin with gtol = 1e-4 to find two
        plane with 2 equilibrium solutions only.
        """
        self.p.m['Valid phases'] = ['x']
        ph_eq, mph_eq, mph_ph = phase_seperation_detection(g_x_test_func,
                                                           self.s, self.p,
                                                           P=101e3, T=300.0,
                                                           n=100,
                                                           gtol=1e-4,
                                                           LLE_only=True,
                                                           Plot_Results=False)
        print ph_eq['x']

        numpy.testing.assert_allclose(ph_eq['x'],
                                      [[numpy.array([ 0.58775493]),
                                        numpy.array([ 0.0045433])],
                                       [numpy.array([ 0.5824251]),
                                        numpy.array([ 0.93304455])
                                        ]],
                                      rtol=9e-02)

    def test_b6(self):
        """
        Phase sep. Mitsos et al. (2007) test 1 bin with gtol = 1e-2 to find
        three equilibrium solutions on the plane.
        """
        self.p.m['Valid phases'] = ['x']
        ph_eq, mph_eq, mph_ph = phase_seperation_detection(g_x_test_func,
                                                           self.s, self.p,
                                                           P=101e3, T=300.0,
                                                           n=100,
                                                           gtol=1e-2,
                                                           LLE_only=True,
                                                           Plot_Results=False)

        numpy.testing.assert_allclose(ph_eq['x'],
                                      numpy.array([[[ 0.591987],
                                                    [ 0.004557],
                                                    [ 0.934769]]])
                                      , rtol=9e-02)


class TestNcompFuncsTern(unittest.TestCase):
    """
    Test ncomp functions
    """
    data = data_handling.ImportData()
    data.comps = ['acetone', 'benzene', 'water']
    data.phases = ['x']
    data.eos = 'DWPM'
    data.model = 'Adachi-Lu'
    data.r = None
    data.s = None
    if len(data.comps) > 1:  # multi component simulation.
        # Load all pure dictionaries data.c[i]
        data.load_pure_data()
        # Load VLE and mixture parameter data
        data.load()

        s, p = n_comp_init(data)
        p.m['r'], p.m['s'] = 1.0, 1.0


    def test_t1(self):
        """
        State and data class definition
        """
        self.assertTrue(type(self.s.c) is list)
        self.assertTrue(type(self.s.c[1]) is dict)
        self.assertTrue(type(self.p.c) is list)
        self.assertTrue(type(self.p.c[1]) is dict)
        self.assertTrue(type(self.p.m['k'][1][2]) is float )

    def test_t2(self):
        """
        Equil. Mitsos et al. (2007) test 2 tern
        """
        Z_0 = numpy.array([0.3, 0.2])
        X_eq, g_eq, phase_eq = phase_equilibrium_calculation(self.s, self.p,
                                                          g_x_test_func2,
                                                          Z_0,
                                                          k=None,
                                                          P=101e3, T=300.0,
                                                          tol=1e-12,
                                                          phase_tol=1e-5,
                                                          n = 100 + 100 * 3,
                                                          Print_Results=False,
                                                          Plot_Results=False)

        # Order phases correctly:
        if X_eq[0][0] < 0.1:
            Ans_X_I = [9.0e-05, 9.99990000e-01]
            Ans_X_II = [0.346864,  0.075812]
        else:
            Ans_X_I = [0.346864,  0.075812]
            Ans_X_II = [9.0e-05, 9.99990000e-01]

        # Note: There is some decimal rounding in the [9.0e-05, 9.9999e-01] ans
        numpy.testing.assert_allclose(X_eq[0],
                                      Ans_X_I,
                                      rtol=1e-02)

        numpy.testing.assert_allclose(X_eq[1],
                                      Ans_X_II,
                                      rtol=1e-02)

def ncomp_suite():
    """
    Gather all the pure tests from this module in a test suite.
    """
    TestNcomp = unittest.TestSuite()
    ncomp_suite1 = unittest.makeSuite(TestNcompFuncsBin)
    ncomp_suite2 = unittest.makeSuite(TestNcompFuncsTern)
    TestNcomp.addTest(ncomp_suite1)
    TestNcomp.addTest(ncomp_suite2)
    return TestNcomp


if __name__ == '__main__':
    TestNcomp = ncomp_suite()
    unittest.TextTestRunner(verbosity=2).run(TestNcomp)

