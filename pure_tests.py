#!/usr/bin/env python
import unittest
import data_handling
import main
import numpy

class TestPureFuncs(unittest.TestCase):
    """
    Test known pure results.
    """
    data = data_handling.ImportData()
    data.comps = ['ethane']
    data.eos = 'DWPM'
    data.load_pure_data()
    data.save_pure = False  # Force test not to save.
    p1_ans1 = 0.533365967206
    p1_ans2 = 6.36225762119e-05
    p2_ans = 0.641931479217
    p3_ans = 0.701094516008
    p4_ans = [1046007.02038, 0.00160041689843, 8.83846522926e-05]

    def test_p1(self):
        """
        Critical parameters
        """
        self.data.c[0]['a_c (Pa m6 mol-2)'][0] = ''
        self.data.c[0]['b_c (m3 mol-1)'][0] = ''
        s, p = main.pure_sim(self.data)
        
        # Redefine for other tests
        self.data.c[0]['a_c (Pa m6 mol-2)'] =[0.533365967206]
        self.data.c[0]['b_c (m3 mol-1)'] = [6.36225762119e-05]

        numpy.testing.assert_allclose([self.p1_ans1, self.p1_ans2],
                                      [p['a_c'], p['b_c']], rtol=1e-03)

    def test_p2(self):
        """
        Soave model optimisation
        """

        self.data.model = 'Soave'
        self.data.c[0]['model'] = 'Soave'
        self.data.c[0]['m (Soave)'][0] = ''
        s, p = main.pure_sim(self.data)
        numpy.testing.assert_allclose([self.p2_ans],[p['m'][0]], rtol=1e-03)

    def test_p3(self):
        """
        Adachi-Lu model optimisation
        """
        self.data.model = 'Adachi-Lu'
        self.data.c[0]['model'] = 'Adachi-Lu'
        self.data.c[0]['m (Adachi-Lu)'][0] = ''
        s, p = main.pure_sim(self.data)
        numpy.testing.assert_allclose([self.p3_ans],[p['m'][0]], rtol=1e-02)

    def test_p4(self):
        """
        Phase equilibrium at fixed temperature
        """
        self.data.T = 243.01
        s, p = main.pure_sim(self.data)
        p4 = [s['P_sat'], s['V_v'], s['V_l']]
        numpy.testing.assert_allclose(self.p4_ans, p4, rtol=1e-02)

    def test_p5(self): #TODO: Add main code and test results
        """
        Phase equilibrium at fixed pressure
        """
        self.data.P = 1046007.02038
        s, p = main.pure_sim(self.data)
        numpy.testing.assert_allclose(1, 1, rtol=1e-02)

class TestNewPure(unittest.TestCase):
    def test_pn1(self):
        self.assertEqual(True, True)

def pure_suite():
    """
    Gather all the pure tests from this module in a test suite.
    """
    TestPure = unittest.TestSuite()
    pure_suite1 = unittest.makeSuite(TestPureFuncs)
    TestPure.addTest(pure_suite1)
    return TestPure


if __name__ == '__main__':
    TestPure = pure_suite()
    unittest.TextTestRunner(verbosity=2).run(TestPure)

