#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test Run examples:

ex.
$ python2 -m unittest -v data_handling_tests.TestDataLoad
$ python2 -m unittest -v data_handling_tests.TestNewData
"""
import data_handling
import unittest
import numpy

def test_data_type(self, data, d_key, data_type=str):
    """
    Test that the loaded data is the expected object type.

    Parameters
    ----------
    data : class
          Contains the loaded data.

    d_key : d_key contains the string to test.

    data_type : object type to test for
    """
    if data_type is list:
        self.assertIs(type(self.p.m[d_key]), data_type)
    if data_type is float:
        self.assertIs(type(self.p.m[d_key][1]), data_type)

def test_container_string(data_keys, expected_string):
    """
    Tests that the correct string for the components are specified.

    Parameters
    ----------
    data : class
           Contains the loaded data.
    """
    return numpy.array_equal(data_keys, expected_string)


class TestDataLoad(unittest.TestCase):
    """
    Known data loading as expected
    """
    import numpy
    data = data_handling.ImportData()
    data.comps = ['carbon_dioxide', 'ethane']
    data.eos = 'DWPM'
    data.model = 'Adachi-Lu'
    data.phases = ['x','y']
    data.r = None
    data.s = None
    data.k_params = None
    data.load_pure_data()
    data.load()
    data.c[0]['a_c (Pa m6 mol-2)'][0] = '' # Test param caculation
    data.c[1]['b_c (m3 mol-1)'][0] = ''



    crit_ans = numpy.array([0.34985866655, # data.c[0]['a_c (Pa m6 mol-2)'][0]
                            6.36225762119e-05 # data.c[1]['b_c (m3 mol-1)'][0]
                            ])

    p = data_handling.MixParameters()
    p.mixture_parameters(data.VLE, data)
    p.m['n'] = len(data.comps)  # Define system size
    for i in range(p.m['n']):  # Set params for all compounds
        p.parameters(data.c[i])  # Defines p.c[i]

    points_ans = numpy.array([216.59,     # p.c[1]['T'][0]
                              304.59,     # p.c[1]['T'][88]
                              517278.8,   # p.c[1]['P'][0]
                              7381661.9,  # p.c[1]['P'][88]
                              7381661.9,  # p.c[1]['P_c']
                              90.356,     # p.c[2]['T'][0]
                              279.56,     # p.c[2]['T'][88]
                              1.1198,     # p.c[2]['P'][0]
                              2781354.5,  # p.c[2]['P'][88]
                              305.3577,   # p.c[2]['T_c']
                              ])

    points_ans_vle = numpy.array([1.0,   # p.m['x'][1][0]
                                  0.926, # p.m['x'][1][1]
                                  1.0,   # p.m['y'][1][0]
                                  0.822  # p.m['y'][1][1]
                                  ])
    def test_d1(self):
        """Pure data load"""
        points_return = numpy.array([self.p.c[1]['T'][0],
                                     self.p.c[1]['T'][88],
                                     self.p.c[1]['P'][0],
                                     self.p.c[1]['P'][88],
                                     self.p.c[1]['P_c'],
                                     self.p.c[2]['T'][0],
                                     self.p.c[2]['T'][88],
                                     self.p.c[2]['P'][0],
                                     self.p.c[2]['P'][88],
                                     self.p.c[2]['T_c']
                                     ])

        print self.points_ans - points_return
        numpy.testing.assert_array_equal(self.points_ans, points_return)


    def test_d2(self):
        """VLE data load"""
        points_return_vle = numpy.array([self.p.m['x'][1][0],
                                         self.p.m['x'][1][1],
                                         self.p.m['y'][1][0],
                                         self.p.m['y'][1][1]
                                         ])

        numpy.testing.assert_array_equal(self.points_ans_vle,
                                         points_return_vle)

    def test_d3(self):
        """Critical parameter calculation"""
        crit_calc = numpy.array([self.p.c[1]['a_c'],
                                 self.p.c[2]['b_c']
                                ])

        numpy.testing.assert_allclose(crit_calc, self.crit_ans, rtol=1e-03)

class TestNewData(unittest.TestCase):
    """
    Test if newly added data is loading as expected.
    TODO: Add simple argparse/config options to test new data
    """
    import numpy
    data = data_handling.ImportData()

    # TODO: Allow for a system to test argument passes
    data.comps = ['carbon_dioxide', 'ethane']
    data.eos = 'DWPM'
    data.model = 'Adachi-Lu'
    data.phases = ['x','y']
    data.r = None
    data.s = None
    data.k_params = None

    data.load_pure_data() # Using data.comps

    def test_nd1(self):
        """
        New dataset key strings
        """
        expected_string = numpy.array(['P_c (Pa)', 'b_c (m3 mol-1)', 'Z_c',
                   'V_c (m3 mol-1)', 'R (m3 Pa K-1 mol-1)', 'm (Soave)',
                   'virialT', 'T (K)','m (Adachi-Lu)', 'P (Pa)', 'w',
                   'virialB', 'T_c (K)', 'model','a_c (Pa m6 mol-2)',
                   'name'])

        passlist = []
        for i in range(len(self.data.c)):
            passlist.append(test_container_string(self.data.c[i].keys(),
                                                  expected_string))

        self.assertTrue(passlist)

    if len(data.comps) > 1:  # pure component simulation.
        # Load VLE and mixture parameter data
        data.load()
        p = data_handling.MixParameters()
        p.mixture_parameters(data.VLE, data)
        p.m['n'] = len(data.comps)  # Define system size
        print p.m.keys()
        for i in range(p.m['n']):  # Set params for all compounds
            p.parameters(data.c[i])  # Defines p.c[i]

        def test_nd3(self):
            """
            New dataset VLE object types
            """
            for d_key in ['T', 'P', 'x']:
                test_data_type(self, self.p.m, d_key, data_type=list)

            for d_key in ['T', 'P']:
                test_data_type(self, self.p.m, d_key, data_type=float)

            for d_key in self.p.m['Valid phases']:
                self.assertIs(type(self.p.m[d_key][1][0]), float)
                self.assertIs(type(self.p.m[d_key][2][0]), float)


def data_handling_suite():
    """
    Gather all the data_handling tests from this module in a test suite.
    """
    TestData = unittest.TestSuite()
    data_handling_suite1 = unittest.makeSuite(TestDataLoad)
    #data_handling_suite2 = unittest.makeSuite(TestNewData)
    TestData.addTest(data_handling_suite1)
    #TestData.addTest(data_handling_suite2)
    return TestData



if __name__ == '__main__':
    TestData = data_handling_suite()
    unittest.TextTestRunner(verbosity=2).run(TestData)




