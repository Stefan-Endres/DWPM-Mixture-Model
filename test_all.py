#!/usr/bin/env python2
"""
#python2 -m unittest -v test_all
#python2 -m unittest -v test_all.AllTestCases.TestFunctions.test_t1
"""

import unittest
import numpy
#from tgo_tests import TestTgo
from tgo_tests import tgo_suite
from data_handling_tests import data_handling_suite


class AllTestCases(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, False)

    #TestFunctions


def test_all_wrap():
    """
    Gather all the data_handling tests from this module in a test suite.
    """
    TestAll = unittest.TestSuite()
    suite1 = unittest.makeSuite(TestTgo)
    suite2 = unittest.makeSuite(DataTests)
    TestAll.addTest(suite1)
    TestAll.addTest(suite2)
    return TestAll


#def suite1():
#    suite1 = unittest.TestSuite()
#    suite1.addTest(MyTestCase)
#    suite1.addTest(TestFunctions)
#    return suite1

if __name__ == '__main__':
    TestTgo = tgo_suite()
    DataTests = data_handling_suite()
    #unittest.TextTestRunner(verbosity=2).run(TestTgo)

    TestAll = data_handling_suite()
    unittest.TextTestRunner(verbosity=2).run(TestAll)



