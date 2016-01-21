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


#def test_all_wrap(TestTgo, DataTests):
def test_all_wrap(test_suite):
    """
    Gather all the data_handling tests from this module in a test suite.
    """
    TestAll = unittest.TestSuite()
    for ts in test_suite:
        TestAll.addTest(ts)
    return TestAll


if __name__ == '__main__':
    TestTgo = tgo_suite()
    DataTests = data_handling_suite()
    TestAll = test_all_wrap((TestTgo, DataTests))
    unittest.TextTestRunner(verbosity=2).run(TestAll)



