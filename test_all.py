#!/usr/bin/env python2
"""
#python2 -m unittest -v test_all
#python2 -m unittest -v test_all.AllTestCases.TestFunctions.test_t1
"""

import unittest
#from tgo_tests import tgo_test_suite
from tgo_tests import tgo_suite


class AllTestCases(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, False)

    #TestFunctions

#def suite1():
#    suite1 = unittest.TestSuite()
#    suite1.addTest(MyTestCase)
#    suite1.addTest(TestFunctions)
#    return suite1

if __name__ == '__main__':

    tgo_test_suite=tgo_suite()
    unittest.TextTestRunner(verbosity=2).run(tgo_test_suite)