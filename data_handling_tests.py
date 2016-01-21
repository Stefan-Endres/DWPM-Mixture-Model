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

def test_data_type(data):
    """
    Test that the loaded data is the expected data

    Parameters
    ----------
    data: class
          Contains the loaded data.

    Returns
    -------

    """


class TestDataLoad(unittest.TestCase):
    import numpy
    def test_d1(self):
        """test"""
        data = data_handling.ImportData()
        numpy.testing.assert_array_equal([1],[1])

    def test_d2(self):
        """test"""
        numpy.testing.assert_array_equal([1],[1])

class TestNewData(unittest.TestCase):
    """
    TODO: Add simple argparse/config options to
    """
    import numpy
    def test_nd1(self):
        """test"""
        numpy.testing.assert_array_equal([1],[1])

def data_handling_suite():
    """
    Gather all the data_handling tests from this module in a test suite.
    """
    TestData = unittest.TestSuite()
    data_handling_suite1 = unittest.makeSuite(TestDataLoad)
    data_handling_suite2 = unittest.makeSuite(TestNewData)
    TestData.addTest(data_handling_suite1)
    TestData.addTest(data_handling_suite2)
    return TestData



if __name__ == '__main__':
    TestData = data_handling_suite()
    unittest.TextTestRunner(verbosity=2).run(TestData)




