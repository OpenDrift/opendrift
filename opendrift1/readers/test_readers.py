#!/usr/bin/env python

import unittest
from readers import Readers

class TestReaders(unittest.TestCase):

    def setUp(self):
        self.Readers = Readers()

    def test_wrong_reader_name(self):
        self.assertRaises(ImportError,
                          self.Readers.set_point_reader,
                          'nonexistent_reader', 'some_parameter')

    def test_missing_parameters(self):
        self.assertRaises(ValueError,
                          self.Readers.set_point_reader,
                          'norkystFake', 'nonexistent_parameter')


if __name__ == '__main__':
    unittest.main()
