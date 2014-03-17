#!/usr/bin/env python

import unittest
from readers import Readers

class TestReaders(unittest.TestCase):

    def setUp(self):
        self.Readers = Readers()
        self.readerName = 'NorKyst800Fake'
        self.Readers.set_point_reader(self.readerName)

    def test_wrong_reader_name(self):
        self.assertRaises(ImportError,
                          self.Readers.set_point_reader('nonexistent_reader'))

    def test_contains_parameter(self):
        self.assertTrue('eastward_sea_water_velocity' in 
        random.shuffle(self.seq)
        self.seq.sort()
        self.assertEqual(self.seq, range(10))

        # should raise an exception for an immutable sequence
        self.assertRaises(TypeError, random.shuffle, (1,2,3))


if __name__ == '__main__':
    unittest.main()
