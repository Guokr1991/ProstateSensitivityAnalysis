import unittest
from prostate27 import Prostate27


class runTest(unittest.TestCase):

    def test_nearest_neighbor(self):
        p = Prostate27()
        self.assertIn('6a', p.nearest_neighbors('5a'))
        self.assertIn('11a', p.nearest_neighbors('15as'))
        self.assertIn('6a', p.nearest_neighbors('15as'))
        self.assertIn('5p', p.nearest_neighbors('11p'))

    def test_location(self):
        p = Prostate27()
        self.assertTrue(p.location('12a') == (2, 1, 3))
        self.assertTrue(p.location('1a') == (0, 1, 1))

    def test_nn_ranges(self):
        p = Prostate27()
        self.assertTrue(p.nn_ranges((1, 2, 3)) ==
                        [[0, 1, 2], [1, 2], [2, 3]])
        self.assertTrue(p.nn_ranges((0, 0, 0)) ==
                        [[0, 1], [0, 1], [0, 1, 2, 3]])

if __name__ == '__main__':
    unittest.main()
