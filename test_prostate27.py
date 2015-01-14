import unittest
from prostate27 import Prostate27


class runTest(unittest.TestCase):

    def test_nearest_neighbor(self):
        p = Prostate27()
        self.assertIn('6a', p.nearest_neighbors('5a'))
        self.assertIn('11a', p.nearest_neighbors('15as'))
        self.assertIn('5p', p.nearest_neighbors('11p'))

if __name__ == '__main__':
    unittest.main()
