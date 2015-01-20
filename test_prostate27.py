import unittest
from prostate27 import Prostate27


class runTest(unittest.TestCase):

    def test_nearest_neighbor(self):
        """
        nearest neighbors are found for individual or list of regions
        """
        p = Prostate27()
        self.assertIn('6a', p.nearest_neighbors(['5a']))
        self.assertIn('11a', p.nearest_neighbors(['15as']))
        self.assertIn('6a', p.nearest_neighbors(['15as']))
        self.assertIn('5p', p.nearest_neighbors(['11p']))
        self.assertIn('12p', p.nearest_neighbors(['5a', '11p']))

    def test_location(self):
        """
        indices corresponding to region
        """
        p = Prostate27()
        self.assertTrue(p.location('12a') == (2, 1, 3))
        self.assertTrue(p.location('1a') == (0, 1, 1))

    def test_nn_ranges(self):
        """
        27 region index bounds not exceeded; AS regions pull all mid-gland
        """
        p = Prostate27()
        self.assertTrue(p.nn_ranges((1, 2, 3)) ==
                        [[0, 1, 2], [1, 2], [2, 3]])
        self.assertTrue(p.nn_ranges((0, 0, 0)) ==
                        [[0, 1], [0, 1], [0, 1, 2, 3]])

    def test_anterior_posterior(self):
        """
        anterior / posterior assignment
        """
        p = Prostate27()
        self.assertTrue(p.anterior_posterior('14as') == 'anterior')
        self.assertTrue(p.anterior_posterior('1p') == 'posterior')
        self.assertTrue(p.anterior_posterior('5a') == 'anterior')

    def test_zone(self):
        """
        peripheral zone / central gland assignment
        """
        p = Prostate27()
        self.assertTrue(p.zone('9a') == 'central gland')
        self.assertTrue(p.zone('11a') == 'central gland')
        self.assertTrue(p.zone('4a') == 'peripheral zone')
        self.assertTrue(p.zone('14as') == 'peripheral zone')


if __name__ == '__main__':
    unittest.main()
