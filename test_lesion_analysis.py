import unittest
from lesion_analysis import LesionAnalysis


class runTest(unittest.TestCase):

    def test_clin_sig(self):
        """
        histology clinical significance
        """
        self.assertTrue(LesionAnalysis.clin_sig(100, 8))
        self.assertFalse(LesionAnalysis.clin_sig(100, 5))
        self.assertFalse(LesionAnalysis.clin_sig(100, 6))
        self.assertTrue(LesionAnalysis.clin_sig(700, 6))

if __name__ == '__main__':
    unittest.main()
