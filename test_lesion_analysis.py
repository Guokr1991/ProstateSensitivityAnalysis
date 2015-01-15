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

    def test_no_index_match(self):
        """
        no histology index lesion
        """
        P = LesionAnalysis(None)
        P.histology['index'] = None
        self.assertFalse(P.check_index_match())

    def test_histology_arfi_index_match(self):
        """
        histology / ARFI index match
        """
        P = LesionAnalysis(None)
        P.histology['index'] = {'region': '5p', 'nn': set()}
        P.arfi['index'] = {'region': '5p'}
        self.assertTrue(P.check_index_match())

if __name__ == '__main__':
    unittest.main()
