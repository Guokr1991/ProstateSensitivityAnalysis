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
        nnset = set(['5p', '11p', '12p'])
        P.histology['index'] = {'region': '11p', 'nn': nnset}

        P.arfi['index'] = {'region': '5p'}
        P.check_index_match()
        self.assertFalse(P.index_match['exact'])
        self.assertTrue(P.index_match['nn'])

        P.arfi['index'] = {'region': '11p'}
        P.check_index_match()
        self.assertTrue(P.index_match['exact'])
        self.assertTrue(P.index_match['nn'])

    def test_valid_dataset(self):
        """
        valid dataset
        """
        P = LesionAnalysis(None)
        P.arfi['index'] = None
        P.histology[None] = None
        P.valid_dataset()
        self.assertFalse(P.valid)
        P.arfi = {'4p': '3'}
        P.valid_dataset()
        self.assertFalse(P.valid)
        P.histology = {'region': '11p'}
        P.valid_dataset()
        self.assertTrue(P.valid)

if __name__ == '__main__':
    unittest.main()
