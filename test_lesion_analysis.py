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
        P11 = LesionAnalysis(11, './testing')
        self.assertFalse(P11.index_match['exact'])
        self.assertTrue(P11.index_match['nn'])

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

    def test_arfi_lesion(self):
        """
        ARFI lesion file parsing
        """
        P = LesionAnalysis(None)
        P.arfi = {'11p': '3', '5a': '1'}
        P.arfi_lesions()
        self.assertTrue(P.arfi['index']['region'] == '11p')
        self.assertTrue(P.arfi['index']['IOS'] == '3')
        # these tests depend on Prostate27, but I'll include them too
        self.assertTrue(P.arfi['index']['location'] == 'posterior')
        self.assertTrue(P.arfi['index']['zone'] == 'peripheral zone')

    def test_no_arfi_lesion(self):
        """
        no ARFI lesions read
        """
        P = LesionAnalysis(None)
        P.arfi = {'read': 'no lesions read'}
        self.assertNotIn('index', P.arfi)

if __name__ == '__main__':
    unittest.main()
