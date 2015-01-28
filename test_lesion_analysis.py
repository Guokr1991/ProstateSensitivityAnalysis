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
        P13 = LesionAnalysis(13, './testing')
        self.assertFalse(P13.check_index_match())

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
        P13 = LesionAnalysis(13, './testing')
        self.assertFalse(P13.valid)

        P11 = LesionAnalysis(11, './testing')
        self.assertTrue(P11.valid)

    def test_arfi_lesion(self):
        """
        ARFI lesion file parsing
        """
        P11 = LesionAnalysis(11, './testing')
        self.assertTrue(P11.arfi['index']['region'] == '11p')
        self.assertTrue(P11.arfi['index']['IOS'] == 3)
        # these tests depend on Prostate27, but I'll include them too
        self.assertTrue(P11.arfi['index']['location'] == 'posterior')
        self.assertTrue(P11.arfi['index']['zone'] == 'peripheral zone')

    def test_no_arfi_lesion(self):
        """
        no ARFI lesions read
        """
        P10 = LesionAnalysis(10, './testing')
        self.assertFalse('index' in P10.arfi)

#    def test_check_benign_match(self):
#        """
#        ARFI read : benign pathology
#        """
#        P14 = LesionAnalysis(14, './testing')
#        self.assertTrue(P14.benign_match['bph'])
#        self.assertFalse(P14.benign_match['atrophy'])
#        P12 = LesionAnalysis(12, './testing')
#        self.assertFalse(P12.benign_match['bph'])
#        self.assertTrue(P12.benign_match['atrophy'])

if __name__ == '__main__':
    unittest.main()
