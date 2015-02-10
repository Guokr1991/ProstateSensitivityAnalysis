import unittest
from lesion_analysis import LesionAnalysis
from mr_analysis import MRAnalysis


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

    def test_check_benign_match(self):
        """
        ARFI read : benign pathology
        """
        P14 = LesionAnalysis(14, './testing')
        self.assertTrue(P14.benign_match['bph'])
        self.assertFalse(P14.benign_match['atrophy'])
        P12 = LesionAnalysis(12, './testing')
        self.assertFalse(P12.benign_match['bph'])
        self.assertTrue(P12.benign_match['atrophy'])
    
#    def test_histology_ECE(self):
#        """
#        Histology JSON ECE parsing
#        """
#        P10 = LesionAnalysis(10, './testing')
#        self.assertTrue(P10.histology['index']['Staging']='T3a')
#        self.assertTrue(P10.histology['index']['ECE_extent']='Established')
#        P10 = MRAnalysis(10, './testing') 
#        self.assertTrue(P10.histology['index']['Staging']='T3a')
#        self.assertTrue(P10.histology['index']['ECE_extent']='Established')  
#        P13 = MRAnalysis(13, './testing') 
#        self.assertTrue(P13.histology['index']['Staging']='T2x')
#        self.assertTrue(P13.histology['index']['ECE_extent']='None')                 
#        
#    def test_valid_MR_dataset(self):
#        """#    def test_valid_MR_dataset(self):
#        """
#        valid MRI dataset
#        """
#        P13 = MRAnalysis(13, './testing')
#        self.assertFalse(P13.valid)
#
#        P11 = MRAnalysis(11, './testing')
#        self.assertTrue(P11.valid)
#
#        P14 = MRAnalysis(14, './testing')
#        self.assertTrue(P14.valid)
#
#    def test_mri_lesion(self):
#        """
#        MRI file parsing for index lesion
#        """
#        P11 = MRAnalysis(11, './testing')
#
#
#        self.assertTrue(P11.mri['index']['region'] == '9p')
#        self.assertTrue(P11.mri['index']['IOS'] == 3)
#        # these tests depend on Prostate27, but I'll include them too
#        self.assertTrue(P11.mri['index']['location'] == 'posterior')
#        self.assertTrue(P11.mri['index']['zone'] == 'peripheral zone')
#        # these tests are only relevant for MRI data
#        self.assertTrue(P11.mri['index']['Gleason'] == 7)
#        self.assertTrue(P11.mri['index']['diameter_mm'] == 8)
#        self.assertTrue(P11.mri['index']['ECE'] == false)
#        self.assertTrue(P11.mri['index']['lesion_number'] == 1)
#        # test for case with multiple IOS3 lesions (top line in text file is index)
#        P14 = MRAnalysis(14, './testing')
#        self.assertTrue(P14.mri['index']['region'] == '10p')
#        self.assertTrue(P14.mri['index']['IOS'] == 3)
#        self.assertFalse(P14.mri['index']['region'] == '3p')
#
#    def test_MR_no_index_match(self):
#        """
#        no histology index lesion using MRAnalysis class
#        """
#        P13 = MRAnalysis(13, './testing')
#        self.assertFalse(P13.check_index_match())
#
#    def test_histology_mri_index_match(self):
#        """
#        histology / MRI index match
#        """
#        P11 = MRAnalysis(11, './testing')
#        self.assertFalse(P11.index_match['exact'])
#        self.assertTrue(P11.index_match['nn'])
#
#        P14 = MRAnalysis(14, './testing')
#        self.assertFalse(P14.index_match['exact'])
#        self.assertTrue(P14.index_match['nn'])
#        valid MRI dataset
#        """
#        P13 = MRAnalysis(13, './testing')
#        self.assertFalse(P13.valid)
#
#        P11 = MRAnalysis(11, './testing')
#        self.assertTrue(P11.valid)
#
#        P14 = MRAnalysis(14, './testing')
#        self.assertTrue(P14.valid)
#
#    def test_mri_lesion(self):
#        """
#        MRI file parsing for index lesion
#        """
#        P11 = MRAnalysis(11, './testing')
#        self.assertTrue(P11.mri['index']['region'] == '9p')
#        self.assertTrue(P11.mri['index']['IOS'] == 3)
#        # these tests depend on Prostate27, but I'll include them too
#        self.assertTrue(P11.mri['index']['location'] == 'posterior')
#        self.assertTrue(P11.mri['index']['zone'] == 'peripheral zone')
#        # these tests are only relevant for MRI data
#        self.assertTrue(P11.mri['index']['Gleason'] == 7)
#        self.assertTrue(P11.mri['index']['diameter_mm'] == 8)
#        self.assertTrue(P11.mri['index']['ECE'] == false)
#        self.assertTrue(P11.mri['index']['lesion_number'] == 1)
#        # test for case with multiple IOS3 lesions (top line in text file is index)
#        P14 = MRAnalysis(14, './testing')
#        self.assertTrue(P14.mri['index']['region'] == '10p')
#        self.assertTrue(P14.mri['index']['IOS'] == 3)
#        self.assertFalse(P14.mri['index']['region'] == '3p')
#
#    def test_MR_no_index_match(self):
#        """
#        no histology index lesion using MRAnalysis class
#        """
#        P13 = MRAnalysis(13, './testing')
#        self.assertFalse(P13.check_index_match())
#
#    def test_histology_mri_index_match(self):
#        """
#        histology / MRI index match
#        """
#        P11 = MRAnalysis(11, './testing')
#        self.assertFalse(P11.index_match['exact'])
#        self.assertTrue(P11.index_match['nn'])
#
#        P14 = MRAnalysis(14, './testing')
#        self.assertFalse(P14.index_match['exact'])
#        self.assertTrue(P14.index_match['nn'])
  
if __name__ == '__main__':
    unittest.main()
