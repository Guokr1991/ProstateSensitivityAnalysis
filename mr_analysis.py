# -*- coding: utf-8 -*-
"""
Created on Mon Jan 26 10:44:57 2015

@author: tjg17
"""
from lesion_analysis import LesionAnalysis
class MRAnalysis(LesionAnalysis):
    """
    class for the prostate lesion analysis b/w MRI and histology
    """
    def __init__(self, pnum, root='/luscinia/ProstateStudy/invivo'):
        self.pnum = pnum
        self.root = '%s/Patient%s' % (root, self.pnum)
        self.histology = self.read_json('%s/Histology/HistologyLesions.json' %
                                        self.root)
        self.arfi = self.read_json('%s/ARFI_Lesions.json' % self.root) # read in arfi
        self.mri = self.read_json('%s/MRI_Images/MRILesions.json' % self.root) # read in mri
        self.valid = self.valid_dataset(self.mri, self.histology)
        if self.valid:
            self.mri_lesions()
            self.arfi_lesions() # keep arfi analysis in MRanalysis class
            self.histology_lesions()
            self.check_index_match()
            self.check_benign_match()
            self.check_clin_sig_match()
            self.check_hist_clin_sig_sensitivity()
                                 
    def mri_lesions(self):
        """
        characterize all MRI lesions and define MRI index lesion
        """
        if self.mri.values()[0] != 'no lesions read':
            try:
                from prostate27 import Prostate27
                p = Prostate27()

                index = {}
                for lesion in self.mri['lesions']:
                    if lesion['index'] is True:
                        index['IOS'] = lesion['IOS']
                        index['region'] = lesion['region']
                        index['location'] = \
                            p.anterior_posterior([lesion['region']])
                        index['zone'] = p.zone(lesion['region'])
                        index['Gleason'] = lesion['Gleason']
                        index['diameter_mm'] = lesion['diameter_mm']
                        index['ECE'] = lesion['extracap']
                        index['lesion_number'] = lesion['lesion_number']
                        self.mri['index'] = index
                        break
            except ValueError:
                self.mri['index'] = None

            # specify location and zone for each lesion
            #for n, lesion in enumerate(self.mri['lesions']):
            #    self.mri['lesions'][n].update({'location': p.anterior_posterior([lesion['region']]),
             #                                   'zone': p.zone(lesion['region'])})                                     'zone': prostate.zone(lesion)}
    def __str__(self):
        """
        print analysis results
        """

        s = ['================= PATIENT %s =================' % self.pnum]
        if self.valid is False:
            s.append('Incomplete dataset; not included in analysis.')
        else:
            s.append('Index lesion EXACT match:\t\t%s' %
                     self.index_match['exact'])
            s.append('Index lesion NEAREST NEIGHBOR match:\t%s' %
                     self.index_match['nn'])
            s.append('MRI:Atrophy Match:\t\t\t%s' %
                     self.benign_match['atrophy'])
            s.append('MRI:BPH Match:\t\t\t%s' % self.benign_match['bph'])

          
        return '\n'.join(s)