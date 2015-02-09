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
        # self.root = '/home/mlp6/Downloads/invivo/Patient%s' % self.pnum
        self.mri_ios = '%s/MRI_Images/MRI_Index_Lesion.txt' % self.root
        self.hist_lesions = '%s/Histology/HistologyLesions.txt' % self.root
        self.read_mri()
        self.arfi = self.mri  # defined so self.arfi points to self.mri for lesion analysis of MRI
        self.read_histology()
        self.valid_dataset()
        if self.valid:
            self.mri_lesions()
            self.histology_lesions()
            self.check_index_match()
            self.check_benign_match()
            self.check_clin_sig_match()
            self.check_hist_clin_sig_sensitivity()
            
    def read_mri(self):
        """
        read MRI lesion IOS data
        """
        self.mri = {}
        read_index ={} #read MRI index here so not out of order
        
        try:
            with open(self.mri_ios, 'r') as f:
                mriios = f.readlines()

            if 'None' not in mriios[0]:
                for lesion in mriios:
                    read_index['region']=lesion.split(', ')[0]
                    break
                for lesion in mriios:
                    lesion = lesion[:-1]
                    self.mri[lesion.split(', ')[0]] = lesion.split(', ')[1:]
                self.mri['index'] = read_index   
            else:
                self.mri['read'] = 'no lesions read'
        except IOError:
            # print "%s does not exist" % self.mri_ios
            self.mri[None] = None
                      
    def mri_lesions(self):
        """
        characterize all MRI lesions and define MRI index lesion dict based on
        highest IOS
        """
        if self.mri.values()[0] != 'no lesions read':
            try:
                from prostate27 import Prostate27
                prostate = Prostate27()

                index = {}
                index['IOS']      = max(self.mri.values())[0] # extract IOS only
                index['gleason']  = max(self.mri.values())[1]
                index['diameter'] = max(self.mri.values())[2]
                index['extracap'] = max(self.mri.values())[3]
                index['region'] = self.mri['index']['region']
                index['location'] = prostate.anterior_posterior(index['region'])
                index['zone'] = prostate.zone(index['region'])
                self.mri['index'] = index
            except ValueError:
                self.mri['index'] = None

            # specify location and zone for each lesion
            for lesion in self.mri:
                if 'index' not in lesion:
                    self.mri[lesion] = {'IOS': self.mri["%s"      % lesion][0],
                                        'gleason': self.mri["%s"  % lesion][1],
                                        'diameter': self.mri["%s" % lesion][2],
                                        'extracap': self.mri["%s" % lesion][3],
                                        'location':
                                         prostate.anterior_posterior(lesion),
                                         'zone': prostate.zone(lesion)}
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