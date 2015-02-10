class LesionAnalysis:
    """
    class for the prostate lesion analysis b/w ARFI and histology
    """

    def __init__(self, pnum, root='/luscinia/ProstateStudy/invivo'):
        self.pnum = pnum
        self.root = '%s/Patient%s' % (root, self.pnum)
        self.histology = self.read_json('%s/Histology/HistologyLesions.json' %
                                        self.root)
        self.arfi = self.read_json('%s/ARFI_Lesions.json' % self.root)
        self.valid = self.valid_dataset(self.arfi, self.histology)
        if self.valid:
            self.arfi_lesions()
            self.histology_lesions()
            self.check_index_match()
            self.check_benign_match()
            self.check_clin_sig_match()
            self.check_hist_clin_sig_sensitivity()

    @staticmethod
    def read_json(json_input):
        """
        read JSON file

        INPUT: json_input (str) - location of JSON input file
        OUTPUT: modality_dict - modality dictionary with lesion information
        """
        import json

        modality_dict = {}
        try:
            with open(json_input, 'r') as f:
                modality_dict = json.load(f)
        except IOError:
            # print "%s does not exist" % json_input
            modality_dict[None] = None

        return modality_dict

    def arfi_lesions(self):
        """
        characterize all ARFI lesions and define ARFI index lesion
        """
        if self.arfi.values()[0] != 'no lesions read':
            try:
                from prostate27 import Prostate27
                p = Prostate27()

                index = {}
                for lesion in self.arfi['lesions']:
                    if lesion['index'] is True:
                        index['IOS'] = lesion['IOS']
                        index['region'] = lesion['region']
                        index['location'] = \
                            p.anterior_posterior([lesion['region']])
                        index['zone'] = p.zone(lesion['region'])
                        self.arfi['index'] = index
                        break
            except ValueError:
                self.arfi['index'] = None

            # specify location and zone for each lesion
            for n, lesion in enumerate(self.arfi['lesions']):
                self.arfi['lesions'][n].update({'location':
                                                p.anterior_posterior(
                                                    [lesion['region']]),
                                                'zone': p.zone(
                                                    lesion['region'])})

    def histology_lesions(self):
        """
        characterize all histology lesions in dict, including index, nearest
        neighbors, clinical significance, location and zone
        """
        from prostate27 import Prostate27
        p = Prostate27()

        try:
            index = {}
            for lesion in self.histology['pca']:
                if lesion['index'] is True:
                    index['region'] = lesion['region']
                    index['Gleason'] = lesion['Gleason']
                    index['nn'] = p.nearest_neighbors([lesion['region']])
                    index['location'] = p.anterior_posterior(
                        [lesion['region']])
                    index['zone'] = p.zone(lesion['region'])
                    self.histology['index'] = index
                    break
            
            
            
            for n, les in enumerate(self.histology['pca']):
                if self.clin_sig(les['volume_cc'], les['Gleason']):
                    self.histology['pca'][n].update(
                        {'Clinically Significant': True})
                else:
                    self.histology['pca'][n].update(
                        {'Clinically Significant': False})
                    self.histology['pca'][n].update(
                        {'location': p.anterior_posterior([les['region']])})
                    self.histology['pca'][n].update(
                        {'zone': p.zone(les['region'])})
                                     
        except ValueError:
            print "No PCA lesion"
            self.histology['index'] = None
        
        if 'index' in self.histology:
            print "Histology Index Exists, P%i" % self.pnum
        else:
            print "No Histology Index Lesion, P%i" % self.pnum
            self.histology['index']=None

    @staticmethod
    def clin_sig(volume, gleason):
        """
        determine lesion clinical significance
        """
        volume = float(volume)
        gleason = int(gleason)
        if (gleason >= 7) or (gleason >= 6 and volume >= 500):
            return True
        else:
            return False

    def check_index_match(self):
        """
        check for an exact and NN matches b/w ARFI and histology index lesions
        """
        self.index_match = {}
        try:
            if self.histology['index'] is None:
                self.index_match['exact'] = False
                self.index_match['nn'] = False
            else:
                # check for an exact match
                if self.arfi['index']['region'] == \
                   self.histology['index']['region']:
                    self.index_match['exact'] = True
                else:
                    self.index_match['exact'] = False

                # check for a NN match
                if self.arfi['index']['region'] in \
                   self.histology['index']['nn']:
                    self.index_match['nn'] = True
                else:
                    self.index_match['nn'] = False
        except KeyError:
            self.index_match['exact'] = False
            self.index_match['nn'] = False

    def check_clin_sig_match(self):
        """
        check if ARFI reads are clinically-significant lesions
        """
        from prostate27 import Prostate27
        p = Prostate27()

        hist_nnset = set()
        histnonsigset = set()
        self.clin_sig_match = []

        # find histology nearest neighbors for all clinically sig lesions
        for i in self.histology['pca']:
            if i['Clinically Significant']:
                hist_nnset.update(p.nearest_neighbors([i['region']]))
            else:
                histnonsigset.update(p.nearest_neighbors([i['region']]))

        self.false_positive = []
        try:
            for lesion in self.arfi['lesions']:
                lesion_region = lesion['region']
                if lesion_region in hist_nnset:
                    self.clin_sig_match.append(
                        [True, lesion['location']])
                else:
                    self.clin_sig_match.append(
                        [False, lesion['location']])
                    # check if something else exists in this region, including:
                    # non-clinically significant PCA, atrophy, BPH
                    try:
                        if lesion_region in self.histology['atrophy']['regions']:
                            self.false_positive.append('atrophy')
                        elif lesion_region in self.histology['bph']['regions']:
                            self.false_positive.append('bph')
                        elif lesion_region in histnonsigset:
                            self.false_postive.append('pca')
                        else:
                            self.false_positive = None
                    except KeyError:
                        self.false_positive = None
        except KeyError:
            self.false_positive = None

    def check_hist_clin_sig_sensitivity(self):
        """
        check if ARFI detected clinically-significant lesions
        """
        from prostate27 import Prostate27
        p = Prostate27()

        self.clin_sig_sensitivity = []
        try:
            for lesion in self.histology['pca']:
                if lesion['Clinically Significant']:
                    nnset = p.nearest_neighbors([lesion['region']])
                    for arfi in self.arfi['lesions']:
                        lesion_region = arfi['region']
                        if ('index' not in lesion_region) and \
                           ('read' not in lesion_region):
                            if lesion_region in nnset:
                                self.clin_sig_sensitivity.append(
                                    [True,
                                     p.anterior_posterior([lesion['region']])])
                            else:
                                self.clin_sig_sensitivity.append(
                                    [False,
                                     p.anterior_posterior([lesion['region']])])
        except KeyError:
            self.clin_sig_sensitivity = None

    def check_benign_match(self):
        """
        check if atrophy or bph is present in the exact or nearest neighbor
        region to an ARFI lesion
        """
        from prostate27 import Prostate27
        p = Prostate27()

        self.benign_match = {}
        for benign in ['atrophy', 'bph']:
            try:
                if self.histology['index'] is None:
                    self.benign_match[benign] = False
                else:
                    benign_regions = self.histology[benign]['regions']
                    arfi_regions_nn = p.nearest_neighbors(
                        [x['region'] for x in self.arfi['lesions']])
                    if (any([x in benign_regions for x in arfi_regions_nn])) \
                        and (self.histology['index']['region'] not in
                             benign_regions):
                        self.benign_match[benign] = True
                    else:
                        self.benign_match[benign] = False
            except KeyError:
                self.benign_match[benign] = False

    @staticmethod
    def valid_dataset(modality, histology):
        """
        check if valid dataset to include in the sensitivity analysis

        INPUT: modality (dict)
               histology (dict)
        OUTPUT: True/False
        """
        if None in modality or None in histology:
            return False
        else:
            return True

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
            s.append('ARFI:Atrophy Match:\t\t\t%s' %
                     self.benign_match['atrophy'])
            s.append('ARFI:BPH Match:\t\t\t\t%s' % self.benign_match['bph'])

        return '\n'.join(s)
