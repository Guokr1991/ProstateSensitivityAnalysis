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
        #    self.histology_lesions()
        #    self.check_index_match()
        #    self.check_benign_match()
        #    self.check_clin_sig_match()
        #    self.check_hist_clin_sig_sensitivity()

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
            print "%s does not exist" % json_input
            modality_dict[None] = None

        return modality_dict

    def arfi_lesions(self):
        """
        characterize all ARFI lesions and define ARFI index lesion dict based
        on highest IOS
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
                                                p.anterior_posterior([lesion['region']]),
                                                'zone': p.zone(lesion['region'])})

    def histology_lesions(self):
        """
        characterize all histology lesions in dict, including index, nearest
        neighbors, clinical significance, location and zone
        """
        from prostate27 import Prostate27
        p = Prostate27()

        try:
            # find max Gleason score, then max volume with that max Gleason
            # preference given to first entry in file
            for lesion in self.histology['pca']:
                # Gleason scores that weren't reported were recorded as '-1'
                # these will be set to 6 for now
                if lesion[2] == '-1':
                    lesion[2] = '6'
                # check if the lesion is clinically significant
                # accept the first row if that is significant since clinical
                # reads of index were entered first
                if lesion[2] >= 7 or (lesion[2] >= 6 and lesion[1] >= 500):
                    index = {}
                    index['region'] = lesion[0]
                    index['Gleason'] = lesion[2]
                    index['nn'] = p.nearest_neighbors([index['region']])
                    index['location'] = \
                        p.anterior_posterior([index['region']])
                    index['zone'] = p.zone(index['region'])
                    self.histology['index'] = index
                    break
            if not self.histology['index']:
                print "No clinically-significant PCA lesion"
                self.histology['index'] = None
            for lcnt, les in enumerate(self.histology['pca']):
                if self.clin_sig(les[1], les[2]):
                    self.histology['pca'][lcnt].append('ClinSig')
                else:
                    self.histology['pca'][lcnt].append('NotClinSig')
                self.histology['pca'][lcnt].append(p.anterior_posterior([les[0]]))
                self.histology['pca'][lcnt].append(p.zone(les[0]))
        except KeyError:
            print "No PCA lesion"
            self.histology['index'] = None

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
            if i[3] == 'ClinSig':
                hist_nnset.update(p.nearest_neighbors([i[0]]))
            if i[3] == 'NotClinSig':
                histnonsigset.update(p.nearest_neighbors([i[0]]))

        self.false_positive = []
        for lesion_region in self.arfi:
            if 'index' not in lesion_region and 'read' not in lesion_region:
                if lesion_region in hist_nnset:
                    self.clin_sig_match.append([True, self.arfi[lesion_region]['location']])
                else:
                    self.clin_sig_match.append([False, self.arfi[lesion_region]['location']])
                    # check if something else exists in this region, including:
                    # non-clinically significant PCA, atrophy, BPH
                    try:
                        if lesion_region in self.histology['atrophy']:
                            self.false_posivite.append('atrophy')
                        elif lesion_region in self.histology['bph']:
                            self.false_posivite.append('bph')
                        elif lesion_region in histnonsigset:
                            self.false_postive.append('pca')
                        else:
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
        for i in self.histology['pca']:
            if i[3] == 'ClinSig':
                nnset = p.nearest_neighbors([i[0]])
                for lesion_region in self.arfi:
                    if 'index' not in lesion_region and 'read' not in lesion_region:
                        if lesion_region in nnset:
                            self.clin_sig_sensitivity.append([True, i[4]])
                        else:
                            self.clin_sig_sensitivity.append([False, i[4]])

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
                    benign_regions = self.histology[benign]
                    arfi_regions_nn = p.nearest_neighbors(self.arfi.keys())
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
