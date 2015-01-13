class LesionAnalysis:
    """
    class for the prostate lesion analysis b/w ARFI and histology
    """

    def __init__(self, pnum):
        self.pnum = pnum
        self.root = '/luscinia/ProstateStudy/invivo/Patient%s' % self.pnum
        # self.root = '/home/mlp6/Downloads/invivo/Patient%s' % self.pnum
        self.arfi_ios = '%s/ARFI_Index_Lesion_IOS.txt' % self.root
        self.hist_lesions = '%s/Histology/HistologyLesions.txt' % self.root
        self.read_arfi()
        self.read_histology()
        self.valid_dataset()
        if self.valid_dataset:
            self.arfi_lesions()
            self.histology_lesions()
            self.check_index_match()
            self.check_benign_match()
            self.check_clin_sig_match()
            self.check_hist_clin_sig_sensitivity()

    def read_arfi(self):
        """
        read ARFI lesion IOS data
        """
        self.arfi = {}
        try:
            with open(self.arfi_ios, 'r') as f:
                arfiios = f.readlines()

            if 'None' not in arfiios:
                for lesion in arfiios:
                    lesion = lesion[:-1]
                    self.arfi[lesion.split(', ')[0]] = lesion.split(', ')[1]
            else:
                self.arfi['read'] = 'no lesions read'
        except IOError:
            # print "%s does not exist" % self.arfi_ios
            self.arfi[None] = None

    def read_histology(self):
        """
        head histology pca, atrophy and bph lesions
        """
        self.histology = {}
        try:
            with open(self.hist_lesions, 'r') as f:
                histread = f.readlines()

            for lesion in histread:
                lesion = lesion[:-1]
                # make sure we hav)e a properly-formatted histology file
                if not any([x in lesion for x in ['pca', 'bph', 'atrophy']]):
                    print "WARNING: Malformed histology lesion file (P%s)." % \
                        self.pnum
                # there can be multiple pca lesions
                if 'pca' in lesion:
                    if 'pca' not in self.histology:
                        self.histology['pca'] = [lesion.split(', ')[1:]]
                    else:
                        self.histology['pca'].append(lesion.split(', ')[1:])
                else:
                    self.histology[lesion.split(', ')[0]] = \
                        lesion.split(', ')[1:]
        except IOError:
            # print "%s does not exist" % self.hist_lesions
            self.histology[None] = None

    def arfi_lesions(self):
        """
        characterize all ARFI lesions and define ARFI index lesion dict based on
        highest IOS
        """
        if self.arfi.values()[0] != 'no lesions read':
            try:
                from prostate27 import Prostate27
                prostate = Prostate27()

                index = {}
                index['IOS'] = max(self.arfi.values())
                index['region'] = [x for x, y in self.arfi.items() if
                                   y == index['IOS']][0]
                index['location'] = prostate.anterior_posterior(index['region'])
                index['zone'] = prostate.zone(index['region'])
                self.arfi['index'] = index
            except ValueError:
                self.arfi['index'] = None

            # specify location and zone for each lesion
            for lesion in self.arfi:
                if 'index' not in lesion:
                    self.arfi[lesion] = {'IOS': self.arfi["%s" % lesion],
                                         'location':
                                         prostate.anterior_posterior(lesion),
                                         'zone': prostate.zone(lesion)}

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
                    index['nn'] = p.nearest_neighbors(index['region'])
                    index['location'] = \
                        p.anterior_posterior(index['region'])
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
                self.histology['pca'][lcnt].append(p.anterior_posterior(les[0]))
                self.histology['pca'][lcnt].append(p.zone(les[0]))
        except KeyError:
            print "No PCA lesion"
            self.histology['index'] = None

    def clin_sig(self, volume, gleason):
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
        self.clin_sig_match = []

        # find histology nearest neighbors for all clinically sig lesions
        for i in self.histology['pca']:
            if i[3] == 'ClinSig':
                hist_nnset.update(p.nearest_neighbors(i[0]))

        for lesion_region in self.arfi:
            if 'index' not in lesion_region and 'read' not in lesion_region:
                if lesion_region in hist_nnset:
                    self.clin_sig_match.append([True, self.arfi[lesion_region]['location']])
                else:
                    self.clin_sig_match.append([False, self.arfi[lesion_region]['location']])

    def check_hist_clin_sig_sensitivity(self):
        """
        check if ARFI detected clinically-significant lesions
        """
        from prostate27 import Prostate27
        p = Prostate27()

        self.clin_sig_sensitivity = []
        for i in self.histology['pca']:
            if i[3] == 'ClinSig':
                nnset = p.nearest_neighbors(i[0])
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
        self.benign_match = {}
        for benign in ['atrophy', 'bph']:
            try:
                if self.histology['index'] is None:
                    self.benign_match[benign] = False
                else:
                    benign_regions = self.histology[benign]
                    if any([x in benign_regions for x in
                       self.arfi.keys()]) and \
                       self.histology['index']['region'] not in \
                       benign_regions:
                        self.benign_match[benign] = True
                    else:
                        self.benign_match[benign] = False
            except KeyError:
                self.benign_match[benign] = False

    def valid_dataset(self):
        """
        check if this is a valid dataset to include in the sensitivity analysis
        """
        if None in self.arfi or None in self.histology:
            self.valid_dataset = False
        else:
            self.valid_dataset = True

    def __str__(self):
        """
        print analysis results
        """

        s = ['================= PATIENT %s =================' % self.pnum]
        if self.valid_dataset is False:
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
