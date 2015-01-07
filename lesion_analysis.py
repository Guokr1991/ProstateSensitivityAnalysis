class LesionAnalysis:
    """
    class for the prostate lesion analysis b/w ARFI and histology
    """

    def __init__(self, pnum):
        self.pnum = pnum
        self.root = '/luscinia/ProstateStudy/invivo/Patient%s' % self.pnum
        self.arfi_ios = '%s/ARFI_Index_Lesion_IOS.txt' % self.root
        self.hist_lesions = '%s/Histology/HistologyLesions.txt' % self.root
        self.read_arfi()
        self.read_histology()
        self.valid_dataset()
        self.check_index_exact_match()
        self.check_index_nn_match()
        self.check_arfi_benign_match('atrophy')
        self.check_arfi_benign_match('bph')

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
            #print "%s does not exist" % self.arfi_ios
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
                    print "WARNING: Malformed histology lesion file (P%s)." % self.pnum
                # there can be multiple pca lesions
                if 'pca' in lesion:
                    if 'pca' not in self.histology:
                        self.histology['pca'] = [lesion.split(', ')[1:]]
                    else:
                        self.histology['pca'].append(lesion.split(', ')[1:])
                else:
                    self.histology[lesion.split(', ')[0]] = lesion.split(', ')[1:]
        except IOError:
            #print "%s does not exist" % self.hist_lesions
            self.histology[None] = None

    def nearest_neighbor(self, region):
        """
        extract the set of nearest neighbor regions; 27 regions based on INSERT
        CITATION HERE

        INPUT: region (string) - find nearest neighbors around this region
        """

        prostate27roi = [[0 for AP in range(3)] for BA in range(3)]
        prostate27roi[0][0] = ['13as' for LAT in range(4)]
        prostate27roi[0][1] = ['2a', '1a', '7a', '8a']
        prostate27roi[0][2] = ['2p', '1p', '7p', '8p']
        prostate27roi[1][0] = ['14as' for LAT in range(4)]
        prostate27roi[1][1] = ['4a', '3a', '9a', '10a']
        prostate27roi[1][2] = ['4p', '3p', '9p', '10p']
        prostate27roi[2][0] = ['15as' for LAT in range(4)]
        prostate27roi[2][1] = ['6a', '5a', '11a', '12a']
        prostate27roi[2][2] = ['6p', '5p', '11p', '12p']

        # find region index
        # TODO: replace with list comprehension
        for i, a in enumerate(prostate27roi):
            for j, b in enumerate(a):
                for k, c in enumerate(b):
                    if c == region:
                        rindices = (i, j, k)

        self.calc_nn_ranges(rindices)

        nn = [[[prostate27roi[i][j][k] for i in self.valid_ranges[0]] for j in self.valid_ranges[1]] for k in self.valid_ranges[2]]

        print nn
        return nn

    def calc_nn_ranges(self, rindices):
        """
        calculate index ranges for nearest neighbor region identification
        """
        self.valid_ranges = [range(x-1, x+2) for x in rindices]

        for i in range(3):
            if min(self.valid_ranges[i]) < 0:
                self.valid_ranges[i] = self.valid_ranges[i][1:]
            if i <= 1:
                if max(self.valid_ranges[i]) > 2:
                    self.valid_ranges[i] = self.valid_ranges[i][:-1]
            else:
                if max(self.valid_ranges[i]) > 3:
                    self.valid_ranges[i] = self.valid_ranges[i][:-1]

        # AS regions, span entire lateral extent
        if rindices[1] == 0:
            self.valid_ranges[2] = range(4)

    def check_index_exact_match(self):
        """
        check for an exact match b/w ARFI and histology index lesions
        """
        arfi_index = self.arfi.keys()[0]
        try:
            hist_index = self.histology['pca'][0][0]
            if arfi_index == hist_index:
                self.index_exact_match = True
            else:
                self.index_exact_match = False
        except KeyError:
            self.index_exact_match = False

    def check_index_nn_match(self):
        """
        check for nearest-neightbor match b/w ARFI and histology index lesions
        """
        self.arfi_index = self.arfi.keys()[0]
        try:
            self.hist_index = self.histology['pca'][0][0]
            self.hist_index_nn = self.nearest_neighbor(self.hist_index)

            if self.arfi_index in self.hist_index_nn:
                self.index_nn_match = True
            else:
                self.index_nn_match = False
        except KeyError:
            self.index_nn_match = False

    def check_arfi_benign_match(self, benign):
        """
        check if atrophy or bph is present in the exact or nearest neighbor
        region to an ARFI lesion

        INPUT: benign expected to be either 'atrophy' or 'bph'

        EXAMPLE check_arfi_benign_match(self, 'atrophy')
        """
        try:
            benign_regions = self.histology[benign][1:]
            if any([x in benign_regions for x in self.arfi.keys()]):
                setattr(self, 'arfi_%s_match' % benign, True)
            else:
                setattr(self, 'arfi_%s_match' % benign, False)
        except KeyError:
            setattr(self, 'arfi_%s_match' % benign, False)

    def valid_dataset(self):
        """
        check if this is a valid dataset to include in the sensitivity analysis
        """
        if None in self.arfi and None in self.histology:
            self.valid = False
        else:
            self.valid = True

    def print_analysis(self):
        """
        print analysis results to the terminal
        """
        print "================= PATIENT %s =================" % self.pnum
        if self.valid is False:
            print "Incomplete dataset; not included in analysis."
        else:
            print "Index lesion EXACT match:\t\t%s" % self.index_exact_match
            print "Index lesion NEAREST NEIGHBOR match:\t%s" % self.index_nn_match
            print "ARFI:Atrophy Match:\t\t\t%s" % self.arfi_atrophy_match
            print "ARFI:BPH Match:\t\t\t\t%s" % self.arfi_atrophy_match
