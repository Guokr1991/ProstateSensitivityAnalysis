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
                self.arfi[None] = None
        except IOError:
            print "%s does not exist" % self.arfi_ios
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
                    print "WARNING: Malformed histology lesion file."
                # there can be multiple pca lesions
                if 'pca' in lesion:
                    if 'pca' not in self.histology:
                        self.histology['pca'] = [lesion.split(', ')[1:]]
                    else:
                        self.histology['pca'].append(lesion.split(', ')[1:])
                else:
                    self.histology[lesion.split(', ')[0]] = lesion.split(', ')[1:]
        except IOError:
            print "%s does not exist" % self.hist_lesions
            self.histology[None] = None

    @staticmethod
    def nearest_neighbor(region):
        """
        extract the set of nearest neighbor regions
        """
        nn27dict = dict(nn1p=set(['2p', '2a', '1a', '7p', '7a', '3p']),
                        nn2p=set(['1p', '1a', '2a', '4p']),
                        nn3p=set(['4p', '3a', '9p', '5p', '1p']),
                        nn4p=set(['4a', '3p', '3a', '2p', '6p']),
                        nn5p=set(['6p', '11p', '5a', '15as', '3p', '6p']),
                        nn6p=set(['5p', '6a', '5a', '4p']),
                        nn7p=set(['1p', '8p', '7a', '8a', '1a', '9p']),
                        nn8p=set(['7p', '8a', '10p', '7a']),
                        nn9p=set(['10p', '3p', '9a', '7a', '11p']),
                        nn10p=set(['9p', '10a', '9a', '8p', '12p']),
                        nn11p=set(['12p', '5p', '11a', '15as', '12a', '9p']),
                        nn12p=set(['11p', '12a', '11a', '10p']),
                        nn1a=set(['13as', '2a', '7a', '3a']),
                        nn2a=set(['1a', '13as', '2p', '4a']),
                        nn3a=set(['5a', '4p', '3p', '14as', '5a', '1p']),
                        nn4a=set(['3a', '14as', '4p', '2a', '6a']),
                        nn5a=set(['15as', '6a', '3a', '5p']),
                        nn6a=set(['5a', '15as', '6p', '4a']),
                        nn7a=set(['1a', '13as', '8a', '7p', '9a']),
                        nn8a=set(['8p', '7a', '13as', '10a']),
                        nn9a=set(['10a', '14as', '7a', '11a']),
                        nn10a=set(['10p', '8a', '12a', '9a']),
                        nn11a=set(['12a', '15as', '12p', '9a']),
                        nn12a=set(['11a', '15as', '12p', '10a']),
                        nn13as=set(['1a', '7a', '8a', '2a', '14as']),
                        nn14as=set(['3a', '9a', '4a', '10a', '13as', '15a']),
                        nn15as=set(['6a', '12a', '5a', '11a', '14as']))

        nearest_neighbors = nn27dict['nn%s' % region]

        return nearest_neighbors

    def check_index_exact_match(self):
        """
        check for an exact match b/w ARFI and histology index lesions
        """
        arfi_index = self.arfi.keys()[0]
        hist_index = self.histology['pca'][0][0]

        if arfi_index == hist_index:
            self.index_exact_match = True
        else:
            self.index_exact_match = False

    def check_index_nn_match(self):
        """
        check for nearest-neightbor match b/w ARFI and histology index lesions
        """
        self.arfi_index = self.arfi.keys()[0]
        self.hist_index = self.histology['pca'][0][0]
        self.hist_index_nn = self.nearest_neighbor(self.hist_index)

        if self.arfi_index in self.hist_index_nn:
            self.index_nn_match = True
        else:
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

    def print_analysis(self):
        """
        print analysis results to the terminal
        """
        print "================= PATIENT %s =================" % self.pnum
        if None in self.arfi or None in self.histology:
            print "Incomplete dataset; not included in analysis."
        else:
            print "Index lesion EXACT match:\t\t%s" % self.index_exact_match
            print "Index lesion NEAREST NEIGHBOR match:\t%s" % self.index_nn_match
            print "ARFI:Atrophy Match:\t\t\t%s" % self.arfi_atrophy_match
            print "ARFI:BPH Match:\t\t\t\t%s" % self.arfi_atrophy_match
