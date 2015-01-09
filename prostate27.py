class Prostate27:
    """
    define 27 regions, anatomic location (anterior/posterior, PZ/CG)
    """
    def __init__(self):
        self.define_regions()

    def define_regions(self):
        """
        define 27 regions based on INSERT REFERENCE HERE
        """
        self.regions = [[0 for AP in range(3)] for BA in range(3)]
        self.regions[0][0] = ['13as' for LAT in range(4)]
        self.regions[0][1] = ['2a', '1a', '7a', '8a']
        self.regions[0][2] = ['2p', '1p', '7p', '8p']
        self.regions[1][0] = ['14as' for LAT in range(4)]
        self.regions[1][1] = ['4a', '3a', '9a', '10a']
        self.regions[1][2] = ['4p', '3p', '9p', '10p']
        self.regions[2][0] = ['15as' for LAT in range(4)]
        self.regions[2][1] = ['6a', '5a', '11a', '12a']
        self.regions[2][2] = ['6p', '5p', '11p', '12p']

    def location(self, lesion_region):
        """
        determine region location indices
        """
        for i, a in enumerate(self.regions):
            for j, b in enumerate(a):
                for k, c in enumerate(b):
                    if c == lesion_region:
                        lesion_region_indices = (i, j, k)

        return lesion_region_indices

    '''
    def nearest_neighbor(self, region):
        """
        extract the set of nearest neighbor for 27 region prostate

        INPUT: region (string) - find nearest neighbors around this region
        """


        self.calc_nn_ranges(rindices)

        nn = [[[prostate27roi[i][j][k] for i in self.valid_ranges[0]]
               for j in self.valid_ranges[1]] for k in self.valid_ranges[2]]

        nnset = set([x for n in nn for m in n for x in m])

        return nnset

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

    def arfi_index(self):
        """
        define ARFI index lesion dict based on highest IOS
        """
        try:
            index = {}
            index['IOS'] = max(self.arfi.values())
            index['region'] = [x for x, y in self.arfi.items() if
                               y == index['IOS']][0]
            self.arfi['index'] = index
        except ValueError:
            self.arfi['index'] = None

    def histology_index(self):
        """
        define histology index lesion dict and nearest neighbor set
        """
        try:
            # find max Gleason score, then max volume with that max Gleason
            maxGleason = 0
            maxVolumeCC = 0
            for lesion in self.histology['pca']:
                # Gleason scores that weren't reported were recorded as '-1'
                # these will be set to 6 for now
                if lesion[2] == '-1':
                    lesion[2] = '6'
                if lesion[2] > maxGleason:
                    maxGleason = lesion[2]
                    maxVolumeCC = lesion[1]
                if lesion[2] == maxGleason:
                    if lesion[1] > maxVolumeCC:
                        maxGleason = lesion[2]
                        maxVolumeCC = lesion[1]

            # make sure the lesion is clinically significant
            if maxGleason == 7 or (maxGleason >= 6 and maxVolumeCC >= 500):
                index = {}
                index['region'] = self.histology['pca'][0][0]
                index['Gleason'] = self.histology['pca'][0][2]
                index['nn'] = \
                    self.nearest_neighbor(index['region'])
                self.histology['index'] = index
            else:
                print "No clinically-significant PCA lesion"
                self.histology['index'] = None
        except KeyError:
            print "No PCA lesion"
            self.histology['index'] = None

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
        if None in self.arfi and None in self.histology:
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
            s.append('Index lesion EXACT match:\t\t%s' % self.index_match['exact'])
            s.append('Index lesion NEAREST NEIGHBOR match:\t%s' %
                     self.index_match['nn'])
            s.append('ARFI:Atrophy Match:\t\t\t%s' % self.benign_match['atrophy'])
            s.append('ARFI:BPH Match:\t\t\t\t%s' % self.benign_match['bph'])

        return '\n'.join(s)
    '''
