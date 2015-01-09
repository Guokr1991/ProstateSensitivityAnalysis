class Prostate27:
    """
    define 27 regions, anatomic location (anterior/posterior, PZ/CG)
    """
    def __init__(self):
        self.define_regions()

    def __str__():
        pass

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

    def nearest_neighbors(self, lesion_region):
        """
        extract the set of nearest neighbors for 27 region prostate
        """
        lesion_region_indices = self.location(lesion_region)

        nnranges = self.nn_ranges(lesion_region_indices)

        nn = [[[self.regions[i][j][k] for i in nnranges[0]]
               for j in nnranges[1]] for k in nnranges[2]]

        nnset = set([x for n in nn for m in n for x in m])

        return nnset

    @staticmethod
    def nn_ranges(lesion_location):
        """
        calculate index ranges for nearest neighbor region identification
        """
        nnranges = [range(x-1, x+2) for x in lesion_location]

        for i in range(3):
            if min(nnranges[i]) < 0:
                nnranges[i] = nnranges[i][1:]
            if i <= 1:
                if max(nnranges[i]) > 2:
                    nnranges[i] = nnranges[i][:-1]
            else:
                if max(nnranges[i]) > 3:
                    nnranges[i] = nnranges[i][:-1]

        # AS regions, span entire lateral extent
        if lesion_location[1] == 0:
            nnranges[2] = range(4)

        return nnranges
