"""
lesion.py - functions for the prostate lesion analysis b/w ARFI, MRI and
histology
"""

__author__ = 'Mark L. Palmeri'
__email__ = 'mlp6@duke.edu'
__date__ = '2015-01-05'


def main():
    print "Nothing happening here"


def nn27(region):
    """
    nn27.py - extract the set of nearest neighbor regions

    INPUT: region (string)
    OUTPUT: set

    EXAMPLE: nearest_neighbors = nn27('5a')
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


if __name__ == "__main__":
    main()
