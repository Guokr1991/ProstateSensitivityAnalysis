#!/bin/env python
"""
analyze_ARFI_Hist_MRI.py

Analyze ARFI, histology and MRI data from individual read files in each patient
directory.
"""

__author__ = "Mark Palmeri"
__email__ = "mlp6@duke.edu"
__date__ = "2014-11-19"

import csv

invivo_root = '/luscinia/ProstateStudy/invivo'

Pnum = range(56,107)

regions27 = ['1p', '2p', '3p', '4p', '5p', '6p', '7p', '8p', '9p', '10p',
             '11p', '12p', '1a', '2a', '3a', '4a', '5a', '6a', '7a', '8a',
             '9a', '10a', '11a', '12a', '13as', '14as', '15as']

arfi = [[0]*27 for i in range(0,Pnum+1)]
hist = [[0]*27 for i in range(0,Pnum+1)]
mri = [[0]*27 for i in range(0,Pnum+1)]


for p in Pnum:
    arfi_file = 'Patient%i/ARFI_Index_Lesion_IOS.txt' % p
    hist_file = 'Patient%i/Histology/HistologyLesions.txt' % p
    mr_file = 'Patient%i/MRI_Images/MRI_Index_Region.txt' % p
    f = fopen('%s/Patient%i/%s' % (invivo_root, p, arfi_ios_file), 'r')
    reader = csv.reader(f)

    for r in reader:


f = open('ARFI_27_Regions.csv')
reader = csv.reader(f)


# skip the header
next(reader, None)

for r in reader:
    # leave off P# and Matlab index
    ios = r[2::]
    lesion_indices = [a for a, b in enumerate(ios) if b != '0']

    arfi_path = '%s/Patient%s' % (invivo_root, r[0])

    try:
        f = open('%s/ARFI_Index_Lesion_IOS.txt' % arfi_path, 'a')
        if len(lesion_indices) == 0:
            f.write('None')
        else:
            for les in lesion_indices:
                f.write('%s, %s\n' % (regions27[les], ios[les]))
                print "%s, %s" % (regions27[les], ios[les])
    except IOError:
        print "Cannot create file to save ARFI IOS data"
    else:
        print "Successfully wrote Patient%s ARFI index lesion IOS data." % r[0]
        f.close()
