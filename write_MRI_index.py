#!/bin/env python
"""
write_MRI_index.py

Write individual MRI index lesion files in each patient directory to do
additional analysis instead of Zach's spreadsheets that challenges plot
generation, etc.

This should be a one-time-use script, with a follow-on script to collect all of
the data for downstream analysis.
"""

__author__ = "Mark Palmeri"
__email__ = "mlp6@duke.edu"
__date__ = "2014-11-18"

import csv

invivo_root = '/luscinia/ProstateStudy/invivo'

f = open('MR_27_Regions_full.csv')
reader = csv.reader(f)

regions27 = ['1p', '2p', '3p', '4p', '5p', '6p', '7p', '8p', '9p', '10p',
             '11p', '12p', '1a', '2a', '3a', '4a', '5a', '6a', '7a', '8a',
             '9a', '10a', '11a', '12a', '13as', '14as', '15as']

# skip the header
next(reader, None)

for r in reader:
    # leave off P# and Matlab index
    ios = r[2::]
    lesion_indices = [a for a, b in enumerate(ios) if b != '']

#    arfi_path = '%s/Patient%s/MRI_Images' % (invivo_root, r[0][:2])
# test in the first patient directory P56
    mri_path = '%s/Patient56/MRI_Images/test' % invivo_root
# change created file to test name
    try:
        f = open('%s/MRI_Index_Region_test7.txt' % mri_path, 'a')
        if len(lesion_indices) == 0:
            f.write('None')
        else:
            for les in lesion_indices[::4]:
                f.write('%s, %s, ,G%s, %smm, %s\n' % (regions27[les/4], ios[les],
                                          ios[les+1],ios[les+2],ios[les+3]))
                print "%s, %s, ,G%s, %smm, %s\n" % (regions27[les/4], ios[les],
                                          ios[les+1],ios[les+2],ios[les+3])
    except IOError:
        print "Cannot create file to save MRI index lesion data"
    else:
        print "Successfully wrote Patient%s MRI index lesion data." % r[0]
        f.close()
