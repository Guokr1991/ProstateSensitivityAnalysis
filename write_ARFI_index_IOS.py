#!/bin/env python
"""
write_ARFI_index_IOS.py

Write individual ARFI index lesion IOS files in each patient directory to do
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

f = open('ARFI_27_Regions.csv')
reader = csv.reader(f)

regions27 = ['1p', '2p', '3p', '4p', '5p', '6p', '7p', '8p', '9p', '10p',
             '11p', '12p', '1a', '2a', '3a', '4a', '5a', '6a', '7a', '8a',
             '9a', '10a', '11a', '12a', '13as', '14as', '15as']

# skip the header
next(reader, None)

for r in reader:
    # leave off P# and Matlab index
    ios = r[2::]
    lesion_indices = [a for a, b in enumerate(ios) if b != '0']

    arfi_path = '%s/Patient%s' % (invivo_root, r[0])

    try:
        f = open('%s/ARFI_Index_Lesion_IOS.txt' % arfi_path, 'w')
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
