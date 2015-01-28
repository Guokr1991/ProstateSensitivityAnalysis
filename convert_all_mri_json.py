import os
from convert_mri_txt_json import mri_txt_to_json

for p in range(56, 107):
    d = '/luscinia/ProstateStudy/invivo/Patient%i/MRI_Images' % p
    if os.path.exists('%s/MRI_Index_Lesion.txt' % d):
        os.chdir(d)
        print "Converting %s" % d
        mri_txt_to_json()
    else:
        print "Not converting %s" % d
