import os
from convert_histology_txt_json import hist_txt_to_json
from convert_arfi_txt_json import arfi_txt_to_json

for p in range(56, 107):
    d = '/luscinia/ProstateStudy/invivo/Patient%i' % p
    if os.path.exists('%s/Histology/HistologyLesions.txt' % d):
        os.chdir(d)
        print "Histology: Converting %s" % d
        hist_txt_to_json()
    else:
        print "Histology: Not converting %s" % d

    if os.path.exists('%s/ARFI_Index_Lesion_IOS.txt' % d):
        os.chdir(d)
        print "ARFI: Converting %s" % d
        arfi_txt_to_json()
    else:
        print "ARFI: Not converting %s" % d
