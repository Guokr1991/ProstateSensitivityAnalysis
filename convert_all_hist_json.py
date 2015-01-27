import os
from convert_histology_txt_json import hist_txt_to_json

for p in range(56, 107):
    d = '/luscinia/ProstateStudy/invivo/Patient%i/Histology' % p
    if os.path.exists('%s/HistologyLesions.txt' % d):
        os.chdir(d)
        print "Converting %s" % d
        hist_txt_to_json()
    else:
        print "Not converting %s" % d
