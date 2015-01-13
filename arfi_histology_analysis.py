from lesion_analysis import LesionAnalysis

Ptotal = []
Pexact = []
Pnn = []
Patrophy = []
Pbph = []
Pmiss = []
Pclinsig = []

for p in range(56, 107):
    P = LesionAnalysis(p)
    if P.valid_dataset:
        Ptotal.append(p)
        if P.index_match['exact']:
            Pexact.append(p)
        if P.index_match['nn']:
            Pnn.append(p)
        else:
            Pmiss.append(p)
        if P.benign_match['atrophy']:
            Patrophy.append(p)
        if P.benign_match['bph']:
            Pbph.append(p)
        Pclinsig.append(P.clin_sig_match)

PexactIOS = []
PexactGleason = []
for p in Pexact:
    P = LesionAnalysis(p)
    PexactIOS.append(P.arfi['index']['IOS'])
    PexactGleason.append(P.histology['index']['Gleason'])

print "ARFI:HISTOLOGY ANALYSIS"
print "======================"
print "Valid Patients (%i): %s" % (len(Ptotal), Ptotal)
print "\nINDEX LESIONS"
print "============="
print "ARFI Sensitivity (Exact) = %i/%i (%.2f)" % (len(Pexact),
                                                   len(Ptotal),
                                                   float(len(Pexact)) /
                                                   float(len(Ptotal)))
print "ARFI Sensitivity (NN) = %i/%i (%.2f)" % (len(Pnn),
                                                len(Ptotal),
                                                float(len(Pnn)) /
                                                float(len(Ptotal)))
print "Exact ARFI:Histology Matches:"
for i, x in enumerate(Pexact):
    print '\t%s (IOS: %s, Gleason: %s)' % (x, PexactIOS[i], PexactGleason[i])
print "NN ARFI:Histology Matches: %s" % Pnn
print "Missed Cases: %s" % Pmiss

print "\nARFI LESIONS"
print "============"
ARFIclinsig = len([j for i in Pclinsig for j in i if j[0]])
ARFIposterior = len([j for i in Pclinsig for j in i if j[1] == 'posterior'])
ARFIanterior = len([j for i in Pclinsig for j in i if j[1] == 'anterior'])
ARFItotal = len([j for i in Pclinsig for j in i])
print "%i/%i (%.2f) were clinically significant lesions" % (ARFIclinsig,
                                                            ARFItotal,
                                                            float(ARFIclinsig) /
                                                            float(ARFItotal))
print "%i/%i (%.2f) read lesions were posterior" % (ARFIposterior,
                                                    ARFItotal,
                                                    float(ARFIposterior) /
                                                    float(ARFItotal))

print "%i/%i (%.2f) read lesions were anterior" % (ARFIanterior,
                                                   ARFItotal,
                                                   float(ARFIanterior) /
                                                   float(ARFItotal))
print "\nBENIGN CONFOUNDERS"
print "=================="
print "Atrophy: %s" % Patrophy
print "BPH: %s" % Pbph
