from lesion_analysis import LesionAnalysis

Ptotal = []
Pexact = []
Pnn = []
Patrophy = []
Pbph = []
Pmiss = []

for p in range(56, 107):
    P = LesionAnalysis(p)
    if P.valid_dataset:
        Ptotal.append(p)
        if P.index_exact_match:
            Pexact.append(p)
        if P.index_nn_match:
            Pnn.append(p)
        else:
            Pmiss.append(p)
        if P.arfi_atrophy_match:
            Patrophy.append(p)
        if P.arfi_bph_match:
            Pbph.append(p)

PexactIOS = []
PexactGleason = []
for p in Pexact:
    P = LesionAnalysis(p)
    PexactIOS.append(P.arfi['index']['IOS'])
    PexactGleason.append(P.hist_index['Gleason'])

print "ARFI:HISTLOGY ANALYSIS"
print "======================"
print "Valid Patients (%i): %s" % (len(Ptotal), Ptotal)
print "INDEX LESIONS"
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
    print '%s (IOS: %s, Gleason: %s)' % (x, PexactIOS[i], PexactGleason[i])
print "NN ARFI:Histology Matches: %s" % Pnn
print "Missed Cases: %s" % Pmiss

print "BENIGN CONFOUNDERS"
print "=================="
print "Atrophy: %s" % Patrophy
print "BPH: %s" % Pbph
