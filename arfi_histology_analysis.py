from lesion_analysis import LesionAnalysis

Ptotal = []
Pexact = []
Pnn = []
Patrophy = []
Pbph = []
Pmiss = []
Pclinsigtotal = []
Pclinsighit = []

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
        Pclinsigtotal.append(len(P.clin_sig_match))
        Pclinsighit.append(P.clin_sig_match.count(True))

PexactIOS = []
PexactGleason = []
for p in Pexact:
    P = LesionAnalysis(p)
    PexactIOS.append(P.arfi['index']['IOS'])
    PexactGleason.append(P.histology['index']['Gleason'])

print "ARFI:HISTOLOGY ANALYSIS"
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

print "CLINICALLY-SIGNIFICANT LESIONS"
print "=============================="
print "Total: %i/%i (%.2f)" % (sum(Pclinsighit), sum(Pclinsigtotal),
                               float(sum(Pclinsighit)) /
                               float(sum(Pclinsigtotal)))
print "Anterior: "
print "Posterior: "

print "BENIGN CONFOUNDERS"
print "=================="
print "Atrophy: %s" % Patrophy
print "BPH: %s" % Pbph
