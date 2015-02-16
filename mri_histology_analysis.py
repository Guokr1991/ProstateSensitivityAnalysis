from mr_analysis import MRAnalysis

Ptotal = []
Pexact = []
Pnn = []
Patrophy = []
Pbph = []
Pmiss = []
Pclinsig = []
Pclinsigsens = []
Pfalsepositive = []

for p in range(56, 107):
    P = MRAnalysis(p)
    if P.valid:
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
        Pclinsigsens.append(P.clin_sig_sensitivity)
        Pfalsepositive.append(P.false_positive)

PexactIOS = []
PexactGleason = []
for p in Pexact:
    P = MRAnalysis(p) # reruns through MRAnalysis class
    PexactIOS.append(P.mri['index']['IOS'])
    PexactGleason.append(P.histology['index']['Gleason'])

print "\n\nMRI:HISTOLOGY ANALYSIS"
print "======================"
print "Valid Patients (%i): %s" % (len(Ptotal), Ptotal)
print "\nINDEX LESIONS"
print "============="
print "MRI Sensitivity (Exact) = %i/%i (%.2f)" % (len(Pexact),
                                                   len(Ptotal),
                                                   float(len(Pexact)) /
                                                   float(len(Ptotal)))
print "MRI Sensitivity (NN) = %i/%i (%.2f)" % (len(Pnn),
                                                len(Ptotal),
                                                float(len(Pnn)) /
                                                float(len(Ptotal)))
print "Exact MRI:Histology Matches:"
for i, x in enumerate(Pexact):
    print '\t%s (IOS: %s, Histology Gleason : %s)' % (x, PexactIOS[i], PexactGleason[i])
print "NN MRI:Histology Matches: %s" % Pnn
print "Missed Cases: %s" % Pmiss

print "\nMRI LESIONS"
print "============"
MRIclinsig = len([j for i in Pclinsig for j in i if j[0]])
MRIposterior = len([j for i in Pclinsig for j in i if j[1] == 'posterior'])
MRIanterior = len([j for i in Pclinsig for j in i if j[1] == 'anterior'])
MRItotal = len([j for i in Pclinsig for j in i])
print "%i/%i (%.2f) were clinically significant lesions" % (MRIclinsig,
                                                            MRItotal,
                                                            float(MRIclinsig) /
                                                            float(MRItotal))
print "\t%i/%i (%.2f) read lesions were posterior" % (MRIposterior,
                                                    MRItotal,
                                                    float(MRIposterior) /
                                                    float(MRItotal))

print "\t%i/%i (%.2f) read lesions were anterior" % (MRIanterior,
                                                   MRItotal,
                                                   float(MRIanterior) /
                                                   float(MRItotal))
print "False MRI reads:"
print "\tNon-clinically-significant PCA: %s" % [x for x in Pfalsepositive if x == 'pca']
print "\tAtrophy: %s" % [x for x in Pfalsepositive if x == 'atrophy']
print "\tBPH: %s" % [x for x in Pfalsepositive if x == 'bph']

print "\nCLINICALLY-SIGNIFICANT HISTOLOGY LESIONS"
print "========================================"
histclinsig = len([j for i in Pclinsigsens for j in i if j[0]])
histposterior = len([j for i in Pclinsigsens for j in i if j[1] == 'posterior'])
histanterior = len([j for i in Pclinsigsens for j in i if j[1] == 'anterior'])
histtotal = len([j for i in Pclinsigsens for j in i])
print "%i/%i (%.2f) of clinically-significant lesions were detected" % (histclinsig,
                                                                        histtotal,
                                                                        float(histclinsig) /
                                                                        float(histtotal))
print "\t%i/%i (%.2f) of these lesions were posterior" % (histposterior,
                                                    histtotal,
                                                    float(histposterior) /
                                                    float(histtotal))

print "\t%i/%i (%.2f) of these lesions were anterior" % (histanterior,
                                                   histtotal,
                                                   float(histanterior) /
                                                   float(histtotal))



print "\nINDEX LESION BENIGN CONFOUNDERS"
print "==============================="
print "Atrophy: %s" % Patrophy
print "BPH: %s" % Pbph
