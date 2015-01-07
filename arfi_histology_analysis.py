from lesion_analysis import LesionAnalysis

Ptotal = []
Pexact = []
Pnn = []

for p in range(56, 107):
    P = LesionAnalysis(p)
    if P.valid:
        Ptotal.append(p)
    if P.index_exact_match:
        Pexact.append(p)
    if P.index_exact_match or P.index_nn_match:
        Pnn.append(p)

print "ARFI Sensitivity (Exact) = %i/%i (%.2f)" % (len(Pexact),
                                                   len(Ptotal),
                                                   float(len(Pexact)) /
                                                   float(len(Ptotal)))
print "ARFI Sensitivity (NN) = %i/%i (%.2f)" % (len(Pnn),
                                                len(Ptotal),
                                                float(len(Pnn)) /
                                                float(len(Ptotal)))


print "=========== Valid Patients ==========="
print Ptotal
print "==== Exact ARFI:Histology Matches ===="
print Pexact
print "=== Exact+NN ARFI:Histology Matches ==="
print Pnn
