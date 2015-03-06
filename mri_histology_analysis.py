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

# check ECE for index lesion matches
Pece_index_estab = []
Pece_index_focal = []
Pece_index_TN = []
Pece_index_FP=[]
Pece_index_miss=[]
# check ECE for patient-level matches
Pece_pm_TP = []
Pece_pm_FN = []
Pece_pm_TN= []
Pece_pm_FP=[]
# check ECE for patient level matches only on Established focal extent lesions
Pece_em_TP = []
Pece_em_FN = []
Pece_em_TN= []
Pece_em_FP=[]



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
        
        # compile ece index match data
        if P.ece_match['index']:
            if P.ece_match['index']['Established']:
                Pece_index_estab.append(p)
            elif P.ece_match['index']['Focal']:
                Pece_index_focal.append(p)
            elif P.ece_match['index']['True_Negative']:
                Pece_index_TN.append(p)
            elif P.ece_match['index']['False_Positive']:
                Pece_index_FP.append(p)
            else:
                Pece_index_miss.append(p)
        # compile ECE patient-level match data
        if P.ece_match['patient_match']:
            if P.ece_match['patient_match']=='True_Positive':
                Pece_pm_TP.append(p)
            elif P.ece_match['patient_match']=='False_Negative':
                Pece_pm_FN.append(p)
            elif P.ece_match['patient_match']=='True_Negative':
                Pece_pm_TN.append(p)
            elif P.ece_match['patient_match']=='False_Positive':
                Pece_pm_FP.append(p)
        # compile ECE patient-level match data for established ECE extent only
        if P.ece_match['patient_match']:
            if P.ece_match['established_match']=='True_Positive':
                Pece_em_TP.append(p)
            elif P.ece_match['established_match']=='False_Negative':
                Pece_em_FN.append(p)
            elif P.ece_match['established_match']=='True_Negative':
                Pece_em_TN.append(p)
            elif P.ece_match['established_match']=='False_Positive':
                Pece_em_FP.append(p)
            
    

PexactIOS = []
PexactGleason = []
# complie ece_match data
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

positive_ece_index_cases = (len(Pece_index_estab + Pece_index_focal + Pece_index_miss))
positive_ece_patient_cases = (len(Pece_pm_TP + Pece_pm_FN))
positive_ece_estab_cases = (len(Pece_em_TP + Pece_em_FN))
print "\nEXTRACAPSULAR EXTENT ANALYSIS"
print "========================================"
print "%i/%i (%.2f) of detected MRI Index lesions had Extracapsular Extent" %  (positive_ece_index_cases,
                                                                    len(Pnn),
                                                                    float (positive_ece_index_cases) /
                                                                    float (len(Pnn)))
                                                                    
print "\t%i/%i (%.2f) of ECE index lesions were False Positive for Extracapsular extent in MRI" %  (len(Pece_index_FP),
                                                            len(Pnn),
                                                            float (len(Pece_index_FP)) / 
                                                            float (len(Pnn)))
                                            
print "\t%i/%i (%.2f) of MRI lesions called matched with Established Extracapsular Extent" %  (len(Pece_index_estab),
                                                            positive_ece_index_cases,
                                                            float (len(Pece_index_estab)) / 
                                                            float (positive_ece_index_cases))

print "\t%i/%i (%.2f) of MRI lesions called matched with Focal Extracapsular Extent" %  (len(Pece_index_focal),
                                                            positive_ece_index_cases,
                                                            float (len(Pece_index_focal)) / 
                                                            float (positive_ece_index_cases))
                                                            
print "\t%i/%i (%.2f) of ECE index lesions did not call Extracapsular extent" %  (len(Pece_index_miss),
                                                            positive_ece_index_cases,
                                                            float (len(Pece_index_miss)) / 
                                                            float (positive_ece_index_cases))

print "\n%i/%i (%.2f) of MRI Patients had any Extracapsular Extent" %  (positive_ece_patient_cases,
                                                                    len(Ptotal),
                                                                    float (positive_ece_patient_cases) /
                                                                    float (len(Ptotal))) 

print "\t%i/%i (%.2f) of MRI patients were False Positive for Extracapsular extent in MRI" %  (len(Pece_pm_FP),
                                                            len(Ptotal),
                                                            float (len(Pece_pm_FP)) / 
                                                            float (len(Ptotal)))
                                                            
print "\t%i/%i (%.2f) of MRI Patients with ECE were called ECE on MRI" %  (len(Pece_pm_TP),
                                                                    positive_ece_patient_cases,
                                                                    float (len(Pece_pm_TP)) /
                                                                    float (positive_ece_patient_cases))

print "\t%i/%i (%.2f) of MRI Patients with ECE were not called ECE on MRI" %  (len(Pece_pm_FN),
                                                                    positive_ece_patient_cases,
                                                                    float (len(Pece_pm_FN)) /
                                                                    float (positive_ece_patient_cases)) 

print "\n%i/%i (%.2f) of MRI Patients had any Established Extracapsular Extent" %  (positive_ece_estab_cases,
                                                                    len(Ptotal),
                                                                    float (positive_ece_estab_cases) /
                                                                    float (len(Ptotal))) 

print "\t%i/%i (%.2f) of MRI patients were False Positive for Established Extracapsular extent in MRI" %  (len(Pece_em_FP),
                                                            len(Ptotal),
                                                            float (len(Pece_em_FP)) / 
                                                            float (len(Ptotal)))
                                                            
print "\t%i/%i (%.2f) of MRI Patients with ECE were called ECE on MRI" %  (len(Pece_em_TP),
                                                                    positive_ece_estab_cases,
                                                                    float (len(Pece_em_TP)) /
                                                                    float (positive_ece_estab_cases))

print "\t%i/%i (%.2f) of MRI Patients with ECE were not called ECE on MRI" %  (len(Pece_em_FN),
                                                                    positive_ece_estab_cases,
                                                                    float (len(Pece_em_FN)) /
                                                                    float (positive_ece_estab_cases)) 


