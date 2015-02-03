def main():
    mri_txt_to_json()

def mri_txt_to_json():
    j = open('MRILesions.json', 'w')
    j.write('{\n')
    j.write('\t"lesions": [\n')

    index = True
    benign = False
    lesion_num = 2
    with open('MRI_Index_Lesion.txt', 'r') as t:
        tfile = t.readlines()
    
    n = 0; # define indexing variable

    for td in tfile:
        td = td[:-1]

        if index:
            j.write('\t\t{\n\t\t\t"region": "%s",\n' % td.split(',')[0])
            j.write('\t\t\t"IOS": %1.f,\n'     % float(td.split(',')[1]))
            j.write('\t\t\t"Gleason":  %1.f,\n'  %     float(td.split(',')[2][-1]))
            j.write('\t\t\t"diameter_mm": %1.f,\n'  %  float(td.split(',')[3][1:len(td.split(',')[3])-2]))
            
            if td.split(',')[4] == ' N':
                j.write('\t\t\t"extracap": false,\n')
            else:
                j.write('\t\t\t"extracap": true,\n')
            
            j.write('\t\t\t"lesion_number":  %1.f,\n'  % 1)
            j.write('\t\t\t"index": true\n\t\t}')
            index = False

            if n < (len(tfile) -1): # space with comma for next lesion or close with ] if only lesion
                j.write(',')
            else:
                j.write(']\n')
            n = n+1
            
        else:
            j.write('\t\t{\n\t\t\t"region": "%s",\n' % td.split(',')[0])
            j.write('\t\t\t"IOS": %1.f,\n'     % float(td.split(',')[1]))
            j.write('\t\t\t"Gleason":  %1.f,\n'  %     float(td.split(',')[2][-1]))
            j.write('\t\t\t"diameter_mm": %1.f,\n'  %  float(td.split(',')[3][1:len(td.split(',')[3])-2]))
            
            if td.split(',')[4] == ' N':
                j.write('\t\t\t"extracap": false,\n')
            else:
                j.write('\t\t\t"extracap": true,\n')
            j.write('\t\t\t"lesion_number":  %1.f,\n'  % lesion_num)
            lesion_num = lesion_num + 1
            j.write('\t\t\t"index": false\n\t\t}')
	
            if n < (len(tfile) -1):
                j.write(',')
            else:
                j.write(']\n')
            n = n+1
            
    j.write('\n}')
            
if __name__ == '__main__':
    main()
