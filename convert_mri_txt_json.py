def main():
    mri_txt_to_json()

def mri_txt_to_json():
    j = open('MRILesions.json', 'w')

    index = True
    benign = False
    lesion_num = 2
    with open('MRI_Index_Lesion.txt', 'r') as t:
        tfile = t.readlines()

    for td in tfile:
        td = td[:-1]

        if index:
            #j.write('\t"pca": [\n')
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
            
if __name__ == '__main__':
    main()
