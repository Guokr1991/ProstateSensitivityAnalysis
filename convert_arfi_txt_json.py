def main():
    arfi_txt_to_json()


def arfi_txt_to_json():
    j = open('ARFI_Lesions.json', 'w')
    j.write('{\n')

    index = True
    with open('ARFI_Index_Lesion_IOS.txt', 'r') as t:
        tfile = t.readlines()

    num_ios = len(tfile)
    if 'None' not in tfile[0]:
        j.write('\t"lesions": [\n')
        for n, td in enumerate(tfile):
            td = td[:-1]
            j.write('\t\t{\n\t\t\t"region": "%s",\n' % td.split(',')[0])
            j.write('\t\t\t"IOS": %i,\n' % int(td.split(',')[1]))
            if index:
                j.write('\t\t\t"index": true\n\t\t}')
                index = False
            else:
                j.write('\t\t\t"index": false\n\t\t}')
            if n < (num_ios - 1):
                j.write(',')
            else:
                j.write(']\n')
    else:
        j.write('\t"read": "no lesions read"')

    j.write('\n}')

if __name__ == '__main__':
    main()
