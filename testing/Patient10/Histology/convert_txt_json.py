j = open('HistologyLesions.json', 'w')
j.write('{\n')

index = True
benign = False
with open('HistologyLesions.txt', 'r') as t:
    tfile = t.readlines()

for td in tfile:
    td = td[:-1]

    if 'pca' in td and index:
        j.write('\t"pca": [\n')
        j.write('\t\t{\n\t\t\t"region": "%s",\n' % td.split(',')[1])
        j.write('\t\t\t"volume_cc": %.1f,\n' % float(td.split(',')[2]))
        j.write('\t\t\t"Gleason": %i,\n' % float(td.split(',')[3]))
        j.write('\t\t\t"index": true\n\t\t}')
        index = False
    elif 'pca' in td and not index:
        j.write(',\n')
        j.write('\t\t{\n\t\t\t"region": "%s",\n' % td.split(',')[1])
        j.write('\t\t\t"volume_cc": %.1f,\n' % float(td.split(',')[2]))
        j.write('\t\t\t"Gleason": %i,\n' % float(td.split(',')[3]))
        j.write('\t\t\t"index": false\n\t\t}')

    if ('atrophy' in td) or ('bph' in td):
        if not benign:
            j.write('],\n')
        else:
            j.write(',\n')
        num_regions = len(td.split(',')[1:])
        j.write('\t"%s": {\n\t\t"regions": [' % td.split(',')[0])
        for n, r in enumerate(td.split(',')[1:]):
            if n < (num_regions-1):
                j.write('"%s", ' % r[1:])
            else:
                j.write('"%s"]\n\t\t}' % r[1:])
        benign = True

j.write('\n}')
