import re

def GmDPC(fastas):
    group = {
        'A': 'A',
        'C': 'CDHNQST',
        'G': 'GPV',
        'E': 'E',
        'F': 'FIMY',
        'K': 'K',
        'L': 'LW',
        'R': 'R'
    } # amino acid groupings

    groupKey = group.keys()
    baseNum = len(groupKey)
    dipeptide = [g1 + '.' + g2 for g1 in groupKey for g2 in groupKey]

    index = {}
    for key in groupKey:
        for aa in group[key]:
            index[aa] = key

    encodings = []
    header = ['#'] + dipeptide
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])

        code = [name]
        myDict = {}
        for t in dipeptide:
            myDict[t] = 0

        sum = 0
        for j in range(len(sequence) - 1):
            myDict[index[sequence[j]] + '.' + index[sequence[j + 1]]] += 1

            sum = sum + 1

        if sum == 0:
            for t in dipeptide:
                code.append(0)
        else:
            for t in dipeptide:
                code.append(myDict[t] / sum)
        encodings.append(code)

    return encodings
