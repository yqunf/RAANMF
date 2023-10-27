import re

def GmAAC(fastas):
    group = {
        'A': 'A',
        'C': 'CDHNQST',
        'G': 'GPV',
        'E': 'E',
        'F': 'FIMY',
        'K': 'K',
        'L': 'LW',
        'R': 'R'
    }       #amino acid grouping

    groupKey = group.keys()

    encodings = []
    header = ['#']
    for key in groupKey:
        header.append(key)
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        count = Counter(sequence)
        myDict = {}
        for key in groupKey:
            for aa in group[key]:
                myDict[key] = myDict.get(key, 0) + count[aa]

        for key in groupKey:
            code.append(myDict[key] / len(sequence))
        encodings.append(code)

    return encodings
