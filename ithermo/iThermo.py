import os
os.environ['KERAS_BACKEND']='tensorflow'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '0'

import warnings
warnings.filterwarnings("ignore")

import csv
import pandas as pd
import numpy as np
from tensorflow import keras
from tensorflow.keras.models import load_model

from collections import Counter
import re, os, sys, platform
import pandas as pd
import math
import fileinput



def read_csv(filename):
    with open(filename, newline='') as f_input:
        return [list(map(float, row)) for row in csv.reader(f_input)]
        
def model_predict(input_sequence_file, model):

    with open(input_sequence_file) as f:
        records = f.read()   
    records = records.split('>')[1:]
    myFasta = []
    myseq = []

    parray=[]
    
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '-', ''.join(array[1:]).upper())
        myFasta.append([name, sequence])
        myseq.append(sequence)
        protein_ID=name
        parray.append(protein_ID)
 
    #Feature extraction
    # 1-Amino acid composition(AAC)
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    encodings = []
    header = ['#']
    for i in AA:
        header.append(i)
    encodings.append(header)
    for i in myFasta:
        name, sequence = i[0], re.sub('-', '', i[1])
        count = Counter(sequence)
        for key in count:
            count[key] = count[key]/len(sequence)
        code = [name]
        for aa in AA:
            code.append(count[aa])
        encodings.append(code)
    #saving extracted features

    file = ('extracted_features/AAC.tsv')
    with open(file, 'w') as f:
            for i in range(len(encodings[0]) - 1):
                f.write(encodings[0][i] + '\t')
            f.write(encodings[0][-1] + '\n')
            for i in encodings[1:]:
                f.write(i[0] + '\t')
                for j in range(1, len(i) - 1):
                    f.write(str(float(i[j])) + '\t')
                f.write(str(float(i[len(i)-1])) + '\n')
    tsv_file=('extracted_features/AAC.tsv')
    csv_table=pd.read_table(tsv_file,sep='\t')
    csv_table.to_csv('extracted_features/AAC.csv',index=False)
    # 2-Pseudo amino acid composition(PAAC)
    pPath = os.path.split('PAAC.txt')[0]
    sys.path.append(pPath)
    dataFile = re.sub('codes$', '', pPath )+ r'other_files/PAAC.txt' if platform.system() == 'Windows' else re.sub('codes$', '', pPath) + 'other_files/PAAC.txt'
    def Rvalue(aa1, aa2, AADict, Matrix):
        return sum([(Matrix[i][AADict[aa1]] - Matrix[i][AADict[aa2]]) ** 2 for i in range(len(Matrix))]) / len(Matrix)
    lambdaValue,w = 30,0.05
    with open(dataFile) as f:
        records = f.readlines()
    AA = ''.join(records[0].rstrip().split()[1:])
    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i
    AAProperty = []
    AAPropertyNames = []
    for i in range(1, len(records)):
        array = records[i].rstrip().split() if records[i].rstrip() != '' else None
        AAProperty.append([float(j) for j in array[1:]])
        AAPropertyNames.append(array[0])
        AAProperty1 = []
    for i in AAProperty:
        meanI = sum(i) / 20
        fenmu = math.sqrt(sum([(j-meanI)**2 for j in i])/20)
        AAProperty1.append([(j-meanI)/fenmu for j in i])
        encodings = []
    header = ['#']
    for aa in AA:
        header.append('Xc1.' + aa)
    for n in range(1, lambdaValue + 1):
        header.append('Xc2.lambda' + str(n))
    encodings.append(header)
    for i in myFasta:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        theta = []
        for n in range(1, lambdaValue + 1):
            theta.append(sum([Rvalue(sequence[j], sequence[j + n], AADict, AAProperty1) for j in range(len(sequence) - n)]) / (len(sequence) - n))
        myDict = {}
        for aa in AA:
            myDict[aa] = sequence.count(aa)
        code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
        code = code + [(w * j) / (1 + w * sum(theta)) for j in theta]
        encodings.append(code)
    #saving extracted features

    file = ('extracted_features/PAAC.tsv')
    with open(file, 'w') as f:
            for i in range(len(encodings[0]) - 1):
                f.write(encodings[0][i] + '\t')
            f.write(encodings[0][-1] + '\n')
            for i in encodings[1:]:
                f.write(i[0] + '\t')
                for j in range(1, len(i) - 1):
                    f.write(str(float(i[j])) + '\t')
                f.write(str(float(i[len(i)-1])) + '\n')
    tsv_file=('extracted_features/PAAC.tsv')
    csv_table=pd.read_table(tsv_file,sep='\t')
    csv_table.to_csv('extracted_features/PAAC.csv',index=False)
    # 3-Amphiphilic Pseudo Amino Acid Composition
    dataFile = re.sub('codes$', '', pPath )+ r'other_files/PAAC.txt' if platform.system() == 'Windows' else re.sub('codes$', '', pPath) + 'other_files/PAAC.txt'
    lambdaValue,w = 30,0.05
    with open(dataFile) as f:
        records = f.readlines()
        AA = ''.join(records[0].rstrip().split()[1:])
        AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i
        AAProperty = []
        AAPropertyNames = []
        for i in range(1, len(records) - 1):
            array = records[i].rstrip().split() if records[i].rstrip() != '' else None
            AAProperty.append([float(j) for j in array[1:]])
            AAPropertyNames.append(array[0])
    AAProperty1 = []
    for i in AAProperty:
        meanI = sum(i) / 20
        fenmu = math.sqrt(sum([(j - meanI) ** 2 for j in i]) / 20)
        AAProperty1.append([(j - meanI) / fenmu for j in i])
    encodings = []
    header = ['#']
    for i in AA:
        header.append('Pc1.' + i)
    for j in range(1, lambdaValue + 1):
        for i in AAPropertyNames:
            header.append('Pc2.' + i + '.' + str(j))
    encodings.append(header)
    for i in myFasta:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        theta = []
        for n in range(1, lambdaValue + 1):
            for j in range(len(AAProperty1)):
                theta.append(sum([AAProperty1[j][AADict[sequence[k]]] * AAProperty1[j][AADict[sequence[k + n]]] for k in
                                  range(len(sequence) - n)]) / (len(sequence) - n))
        myDict = {}
        for aa in AA:
            myDict[aa] = sequence.count(aa)

        code = code + [myDict[aa] / (1 + w * sum(theta)) for aa in AA]
        code = code + [w * value / (1 + w * sum(theta)) for value in theta]
        encodings.append(code)
    file = ('extracted_features/APAAC.tsv')
    with open(file, 'w') as f:
        
            for i in range(len(encodings[0]) - 1):
                f.write(encodings[0][i] + '\t')
            f.write(encodings[0][-1] + '\n')
            for i in encodings[1:]:
                f.write(i[0] + '\t')
                for j in range(1, len(i) - 1):
                    f.write(str(float(i[j])) + '\t')
                f.write(str(float(i[len(i)-1])) + '\n')
    tsv_file=('extracted_features/APAAC.tsv')
    csv_table=pd.read_table(tsv_file,sep='\t')
    csv_table.to_csv('extracted_features/APAAC.csv',index=False)
    # 4- Dipeptide composition(DPC)
    encodings = []
    diPeptides = [aa1 + aa2 for aa1 in AA for aa2 in AA]
    header = ['#'] + diPeptides
    encodings.append(header)

    AADict = {}
    for i in range(len(AA)):
            AADict[AA[i]] = i

    for i in myFasta:
            name, sequence = i[0], re.sub('-', '', i[1])
            code = [name]
            tmpCode = [0] * 400
            for j in range(len(sequence) - 2 + 1):
                tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j+1]]] = tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j+1]]] +1
            if sum(tmpCode) != 0:
                tmpCode = [i/sum(tmpCode) for i in tmpCode]
            code = code + tmpCode
            encodings.append(code)
    #saving extracted features

    file = ('extracted_features/DPC.tsv')
    with open(file, 'w') as f:
            for i in range(len(encodings[0]) - 1):
                f.write(encodings[0][i] + '\t')
            f.write(encodings[0][-1] + '\n')
            for i in encodings[1:]:
                f.write(i[0] + '\t')
                for j in range(1, len(i) - 1):
                    f.write(str(float(i[j])) + '\t')
                f.write(str(float(i[len(i)-1])) + '\n')
    tsv_file=('extracted_features/DPC.tsv')
    csv_table=pd.read_table(tsv_file,sep='\t')
    csv_table.to_csv('extracted_features/DPC.csv',index=False)
    # 5 Dipeptide deviation from expected mean (DDE)
    myCodons = {
        'A': 4,
        'C': 2,
        'D': 2,
        'E': 2,
        'F': 2,
        'G': 4,
        'H': 2,
        'I': 3,
        'K': 2,
        'L': 6,
        'M': 1,
        'N': 2,
        'P': 4,
        'Q': 2,
        'R': 6,
        'S': 6,
        'T': 4,
        'V': 4,
        'W': 1,
        'Y': 2
    }
    encodings = []
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    diPeptides = [aa1 + aa2 for aa1 in AA for aa2 in AA]
    header = ['#'] + diPeptides
    encodings.append(header)

    myTM = []
    for pair in diPeptides:
        myTM.append((myCodons[pair[0]] / 61) * (myCodons[pair[1]] / 61))

    AADict = {}
    for i in range(len(AA)):
        AADict[AA[i]] = i

    for i in myFasta:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        tmpCode = [0] * 400
        for j in range(len(sequence) - 2 + 1):
            tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j+1]]] = tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j+1]]] +1
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]

        myTV = []
        for j in range(len(myTM)):
            myTV.append(myTM[j] * (1-myTM[j]) / (len(sequence) - 1))

        for j in range(len(tmpCode)):
            tmpCode[j] = (tmpCode[j] - myTM[j]) / math.sqrt(myTV[j])

        code = code + tmpCode
        encodings.append(code)
    #saving extracted features

    file = ('extracted_features/DDE.tsv')
    with open(file, 'w') as f:
            for i in range(len(encodings[0]) - 1):
                f.write(encodings[0][i] + '\t')
            f.write(encodings[0][-1] + '\n')
            for i in encodings[1:]:
                f.write(i[0] + '\t')
                for j in range(1, len(i) - 1):
                    f.write(str(float(i[j])) + '\t')
                f.write(str(float(i[len(i)-1])) + '\n')
    tsv_file=('extracted_features/DDE.tsv')
    csv_table=pd.read_table(tsv_file,sep='\t')
    csv_table.to_csv('extracted_features/DDE.csv',index=False)
    for line in fileinput.input('extracted_features/DDE.csv', inplace=1):
        print(line.lower(), end='') 
    # 6 Composition of k-spaced Amino Acid Pairs (CKSAAP)
    def CKSAAP(myFasta, gap=5, **kw):
        AA = kw['order'] if kw['order'] != None else 'ACDEFGHIKLMNPQRSTVWY'
        encodings = []
        aaPairs = []
        for aa1 in AA:
            for aa2 in AA:
                aaPairs.append(aa1 + aa2)
        header = ['#']
        for g in range(gap+1):
            for aa in aaPairs:
                header.append(aa + '.gap' + str(g))
        encodings.append(header)
        for i in myFasta:
            name, sequence = i[0], i[1]
            code = [name]
            for g in range(gap+1):
                myDict = {}
                for pair in aaPairs:
                    myDict[pair] = 0
                sum = 0
                for index1 in range(len(sequence)):
                    index2 = index1 + g + 1
                    if index1 < len(sequence) and index2 < len(sequence) and sequence[index1] in AA and sequence[index2] in AA:
                        myDict[sequence[index1] + sequence[index2]] = myDict[sequence[index1] + sequence[index2]] + 1
                        sum = sum + 1
                for pair in aaPairs:
                    code.append(myDict[pair] / sum)
            encodings.append(code)


        return encodings


    if __name__ == '__main__':
        
        kw = {'order': 'ACDEFGHIKLMNPQRSTVWY'}
        encodings = CKSAAP(myFasta, gap = 5, **kw)
    #saving extracted features

    file = ('extracted_features/CKSAAP.tsv')
    with open(file, 'w') as f:
            for i in range(len(encodings[0]) - 1):
                f.write(encodings[0][i] + '\t')
            f.write(encodings[0][-1] + '\n')
            for i in encodings[1:]:
                f.write(i[0] + '\t')
                for j in range(1, len(i) - 1):
                    f.write(str(float(i[j])) + '\t')
                f.write(str(float(i[len(i)-1])) + '\n')
    tsv_file=('extracted_features/CKSAAP.tsv')
    csv_table=pd.read_table(tsv_file,sep='\t')
    csv_table.to_csv('extracted_features/CKSAAP.csv',index=False)
    # 7-CTD composition
    def Count(seq1, seq2):
        sum = 0
        for aa in seq1:
            sum = sum + seq2.count(aa)
        return sum
    group1 = {
        'hydrophobicity_PRAM900101': 'RKEDQN',
        'hydrophobicity_ARGP820101': 'QSTNGDE',
        'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
        'hydrophobicity_PONP930101': 'KPDESNQT',
        'hydrophobicity_CASG920101': 'KDEQPSRNTG',
        'hydrophobicity_ENGD860101': 'RDKENQHYP',
        'hydrophobicity_FASG890101': 'KERSQD',
        'normwaalsvolume': 'GASTPDC',
        'polarity':        'LIFWCMVY',
        'polarizability':  'GASDT',
        'charge':          'KR',
        'secondarystruct': 'EALMQKRH',
        'solventaccess':   'ALFCGIVW'
    }
    group2 = {
        'hydrophobicity_PRAM900101': 'GASTPHY',
        'hydrophobicity_ARGP820101': 'RAHCKMV',
        'hydrophobicity_ZIMJ680101': 'HMCKV',
        'hydrophobicity_PONP930101': 'GRHA',
        'hydrophobicity_CASG920101': 'AHYMLV',
        'hydrophobicity_ENGD860101': 'SGTAW',
        'hydrophobicity_FASG890101': 'NTPG',
        'normwaalsvolume': 'NVEQIL',
        'polarity':        'PATGS',
        'polarizability':  'CPNVEQIL',
        'charge':          'ANCQGHILMFPSTWYV',
        'secondarystruct': 'VIYCWFT',
        'solventaccess':   'RKQEND'
    }
    group3 = {
        'hydrophobicity_PRAM900101': 'CLVIMFW',
        'hydrophobicity_ARGP820101': 'LYPFIW',
        'hydrophobicity_ZIMJ680101': 'LPFYI',
        'hydrophobicity_PONP930101': 'YMFWLCVI',
        'hydrophobicity_CASG920101': 'FIWC',
        'hydrophobicity_ENGD860101': 'CVLIMF',
        'hydrophobicity_FASG890101': 'AYHWVMFLIC',
        'normwaalsvolume': 'MHKFRYW',
        'polarity':        'HQRKNED',
        'polarizability':  'KMHFRYW',
        'charge':          'DE',
        'secondarystruct': 'GNPSD',
        'solventaccess':   'MSPTHY'
    }
    groups = [group1, group2, group3]
    property = (
    'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
    'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
    'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')
    encodings = []
    header = ['#']
    for p in property:
        for g in range(1, len(groups) + 1):
            header.append(p + '.G' + str(g))
    encodings.append(header)
    for i in myFasta:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [name]
        for p in property:
            c1 = Count(group1[p], sequence) / len(sequence)
            c2 = Count(group2[p], sequence) / len(sequence)
            c3 = 1 - c1 - c2
            code = code + [c1, c2, c3]
        encodings.append(code)
    #saving extracted features

    file = ('extracted_features/CTDC.tsv')
    with open(file, 'w') as f:
            for i in range(len(encodings[0]) - 1):
                f.write(encodings[0][i] + '\t')
            f.write(encodings[0][-1] + '\n')
            for i in encodings[1:]:
                f.write(i[0] + '\t')
                for j in range(1, len(i) - 1):
                    f.write(str(float(i[j])) + '\t')
                f.write(str(float(i[len(i)-1])) + '\n')
    tsv_file=('extracted_features/CTDC.tsv')
    csv_table=pd.read_table(tsv_file,sep='\t')
    csv_table.to_csv('extracted_features/CTDC.csv',index=False)
    #feature selection
    aac =pd.read_csv('extracted_features/AAC.csv')[['A', 'Q', 'I', 'K', 'T', 'E','S','Y', 'V','D']]
    paac = pd.read_csv('extracted_features/PAAC.csv')[['Xc1.A','Xc1.E','Xc1.T','Xc1.D', 'Xc1.S', 'Xc1.Q','Xc1.G', 'Xc1.K']]
    apaac = pd.read_csv('extracted_features/APAAC.csv')[['Pc1.A', 'Pc1.Q','Pc1.K', 'Pc1.G', 'Pc1.S', 'Pc1.E', 'Pc1.T', 'Pc1.I']]
    dpc = pd.read_csv('extracted_features/DPC.csv')[['AA','EE', 'KE', 'EK', 'KK', 'LK', 'EI', 'KI','IK', 'LA']]
    dpde = pd.read_csv('extracted_features/DDE.csv')[['ke','aa', 'ee']] 
    cksaap = pd.read_csv('extracted_features/CKSAAP.csv')[['EK.gap2', 'EK.gap3','AA.gap2', 'AA.gap1']] 
    ctdc=pd.read_csv('extracted_features/CTDC.csv').drop(labels=['#'], axis=1)
    #feature fusion
    cf1 =pd.merge(paac.reset_index(),apaac.reset_index(),how='inner').drop(labels=['index'], axis=1)
    cf2 =pd.merge(cf1.reset_index(), ctdc.reset_index(), how='inner').drop(labels=['index'], axis=1)
    cf3 =pd.merge(cf2.reset_index(), cksaap.reset_index(),how='inner').drop(labels=['index'], axis=1)
    cf4 =pd.merge(cf3.reset_index(), dpde.reset_index(), how='inner').drop(labels=['index'], axis=1)
    cf5 =pd.merge(cf4.reset_index(), aac.reset_index(), how='inner').drop(labels=['index'], axis=1)
    cf6 =pd.merge(cf5.reset_index(), dpc.reset_index(), how='inner').drop(labels=['index'], axis=1)
    # print the fused features under same label
    cf6.to_csv('extracted_features/fused_features.csv', index=False)
    cf7 = pd.read_csv('extracted_features/fused_features.csv')
    
    
    #standard scalling of data
    from sklearn.preprocessing import StandardScaler
    import pickle
    with open('other_files/sc.pickle', 'rb') as inputfile:
        sc = pickle.load(inputfile)
    sample = sc.transform(cf7)
    

    #Calculationg probability and prediction
    prob=model.predict(sample)
    preds=(model.predict(sample) > 0.5).astype("int32")
    print("preds", preds)

    return preds,prob,parray


# Model saved with Keras model.save()
MODEL_PATH = 'models/model.h5'

# Loading trained model
model = load_model(MODEL_PATH)
model.summary()


input_sequence_file = r"D:\Mystudy\smallpaper\new\protein stablility\ithermo\input\example.txt"
outfile= r"D:\Mystudy\smallpaper\new\protein stablility\ithermo\output\result.txt"
preds1,prob1,ID = model_predict(input_sequence_file, model)
gg = open(outfile,'w')
gg.write('Protein_ID\tProbability\tType\n')
for i in range(len(ID)):
    if preds1[i] == 0:
        result_txt="Mesophilic"
        #print(ID[i]+'\t'+str(prob1[i][0])+'\t'+str(preds1[i][0])+'\n')
        gg.write(ID[i]+'\t'+str(prob1[i][0])+'\t'+result_txt+'\n')
    else:
        result_txt="Thermophilic"
        gg.write(ID[i]+'\t'+str(prob1[i][0])+'\t'+result_txt+'\n')

gg.close()

print("Please check output folder for predicted results")


