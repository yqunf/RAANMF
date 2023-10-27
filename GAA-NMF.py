import numpy as np
import re, os, sys
from sklearn.decomposition import NMF
from numpy import *
import readfast
from collections import Counter
from scipy.cluster.hierarchy import linkage, cophenet


Amino_Acid_type = 'ACDEFGHIKLMNPQRSTVWY'

def AAC(fastas):
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    encodings = []
    header = ['#']
    for i in AA:
        header.append(i)
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        count = Counter(sequence)
        for key in count:
            count[key] = count[key] / len(sequence)
        code = [name]
        for aa in AA:
            code.append(count[aa])
        encodings.append(code)
    return encodings

def data(file):
    X = AAC(readFasta(file))
    del X[0]
    n = []
    for i in X:
        del i[0]
        n.append((i))
    N = np.array(n)
    return N

V = data("D:\研究生学习\算法\数据集\嗜热蛋白数据\最新2021数据集\最新训练.txt")

X = np.zeros((20, 20))
jieguo = []
for i in range(100):
    ac = {}
    result = []
    _dict = {}
    max_components = m  # m=2-19
    nmf = NMF(n_components=max_components, init='random', max_iter=1000, alpha=0.08)
    W = np.round(nmf.fit_transform(V), 8)
    H = np.round(nmf.components_, 8)
    for j in range(20):
        ac[Amino_Acid_type[j]] = H.argmax(axis=0).tolist()[j]
    for key, value in ac.items():
        if value not in _dict.keys():
            _dict[value] = []
        _dict[value].append(key)
    Q=np.zeros((20,20))
    for key, value in _dict.items():
        if len(value) > 1:
            _str = ",".join([str(x) for x in value])
            print("分在{key}组的氨基酸:{values}".format(key=key, values=_str))

            for x in Amino_Acid_type:
                for y in Amino_Acid_type:
                    if x==y:
                        Q[Amino_Acid_type.index(x),Amino_Acid_type.index(y)]=0
                    if x !=y:
                        if x in _str and y in _str:
                            Q[Amino_Acid_type.index(x), Amino_Acid_type.index(y)] = 1
    X+=Q
    print('-' * 100)
Y=X/100
C=np.array([[]])
for i in range(20):
    B=np.array([Y[i,i+1:20]])
    C=np.column_stack((C,B))
Z=linkage(1-C, 'average')
print(cophenet(Z,1-C))
