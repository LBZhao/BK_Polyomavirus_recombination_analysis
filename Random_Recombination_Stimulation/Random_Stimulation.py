'''
This script is used to stimulate homologylength during random recombination.
'''

import sys
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def Extension(path,n):
    input=SeqIO.read(path,"fasta")
    inputextended=input.seq.upper()+input.seq[:n].upper()
    return(inputextended)

def ListBuild(sequence,lenth):
    sequenced_list=[]
    for i in range(len(sequence)-(lenth-1)):
        sequenced_list.append(sequence[i:i+lenth])
    return(sequenced_list)

def Pair(a,b):
    if a == "A" and b=="T" or a=="T" and b=="A":
        return(True)
    elif a== "G" and b== "C" or a=="C" and b=="G":
        return(True)
    else: return(False)

extension = 26
def Homologylength(a,b):
    b=b.complement()
    n=0
    if Pair(a[12],b[12]):
        for i in range(14):
            if Pair(a[12+i],b[12+i]):
                n += 1
            else: break
        for i in range(12):
            if Pair(a[11-i],b[11-i]):
                n += 1
            else: break
        return(n)
    elif Pair(a[13],b[13]):
        for i in range(13):
            if Pair(a[13+i],b[13+i]):
                n += 1
            else: break
        for i in range(13):
            if Pair(a[12-i],b[12-i]):
                n += 1
            else: break
        return(n)
    else: return(n)


#Build a dictionary for BKPyV DIK
temp=Extension(r"./BK_Dik.fasta",extension-1)
DIK_set=ListBuild(temp,extension)
del temp

lenth=len(DIK_set)

Homologylength_result={}
for i in range(27):
    Homologylength_result[i]=0
for i in range(int(sys.argv[1])):
    left=random.randint(0,lenth-1)
    right=random.randint(0,lenth-1)
    if left == right:
        continue
    Homologylength_result[Homologylength(DIK_set[left],DIK_set[right])]+=1
for i in range(27):
    print ('{:<10}{:<10}'.format(str(i),str(float(Homologylength_result[i])/float(sys.argv[1])*100)))
