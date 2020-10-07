'''
This script is used to stimulate homology length during random recombination.
The hypothesis is that double-stranded DNA damage will generate various ends at BK polyomavirus replication foci.
Damage ends will randomly anneal at the joining site. If no base-pairs form during annealing, DNA ends will be joined directly.
Test number and repeat times are passed to this scrpit as variable. Default test times is 2000. Default repeat times is 1.

Options:
In order to achieve multiple functions with one script, the following options are created:
-e  ends slide simulation.
        During this simulation: DNA damages ends slide into each other step by step, and the longgest annealing length is recorded.

            5'NNNNNNNNNNNNNNNNNNNNN3'

                                           3'NNNNNNNNNNNNNNNNNNNNNNNNNNNNN5'

            ------------------------------------------------------------

            5'NNNNNNNNNNNNNNNNNNNNN3'
                                  |
                                3'NNNNNNNNNNNNNNNNNNNN5'

            ------------------------------------------------------------

            5'NNNNNNNNNNNNNNNNNNNNN3'
                                 |
                               3'NNNNNNNNNNNNNNNNNNNN5'

            ------------------------------------------------------------

            5'NNNNNNNNNNNNNNNNNNNNN3'
                                | |
                              3'NNNNNNNNNNNNNNNNNNNN5'

            ------------------------------------------------------------

            5'NNNNNNNNNNNNNNNNNNNNN3'
                         ||||  |  |
                    3'NNNNNNNNNNNNNNNNNNNN5'

-i  ends invasion simulation.
            During this simulation: DNA damages ends invade into double stranded DNA, and the longgest annealing length is recorded.

                5'NNNNNNNNNNNNNNNNNNNNN3'      5'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN3'
                                                 ||||||||||||||||||||||||||||||||||||||
                                               3'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN5'

                ------------------------------------------------------------

                 5'NNNNNNNNNNNNNNNNNNNN
                                      N3'
                                      |
                   3'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN5'

                ------------------------------------------------------------

                 5'NNNNNNNNNNNNNNNNNNN
                                     NN3'
                                     ||
                   3'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN5'

                ------------------------------------------------------------

                 5'NNNNNNNNNNNNNNNNNN
                                    NNN3'
                                    |||
                   3'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN5'

                ------------------------------------------------------------

                 5'NNNNNNNNNNNNNNNNN
                                   NNNN3'
                                   ||||
                   3'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN5'

-m  middle annealing. Please refer to the suppliment data for details.

Example Bash command:
To simulate 3000 random annealing for 10 times:
python Random_stimulation.py 3000 10
To simulate 5000 random annealing for 1 time:
python Random_stimulation.py 5000
'''

import sys
import random
from Bio import SeqIO
from Bio.Seq import Seq

end_slide_model         = "-e" in sys.argv
end_invasion_model      = "-i" in sys.argv
middle_annealing_model  = "-m" in sys.argv
if not end_slide_model and not end_invasion_model:
    middle_annealing_model = True

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

def Homologylength_end(a,b):
    b=b.complement()
    max=0
    for i in range(15):
        walkable = True
        for j in range(i+1):
            if Pair(a[25-i+j],b[j]):
                pass
            else:
                walkable = False
                break
        if walkable:
            max=i+1
    return(max)

def Homologylength_end_invasion(a,b):
    b=b.complement()
    n=0
    if Pair(a[25],b[12]):
        for i in range(12):
            if Pair(a[25-i],b[12-i]):
                n += 1
            else: break
        return(n)
    elif Pair(a[25],b[13]):
        for i in range(13):
            if Pair(a[15-i],b[13-i]):
                n += 1
            else: break
        return(n)
    else: return(n)






#Build a dictionary for BKPyV DIK
temp=Extension(r"./BK_Dik.fasta",extension-1)
DIK_set=ListBuild(temp,extension)
del temp

lenth=len(DIK_set)

try: test_times = int(sys.argv[1])
except: test_times = 2000
try: test_repeats = int(sys.argv[2])
except: test_repeats = 1
combined_result=[]
Homologylength_result={}
for i in range(test_repeats):
    for j in range(27):
        Homologylength_result[j]=0
    for k in range(test_times):
        left=random.randint(0,lenth-1)
        right=random.randint(0,lenth-1)
        if left == right:
            continue
        if middle_annealing_model:
            Homologylength_result[Homologylength(DIK_set[left],DIK_set[right])]+=1
            continue
        if end_slide_model:
            Homologylength_result[Homologylength_end(DIK_set[left],DIK_set[right])]+=1
            continue
        if end_invasion_model:
            Homologylength_result[Homologylength_end_invasion(DIK_set[left],DIK_set[right])]+=1
            continue
    temp=[]
    for l in range(27):
        #print ('{:<10}{:<10}'.format(str(i),str(float(Homologylength_result[i])/float(sys.argv[1])*100)))
        temp.append(float(Homologylength_result[l])/float(sys.argv[1])*100)
    combined_result.append(temp)
for i in range(27):
    for j in range(test_repeats):
        if j < test_repeats-1:
            print('{:.2f}'.format(float(combined_result[j][i])), end = ',')
        else: print('{:.2f}'.format(float(combined_result[j][i])))
