import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def Extension(path,n):
    input=SeqIO.read(path,"fasta")
    inputextended=input.seq.upper()+input.seq[:n].upper()
    inputextendedreverse=inputextended.reverse_complement().upper()
    return(str(inputextended),str(inputextendedreverse))


#Circularize BKPyV Dunlop Genome by adding 250bp from the beginning
DUNextended,DUNextendedreverse=Extension(r"/mnt/d/SynologyDrive/Python_Script/Database/Dunlop.fasta",250)
#Circularize BKPyV DIK Genome by adding 250bp from the beginning
DIKextended,DIKextendedreverse=Extension(r"/mnt/d/SynologyDrive/Python_Script/Database/DIK.fasta",250)
#Circularize pGEM7-DIK plasmid by adding 250bp from the beginning
pGEMDIKextended,pGEMDIKextendedreverse=Extension(r"/mnt/d/SynologyDrive/Python_Script/Database/pGEM-DIK.fas",250)
#print((DUNextended==DUNextendedreverse),(DIKextended==DIKextendedreverse),(pGEMDIKextended==pGEMDIKextendedreverse))


def DicBuild(sequence,reversecomplementsequence,lenth):
    sequencedic={}
    for i in range(len(sequence)-(lenth-1)):
        if sequence[i:i+lenth] not in sequencedic:
            sequencedic[sequence[i:i+lenth]]=[i]
        else:sequencedic[sequence[i:i+lenth]].append(i)
    for i in range(len(reversecomplementsequence)-(lenth-1)):
        if reversecomplementsequence[i:i+lenth] not in sequencedic:
            sequencedic[reversecomplementsequence[i:i+lenth]]=[-i-1]
        else:sequencedic[reversecomplementsequence[i:i+lenth]].append(-i-1)
    return(sequencedic)


#Build a dictionary for BKPyV Dunlop
DUNDIC=DicBuild(DUNextended,DUNextendedreverse,251)
#Build a dictionary for BKPyV DIK
DIKDIC=DicBuild(DIKextended,DIKextendedreverse,251)
#Build a dictionary for pGEM7-DIK
pGEMDIKDIC=DicBuild(pGEMDIKextended,pGEMDIKextendedreverse,251)


#Build a dictionary for BKPyV Dunlop
a,b=Extension(r"/mnt/d/SynologyDrive/Python_Script/Database/Dunlop.fasta",24)
DUNDIC25=DicBuild(a,b,25)
#Build a dictionary for BKPyV DIK
a,b=Extension(r"/mnt/d/SynologyDrive/Python_Script/Database/DIK.fasta",24)
DIKDIC25=DicBuild(a,b,25)
#Build a dictionary for pGEM7-DIK
a,b=Extension(r"/mnt/d/SynologyDrive/Python_Script/Database/pGEM-DIK.fas",24)
pGEMDIKDIC25=DicBuild(a,b,25)
del a,b

def HammingDistance(seq1,seq2,Giveupthreshold):
    '''
    This function takes two DNA sequence that have the same lenth as string.
    If seq1 is longer than seq 2. It will generate a error.
    Giveupthreshold is used to increse the speed of comparason.
    Because this function takes tremendous time to compare each nucleitide, once the Hamming Distance reaches the threshhold, further comparison is usless.
    Therefore, once the HammingDistance reach Giveupthreshold, this function will just return a huge number indicating that the humming distance is useless.
    '''
    n=0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            n += 1
            if n >= Giveupthreshold:
                n= 9999
                break
    return n

def ismatch(seq,sequencedic):
    global pGEMDIKextended,pGEMDIKextendedreverse
    location=[]
    lookback={}
    not_the_same_at_all_check=0
    for i in range(10):
        if seq[i*25:(i*25+25)] in sequencedic:
            lookback[i]=sequencedic[seq[i*25:(i*25+25)]][0]
            location.append(sequencedic[seq[i*25:(i*25+25)]][0])
        else:
            lookback[i]=False
            not_the_same_at_all_check+=1
    if not_the_same_at_all_check==10:
        return False
    for i in lookback:
        if type(lookback[i])==type(0):
            index=i
            break
    print(location)
    for i in lookback:
        if type(lookback[i])==type(False):
            print(seq[i*25:i*25+25],pGEMDIKextended[(lookback[index]-(index-i)*25):((lookback[index]-(index-i)*25)+25)])
                #if sequencedic[seq[index*25:(index*25+25)][0] > 0:
                #a=HammingDistance(seq[item*25:item*25+25],pGEMDIKextended[(lookback[index]-(index-item)*25):((lookback[index]-(index-item)*25)+25)],3)
                #if  a<=3:
                    #location.append((lookback[index]-(index-item)*25))
                #else: continue
            #else:
                #a=HammingDistance(seq[item*25:(item*25+25)],pGEMDIKextendedreverse[(lookback[index]+(index-item)*25):((lookback[index]+(index-item)*25)+25)],3)
                #if a <=3:
                    #location.append((lookback[index]+(index-item)*25))
                #else: continue

with gzip.open(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118033/118033_GGACTCCT-TATCCTCT_S1_L001_R1_001.fastq.gz', 'rt') as reads, open(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118033/collect1', "w") as collect:
    for record in SeqIO.parse(reads,"fastq"):
        if record.seq in pGEMDIKDIC:
            pass
        else:
            ismatch(record.seq,pGEMDIKDIC25)
#            SeqIO.write(record,collect,"fastq")
