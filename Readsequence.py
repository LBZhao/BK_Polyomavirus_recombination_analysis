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
    for i in range(25):
        if seq1[i] != seq2[i]:
            n += 1
            if n >= Giveupthreshold:
                n= 9999
                break
    return n
'''
#This function will not be able to solve repetitive sequence problem. Need to modify later.
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
    for i in lookback:
        if type(lookback[i])==type(False):
            if sequencedic[seq[index*25:(index*25+25)]][0] > 0:
                start=lookback[index]-(index-i)*25
                end=(lookback[index]-(index-i)*25)+25
                """Important: this number(8138, the length of pGEM7-DIK)is used to solve the cirular problem."""
                if start < 0:
                    start+=8138
                    end+=8138
                if HammingDistance(seq[i*25:i*25+25],pGEMDIKextended[start:end],3)<3:
                    location.append(lookback[index]-(index-i)*25)
                else: return False
            else:
                start=-(lookback[index]+(index-i)*25)+226
                end=-((lookback[index]+(index-i)*25)-26)+226
                if start < 0:
                    start+=8138
                    end+=8138
                print(seq[i*25:(i*25+25)],pGEMDIKextendedreverse[start:end],HammingDistance(seq[i*25:(i*25+25)],pGEMDIKextendedreverse[start:end],3))
                #a=HammingDistance(seq[item*25:(item*25+25)],pGEMDIKextendedreverse[(lookback[index]+(index-item)*25):((lookback[index]+(index-item)*25)+25)],3)
                #if a <=3:
                    #location.append((lookback[index]+(index-item)*25))
                #else: continue
'''

filelist=(
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118033/118033_GGACTCCT-TATCCTCT_S1_L001_R1_001.fastq.gz',r'/mnt/d/01_pGEM7-DIK.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118034/118034_TCGACTAG-TATCCTCT_S2_L001_R1_001.fastq.gz',r'/mnt/d/02_DIK_Mock.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118035/118035_CCTAGAGT-TATCCTCT_S3_L001_R1_001.fastq.gz',r'/mnt/d/03_DIK_01DPI.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118036/118036_GCGTAAGA-TATCCTCT_S4_L001_R1_001.fastq.gz',r'/mnt/d/04_DIK_03DPI.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118037/118037_TTATGCGA-TATCCTCT_S5_L001_R1_001.fastq.gz',r'/mnt/d/05_DIK_05DPI.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118038/118038_TCGCCTTA-TATCCTCT_S6_L001_R1_001.fastq.gz',r'/mnt/d/06_DIK_07DPI.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118039/118039_CGTCTAAT-GTAAGGAG_S7_L001_R1_001.fastq.gz',r'/mnt/d/07_DIK_10DPI.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118040/118040_TCGACTAG-GTAAGGAG_S8_L001_R1_001.fastq.gz',r'/mnt/d/08_DIK_15DPI.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118041/118041_CCTAGAGT-GTAAGGAG_S9_L001_R1_001.fastq.gz',r'/mnt/d/09_DIK_20DPI.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118042/118042_GCGTAAGA-GTAAGGAG_S10_L001_R1_001.fastq.gz',r'/mnt/d/10_DIK_30DPI.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118043/118043_TTATGCGA-GTAAGGAG_S11_L001_R1_001.fastq.gz',r'/mnt/d/11_DIK_40DPI.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118044/118044_TCGCCTTA-GTAAGGAG_S12_L001_R1_001.fastq.gz',r'/mnt/d/12_DIK_60DPI.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118045/118045_CGTCTAAT-ACTGCATA_S13_L001_R1_001.fastq.gz',r'/mnt/d/13_DIK_DIK_DUN_60DPI.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118046/118046_TCGACTAG-ACTGCATA_S14_L001_R1_001.fastq.gz',r'/mnt/d/14_D3_Mock.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118047/118047_CCTAGAGT-ACTGCATA_S15_L001_R1_001.fastq.gz',r'/mnt/d/15_D3_DIK.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118048/118048_GCGTAAGA-ACTGCATA_S16_L001_R1_001.fastq.gz',r'/mnt/d/16_D3_DIK_Mutant.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118049/118049_TTATGCGA-ACTGCATA_S17_L001_R1_001.fastq.gz',r'/mnt/d/17_D3_Dun.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118050/118050_TCGCCTTA-ACTGCATA_S18_L001_R1_001.fastq.gz',r'/mnt/d/18_D5_Mock.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118051/118051_CGTCTAAT-GCGTAAGA_S19_L001_R1_001.fastq.gz',r'/mnt/d/19_D5_DIK_Mutant.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118052/118052_TCGACTAG-GCGTAAGA_S20_L001_R1_001.fastq.gz',r'/mnt/d/20_D5_DIK.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118053/118053_CCTAGAGT-GCGTAAGA_S21_L001_R1_001.fastq.gz',r'/mnt/d/21_D5_Dunlop.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118054/118054_GCGTAAGA-GCGTAAGA_S22_L001_R1_001.fastq.gz',r'/mnt/d/22_D7_Mock.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118055/118055_TTATGCGA-GCGTAAGA_S23_L001_R1_001.fastq.gz',r'/mnt/d/23_D7_DIK.fastq'),
(r'/mnt/d/Box Sync/Imperiale Lab Notebooks/Linbo Zhao/Projects/BKPyV recombination/20181011 Sequencing Results/Run_2500/imperiale/Sample_118056/118056_TCGCCTTA-GCGTAAGA_S24_L001_R1_001.fastq.gz',r'/mnt/d/24_D7_DIK_Mutant.fastq'))


for (path,outputpath) in filelist:
    totalcounter=0
    dumpcounter=0
    collectcounter=0
    with gzip.open(path, 'rt') as reads, open(outputpath, "w+") as collect:
            for record in SeqIO.parse(reads,"fastq"):
                totalcounter+=1
                if record.seq in DIKDIC:
                    dumpcounter+=1
                elif record.seq in DUNDIC:
                    dumpcounter+=1
                else:
                    for j in range(10):
                        if record.seq[j*25:(j*25+25)] in DIKDIC25:
                            SeqIO.write(record,collect,"fastq")
                            collectcounter+=1
                            break
                        elif record.seq[j*25:(j*25+25)] in DUNDIC25:
                            SeqIO.write(record,collect,"fastq")
                            collectcounter+=1
                            break
    print("%.3f%% dumped in the first step. %.3f%% collect in the final step. %.3f%% items do not match at all." %((dumpcounter*100/totalcounter),(collectcounter*100/totalcounter),((totalcounter-dumpcounter-collectcounter)*100/totalcounter)))
    reads.close()
    collect.close()
