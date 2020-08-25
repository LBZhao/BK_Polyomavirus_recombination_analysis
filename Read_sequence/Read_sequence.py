import gzip
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def Extension(path,n):
    input=SeqIO.read(path,"fasta")
    inputextended=input.seq.upper()+input.seq[:n].upper()
    inputextendedreverse=inputextended.reverse_complement().upper()
    return(str(inputextended),str(inputextendedreverse))

#Circularize BKPyV DIK Genome by adding 250bp from the beginning
DIKextended,DIKextendedreverse=Extension(r"./BK_Dik.fasta",250)
#Circularize pGEM7-DIK plasmid by adding 250bp from the beginning
pGEMDIKextended,pGEMDIKextendedreverse=Extension(r"./pGEM7-DIK.fasta",250)

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


#Build a dictionary for BKPyV DIK
DIKDIC=DicBuild(DIKextended,DIKextendedreverse,251)
#Build a dictionary for pGEM7-DIK
pGEMDIKDIC=DicBuild(pGEMDIKextended,pGEMDIKextendedreverse,251)


#Build a dictionary for BKPyV DIK
a,b=Extension(r"./BK_Dik.fasta",24)
DIKDIC25=DicBuild(a,b,25)
#Build a dictionary for pGEM7-DIK
a,b=Extension(r"./pGEM7-DIK.FASTA",24)
pGEMDIKDIC25=DicBuild(a,b,25)
del a,b

for j in range(1,len(sys.argv)):
    with gzip.open(sys.argv[j], 'rt') as reads, open(sys.argv[j][:-9]+".fas", "w") as collect:
        for record in SeqIO.parse(reads,"fastq"):
            if record.seq in pGEMDIKDIC:
                pass
            elif record.seq in DIKDIC:
                pass
            else:
                for i in range(10):
                    if record.seq[i*25:(i*25+25)] in pGEMDIKDIC25:
                        SeqIO.write(record,collect,"fasta")
                        break
                    elif record.seq[i*25:(i*25+25)] in pGEMDIKDIC25:
                        SeqIO.write(record,collect,"fasta")
                        break
                    elif record.seq[i*25:(i*25+25)] in pGEMDIKDIC25:
                        SeqIO.write(record,collect,"fasta")
                        break
