'''
1. cat all sequence files with bash command:
for folder in $folderlist; do cd $folder; cat * > ../"${folder:2}.fastq"; cd ..; done
'''


import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


for j in range(1,len(sys.argv)):
    with open(sys.argv[j], 'rt') as reads, open(sys.argv[j][:-6]+".fas", "w") as collect:
        for record in SeqIO.parse(reads,"fastq"):
            if len(record.seq) >1000 and len(record.seq) <8500:
                SeqIO.write(record,collect,"fasta")
