'''
This file is used to convert Fastq file into Fasta file, the total number of record would be given after convertion.
'''
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
filepath=sys.argv[1]

with gzip.open(filepath, "rt") as inp, open(filepath[:-5], "w+") as output:
    for record in SeqIO.parse(inp, "fastq"):
        SeqIO.write(record,output,"fasta")
inp.close()
output.close()
