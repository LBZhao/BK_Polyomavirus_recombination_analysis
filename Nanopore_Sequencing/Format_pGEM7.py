'''
To use this script. Run BLAST with command

Examples:
makeblastdb -in './BK_Dik.fas' -title 'DIK' -parse_seqids -out DIK -dbtype nucl
blastn -query 'input FASTA file path' -db nt -out 'output file path' -outfmt "10 qseqid qstart qend sstart send frames sseq qseq stitle" -gilist 'GI mask file' -num_threads 3

Failed to use the -outfmt "10 qseqid qstart qend sstart send frames sseq qseq stitle" will cause an error.
The idea of this script is to format the blast file to reflect recombination.
To achieve multiple functions with one script, the following options are created:
-m  output microhomology length in CSV file
    -n  non-gap, a gap will not be considered as a match in this mode
        for example, AATAA and AA-AA are considered as a match without this mode.
    -d  duplication will be removed from final results
    -b  output NHEJ and MMEJ break sites, break sites will be output in pairs into CSV file
        -d  duplication will be removed from final results
        -s  output break sites in one column into CSV file
-a  output genome sequence and reads alignment to screen
    -d  duplication will be removed from final results
-b  output break sites in pairs into CSV file
    -d  duplication will be removed from final results
    -s  output break sites in one column into CSV file
    -r  Circular genome rotation. This option is used to rotate the genome to facilitate circular diagram drawing. To use this function, provide an integer following this option.
        For example, -r 2000 will make the 2001st nucleotide as the first base. Rotate the circular genome by 2000 bp.
Example Bash command:
    faslist=`ls ./*.fas`
    for eachfile in $faslist
    do
        python ./Format_Blast_result.py  $eachfile
    done
'''


import sys
from Bio import SeqIO
from Bio.Seq import Seq
from operator import itemgetter, attrgetter


#Read Blast result
inp=open(sys.argv[1][:-4]+"_blast",'r')
oneline=inp.readline()
updated=False

breakcatch=[]
breakpoint_file=open('./%s_Breakpoints.csv' % sys.argv[1][:-4],'w+')
catched_sequence=open('./%s_catched_sequence.fas' % sys.argv[1][:-4],'w+')

#Read FASTA file for alignment. !!Not original BLAST File.
for record in SeqIO.parse(sys.argv[1], "fasta"):
    recombination_result=[]
    temp=oneline.split(',')
    if temp != ['']:
        for i in range(1,5):
            temp[i]=int(temp[i])
    '''
    Check if reading ID and BLAST ID are consistant.
    oneline is a varible that read BLSAT result seperated by ','.
    After split oneline by ',',
    oneline[0] is record.ID (M02127:61:000000000-A6YMP:1:1101:17150:1605)
    oneline[1] is start position on Read
    oneline[2] is end position on Read
    oneline[3] is start position on reference Database
    oneline[4] is end position on reference Database
    oneline[5] is aliment direction. (1/-1)
    oneline[6] is subject sequence
    oneline[7] is query sequencepyt
    oneline[8] is Subject Name
    '''
    while temp[0]==record.id:
        #Read all blast result into a 2D list named recombination_result
        recombination_result.append(temp)
        oneline=inp.readline()
        temp=oneline.split(',')
        if temp != ['']:
            for i in range(1,5):
                temp[i]=int(temp[i])
    #In rare cases, BLAST will fail to provide any alignment. If it happends, the following 2 lines will make sure this Read is dropped.
    if len(recombination_result)==0 or len(recombination_result)==1:
        continue

    recombination_result=sorted(recombination_result, key=itemgetter(1,2))
#    for item in recombination_result:
#        print(item[0:6])
    temp=[]
    temp1=[]
    temp2=[""]
    recombination_index=0
    recombination_index_count=0
    for item in recombination_result:
        temp1.append(item[1])
        temp1.append(item[2])
        temp.append(item[3])
        temp.append(item[4])

    if item[5]=="1/-1":
        temp1.reverse()
        temp.reverse()
        for i in range(2,len(temp1),2):
            temp2.append(","+str(temp1[i-1]-temp1[i]))
            recombination_index+=abs(temp1[i-1]-temp1[i])
            recombination_index_count+=1

    else:
        for i in range(2,len(temp1),2):
            temp2.append(","+str(temp1[i]-temp1[i-1]))
            recombination_index+=abs(temp1[i]-temp1[i-1])
            recombination_index_count+=1

    if recombination_index/recombination_index_count <=10:
        breakcatch.append(temp2)
        breakcatch.append(temp)
        SeqIO.write(record,catched_sequence,"fasta")

for item in breakcatch:
    breakpoint_file.write(','.join([str(i) for i in item]))
    breakpoint_file.write("\n")

breakpoint_file.close()
catched_sequence.close()
