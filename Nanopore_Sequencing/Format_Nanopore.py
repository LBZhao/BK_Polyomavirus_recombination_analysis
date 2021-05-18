'''
To use this script. Run BLAST with command

Examples:
makeblastdb -in './BK_Dik.fas' -title 'DIK' -parse_seqids -out DIK -dbtype nucl
blastn -query 'input FASTA file path' -db nt -out 'output file path' -outfmt "10 qseqid qstart qend sstart send frames sseq qseq stitle" -gilist 'GI mask file' -num_threads 3

Failed to use the -outfmt "10 qseqid qstart qend sstart send frames sseq qseq stitle" will cause an error.
The idea of this script is to format the blast file to reflect recombination.
To achieve multiple functions with one script, the following options are created:
-a  output assembled maps
-b  output break sites in pairs into CSV file
    -r  Circular genome rotation. This option is used to rotate the genome to facilitate circular diagram drawing. To use this function, provide an integer following this option.
        For example, -r 2000 will make the 2001st nucleotide as the first base. Rotate the circular genome by 2000 bp.

'''


import sys
from Bio import SeqIO
from Bio.Seq import Seq
from operator import itemgetter, attrgetter

assemble                = "-a" in sys.argv
break_function          = "-b" in sys.argv
rotation_function       = "-r" in sys.argv

if rotation_function:
    genomelength=5141
    if sys.argv[sys.argv.index("-r")+1].isnumeric():
        rotation_distance=int(sys.argv[sys.argv.index("-r")+1])
    else: raise ValueError("No rotation parameter!")
#Read Blast result
inp=open(sys.argv[1][:-4]+"_blast",'r')
oneline=inp.readline()
updated=False

breakcatch=[]
if assemble:
    breakpoint_file=open('./%s_AssembleBreakpoints.csv' % sys.argv[1][:-4],'w+')
    catched_sequence=open('./%s_Specific_catched_sequence.fas' % sys.argv[1][:-4],'w+')
if break_function:
    breakpoint_file=open('./%s_Paired_Breakpoints.csv' % sys.argv[1][:-4],'w+')
    pass
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

    if assemble:
        recombination_result=sorted(recombination_result, key=itemgetter(1,2))
    #    for item in recombination_result:
    #        print(item[0:6])
        reference_location=[]
        read_location=[]
        distance=[""]
        recombination_index=0
        recombination_index_count=0
        for item in recombination_result:
            read_location.append(item[1])
            read_location.append(item[2])
            reference_location.append(item[3])
            reference_location.append(item[4])

        if item[5]=="1/-1":
            read_location.reverse()
            reference_location.reverse()
            for i in range(2,len(read_location),2):
                distance.append(","+str(read_location[i-1]-read_location[i]))
                recombination_index+=abs(read_location[i-1]-read_location[i])
                recombination_index_count+=1

        else:
            for i in range(2,len(read_location),2):
                distance.append(","+str(read_location[i]-read_location[i-1]))
                recombination_index+=abs(read_location[i]-read_location[i-1])
                recombination_index_count+=1

        if recombination_index/recombination_index_count <=10:
            breakcatch.append(distance)
            breakcatch.append(reference_location)

#           if reference_location==[1,3677,3575,5141]:
#               SeqIO.write(record,catched_sequence,"fasta")
    elif break_function:
        recombination_result=sorted(recombination_result, key=itemgetter(1,2))
    #    for item in recombination_result:
    #        print(item[0:6])
        reference_location=[]
        read_location=[]
        recombination_index=0
        recombination_index_count=0
        for item in recombination_result:
            read_location.append(item[1])
            read_location.append(item[2])
            reference_location.append(item[3])
            reference_location.append(item[4])

        if item[5]=="1/-1":
            read_location.reverse()
            reference_location.reverse()
            for i in range(2,len(read_location),2):
                recombination_index+=abs(read_location[i-1]-read_location[i])
                recombination_index_count+=1

        else:
            for i in range(2,len(read_location),2):
                recombination_index+=abs(read_location[i]-read_location[i-1])
                recombination_index_count+=1

        if recombination_index/recombination_index_count <=10:
            for location in range(1,len(reference_location)-1,2):
                breakcatch.append([reference_location[location],reference_location[location+1]])


if rotation_function:
    for item in breakcatch:
        (a,b)=item
        if a > rotation_distance:
            a -= rotation_distance
        else: a += (genomelength - rotation_distance)
        if b > rotation_distance:
            b -= rotation_distance
        else: b += (genomelength - rotation_distance)
        breakpoint_file.write(str(a)+','+str(b)+'\n')
else:
    breakpoint_file.write(','.join([str(i) for i in pairs]))
    breakpoint_file.write("\n")


breakpoint_file.close()
if assemble:
    catched_sequence.close()
