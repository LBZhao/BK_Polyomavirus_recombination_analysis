'''
This script is used to analyze recombination between BK polyomavirus and host cell genome (whole human genome).
To use this script. Run blast with command
blastn -query 'input FASTA file path' -db nt -out 'output file path' -outfmt "10 qseqid qstart qend sstart send frames sseq qseq stitle" -gilist 'GI mask file' -num_threads 3
Failed to use the -outfmt "10 qseqid qstart qend sstart send frames sseq qseq stitle" will cause an error.
The idea of this script is to format the blast file to reflect recombination.
In order to achieve multiple functions with one script, the following options are created:
-m  output microhomology length in CSV file
    -n  non-gap, a gap will not be considered as a match in this mode
        for example, AATAA and AA-AA is considered as a match without this mode.
    -d  duplication will be removed from final results
    -b  output NHEJ and MMEJ break sites, break sites will be output in pairs into CSV file
        -d  duplication will be removed from final results
        -s  output break sites in one column into CSV file
-a  output genome sequence and reads alignment to screen
    -d  duplication will be removed from final results
-b  output break sites in pairs into CSV file
    -d  duplication will be removed from final results
    -s  output break sites in one column into CSV file

Example Bash command:
    faslist=`ls ./*.fas`
    for eachfile in $faslist
    do
        python ./Format_Blast_result.py  $eachfile
    done

Example BLAST command:
    blastdbcmd -db nt -taxids 9606 -outfmt %g > taxid9606.gi
    blastn -query ./input.fas -db nt -out ./output -gilist taxi9606.gi -outfmt "10 qseqid qstart qend sstart send frames sseq qseq stitle" -num_threads 4
'''


import sys
from Bio import SeqIO
from Bio.Seq import Seq
from operator import itemgetter, attrgetter

counter=0
#Define a function that formats and prints microhomology
def printmicrohomology(output1,output2,output3):
    #print section 1
    print('%8s' % leftsubjectstart + '  ' + output1[:100] + ' '*3 + "DIK")
    print('%8s' % '1' + '  ' + output2[:100] + ' ' * 3 + 'Read')
    print('%8s' % rightsubjectstart + '  ' + output3[:100]+ ' '*3 + "DIK")
    print('\n')

    #print section 2
    if leftdirection=='1/1':
        print('%8s' % (leftsubjectstart+100) + '  '+output1[100:200]+' '*3+"DIK")
    else:
        print('%8s' % (leftsubjectstart-100) + '  '+output1[100:200]+' '*3+"DIK")
    print('%8s' % '101'+'  '+output2[100:200]+' '*3+'Read')
    if rightdirection=='1/1':
        print('%8s' % (rightsubjectstart+ 100) +'  '+output3[100:200]+' '*3+"DIK")
    else:
        print('%8s' % (rightsubjectstart- 100) +'  '+output3[100:200]+' '*3+"DIK")
    print('\n')

    #print section 3
    if leftdirection=='1/1':
        print('%8s' % (leftsubjectstart+200)+'  '+output1[200:]+ '  ' + '%5s' % leftdirection + '  '+'3\'end-'+ '%-8s' % leftsubjectend + ' '*3+"DIK")
    else:
        print('%8s' % (leftsubjectstart-200)+'  '+output1[200:]+ '  ' + '%5s' % leftdirection + '  '+ '3\'end-'+'%-8s' % leftsubjectend + ' '*3+"DIK")
    print('%8s' % '201'+'  '+output2[200:]+ '  ' + '  N/A'+ '  '+ '3\'end-'+'%-8s' % '251' + ' '*3+ 'Read')
    if rightdirection=='1/1':
        print('%8s' % (rightsubjectstart+200)+'  '+output3[200:]+ '  ' + '%5s' % rightdirection + '  '+ '3\'end-'+'%-8s' % rightsubjectend + ' '*3+"DIK")
    else:
        print('%8s' % (rightsubjectstart-200)+'  '+output3[200:]+ '  ' + '%5s' % rightdirection + '  '+ '3\'end-'+'%-8s' % rightsubjectend + ' '*3+"DIK")
    #Print Space between formated results
    print('\n')
    print('\n')
    print('\n')


for record in SeqIO.parse('./BK_Dik.fasta', "fasta"):
    BKloop=record.seq.upper()
    BKloopr=BKloop.reverse_complement()
genomelength=len(BKloop)


#Read Blast result
inp=open(sys.argv[1][:-4],'r')
oneline=inp.readline()
updated=False

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

    #Starting from this line. All information that is necesssory to build the graph will be acquired.
    rightend=0
    n=len(recombination_result)
    #print(n,recombination_result)
    for j in range(n-1):
        leftstart=recombination_result[j][1]
        leftsubject=recombination_result[j][6]
        leftquery=recombination_result[j][7]
        leftend=recombination_result[j][2]
        leftname=recombination_result[j][8]
        leftsubjectstart=recombination_result[j][3]
        leftsubjectend=recombination_result[j][4]
        leftdirection=recombination_result[j][5]
        rightsubject=recombination_result[j+1][6]
        rightquery=recombination_result[j+1][7]
        rightstart=recombination_result[j+1][1]
        rightend=recombination_result[j+1][2]
        rightname=recombination_result[j+1][8]
        rightsubjectstart=recombination_result[j+1][3]
        rightsubjectend=recombination_result[j+1][4]
        rightdirection=recombination_result[j+1][5]

        '''
        This part align readings into 3 lines.
        output1     the left part BLAST result.
        output2     the READ.
        output3     the right part BLAST result.
        '''

        if leftname != rightname:
            if leftname == "BK polyomavirus DNA" or rightname == "BK polyomavirus DNA":
                if leftstart <= rightstart and leftend >= rightend:
                    pass
                elif rightstart <= leftstart and rightend >= leftend:
                    pass
                elif leftname=="Homo sapiens SINE Alu" or rightname=="Homo sapiens SINE Alu":
                    pass
                else:
                    counter+=1
                    print("Leftstart: {:7}  Leftend:   {:7}  {}\nRightstart:{:7}  Rightend:  {:7}  {}\n".format(leftstart,leftend,leftname,rightstart,rightend,rightname))
print(counter)
