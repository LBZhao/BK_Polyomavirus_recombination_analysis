'''
To use this script. Run blast with command
blastn -query 'input FASTA file path' -db nt -out 'output file path' -outfmt "10 qseqid qstart qend sstart send frames sseq qseq stitle" -gilist 'GI mask file' -num_threads 3
Failed to use the -outfmt -outfmt "10 qseqid qstart qend sstart send frames sseq qseq stitle" will cause error.
The idea of this script is format the blast file to reflect recombination.
'''
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from operator import itemgetter, attrgetter

#recombi_readings=open('./specific_recombination.fas','w+')

m = 0
for record in SeqIO.parse('./BK_DIK_NCCR2000.fas', "fasta"):
    BKloop=record.seq.upper()
    BKloopr=BKloop.reverse_complement()
genomelength=len(BKloop)
#Read BLAST result.
leftbreakpoint=[]
rightbreakpoint=[]
totalbreakpoint=[]
#Read Blast result
inp=open(sys.argv[1][:-4],'r')
oneline=inp.readline()
updated=False

#Define a function to isolate the joint.
def jointfinder(left,read,right):
    result=[[],[],[]]
    for i in range(251):
        if len(left)<251 or len(read)<251 or len(right)<251:
            break
        if left[i] != " ":
            if read[i] != " ":
                if right[i] != " ":
                    result[0].append(left[i])
                    result[1].append(read[i])
                    result[2].append(right[i])
    return result

#Define a function to identify microhomology
def microhomologyfinder(inputdic):
    def match(a,b):
        if a==b:
            return True
        elif a=="-":
            return True
        elif b=="-":
            return True
        else: return False
    homologylength=0
    walk=[[],[]]
    for i in range(len(inputdic[0])):
        if match(inputdic[0][i],inputdic[1][i]):
            walk[0].append(True)
        else: walk[0].append(False)
        if match(inputdic[1][i],inputdic[2][i]):
            walk[1].append(True)
        else: walk[1].append(False)
    microhomology=False
    for o in range(len(walk[0])):
        walkable=True
        if walk[0][o] and walk[1][o]:
            for p in range(o):
                if walk[0][p]==False:
                    walkable=False
                    break
            for q in range(o,len(walk[0])):
                if walk[1][q]==False:
                    walkable=False
                    break
            if walkable:
                homologylength+=1
                microhomology=True
    #return(microhomology,homologylength)
    print(homologylength)

#Define a function that format and print microhomology
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
    print('\n')
    print('\n')


#Read FASTA file for alignment. !!Not original BLAST File.
for record in SeqIO.parse(sys.argv[1], "fasta"):
    recombination_result=[]
    temp=oneline.split(',')
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
    if len(recombination_result)==0:
        continue

    recombination_result=sorted(recombination_result, key=itemgetter(1,2))

    print (recombination_result)

    #Solve Circular issue by connect ends.
    for i in range(len(recombination_result)-1):
        #Solve Circular issue by connect ends.
        if recombination_result[i][4] == 1 and recombination_result[i+1][3] == len(BKloop):
            print ("solution 1 ", recombination_result)
            recombination_result.insert(i,[recombination_result[i][0],recombination_result[i][1],recombination_result[i+1][2],recombination_result[i][3],recombination_result[i+1][4],recombination_result[i][5],recombination_result[i][6]+recombination_result[i+1][6],recombination_result[i][7]+recombination_result[i+1][7],recombination_result[i][8]])
            recombination_result.pop(i+1)
            recombination_result.pop(i+1)
            print ("solution 1 solved", recombination_result)
            #print("alart", [i[1:5] for i in recombination_result])
            break
        if recombination_result[i][4] == len(BKloop) and recombination_result[i+1][3] == 1:
            print ("solution 2 ",  recombination_result)
            recombination_result.insert(i,[recombination_result[i][0],recombination_result[i][1],recombination_result[i+1][2],recombination_result[i][3],recombination_result[i+1][4],recombination_result[i][5],recombination_result[i][6]+recombination_result[i+1][6],recombination_result[i][7]+recombination_result[i+1][7],recombination_result[i][8]])
            recombination_result.pop(i+1)
            recombination_result.pop(i+1)
            print ("solution 2 solved ",  recombination_result)
            #print("alart1", [i[1:5] for i in recombination_result])
            break









    #Starting from this line. All information that is necesssory to build the graph will be acquired.
    rightend=0
    n=len(recombination_result)
    if 0<int(recombination_result[0][1])<=50:
        leftstart=int(recombination_result[0][1])
        for j in range(n):
            if leftstart == int(recombination_result[j][1]):
                leftsubject=recombination_result[j][6]
                leftquery=recombination_result[j][7]
                leftend=int(recombination_result[j][2])
                leftname=recombination_result[j][8]
                leftsubjectstart=int(recombination_result[j][3])
                leftsubjectend=int(recombination_result[j][4])
                leftdirection=recombination_result[j][5]
                for i in range(1,n):
                    updated=False
                    '''Note: The following if statement will eliminate all BLAST result that overlap more than 10 bp.
                    (For example, the left part align to human genome and end at 150th bp of the reading. The right part of the Read align to BK genome,
                    and it start at 135th bp. This BLAST alignment will be discarded.)
                    The following if statement will also eliminate all BLAST result that is over 240 bp longer.
                    Once the first alignment is over 240 bp in lenth, the rest of the Read will not be able to provide a reliable BLAST result.
                    Because of the reason stated above, I wrote the next if statement.'''
                    if 200<int(recombination_result[i][2])<252 and int(recombination_result[i][1])>(leftend-50):
                        updated=True
                        rightsubject=recombination_result[i][6]
                        rightquery=recombination_result[i][7]
                        rightstart=int(recombination_result[i][1])
                        rightend=int(recombination_result[i][2])
                        rightname=recombination_result[i][8]
                        rightsubjectstart=int(recombination_result[i][3])
                        rightsubjectend=int(recombination_result[i][4])
                        rightdirection=recombination_result[i][5]
                        break
    else:
        continue
    if rightend < leftend-26:
        continue
    if leftsubjectend==1 and rightsubjectstart==genomelength:
        continue
    if leftsubjectend==genomelength and rightsubjectstart==1:
        continue

    if updated:
        '''
        leftbreakpoint.append(leftsubjectend)
        totalbreakpoint.append(leftsubjectend)
        rightbreakpoint.append(rightsubjectstart)

        if leftsubjectend in (2191,2235,2240,2247) and rightsubjectstart in (2191,2235,2240,2247):
            #SeqIO.write(record, recombi_readings, "fasta")
            pass
        totalbreakpoint.append((leftsubjectend,rightsubjectstart))

        m += 1


        This part align readings into 3 lines.
        output1     the left part BLAST result.
        output2     the READ.
        output3     the right part BLAST result.
        This part is complicated because blast will put insert or deletion into the output sequence.
        '''
        genomelength=len(BKloop)

        if leftdirection=='1/1':
            output1=' '*(leftstart-1) + leftsubject + BKloop[leftsubjectend:leftsubjectend+10]+' '*(241-leftend)
        else:
            output1=' '*(leftstart-1) + leftsubject + BKloopr[genomelength-leftsubjectend+1:genomelength-leftsubjectend+11]+' '*(241-leftend)
        output2=record.seq
        if rightdirection=='1/1':
            output3=' '*(rightstart-11) +BKloop[rightsubjectstart-11:rightsubjectstart-1]+ rightsubject+' '*(251-rightend)
        else:
            output3=' '*(rightstart-11) +BKloopr[genomelength-rightsubjectstart-10:genomelength-rightsubjectstart]+ rightsubject + ' '*(251-rightend)
        microhomologyfinder(jointfinder(output1,output2,output3))

        print("leftstart:"+str(leftstart)+"  rightstart:"+str(rightstart)+"    leftsubjectend:"+str(leftsubjectend)+"    rightsubjectstart:"+str(rightsubjectstart))
        printmicrohomology(output1,output2,output3)


'''
        #print('\n')
        #for i in jointfinder(output1,output2,output3):
        #    print("".join(i))

        #print('\n')
'''
print(m)
#leftbreak=open('./leftbreak.csv','w+')
#rightbreak=open('./rightbreak.csv','w+')
totalbreak=open('./%s.csv' % sys.argv[1][:-4],'w+')

for item in leftbreakpoint:
    leftbreak.write(str(item)+'\n')
    totalbreak.write(str(item)+'\n')
for item in rightbreakpoint:
    rightbreak.write(str(item)+'\n')
    totalbreak.write(str(item)+'\n')

for item in totalbreakpoint:
    (a,b)=item
    totalbreak.write(str(a)+','+str(b)+'\n')

#leftbreak.close()
#rightbreak.close()
totalbreak.close()
#recombi_readings.close()
