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
homologyset=set()

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
    global homologyset
    temp=''.join(inputdic[1])
    if temp in homologyset:
        return(-1)
    else: homologyset.add(temp)
    for i in range(len(inputdic[0])):
        if match(inputdic[0][i],inputdic[1][i]):
            walk[0].append(True)
        else: walk[0].append(False)
        if match(inputdic[1][i],inputdic[2][i]):
            walk[1].append(True)
        else: walk[1].append(False)
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
    return(homologylength)

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
    return(0)
    #Print Space between formated results
    print('\n')
    print('\n')




#recombi_readings=open('./specific_recombination.fas','w+')

for record in SeqIO.parse('./BK_DIK_NCCR2000.fas', "fasta"):
    BKloop=record.seq.upper()
    BKloopr=BKloop.reverse_complement()
genomelength=len(BKloop)
#Read BLAST result.
nonhomocatch=[]
yeshomocatch=[]
#Read Blast result
inp=open(sys.argv[1][:-4],'r')
oneline=inp.readline()
updated=False

nonhomo=open('./%sNonhomology01.csv' % sys.argv[1][:-4],'w+')
yeshomo=open('./%sYeshomo01.csv' % sys.argv[1][:-4],'w+')
#Define a function to isolate the joint.



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

    #Solve Circular issue by connect ends.
    for i in range(len(recombination_result)-1):
        #Solve Circular issue by connect ends.
        if recombination_result[i][4] == 1 and recombination_result[i+1][3] == len(BKloop):
            recombination_result.insert(i,[recombination_result[i][0],recombination_result[i][1],recombination_result[i+1][2],recombination_result[i][3],recombination_result[i+1][4],recombination_result[i][5],recombination_result[i][6]+recombination_result[i+1][6],recombination_result[i][7]+recombination_result[i+1][7],recombination_result[i][8]])
            recombination_result.pop(i+1)
            recombination_result.pop(i+1)
            #print("alart", [i[1:5] for i in recombination_result])
            break
        if recombination_result[i][4] == len(BKloop) and recombination_result[i+1][3] == 1:
            recombination_result.insert(i,[recombination_result[i][0],recombination_result[i][1],recombination_result[i+1][2],recombination_result[i][3],recombination_result[i+1][4],recombination_result[i][5],recombination_result[i][6]+recombination_result[i+1][6],recombination_result[i][7]+recombination_result[i+1][7],recombination_result[i][8]])
            recombination_result.pop(i+1)
            recombination_result.pop(i+1)
            #print("alart1", [i[1:5] for i in recombination_result])
            break

    if len(recombination_result)==1:
        continue

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




        if leftdirection=='1/1':
            output1=' '*(leftstart-1) + leftsubject + BKloop[leftsubjectend:leftsubjectend+10]+' '*(241-leftend)
        else:
            output1=' '*(leftstart-1) + leftsubject + BKloopr[genomelength-leftsubjectend+1:genomelength-leftsubjectend+11]+' '*(241-leftend)
        output2=record.seq
        if rightdirection=='1/1':
            output3=' '*(rightstart-11) +BKloop[rightsubjectstart-11:rightsubjectstart-1]+ rightsubject+' '*(251-rightend)
        else:
            output3=' '*(rightstart-11) +BKloopr[genomelength-rightsubjectstart-10:genomelength-rightsubjectstart]+ rightsubject + ' '*(251-rightend)

        temp=microhomologyfinder(jointfinder(output1,output2,output3))
        if temp == 0 or temp ==1:
            nonhomocatch.append((leftsubjectend,rightsubjectstart))
        if temp > 1:
            yeshomocatch.append((leftsubjectend,rightsubjectstart))

        #print("leftstart:"+str(leftstart)+"  rightstart:"+str(rightstart)+"    leftsubjectend:"+str(leftsubjectend)+"    rightsubjectstart:"+str(rightsubjectstart))
        #print(microhomologyfinder(jointfinder(output1,output2,output3)))
        #print(printmicrohomology(output1,output2,output3))


for item in nonhomocatch:
    (a,b)=item
    nonhomo.write(str(a)+','+str(b)+'\n')
for item in yeshomocatch:
    (a,b)=item
    yeshomo.write(str(a)+','+str(b)+'\n')
nonhomo.close()
yeshomo.close()


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



    print('\n')
    for i in jointfinder(output1,output2,output3):
        print("".join(i))

    print('\n')


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
'''
