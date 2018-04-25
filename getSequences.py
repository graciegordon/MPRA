'''
This script will take in the log2foldchange result file and get the gene sequence from the reference folder and then make the mutations in the sequence. This script assumes that there are 20 bp added to the begining and end of the gene, causing the mutation location to be off by 20bp. This is refered to as the template. 


Author: Gracie Gordon 2018
'''

def readFASTA(filename):
    seq=''
    with open(filename,'r') as fa:
        for line in fa:
            line=line.strip()
            if line[0]!='>':
                seq+=line
    return seq

def readResults(filename):
    labels=[]
    logs=[]
    #resultsdict={}
    with open(filename,'r') as r:
        for line in r:
            #print(line)
            line=line.strip()
            line=line.split('\t')
            #print(line)
            labels.append(line[0])
            logs.append(line[1])
    return labels,logs

def getMuts(labels,refseq):
    seqs=[]
    clipref=list(refseq)
    clipref=refseq[20:(len(refseq)-20)]
    clipref="".join(clipref)

    for l in labels:
        if l =='ref':
            seqs.append(clipref)
        else:
            curlist=l.split(',')
            newseq=list(refseq)
            #print(newseq)
            #print(curlist)
            for mut in curlist:
                curmut=mut.split(':')
                pos=int(curmut[1])-1
                ref=curmut[2][0]
                var=curmut[2][2]
                #print(pos,ref,var)
                #print(newseq[pos])
                newseq[pos]=var
                #print(newseq[pos])
            #clip off 20 bp in front and end
            #print(len(newseq))
            newseq=newseq[20:(len(newseq)-20)]
            #print(len(newseq))
            newseq="".join(newseq)
            
            seqs.append(newseq)
    return seqs


def sliceBuffer():
    assert 1==1
#get refseq
reffile="/Users/student/Documents/AhituvRotation/data/members/mkircher/access/regulatory_tests/sat_mutagenesis/data_release/processed_data/references/templates/IRF6.fa"

ref_seq=readFASTA(reffile)
#print(ref_seq)

#get Result file
oldfile="/Users/student/Documents/AhituvRotation/scripts/toyResults.txt"
labels,logs=readResults(oldfile)
#print(labels)

seqs=getMuts(labels,ref_seq)
#print(seqs)

for i,j,k in zip(labels,seqs,logs):
    print(i,j,k)































