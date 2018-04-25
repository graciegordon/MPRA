'''
Script that will process an MPRA. It will take in the UMI counts and the assignments and calculate the log2FoldChange for each sequence. This will operate on one gene of interest at a time

Author: Gracie Gordon 2018
'''

import numpy as np

def MakeAssignments(filename):
    #create dictionary mapping the mutation to the barcode
    with open(filename, 'r') as assignments:
        barcode_key={}
        for line in assignments:
            curline=line.split()
            if curline[1:] == []:
                barcode_key[curline[0]]=['ref']
            else:
                barcode_key[curline[0]]=curline[1:]
    return barcode_key

def GetCounts(filename):
    #make dictionary of dictionaries where the mutationbarcode:otherbarcode:count
    with open(filename, 'r') as f:
        counts={}
        for line in f:
            curline=line.split()
            if curline[0] in counts:
                #print('org')
                #print(counts[curline[0]])
                
                #print('org',counts[curline[0]])
                temp=counts[curline[0]]
                temp[curline[1]]=curline[2]
                counts[curline[0]]=temp
                #print('new',counts[curline[0]])
        
                #temp=counts[curline[0]]
                #temp.append((curline[1],curline[2]))
                #counts[curline[0]]=temp
                #print('new')
                #print(temp)
            else:
                curdict={}
                curdict[curline[1]]=curline[2]
                counts[curline[0]]=curdict
                #print(curdict)
                #counts[curline[0]]=[(curline[1],curline[2])]

            #counts.append((curline[0],curline[1],curline[2]))

    return counts

def log2Calc(RNA_counts,DNA_counts):
    #calculate the log2fold change for each mutation barcode, taking the log2(RNA/DNA) 
    #for each barcode then the average for mutation barcodes
    log2foldchange={}
    alllog2=[]
    for key in RNA_counts:
        #print('key', key)
        if key in DNA_counts:
            RNA=RNA_counts[key]
            DNA=DNA_counts[key]
            #print('rna', RNA)
            #print('dna', DNA)
            barcode=[]
            
            for key2 in RNA:
                if key2 in DNA:
                    countDNA=int(DNA[key2])   
                    countRNA=int(RNA[key2])
                    #print('dna','rna')
                    #print(countDNA,countRNA)
                    fold=np.log2(countRNA/countDNA)
                    barcode.append(fold)
            
            if barcode != []:
                #print(barcode)
                total=sum(barcode)/float(len(barcode))
                #print(total)
                alllog2.append(total)
                log2foldchange[key]=total
    #print(log2foldchange)
    return log2foldchange,alllog2

def MutAssign(barcode_dict,logfolddict):
    #prepare output of mutations to logfold
    output=[]
    for key in logfolddict:
        #print(key)
        if key in barcode_dict:
            mut=barcode_dict[key]
            fold=logfolddict[key]
            output.append((mut,fold))
        #print(output)
    return output

def calcPvalue(allog2,logfolddict):
    #fix this, not right
    null=0
    log_mean=np.mean(alllog2)
    log_std=np.std(alllog2)
    obs=len(alllog2)
    

def reverseComplement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

#define variables
filenameVars="/Users/student/Documents/AhituvRotation/data/members/mkircher/access/regulatory_tests/sat_mutagenesis/data_release/processed_data/assignments/IRF6.variants.txt"
#filenameCountsDNA="/Users/student/Documents/AhituvRotation/data/members/mkircher/access/regulatory_tests/sat_mutagenesis/data_release/processed_data/UMI_tag_counts/IRF6-DNA-1.tsv"
#filenameCountsRNA="/Users/student/Documents/AhituvRotation/data/members/mkircher/access/regulatory_tests/sat_mutagenesis/data_release/processed_data/UMI_tag_counts/IRF6-RNA-1.tsv"

filenameCountsDNA="/Users/student/Documents/AhituvRotation/data/members/mkircher/access/regulatory_tests/sat_mutagenesis/data_release/processed_data/UMI_tag_counts/IRF6_DNA_all.tsv"
filenameCountsRNA="/Users/student/Documents/AhituvRotation/data/members/mkircher/access/regulatory_tests/sat_mutagenesis/data_release/processed_data/UMI_tag_counts/IRF6_RNA_all.tsv"
#filenameCountsDNA="toyIRF6_DNA.tsv"
#filenameCountsRNA="toyIRF6_RNA.tsv"

barcode_dict=MakeAssignments(filenameVars)
#print(barcode_dict)
RNA_counts=GetCounts(filenameCountsRNA)
#print(RNA_counts)
DNA_counts=GetCounts(filenameCountsDNA)

logfolddict,alllog2=log2Calc(RNA_counts, DNA_counts)



final=MutAssign(barcode_dict,logfolddict)
#print(final)
with open('IRF6_log2fold.txt','w') as out:
    #print('final')
    
    for i in final:
        muts=''
        for x in i[0]:
            muts+=x
            muts+=','
        muts=muts[:-1]
        print(str(muts)+"\t"+str(i[1])+"\n")
        out.write(str(muts)+"\t"+str(i[1])+"\n")

out.close()

#print(RNA_counts['AAAAAACTTTTTAAG'])
#print(DNA_counts['AAAAAACTTTTTAAG'])






