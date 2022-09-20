#'This script is used to extract the gene matrix based on the single cell barcode list
#'The barcode list and the gene matrix is required,
#Barcode list: 
#  "AAACCCAAGAAGGTAG"
#  "AAACCCAAGACCTCCG"
#  "AAACCCACAACTTGGT"
#  "AAACCCAGTTATCTTC"
#  "AAACGAACACACTTAG"
#  "AAACGAAGTAGGCTGA"
#  "AAACGCTTCACGATAC"
#  "AAAGGATAGACCTTTG"
#  "AAAGGATAGTCTAGAA"
#   Gene matrix data
#  "AAACCCAAGAAGGTAG" "AAACCCAAGACCTCCG" "AAACCCACAACTTGGT" "AACCCAGTTATCTTC"
#  "RP11-34P13.7" 1 0 0 0
#  "FO538757.2" 2 0 0 0
#  "AP006222.2" 0 0 0 0
#  "RP4-669L17.10" 0 0 0 1
#  "RP5-857K21.4" 0 0 1 1
#  "RP11-206L10.9" 5 0 1 1
#  "RP11-54O7.16" 0 0 0 0
#  "RP11-54O7.3" 0 0 0 0
#  "SAMD11" 0 0 0 0
#  "NOC2L" 2 0 0 0
#  "KLHL17" 0 0 0 0


import os
import re
import sys


feature=open(sys.argv[1],'r')#Single cell barcode list
TrEx1=open(sys.argv[2],'r')#scRNA-Seq matrix
TrEx2=open(sys.argv[2],'r')#scRNA-Seq matrix
subEx=open(sys.argv[3],'w')#ouput of matrix data


mkl=[]
mkls={}


for line in feature:
    line=line.strip()
    line=line.split()
    mkl.append(line[0])
    mkls[line[0]]=line[0]


nte=[]

for line in TrEx1:
    line=line.strip()
    lineE=line.split()
    a=0
    if '"' in lineE[1]:
        for n in range(0,int(len(lineE))):
                a=a+1
                if lineE[n] in mkls:
                    nte.append(a)


#getting the expression of features
cellID=''
exp=''

for line in TrEx2:
    line=line.strip()
    lineE=line.split()
    if '"' in lineE[1]:
        for id in range(0,int(len(lineE))):
            if lineE[id] in mkls.keys():
                cellID=cellID+lineE[id]+'       '
        print >>subEx,cellID
        cellID

    if '"' not in lineE[1]:
        
        for value in nte:
            #print value,len(lineE)
            exp=exp+lineE[int(value)]+' '
        
        print >>subEx,lineE[0]+'        '+exp
        exp=''
        

feature.close()
TrEx1.close()
TrEx2.close()
subEx.close()

