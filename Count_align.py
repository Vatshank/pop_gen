from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import AlignIO
import numpy as np

alignment = AlignIO.read("HIV1_REF_2010_env_DNA.fasta", "fasta")
matrix= np.array(alignment)

b= matrix.T

for i in range(len(b)) :
    row = b[i,:]
    count_A = 0
    count_T = 0
    count_G = 0
    count_C = 0
    count_X = 0
    	
    for j in range(len(row)):
        if row[j] == 'A':
            count_A +=1
        elif row[j] == 'T':
            count_T +=1
        elif row[j] == 'G':
            count_G +=1     
        elif row[j] == 'C':
            count_C +=1
        elif row[j] == '-':
            count_X +=1
    print "For position", i, "in the alignment, counts are"         
    print "A=",count_A, "\tT=",count_T, "\tG=",count_G, "\tC=",count_C, "\t-=",count_X
         
        