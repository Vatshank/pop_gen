import numpy as np
from Bio import AlignIO
import pickle

with open('new_protein_list.pickle','r') as f:
    new_protein_list = pickle.load(f)

list_rebels = []

genes_not_in_reference = []

for count in range(len(new_protein_list)):
    print new_protein_list[count]
    
    alignment_protein = AlignIO.read(new_protein_list[count][0],"fasta")
    
    if alignment_protein[7].id!='Pristionchus_Hybrid_assembly':
        print "Gene not in reference"
        genes_not_in_reference.append([new_protein_list[count][0].split('/')[-1].split('.fa')[0],len(alignment_protein)])
    #if ('*' in str(alignment_protein[7][0][:-1])) and alignment_protein[7].id=='Pristionchus_Hybrid_assembly':
            #print 'stop codon found midaway'
            #
            #list_rebels.append(new_protein_list[count][0].split('/')[2])
    
    
