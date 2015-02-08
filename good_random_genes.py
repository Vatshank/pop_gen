import numpy as np
from Bio import AlignIO
import pickle

good_genes = ['Contig100-snap.103',
 'Contig194-snap.4',
 'Contig111-snap.38',
 'Contig162-snap.15',
 'Contig116-snap.35',
 'Contig304-snap.5',
 'Contig125-snap.25',
 'Contig176-snap.11',
 'Contig106-snap.58',
 'Contig106-snap.57',
 'Contig134-snap.33',
 'Contig146-snap.10']

#good_genes = ['Contig304-snap.5']

list_rebels =[] 
   
dict_random = {}

for count in range(len(good_genes)): #range(2000)            ## looping over files in the transcript list
    print good_genes[count]
    
    dict_random[good_genes[count]]= {}
        
    #gene_id = good_genes[count]

    alignment_trans = AlignIO.read('data_contig/'+good_genes[count].split('-')[0]+'_transcripts/'+good_genes[count]+'.fa',"fasta")
    
    alignment_protein = AlignIO.read('data_contig/'+good_genes[count].split('-')[0]+'_protein/'+good_genes[count]+'.fa',"fasta") ## reading the alignment
    
                
    if '*' in str(alignment_protein[0][:-1]):
        print 'stop codon found midway'
        list_rebels.append(good_genes[count])
        
    codons=[]
    
    for codon in range(alignment_trans.get_alignment_length()/3):                        ## getting codons from the transcript
        
        codons.append(map(str,[X.seq for X in alignment_trans[:,(3*codon):(codon+1)*3]]))

    Pp_seqs = np.array([ seq.id not in ['outgroup.fa', 'RS5522B'] for seq in alignment_trans], dtype=bool)
    
    codons = (np.array(codons).T)[Pp_seqs,:]    
                
    matrix_trans = codons                                ## using transpose of the matrix
    matrix_protein  = np.array(alignment_protein)[Pp_seqs,:]  
    print matrix_trans.shape, matrix_protein.shape
    
    #length = matrix_trans.shape[1]
    #length_sequence.append(length)
    #
    #gc = 1.0*(str(alignment_trans[0].seq).count('C')+str(alignment_trans[0].seq).count('G'))/(matrix_trans.shape[1]*3.0)
    #gc_sequence.append(gc)
                                    
    
    uni_trans= []
    uni_protein = []                                
    uni_trans_no_singleton = []
    uni_protein_no_singleton = []
                    
    for n in range(matrix_trans.shape[1]):                            ## getting uni, which has the number of different types codons at each site	        
        number_poly_trans = np.unique(matrix_trans[:,n]).shape[0]
        
        uni_trans.append(number_poly_trans)
        
        if (number_poly_trans ==1):                                 ##removing singletons
            uni_trans_no_singleton.append(number_poly_trans)
        else:
            uni_column = np.unique(matrix_trans[:,n])
            sum_uni_trans = []
            for nucleotide in uni_column:
                sum_uni_trans.append(sum(matrix_trans[:,n]==nucleotide))
            
            uni_trans_no_singleton.append(number_poly_trans - sum_uni_trans.count(1))        
            
    for n in range(matrix_protein.shape[1]):
        number_poly_protein = np.unique(matrix_protein[:,n]).shape[0]
        
        uni_protein.append(number_poly_protein)
        
        if (number_poly_protein ==1):
            uni_protein_no_singleton.append(number_poly_protein)
        else:
            uni_column = np.unique(matrix_protein[:,n])
            sum_uni_protein = []
            for amino_acid in uni_column:
                sum_uni_protein.append(sum(matrix_protein[:,n]==amino_acid))
            
            uni_protein_no_singleton.append(number_poly_protein - sum_uni_protein.count(1))
    
    arr_uni_trans = np.array(uni_trans)>1                                ## same as->  arr_uni[arr_uni==1] = 0 #arr_uni[arr_uni>1] = 1  
    arr_uni_protein = np.array(uni_protein)>1
    arr_uni_trans_no_singleton = np.array(uni_trans_no_singleton)>1                                ## same as->  arr_uni[arr_uni==1] = 0 #arr_uni[arr_uni>1] = 1  
    arr_uni_protein_no_singleton  = np.array(uni_protein_no_singleton)>1
    
    arr_uni_nonsyn = arr_uni_trans * arr_uni_protein
    arr_uni_syn = arr_uni_trans * (1-arr_uni_protein)
    arr_uni_nonsyn_no_singleton = arr_uni_trans_no_singleton * arr_uni_protein_no_singleton
    arr_uni_syn_no_singleton = arr_uni_trans_no_singleton * (1-arr_uni_protein_no_singleton)
    
    dict_random[good_genes[count]] = {'syn': arr_uni_syn, 'nonsyn': arr_uni_nonsyn, 'syn_no_singleton' : arr_uni_syn_no_singleton,
                                       'nonsyn_no_singleton': arr_uni_nonsyn_no_singleton}