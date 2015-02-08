#### Calculates autocorr for elegans genes which have gc > 0.5, length > 200 and mean_nonsyn_polymorphism_density < 0.15 ####

import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
import glob
import pickle
import datetime
from Bio.Align import MultipleSeqAlignment


def autocorrelation(arr, dmax):                                ## defining the autocorrelation function
    acorr = np.zeros(dmax)
    for d in range(1,min(dmax,arr.shape[0])):
        acorr[d]+=np.sum(arr[d:] * arr[:-d])
     
    return acorr  
    


transcript_list = glob.glob('data_contig/*trans*/*fa')

with open('new_protein_list.pickle','r') as f:
    new_protein_list = pickle.load(f)
        
            
nbins = 150
acorr_syn = np.zeros(nbins)
acorr_nonsyn = np.zeros(nbins)
normalizing_vector = np.zeros(nbins)


mean_poly_density_syn = []
mean_poly_density_nonsyn = []

length_sequence = []
gc_sequence = []


for count in range(len(transcript_list)): #range(2000)            ## looping over files in the transcript list
    print transcript_list[count]
    
    if (transcript_list[count].split('/')[2] == new_protein_list[count][0].split('/')[2]) :
    
        alignment_trans = AlignIO.read(transcript_list[count],"fasta")
        
        alignment_protein = AlignIO.read(new_protein_list[count][0],"fasta") ## reading the alignment
        
                    
        print '*' in str(alignment_protein[0][:-1])
        codons=[]
        
        for codon in range(alignment_trans.get_alignment_length()/3):                        ## getting codons from the transcript
            
            codons.append(map(str,[X.seq for X in alignment_trans[:,(3*codon):(codon+1)*3]]))
    
        Pp_seqs = np.array([ seq.id not in ['outgroup.fa', 'RS5522B'] for seq in alignment_trans], dtype=bool)
        
        codons = (np.array(codons).T)[Pp_seqs,:]    
                    
        matrix_trans = codons                                ## using transpose of the matrix
        matrix_protein  = np.array(alignment_protein)[Pp_seqs,:]  
        print matrix_trans.shape, matrix_protein.shape
        
        length = matrix_trans.shape[1]
        length_sequence.append(length)
        
        gc = 1.0*(str(alignment_trans[0].seq).count('C')+str(alignment_trans[0].seq).count('G'))/(matrix_trans.shape[1]*3.0)
        gc_sequence.append(gc)
                                        
        
        uni_trans= []
        uni_protein = []        

                        
        for n in range(matrix_trans.shape[1]):                            ## getting uni, which has the number of different types codons at each site
            uni_trans.append(np.unique(matrix_trans[:,n]).shape[0])	        
        
        for n in range(matrix_protein.shape[1]):
            uni_protein.append(np.unique(matrix_protein[:,n]).shape[0])
        
        arr_uni_trans = np.array(uni_trans)>1                                ## same as->  arr_uni[arr_uni==1] = 0 #arr_uni[arr_uni>1] = 1  
        arr_uni_protein  = np.array(uni_protein)>1
                                                                                                                
        
        arr_uni_nonsyn = arr_uni_trans * arr_uni_protein
        arr_uni_syn = arr_uni_trans * (1-arr_uni_protein)
        
        mean_poly_density_syn.append(np.mean(arr_uni_syn))
        mean_poly_density_nonsyn.append(np.mean(arr_uni_nonsyn))
        
        mean_syn = np.mean(arr_uni_syn)
        mean_nonsyn = np.mean(arr_uni_nonsyn)
        
        if(gc>0.50) and (length>200) and (mean_nonsyn<0.15) and (new_protein_list[count][1]==1): 
            acorr_syn+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn+=autocorrelation(arr_uni_nonsyn, nbins)                             
            normalizing_vector+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))

    else:
        print 'Well, this was unexpected'    
                

acorr_syn/=normalizing_vector
acorr_nonsyn/=normalizing_vector   

with open(str(datetime.date.today())+'_results_filtered_genes.pickle','w') as f:
    pickle.dump((acorr_syn, acorr_nonsyn, normalizing_vector),f)               

plt.figure()
plt.plot(acorr_syn[1:]/np.mean(acorr_syn[-10:]), label = 'syn'+' asym: '+ str(round(np.mean(acorr_syn[-10:]), 4)))
plt.plot(acorr_nonsyn[1:]/np.mean(acorr_nonsyn[-10:]), label = 'nonsyn'+' asym: '+ str(round(np.mean(acorr_nonsyn[-10:]), 4)))           
plt.legend()
plt.show()                                                            
                                                                                                                                                                                