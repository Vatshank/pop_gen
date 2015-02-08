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
    


transcript_list = glob.glob('data_contig_called_in_all/*trans*/*fa')

protein_list =  glob.glob('data_contig_called_in_all/*prot*/*fa')

#with open('new_protein_list.pickle','r') as f:
#    new_protein_list = pickle.load(f)
    
with open('new_protein_dict.pickle','r') as f:
    new_protein_dict = pickle.load(f)
    
#with open('mean_coverage_values.pickle','r') as f:
#    mean_coverage_values = pickle.load(f)
    
            
nbins = 150
acorr_syn = np.zeros(nbins)
acorr_nonsyn = np.zeros(nbins)
normalizing_vector = np.zeros(nbins)




#acorr_syn ={}
#acorr_nonsyn = {}
#normalizing_vector = {}

#acorr_syn[('coverage',0,8)] = np.zeros(nbins)
#acorr_syn[('coverage',8,12)] = np.zeros(nbins)
#acorr_syn[('coverage',12,0)] = np.zeros(nbins)
#
#for key in acorr_syn:
#    acorr_nonsyn[key] = np.zeros(nbins)
#    normalizing_vector[key] = np.zeros(nbins)


mean_poly_density_syn = []
mean_poly_density_nonsyn = []

length_sequence = []
gc_sequence = []


for count in range(200):#range(len(transcript_list)):              ## looping over files in the transcript list
    print transcript_list[count]
    
    gene_id = transcript_list[count].split('/')[2]
    
    if (transcript_list[count].split('/')[2] == protein_list[count].split('/')[2]) :
    
        alignment_trans = AlignIO.read(transcript_list[count],"fasta")

        alignment_protein = AlignIO.read(protein_list[count],"fasta") ## reading the alignment
        
        Pp_seqs = np.array([ seq.id not in ['outgroup', 'RS5522B'] for seq in alignment_trans], dtype=bool)
        
        print '*' in str(alignment_protein[0][:-1])
        codons=[]
        
        for codon in range(alignment_trans.get_alignment_length()/3):                        ## getting codons from the transcript
            
            codons.append(map(str,[X.seq for X in alignment_trans[:,(3*codon):(codon+1)*3]]))
    
        codons = (np.array(codons).T)[Pp_seqs,:]
        matrix_trans = codons                               ## using transpose of the matrix
        matrix_protein  = np.array(alignment_protein)[Pp_seqs,:]
        print matrix_trans.shape, matrix_protein.shape
        
        length = matrix_trans.shape[1]
        length_sequence.append(length)
        
        gc = 1.0*(str(alignment_trans[0].seq).count('C')+str(alignment_trans[0].seq).count('G'))/(matrix_trans.shape[1]*3.0)
        gc_sequence.append(gc)
                                        
        
        uni_trans= []
        uni_protein = []        
        uni_trans_no_singleton = []
        uni_protein_no_singleton = []
                        
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
        
        
        if (new_protein_dict['data_contig/'+gene_id.split('-')[0]+'_protein/'+gene_id]==1):                ##filtering proteins in C.elegans
            acorr_syn+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn+=autocorrelation(arr_uni_nonsyn, nbins)                             
            normalizing_vector+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))
            
            
    else:
        print 'Well, this was unexpected'    
                

acorr_syn/=normalizing_vector
acorr_nonsyn/=normalizing_vector                  

plt.figure()
plt.plot(acorr_syn[1:]/np.mean(acorr_syn[-10:]), label = 'syn'+' asym: '+ str(round(np.mean(acorr_syn[-10:]), 4)))
plt.plot(acorr_nonsyn[1:]/np.mean(acorr_nonsyn[-10:]), label = 'nonsyn'+' asym: '+ str(round(np.mean(acorr_nonsyn[-10:]), 4)))           
plt.legend()
plt.show()                                                            
                                                                                                                                                                                    