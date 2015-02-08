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
    
def crosscorrelation(arr_syn,arr_nonsyn,dmax):
    crosscorr = np.zeros(2*dmax)
    for d in range(1,min(dmax,arr_syn.shape[0])):
        crosscorr[d]+=np.sum(arr_syn[d:] * arr_nonsyn[:-d])
        crosscorr[d+dmax]+=np.sum(arr_nonsyn[d:] * arr_syn[:-d])
     
    return crosscorr

transcript_list = glob.glob('data_contig/*trans*/*fa')

with open('new_protein_list.pickle','r') as f:
    new_protein_list = pickle.load(f)
        
            
nbins = 100
acorr_syn = np.zeros(nbins)
acorr_nonsyn = np.zeros(nbins)
acorr_syn_no_singleton= np.zeros(nbins)
acorr_nonsyn_no_singleton = np.zeros(nbins)

crosscorr = np.zeros(2*nbins)
crosscorr_no_singleton = np.zeros(2*nbins)
normalizing_vector = np.zeros(nbins)
normalizing_vector_no_singleton = np.zeros(nbins)

mean_poly_density_syn = []
mean_poly_density_nonsyn = []

length_sequence = []
gc_sequence = []

list_rebels = []

#gene_scores = {}

for count in range(len(transcript_list[:5000])): #range(2000)            ## looping over files in the transcript list
    print transcript_list[count]
    
    if (transcript_list[count].split('/')[2] == new_protein_list[count][0].split('/')[2]) :
        
        gene_id = transcript_list[count].split('/')[2].split('.fa')[0]
    
        alignment_trans = AlignIO.read(transcript_list[count],"fasta")
        
        alignment_protein = AlignIO.read(new_protein_list[count][0],"fasta") ## reading the alignment
        
                    
        if '*' in str(alignment_protein[0][:-1]):
            print 'stop codon found midaway'
            list_rebels.append(transcript_list[count].split('/')[2])
            
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
        arr_uni_protein  = np.array(uni_protein)>1
        arr_uni_trans_no_singleton = np.array(uni_trans_no_singleton)>1                                ## same as->  arr_uni[arr_uni==1] = 0 #arr_uni[arr_uni>1] = 1  
        arr_uni_protein_no_singleton  = np.array(uni_protein_no_singleton)>1
        
        arr_uni_nonsyn = arr_uni_trans * arr_uni_protein
        arr_uni_syn = arr_uni_trans * (1-arr_uni_protein)
        arr_uni_nonsyn_no_singleton = arr_uni_trans_no_singleton * arr_uni_protein_no_singleton
        arr_uni_syn_no_singleton = arr_uni_trans_no_singleton * (1-arr_uni_protein_no_singleton)                                                                                                                                                                                                                 
        
    #    arr_uni_nonsyn_a = arr_uni_trans * arr_uni_protein
    #    arr_uni_syn_a = arr_uni_trans * (~arr_uni_protein)	
    #    
    #    scramble_tmp = np.random.rand(arr_uni_nonsyn_a.shape[0])>0.5
    #
    #    arr_uni_nonsyn = scramble_tmp*arr_uni_nonsyn_a + (~scramble_tmp)*arr_uni_syn_a
    #    arr_uni_syn = (~scramble_tmp)*arr_uni_nonsyn_a + (scramble_tmp)*arr_uni_syn_a      
                        
        mean_poly_density_syn.append(np.mean(arr_uni_syn))
        mean_poly_density_nonsyn.append(np.mean(arr_uni_nonsyn))
        
        mean_syn = np.mean(arr_uni_syn)
        mean_nonsyn = np.mean(arr_uni_nonsyn)
        
        mean_syn_no_singleton = np.mean(arr_uni_syn_no_singleton)
        mean_nonsyn_no_singleton = np.mean(arr_uni_nonsyn_no_singleton)
        
        if(gc>0.50) and (length>200) and (mean_nonsyn<0.15) and (new_protein_list[count][1]==1): 
            tmp1 = autocorrelation(arr_uni_syn, nbins)
            acorr_syn+=tmp1
            tmp2 = autocorrelation(arr_uni_nonsyn, nbins)
            acorr_nonsyn+=tmp2
            tmp3 = crosscorrelation(arr_uni_syn,arr_uni_nonsyn, nbins)                 
            crosscorr+=tmp3        
            tmp4=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))
            normalizing_vector+=tmp4
            #gene_scores[gene_id] = [tmp1, tmp2, tmp3, tmp4, gc, length, mean_nonsyn, mean_syn]
            
        if (new_protein_list[count][1]==1):
            tmp1 = autocorrelation(arr_uni_syn_no_singleton, nbins)
            acorr_syn_no_singleton+=tmp1
            tmp2 = autocorrelation(arr_uni_nonsyn_no_singleton, nbins)
            acorr_nonsyn_no_singleton+=tmp2
            tmp3 = crosscorrelation(arr_uni_syn_no_singleton,arr_uni_nonsyn_no_singleton, nbins)                 
            crosscorr_no_singleton+=tmp3        
            tmp4=np.maximum(len(arr_uni_nonsyn_no_singleton)-np.arange(0,nbins), np.zeros(nbins))
            normalizing_vector_no_singleton+=tmp4
                
         
    else:
        print 'Well, this was unexpected'    
                
#with open('corr_vectors.pickle', 'w') as outfile:
#    pickle.dump(gene_scores, outfile)

acorr_syn/=normalizing_vector
acorr_nonsyn/=normalizing_vector
crosscorr[:nbins]/=normalizing_vector  
crosscorr[nbins:2*nbins]/=normalizing_vector 

acorr_syn_no_singleton/=normalizing_vector_no_singleton
acorr_nonsyn_no_singleton/=normalizing_vector_no_singleton
crosscorr_no_singleton[:nbins]/=normalizing_vector_no_singleton  
crosscorr_no_singleton[nbins:2*nbins]/=normalizing_vector_no_singleton

#with open(str(datetime.date.today())+'_corr_results_filtered_genes.pickle','w') as f:
#    pickle.dump((acorr_syn, acorr_nonsyn, crosscorr, normalizing_vector),f)               

plt.figure()
plt.plot(acorr_syn[1:]/np.mean(acorr_syn[-10:]), label = 'syn'+' asym: '+ str(round(np.mean(acorr_syn[-10:]), 4)))
plt.plot(acorr_nonsyn[1:]/np.mean(acorr_nonsyn[-10:]), label = 'nonsyn'+' asym: '+ str(round(np.mean(acorr_nonsyn[-10:]), 4)))
plt.plot(crosscorr[1:nbins]/np.mean(crosscorr[-10+nbins:nbins]), label = 'crosscorr'+' asym: '+ str(round(np.mean(crosscorr[-10:]), 4)))           
plt.legend()
plt.show() 

plt.figure()
plt.plot(acorr_syn_no_singleton[1:]/np.mean(acorr_syn_no_singleton[-10:]), label = 'syn_no_singleton'+' asym: '+ str(round(np.mean(acorr_syn_no_singleton[-10:]), 4)))
plt.plot(acorr_nonsyn_no_singleton[1:]/np.mean(acorr_nonsyn_no_singleton[-10:]), label = 'nonsyn_no_singleton'+' asym: '+ str(round(np.mean(acorr_nonsyn_no_singleton[-10:]), 4)))
plt.plot(crosscorr_no_singleton[1:nbins]/np.mean(crosscorr_no_singleton[-10+nbins:nbins]), label = 'crosscorr_no_singleton'+' asym: '+ str(round(np.mean(crosscorr_no_singleton[-10:]), 4)))           
plt.legend()
plt.show()