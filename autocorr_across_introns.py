import numpy as np
import glob
import pickle
from Bio import AlignIO
import matplotlib.pyplot as plt



def autocorrelation_across_introns(arr, interface, dmax, len_exon_back, len_exon_front):                                ## defining the autocorrelation function
    acorr = np.zeros(dmax)
    
    #exon_back = arr[:interface]
    #exon_front = arr[interface:]    
    for d in range(1,int(min(dmax,len_exon_back-1, len_exon_front-1))):
        #print interface, len(arr),d, arr[interface-d:interface+1].shape, arr[interface+1:interface+1+d+1].shape, len_exon_back, len_exon_front 
        acorr[d]+=np.sum(arr[interface-d:interface+1] * arr[interface+1:interface+1+d+1])
     
    return acorr  
    
with open('new_protein_dict.pickle','r') as f:
    new_protein_dict = pickle.load(f)
    

nbins = 150    

acorr_syn = np.zeros(nbins)
acorr_nonsyn = np.zeros(nbins)
normalizing_vector = np.zeros(nbins)
#normalizing_vector_earlier = np.zeros(nbins)

mRNA_list   = glob.glob('data_contig/*mRNA/*fa')
contig_list = glob.glob('annotations/*pickle')

length_sequence = []
gc_sequence = []

mean_poly_density_syn = []
mean_poly_density_nonsyn = []

for contig in contig_list[0:1]:
    
    contig_name =  contig.split('/')[1].split('_')[0]
    
    with open(contig,'r') as f:
        annotations_list=pickle.load(f)
        
    for gene_id,anno in annotations_list[:500]:
        
        print gene_id
               
        alignment_trans= AlignIO.read('data_contig/'+gene_id.split('-')[0]+'_transcripts/'+gene_id+'.fa',"fasta")
                
        alignment_protein = AlignIO.read('data_contig/'+gene_id.split('-')[0]+'_protein/'+gene_id+'.fa',"fasta") ## reading the alignment

        
        #print '*' in str(alignment_protein[0][:-1])
        codons=[]
        
        for codon in range(alignment_trans.get_alignment_length()/3):                        ## getting codons from the transcript
            
            codons.append(map(str,[X.seq for X in alignment_trans[:,(3*codon):(codon+1)*3]]))
    
            
        matrix_trans = np.array(codons).T                                ## using transpose of the matrix
        matrix_protein  = np.array(alignment_protein)
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
        
        
        
    ###############
    
        if (new_protein_dict['data_contig/'+gene_id.split('-')[0]+'_protein/'+gene_id+'.fa'] ==1) : 
            
            interface = 0
            
            if anno['strand']=='+':
                
                
                
                for i in range(0,len(anno['exons'])-1):
                    interface+=(anno['exons'][i][1] - anno['exons'][i][0]+1)/3.0
                    
                    len_exon_back = ((anno['exons'][i][1] - anno['exons'][i][0]+1)/3)
                    len_exon_front = ((anno['exons'][i+1][1] - anno['exons'][i+1][0]+1)/3)
                    #print interface, len_exon_front, interface+len_exon_front, len(arr_uni_syn)
                    
                    acorr_syn+=autocorrelation_across_introns(arr_uni_syn, int(interface), nbins, len_exon_back, len_exon_front)
                    acorr_nonsyn+=autocorrelation_across_introns(arr_uni_nonsyn, int(interface), nbins, len_exon_back, len_exon_front)
                    normalizing_vector[:int(min(nbins,len_exon_back-1,len_exon_front-1))]+=np.arange(1,int(min(nbins,len_exon_back-1,len_exon_front-1))+1)
                    
                    #normalizing_vector_earlier[:int(min(nbins,len_exon_back-1,len_exon_front-1))]+=np.arange(1,int(min(nbins,len_exon_back-1,len_exon_front-1))+1)
            else:
                for i in range(0,len(anno['exons'])-1):
                    interface+= (anno['exons'][::-1][i][1] - anno['exons'][::-1][i][0]+1)/3.0 
                    len_exon_back = ((anno['exons'][::-1][i][1] - anno['exons'][::-1][i][0]+1)/3)
                    len_exon_front = ((anno['exons'][::-1][i+1][1] - anno['exons'][::-1][i+1][0]+1)/3)
                    
                    #print interface, len_exon_front, interface+len_exon_front, len(arr_uni_syn) 
                    
                    acorr_syn+=autocorrelation_across_introns(arr_uni_syn, int(interface), nbins, len_exon_back, len_exon_front)
                    acorr_nonsyn+=autocorrelation_across_introns(arr_uni_nonsyn, int(interface), nbins, len_exon_back, len_exon_front)
                    normalizing_vector[:int(min(nbins,len_exon_back-1,len_exon_front-1))]+=np.arange(1,int(min(nbins,len_exon_back-1,len_exon_front-1))+1)
                    
                    #normalizing_vector_earlier[:int(min(nbins,len_exon_back,len_exon_front))]+=np.arange(1,int(min(nbins,len_exon_back,len_exon_front))+1)
                 
acorr_syn/=normalizing_vector
acorr_nonsyn/=normalizing_vector

            
            
    ################
    
plt.figure()    
plt.plot(acorr_syn,label = 'syn')
plt.plot(acorr_nonsyn,label = 'nonsyn')
plt.legend()
plt.show()    
        
plt.figure()    
plt.plot(acorr_syn/np.mean(acorr_syn[-10:]), label = 'syn'+' asym: '+ str(round(np.mean(acorr_syn[-10:]), 4)))
plt.plot(acorr_nonsyn/np.mean(acorr_nonsyn[-10:]), label = 'nonsyn'+' asym: '+ str(round(np.mean(acorr_nonsyn[-10:]), 4)))
plt.legend()
plt.show()    