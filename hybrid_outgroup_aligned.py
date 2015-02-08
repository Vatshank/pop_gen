import numpy as np
import glob
import pickle
from Bio import AlignIO
import matplotlib.pyplot as plt
from Bio.Align import MultipleSeqAlignment



def autocorrelation(arr, dmax):                                ## defining the autocorrelation function
    acorr = np.zeros(dmax)
    for d in range(1,min(dmax,arr.shape[0])):
        acorr[d]+=np.sum(arr[d:] * arr[:-d])
     
    return acorr   
    
with open('new_protein_dict.pickle','r') as f:
    new_protein_dict = pickle.load(f)
    
with open('alignment_to_outgroup.pickle','r') as f:
    alignment_dict = pickle.load(f)

nbins = 150    

acorr_syn = np.zeros(nbins)
acorr_nonsyn = np.zeros(nbins)
normalizing_vector = np.zeros(nbins)

mRNA_list   = glob.glob('data_contig/*mRNA/*fa')
contig_list = glob.glob('annotations/*pickle')

length_sequence = []
gc_sequence = []

mean_poly_density_syn = []
mean_poly_density_nonsyn = []

for contig in contig_list:
    
    contig_name =  contig.split('/')[1].split('_')[0]
    
    with open(contig,'r') as f:
        annotations_list=pickle.load(f)
        
    for gene_id,anno in annotations_list:
        
        print gene_id
               
        alignment_trans= AlignIO.read('data_contig/'+gene_id.split('-')[0]+'_transcripts/'+gene_id+'.fa',"fasta")
                
        alignment_protein = AlignIO.read('data_contig/'+gene_id.split('-')[0]+'_protein/'+gene_id+'.fa',"fasta") ## reading the alignment

        if (len(alignment_trans)==106) and (len(alignment_protein)==106) and (alignment_trans[7].id == 'Pristionchus_Hybrid_assembly') and (alignment_trans[105].id == 'outgroup.fa') and (alignment_protein[7].id == 'Pristionchus_Hybrid_assembly') and (alignment_protein[105].id == 'outgroup.fa') :
            
            trans_hybrid_assembly = alignment_trans[7]
            trans_outgroup = alignment_trans[105]
        
            new_alignment_trans = MultipleSeqAlignment([trans_hybrid_assembly,trans_outgroup]) ##transcript alignment of hybrid assembly and outgroup
            
            #alignment_protein = AlignIO.read(new_protein_list[count][0],"fasta") ## reading the alignment
            protein_hybrid_assembly = alignment_protein[7]
            protein_outgroup = alignment_protein[105]
            
            new_alignment_protein = MultipleSeqAlignment([protein_hybrid_assembly,protein_outgroup])    ##protein alignment of hybrid assembly and outgroup
        
                
            print '*' in str(alignment_protein[0][:-1])
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
                
             ####### Compare the gene position with the values in the alignment dictionary#########       
                    
                    
                    acorr_syn+=autocorrelation(arr_uni_syn, nbins)
                    acorr_nonsyn+=autocorrelation(arr_uni_nonsyn, nbins)                             
                    normalizing_vector+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))
                    
        ################
                    
        else:
            continue
                    
                
        
acorr_syn/=normalizing_vector
acorr_nonsyn/=normalizing_vector

        
plt.figure()    
plt.plot(acorr_syn/np.mean(acorr_syn[-10:]), label = 'syn'+' asym: '+ str(round(np.mean(acorr_syn[-10:]), 4)))
plt.plot(acorr_nonsyn/np.mean(acorr_nonsyn[-10:]), label = 'nonsyn'+' asym: '+ str(round(np.mean(acorr_nonsyn[-10:]), 4)))
plt.legend()
plt.show()    