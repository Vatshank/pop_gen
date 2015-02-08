import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
import glob

def autocorrelation(arr, dmax):                                ## defining the autocorrelation function
    acorr = np.zeros(dmax)
    for d in range(1,min(dmax,arr.shape[0])):
        acorr[d]+=np.sum(arr[d:] * arr[:-d])
     
    return acorr               
            

    
transcript_list = glob.glob('data_contig/*trans*/*fa')        ## making a list of all files in the respective directories
protein_list = glob.glob('data_contig/*pro*/*fa')



codons = {}


acorr_syn = np.zeros(150)
acorr_nonsyn = np.zeros(150)

normalizing_vector_syn = np.zeros_like(acorr_syn)
normalizing_vector_nonsyn = np.zeros_like(acorr_nonsyn)

mean_poly_density_syn = []
mean_poly_density_nonsyn = []

for count in range(100, 200): #range(len(transcript_list)): #len(transcript_list)):                    ## looping over files in the transcript list
    print transcript_list[count]
    
    if (transcript_list[count].split('/')[2] == protein_list[count].split('/')[2]) :
    
        alignment_trans = AlignIO.read(transcript_list[count],"fasta")
        alignment_protein = AlignIO.read(protein_list[count],"fasta") ## reading the alignment
        print '*' in str(alignment_protein[0][:-1])
        codons=[]
        
        for codon in range(alignment_trans.get_alignment_length()/3):                        ## getting codons from the transcript
            
            codons.append(map(str,[X.seq for X in alignment_trans[:,(3*codon):(codon+1)*3]]))
    
            
        matrix_trans = np.array(codons).T                                ## using transpose of the matrix
        matrix_protein  = np.array(alignment_protein)
        print matrix_trans.shape, matrix_protein.shape
            
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
        
    
        acorr_syn+=autocorrelation(arr_uni_syn, len(acorr_syn))                ## adding the autocorrelation values for the transripts
        acorr_nonsyn+=autocorrelation(arr_uni_nonsyn, len(acorr_nonsyn))
        
        normalizing_vector_syn+=len(arr_uni_syn)-np.arange(0,len(acorr_syn))
        normalizing_vector_nonsyn+=len(arr_uni_nonsyn)-np.arange(0,len(acorr_nonsyn))
        
    else:
        print "Well, this was unexpected"
        break   
    
acorr_syn/=normalizing_vector_syn
acorr_nonsyn/=normalizing_vector_nonsyn   
    
plt.figure()
plt.plot(acorr_syn)
plt.plot(acorr_nonsyn)
plt.plot(np.arange(len(acorr_syn)), np.mean(acorr_syn[-10:])*np.ones_like(acorr_syn))
plt.plot(np.arange(len(acorr_nonsyn)), np.mean(acorr_nonsyn[-10:])*np.ones_like(acorr_nonsyn))
