import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
import glob
import pickle
import datetime

def autocorrelation(arr, dmax):                                ## defining the autocorrelation function
    acorr = np.zeros(dmax)
    for d in range(1,min(dmax,arr.shape[0])):
        acorr[d]+=np.sum(arr[d:] * arr[:-d])
     
    return acorr               
            

    
transcript_list = glob.glob('data_contig/*trans*/*fa')        ## making a list of all files in the respective directories
protein_list = glob.glob('data_contig/*pro*/*fa')

codons = {}

acorr_syn = {}
acorr_nonsyn = {}
normalizing_vector = {}

acorr_syn[('syn', 0, 0.05)] = []
acorr_syn[('syn', 0.05, 0.10)] = []
acorr_syn[('syn', 0.10, 0.15)] = []
acorr_syn[('syn', 0.15, 0.20)] = []
acorr_syn[('syn', 0.20, 1)] = []

acorr_syn[('nonsyn', 0, 0.05)] = []
acorr_syn[('nonsyn', 0.05, 0.10)] = []
acorr_syn[('nonsyn', 0.10, 0.15)] = []
acorr_syn[('nonsyn', 0.15, 0.20)] = []
acorr_syn[('nonsyn', 0.20, 1)] = []

acorr_syn[('length', 0, 113)] = []
acorr_syn[('length', 113, 188)] = []
acorr_syn[('length', 188, 291)] = []
acorr_syn[('length', 291, 431)] = []
acorr_syn[('length', 431, 00)] = []

acorr_syn[('GC', 0, 0.48)] = []
acorr_syn[('GC', 0.48, 0.51)] = []
acorr_syn[('GC', 0.51, 0.54)] = []
acorr_syn[('GC', 0.54, 0.58)] = []
acorr_syn[('GC', 0.58, 1)] = []


nbins = 150
for key in acorr_syn:
    acorr_syn[key]= np.zeros(nbins)
    acorr_nonsyn[key] = np.zeros(nbins)
    normalizing_vector[key] = np.zeros(nbins)
    
mean_poly_density_syn = []
mean_poly_density_nonsyn = []

length_sequence = []
gc_sequence = []


for count in range(len(transcript_list)): #range(0,2000)                   ## looping over files in the transcript list
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
        
        if (mean_syn<0.05):
            acorr_syn[('syn', 0, 0.05)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('syn', 0, 0.05)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('syn', 0, 0.05)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))  
        elif (mean_syn<0.10):
            acorr_syn[('syn', 0.05, 0.10)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('syn', 0.05, 0.10)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('syn', 0.05, 0.10)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))  
        elif (mean_syn<0.15):
            acorr_syn[('syn', 0.10, 0.15)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('syn', 0.10, 0.15)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('syn', 0.10, 0.15)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))  
        elif (mean_syn<0.20):
            acorr_syn[('syn', 0.15, 0.20)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('syn', 0.15, 0.20)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('syn', 0.15, 0.20)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))  
        else:
            acorr_syn[('syn', 0.20, 1)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('syn', 0.20, 1)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('syn', 0.20, 1)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))  
        
        
        
        if (mean_nonsyn<0.05):
            acorr_syn[('nonsyn', 0, 0.05)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('nonsyn', 0, 0.05)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('nonsyn', 0, 0.05)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))  
        elif (mean_nonsyn<0.10):
            acorr_syn[('nonsyn', 0.05, 0.10)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('nonsyn', 0.05, 0.10)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('nonsyn', 0.05, 0.10)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))  
        elif (mean_nonsyn<0.15):
            acorr_syn[('nonsyn', 0.10, 0.15)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('nonsyn', 0.10, 0.15)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('nonsyn', 0.10, 0.15)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))  
        elif (mean_nonsyn<0.20):
            acorr_syn[('nonsyn', 0.15, 0.20)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('nonsyn', 0.15, 0.20)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('nonsyn', 0.15, 0.20)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))  
        else:
            acorr_syn[('nonsyn', 0.20, 1)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('nonsyn', 0.20, 1)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('nonsyn', 0.20, 1)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))     
        
        
        if (gc < 0.48):
            acorr_syn[('GC', 0, 0.48)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('GC', 0, 0.48)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('GC', 0, 0.48)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))  
        elif (gc<0.51):
            acorr_syn[('GC', 0.48, 0.51)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('GC', 0.48, 0.51)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('GC', 0.48, 0.51)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))             
        elif (gc<0.54):
            acorr_syn[('GC', 0.51, 0.54)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('GC', 0.51, 0.54)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('GC', 0.51, 0.54)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))             
        elif (gc<0.58):
            acorr_syn[('GC', 0.54, 0.58)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('GC', 0.54, 0.58)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('GC', 0.54, 0.58)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))             
        else:
            acorr_syn[('GC', 0.58, 1)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('GC', 0.58, 1)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('GC', 0.58, 1)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))  
    

        if (length<113):
            acorr_syn[('length', 0, 113)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('length', 0, 113)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('length', 0, 113)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))  
        elif(length<188):
            acorr_syn[('length', 113, 188)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('length', 113, 188)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('length', 113, 188)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))             
        elif(length<291):
            acorr_syn[('length', 188, 291)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('length', 188, 291)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('length', 188, 291)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))        
        elif(length<431):
            acorr_syn[('length', 291, 431)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('length', 291, 431)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('length', 291, 431)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))        
        else:
            acorr_syn[('length', 431, 00)]+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn[('length', 431, 00)]+=autocorrelation(arr_uni_nonsyn, nbins)
            normalizing_vector[('length', 431, 00)]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))            
        
        
    else:
        print "Well, this was unexpected"
        break   

for key in acorr_syn:    
    acorr_syn[key]/=normalizing_vector[key]
    acorr_nonsyn[key]/=normalizing_vector[key]    

with open(str(datetime.date.today())+'_polymorphism_correlation.pickle', 'w') as f:
   pickle.dump((acorr_syn, acorr_nonsyn, normalizing_vector), f)     
                        
#for key in acorr_syn:    
#    plt.figure()
#    plt.plot(acorr_syn[key][1:], label = 'syn'+str(key))
#    plt.plot(acorr_nonsyn[key][1:], label = 'nonsyn'+str(key))
#    plt.plot(np.arange(len(acorr_syn[key])), np.mean(acorr_syn[key][-10:])*np.ones_like(acorr_syn[key]))
#    plt.plot(np.arange(len(acorr_nonsyn[key])), np.mean(acorr_nonsyn[key][-10:])*np.ones_like(acorr_nonsyn[key]))
#    plt.legend()
#    plt.show()

for key in acorr_syn:    
    plt.figure()
    plt.plot(acorr_syn[key][1:]/np.mean(acorr_syn[key][-10:]), label = 'syn'+str(key)+' asym: '+ str(round(np.mean(acorr_syn[key][-10:]), 4)))
    plt.plot(acorr_nonsyn[key][1:]/np.mean(acorr_nonsyn[key][-10:]), label = 'nonsyn'+str(key)+' asym: '+ str(round(np.mean(acorr_nonsyn[key][-10:]), 4)))
    plt.legend()
    plt.show()
