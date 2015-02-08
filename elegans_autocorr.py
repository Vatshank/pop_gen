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
#protein_list = glob.glob('data_contig/*pro*/*fa')

with open('new_protein_list.pickle','r') as f:
    new_protein_list = pickle.load(f)

#with open('real_proteins.pickle','r') as f:
#    elegans_list = pickle.load(f)

categories = { 'length': [ [0, 113], [113,188], [188,291], [291,431], [431,99999]],
                'GC': [[0,0.48], [0.48,0.51], [0.51,0.54], [0.54,0.58], [0.58,1]],
                'syn': [[0,0.05], [0.05,0.10],[0.10,0.15],[0.15,0.20],[0.20,1]],
                'nonsyn':[[0,0.05], [0.05,0.10],[0.10,0.15],[0.15,0.20],[0.20,1]]}
                

codons = {}                                                    ##defining dictionaries

acorr_syn = {}
acorr_nonsyn = {}
acorr_syn_no_singleton = {}
acorr_nonsyn_no_singleton = {}
normalizing_vector = {}
normalizing_vector_no_singleton = {}

acorr_syn_elegans = {}
acorr_nonsyn_elegans = {}
normalizing_vector_elegans = {}
acorr_syn_no_elegans = {}
acorr_nonsyn_no_elegans = {}
normalizing_vector_no_elegans = {}



nbins = 150



acorr_syn[('syn', 0, 0.05)] = []                                ##choosing keys for conditioning on syn polymorphism density
acorr_syn[('syn', 0.05, 0.10)] = []
acorr_syn[('syn', 0.10, 0.15)] = []
acorr_syn[('syn', 0.15, 0.20)] = []
acorr_syn[('syn', 0.20, 1)] = []
acorr_syn[('nonsyn', 0, 0.05)] = []                            ##choosing keys for conditioning on non-syn polymorphism density
acorr_syn[('nonsyn', 0.05, 0.10)] = []
acorr_syn[('nonsyn', 0.10, 0.15)] = []
acorr_syn[('nonsyn', 0.15, 0.20)] = []
acorr_syn[('nonsyn', 0.20, 1)] = []
acorr_syn[('length', 0, 113)] = []                            ##choosing keys for conditioning on protein length
acorr_syn[('length', 113, 188)] = []
acorr_syn[('length', 188, 291)] = []
acorr_syn[('length', 291, 431)] = []
acorr_syn[('length', 431, 99999)] = []
acorr_syn[('GC', 0, 0.48)] = []                                ##choosing keys for conditioning on gc content
acorr_syn[('GC', 0.48, 0.51)] = []
acorr_syn[('GC', 0.51, 0.54)] = []
acorr_syn[('GC', 0.54, 0.58)] = []
acorr_syn[('GC', 0.58, 1)] = []

acorr_syn_total = np.zeros(nbins)
acorr_nonsyn_total = np.zeros(nbins)
normalizing_vector_total = np.zeros(nbins)


acorr_syn_no_singleton_total = np.zeros(nbins)
acorr_nonsyn_no_singleton_total = np.zeros(nbins)
normalizing_vector_no_singleton_total = np.zeros(nbins)

acorr_syn_elegans_total = np.zeros(nbins)
acorr_nonsyn_elegans_total = np.zeros(nbins)
normalizing_vector_elegans_total = np.zeros(nbins)

acorr_syn_no_elegans_total = np.zeros(nbins)
acorr_nonsyn_no_elegans_total = np.zeros(nbins)
normalizing_vector_no_elegans_total = np.zeros(nbins)

for key in acorr_syn:	                                        ##setting initial lists to zero
    acorr_syn[key]= np.zeros(nbins)
    acorr_nonsyn[key] = np.zeros(nbins)
    acorr_syn_no_singleton[key]= np.zeros(nbins)
    acorr_nonsyn_no_singleton[key] = np.zeros(nbins)
    normalizing_vector[key] = np.zeros(nbins)
    normalizing_vector_no_singleton[key] = np.zeros(nbins)

    acorr_syn_elegans[key]= np.zeros(nbins)
    acorr_nonsyn_elegans[key] = np.zeros(nbins)
    normalizing_vector_elegans[key] = np.zeros(nbins)
    acorr_syn_no_elegans[key]= np.zeros(nbins)
    acorr_nonsyn_no_elegans[key] = np.zeros(nbins)
    normalizing_vector_no_elegans[key] = np.zeros(nbins)

mean_poly_density_syn = []
mean_poly_density_nonsyn = []

length_sequence = []
gc_sequence = []


for count in range(len(transcript_list)): #range(0,3000)             ## looping over files in the transcript list
    print transcript_list[count]
    
    if (transcript_list[count].split('/')[2] == new_protein_list[count][0].split('/')[2]) :
    
        alignment_trans = AlignIO.read(transcript_list[count],"fasta")
        alignment_protein = AlignIO.read(new_protein_list[count][0],"fasta") ## reading the alignment
        
        Pp_seqs = np.array([ seq.id not in ['outgroup.fa', 'RS5522B'] for seq in alignment_trans], dtype=bool)
        
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
        
        mean_poly_density_syn.append(np.mean(arr_uni_syn))
        mean_poly_density_nonsyn.append(np.mean(arr_uni_nonsyn))
        
        mean_syn = np.mean(arr_uni_syn)
        mean_nonsyn = np.mean(arr_uni_nonsyn)
        
        mean_syn_no_singleton = np.mean(arr_uni_syn_no_singleton)
        mean_nonsyn_no_singleton = np.mean(arr_uni_nonsyn_no_singleton)
        
        gene_properties = {'length':length, 'GC':gc, 'syn':mean_syn, 'nonsyn':mean_nonsyn}        
        
        gene_properties_no_singleton = {'syn': mean_syn_no_singleton, 'nonsyn': mean_nonsyn_no_singleton}
        
        
        acorr_syn_total+= autocorrelation(arr_uni_syn, nbins)
        acorr_nonsyn_total+= autocorrelation(arr_uni_nonsyn, nbins)
        normalizing_vector_total+= np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))
        
        acorr_syn_no_singleton_total+= autocorrelation(arr_uni_syn_no_singleton, nbins)
        acorr_nonsyn_no_singleton_total+= autocorrelation(arr_uni_nonsyn_no_singleton, nbins)
        normalizing_vector_no_singleton_total+= np.maximum(len(arr_uni_nonsyn_no_singleton)-np.arange(0,nbins), np.zeros(nbins))
        
        if (new_protein_list[count][1] == 1):
            acorr_syn_elegans_total+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn_elegans_total+=autocorrelation(arr_uni_nonsyn, nbins)                             
            normalizing_vector_elegans_total+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))
        if (new_protein_list[count][1] == 0):       
            acorr_syn_no_elegans_total+=autocorrelation(arr_uni_syn, nbins)
            acorr_nonsyn_no_elegans_total+=autocorrelation(arr_uni_nonsyn, nbins)                             
            normalizing_vector_no_elegans_total+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))
        
                
        for categ in categories:
            
            for interval in categories[categ]:
                                                
                if categ == 'GC' or categ == 'length':                
                    if gene_properties[categ]>interval[0] and gene_properties[categ]<=interval[1]:
                        acorr_syn[(categ,interval[0],interval[1])]+=autocorrelation(arr_uni_syn, nbins)
                        acorr_nonsyn[(categ,interval[0],interval[1])]+=autocorrelation(arr_uni_nonsyn, nbins)                             
                        normalizing_vector[(categ,interval[0],interval[1])]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))    
                        
                        acorr_syn_no_singleton[(categ,interval[0],interval[1])]+=autocorrelation(arr_uni_syn_no_singleton, nbins)
                        acorr_nonsyn_no_singleton[(categ,interval[0],interval[1])]+=autocorrelation(arr_uni_nonsyn_no_singleton, nbins)                             
                        normalizing_vector_no_singleton[(categ,interval[0],interval[1])]+=np.maximum(len(arr_uni_nonsyn_no_singleton)-np.arange(0,nbins), np.zeros(nbins)) 
                    
                        if (new_protein_list[count][1] == 1):
                            acorr_syn_elegans[(categ,interval[0],interval[1])]+=autocorrelation(arr_uni_syn, nbins)
                            acorr_nonsyn_elegans[(categ,interval[0],interval[1])]+=autocorrelation(arr_uni_nonsyn, nbins)                             
                            normalizing_vector_elegans[(categ,interval[0],interval[1])]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))
                        if (new_protein_list[count][1] == 0):       
                            acorr_syn_no_elegans[(categ,interval[0],interval[1])]+=autocorrelation(arr_uni_syn, nbins)
                            acorr_nonsyn_no_elegans[(categ,interval[0],interval[1])]+=autocorrelation(arr_uni_nonsyn, nbins)                             
                            normalizing_vector_no_elegans[(categ,interval[0],interval[1])]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))                    
                
                else:
                    
                    if gene_properties[categ]>interval[0] and gene_properties[categ]<=interval[1]:
                        acorr_syn[(categ,interval[0],interval[1])]+=autocorrelation(arr_uni_syn, nbins)
                        acorr_nonsyn[(categ,interval[0],interval[1])]+=autocorrelation(arr_uni_nonsyn, nbins)                             
                        normalizing_vector[(categ,interval[0],interval[1])]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))    
                        
                        if (new_protein_list[count][1] == 1):
                            acorr_syn_elegans[(categ,interval[0],interval[1])]+=autocorrelation(arr_uni_syn, nbins)
                            acorr_nonsyn_elegans[(categ,interval[0],interval[1])]+=autocorrelation(arr_uni_nonsyn, nbins)                             
                            normalizing_vector_elegans[(categ,interval[0],interval[1])]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))
                        if (new_protein_list[count][1] == 0):       
                            acorr_syn_no_elegans[(categ,interval[0],interval[1])]+=autocorrelation(arr_uni_syn, nbins)
                            acorr_nonsyn_no_elegans[(categ,interval[0],interval[1])]+=autocorrelation(arr_uni_nonsyn, nbins)                             
                            normalizing_vector_no_elegans[(categ,interval[0],interval[1])]+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))  
                        
                                        
                                                                        
                    if gene_properties_no_singleton[categ]>interval[0] and gene_properties_no_singleton[categ]<=interval[1]:
                        acorr_syn_no_singleton[(categ,interval[0],interval[1])]+=autocorrelation(arr_uni_syn_no_singleton, nbins)
                        acorr_nonsyn_no_singleton[(categ,interval[0],interval[1])]+=autocorrelation(arr_uni_nonsyn_no_singleton, nbins)                             
                        normalizing_vector_no_singleton[(categ,interval[0],interval[1])]+=np.maximum(len(arr_uni_nonsyn_no_singleton)-np.arange(0,nbins), np.zeros(nbins))     
        
           
        
        
    else:
        print "Well, this was unexpected"
        break   

for key in acorr_syn:                                                   ##Normalizing the acorr list   
    acorr_syn[key]/=normalizing_vector[key]
    acorr_nonsyn[key]/=normalizing_vector[key]    
    acorr_syn_no_singleton[key]/=normalizing_vector_no_singleton[key]
    acorr_nonsyn_no_singleton[key]/=normalizing_vector_no_singleton[key]
    acorr_syn_elegans[key]/=normalizing_vector_elegans[key]
    acorr_nonsyn_elegans[key]/=normalizing_vector_elegans[key]
    acorr_syn_no_elegans[key]/=normalizing_vector_no_elegans[key]
    acorr_nonsyn_no_elegans[key]/=normalizing_vector_no_elegans[key]

acorr_syn_total/=normalizing_vector_total
acorr_nonsyn_total/=normalizing_vector_total    
acorr_syn_no_singleton_total/=normalizing_vector_no_singleton_total
acorr_nonsyn_no_singleton_total/=normalizing_vector_no_singleton_total
acorr_syn_elegans_total/=normalizing_vector_elegans_total
acorr_nonsyn_elegans_total/=normalizing_vector_elegans_total
acorr_syn_no_elegans_total/=normalizing_vector_no_elegans_total
acorr_nonsyn_no_elegans_total/=normalizing_vector_no_elegans_total


with open(str(datetime.date.today())+'_polymorphism_correlation.pickle', 'w') as f:
   pickle.dump((acorr_syn, acorr_nonsyn, normalizing_vector, acorr_syn_no_singleton, acorr_nonsyn_no_singleton, normalizing_vector_no_singleton, acorr_syn_elegans, acorr_nonsyn_elegans, normalizing_vector_elegans, acorr_syn_no_elegans, acorr_nonsyn_no_elegans, normalizing_vector_no_elegans,
               acorr_syn_total, acorr_nonsyn_total, normalizing_vector_total, acorr_syn_no_singleton_total, acorr_nonsyn_no_singleton_total, normalizing_vector_no_singleton_total, acorr_syn_elegans_total, acorr_nonsyn_elegans_total, normalizing_vector_elegans_total, acorr_syn_no_elegans_total, acorr_nonsyn_no_elegans_total, normalizing_vector_no_elegans_total), f)     
                        

#for key in acorr_syn:                                                        ## plotting   
#    plt.figure()
#    #plt.plot(acorr_syn[key][1:]/np.mean(acorr_syn[key][-10:]), label = 'syn'+str(key)+' asym: '+ str(round(np.mean(acorr_syn[key][-10:]), 4)))
#    #plt.plot(acorr_nonsyn[key][1:]/np.mean(acorr_nonsyn[key][-10:]), label = 'nonsyn'+str(key)+' asym: '+ str(round(np.mean(acorr_nonsyn[key][-10:]), 4)))
#    #
#    #plt.plot(acorr_syn_no_singleton[key][1:]/np.mean(acorr_syn_no_singleton[key][-10:]), label = 'syn_no_singleton'+str(key)+' asym: '+ str(round(np.mean(acorr_syn_no_singleton[key][-10:]), 4)))
#    #plt.plot(acorr_nonsyn_no_singleton[key][1:]/np.mean(acorr_nonsyn_no_singleton[key][-10:]), label = 'nonsyn_no_singleton'+str(key)+' asym: '+ str(round(np.mean(acorr_nonsyn_no_singleton[key][-10:]), 4)))
#    #
#    plt.plot(acorr_syn_elegans[key][1:]/np.mean(acorr_syn_elegans[key][-10:]), label = 'syn_elegans'+str(key)+' asym: '+ str(round(np.mean(acorr_syn_elegans[key][-10:]), 4)))
#    plt.plot(acorr_nonsyn_elegans[key][1:]/np.mean(acorr_nonsyn_elegans[key][-10:]), label = 'nonsyn_elegans'+str(key)+' asym: '+ str(round(np.mean(acorr_nonsyn_elegans[key][-10:]), 4)))
#    
#    plt.plot(acorr_syn_no_elegans[key][1:]/np.mean(acorr_syn_no_elegans[key][-10:]), label = 'syn_no_elegans'+str(key)+' asym: '+ str(round(np.mean(acorr_syn_no_elegans[key][-10:]), 4)))
#    plt.plot(acorr_nonsyn_no_elegans[key][1:]/np.mean(acorr_nonsyn_no_elegans[key][-10:]), label = 'nonsyn_no_elegans'+str(key)+' asym: '+ str(round(np.mean(acorr_nonsyn_no_elegans[key][-10:]), 4)))
#    plt.legend()
#    plt.show()
