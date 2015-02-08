#### takes as input the tree and ancestral_sequences(pickle) files ####
#### record the mutations as a dictionary and outputs the autocorrelation over the ancestors ####


import glob
import numpy as np
from cogent.phylo.maximum_likelihood import ML
from cogent.evolve.models import F81,H04G
from cogent import LoadTree, LoadSeqs
from cogent.phylo import nj
import pickle
import datetime
from cogent.core.genetic_code import DEFAULT as standard_code
import matplotlib.pyplot as plt

## calculates the autocorrelation given a uni_array ##
def autocorrelation(arr, dmax):                                
    acorr = np.zeros(dmax)
    for d in range(1,min(dmax,arr.shape[0])):
        acorr[d]+=np.sum(arr[d:] * arr[:-d])
     
    return acorr 
    
## makes the uni_arrays like before and calls the autocorrelation function ## 
def make_arr_uni(dict_gene_id):                
    length_arr = dict_gene_id['length']
    
    nbins = 150
    
    acorr_edge_syn = np.zeros(nbins)
    acorr_edge_nonsyn = np.zeros(nbins)
    acorr_syn = np.zeros(nbins)
    acorr_nonsyn = np.zeros(nbins)
    normalizing_vector = np.zeros(nbins)
    
    arr_uni_syn = np.zeros(length_arr)
    arr_uni_nonsyn = np.zeros(length_arr)
    for edge in dict_gene_id['edges']:
        #print edge
        arr_edge_uni_syn = np.zeros(length_arr)
        arr_edge_uni_nonsyn = np.zeros(length_arr)
        for subs in dict_gene_id['edges'][edge]:
            #print subs
            if str(subs[2])!=str(subs[3]):                    ## Checking for non-synonymous substitutions
                arr_edge_uni_nonsyn[int(subs[4])] = 1
                arr_uni_nonsyn[int(subs[4])] = 1
            else:
                arr_edge_uni_syn[int(subs[4])] = 1            ## Synonymous substitutions
                arr_uni_syn[int(subs[4])] = 1
        
        acorr_edge_syn+=autocorrelation(arr_edge_uni_syn,nbins)
        acorr_edge_nonsyn+=autocorrelation(arr_edge_uni_nonsyn,nbins)
    acorr_syn+=autocorrelation(arr_uni_syn,nbins)
    acorr_nonsyn+=autocorrelation(arr_uni_nonsyn,nbins)
    print np.sum(acorr_edge_syn), np.sum(acorr_edge_nonsyn)
    print np.sum(acorr_syn), np.sum(acorr_nonsyn)
    normalizing_vector+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))

    return acorr_syn,acorr_nonsyn,acorr_edge_syn,acorr_edge_nonsyn,normalizing_vector


## recursive function that makes a dictionary of the mutations in the ancestors ##
def rec_mutations_tree(T,ancestors):
    sequence_root = ancestors.getSeq(T.Name)
    
    codons_root = []
    for codon in range(len(sequence_root)/3):                        ##transcript to codons
        codons_root.append(str(sequence_root[(3*codon):(codon+1)*3]))
    
    codons_root = np.array(codons_root)
        
    length_aln = len(codons_root)
    mutations[gene_id]['length'] = length_aln    
    
    
    
    for child in T.Children:
        if (child.istip()!=True):
            #print child.Name            
            sequence_child = ancestors.getSeq(child.Name)
            codons_child = []
            for codon in range(len(sequence_child)/3):
                codons_child.append(str(sequence_child[(3*codon):(codon+1)*3]))
            codons_child = np.array(codons_child)
                
            #print np.where(sequence_root!=sequence_child)            
            #mutations = np.where(sequence_root!=sequence_child)
            mutations[gene_id]['edges'][child.Name]= []
            
            for where in np.where(codons_root!=codons_child)[0]:                    ##looking for differences in the child and root sequences 
                mutations[gene_id]['edges'][child.Name].append([str(codons_root[where]),str(codons_child[where]), standard_code[str(codons_root[where])], standard_code[str(codons_child[where])],str(where)])
                
            rec_mutations_tree(child,ancestors)
         

nbins = 150

acorr_syn = np.zeros(nbins)
acorr_nonsyn = np.zeros(nbins)
acorr_edge_syn = np.zeros(nbins)
acorr_edge_nonsyn = np.zeros(nbins)
normalizing_vector = np.zeros(nbins)

mutations = {}

tree_list = glob.glob('njtrees_transcript/*tree')

ancestor_list = glob.glob('njtrees_transcript/*sequences')

with open('new_protein_dict.pickle','r') as f:
    new_protein_dict = pickle.load(f)
    

for count in range(len(tree_list[:10000])):
    
    print tree_list[count]
    
    if (new_protein_dict['data_contig/'+tree_list[count].split('/')[1].split('-')[0]+'_protein/'+tree_list[count].split('/')[1].split('_')[0]+'.fa']==1): ## filtering elegans  proteins
    
        if tree_list[count].split('/')[-1].split('_')[0]==ancestor_list[count].split('/')[-1].split('_')[0]:
        
            gene_id = tree_list[count].split('/')[-1].split('_')[0]
            
            njtree = LoadTree(tree_list[count])
            
            with open(ancestor_list[count],'r') as f:
                ancestors = pickle.load(f)
            
            mutations[gene_id] = {}            
            mutations[gene_id]['edges'] = {}
                
            rec_mutations_tree(njtree,ancestors)   ## Creating the muatations dictionary
            
            new_acorr_syn,new_acorr_nonsyn,new_edge_acorr_syn,new_edge_acorr_nonsyn, new_normalizing_vector=make_arr_uni(mutations[gene_id]) ## Making uni_array and calculating the autocorrelation
            
            acorr_syn+=new_acorr_syn
            acorr_nonsyn+=new_acorr_nonsyn
            acorr_edge_syn+=new_edge_acorr_syn
            acorr_edge_nonsyn+=new_edge_acorr_nonsyn
            normalizing_vector+=new_normalizing_vector
        
        else:
            print 'Well, this was unexpected'
            break
        

acorr_syn/=normalizing_vector
acorr_nonsyn/=normalizing_vector


## plotting ##
plt.figure()
plt.plot(acorr_syn[1:]/np.mean(acorr_syn[-10:]), label = 'syn'+' asym: '+ str(round(np.mean(acorr_syn[-10:]), 4)))
plt.plot(acorr_nonsyn[1:]/np.mean(acorr_nonsyn[-10:]), label = 'nonsyn'+' asym: '+ str(round(np.mean(acorr_nonsyn[-10:]), 4)))
plt.plot(acorr_edge_syn[1:]/np.mean(acorr_edge_syn[-10:]), label = 'syn'+' asym_edge: '+ str(round(np.mean(acorr_edge_syn[-10:]), 4)))
plt.plot(acorr_edge_nonsyn[1:]/np.mean(acorr_edge_nonsyn[-10:]), label = 'nonsyn'+' asym_edge: '+ str(round(np.mean(acorr_edge_nonsyn[-10:]), 4)))
plt.legend()
plt.show()

plt.figure()
plt.plot(acorr_edge_nonsyn[1:]/np.mean(acorr_edge_nonsyn[-10:])-acorr_edge_syn[1:]/np.mean(acorr_edge_syn[-10:]), label = 'syn'+' asym_edge: '+ str(round(np.mean(acorr_edge_syn[-10:]), 4)))
