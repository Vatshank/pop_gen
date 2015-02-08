
import glob
import random
import numpy as np
from cogent.phylo.maximum_likelihood import ML
from cogent.evolve.models import F81,H04G
from cogent import LoadTree, LoadSeqs
from cogent.phylo import nj
import pickle
import datetime
from cogent.core.genetic_code import DEFAULT as standard_code
import matplotlib.pyplot as plt


## calculates the charge correlation given a delta charge array ##
def charge_corr(arr_delta_charge,dmax):
    corr_rein= np.zeros(dmax)
    corr_comp= np.zeros(dmax)

    for d in range(1,min(dmax,arr_delta_charge.shape[0])):
        corr_rein[d]+=np.sum((arr_delta_charge[d:]*arr_delta_charge[:-d])>0)
        corr_comp[d]+=np.sum((arr_delta_charge[d:]*arr_delta_charge[:-d])<0)
        
    return corr_rein, corr_comp


def scramble(some_dict):

    list_muts = []
    for edge in B['edges']:
        if len(B['edges'][edge])>0:
            for sub in B['edges'][edge]:
                list_muts.append(sub)
                
    random.shuffle(list_muts)
    
    for edge in B['edges']:
        if len(B['edges'][edge])>0:
            for count,sub in enumerate(B['edges'][edge]):
                B['edges'][edge][count] = list_muts[0]
                list_muts.pop(0)
    
    if len(list_muts!=0):
        print 'you screwed up, again'
    return

## calculates the crosscorrelation given both uni_syn and uni_nonsyn arrays ##
def crosscorrelation(arr_syn,arr_nonsyn,dmax):
    crosscorr = np.zeros(2*dmax)
    for d in range(1,min(dmax,arr_syn.shape[0])):
        crosscorr[d]+=np.sum(arr_syn[d:] * arr_nonsyn[:-d])
        crosscorr[d+dmax]+=np.sum(arr_nonsyn[d:] * arr_syn[:-d])
     
    return crosscorr


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
    
    corr_charge_reinforced_scram = np.zeros(nbins)
    corr_charge_compensated_scram = np.zeros(nbins)    
    crosscorr_scram = np.zeros(2*nbins)
    crosscorr_edge_scram = np.zeros(2*nbins)
    acorr_edge_syn_scram = np.zeros(nbins)
    acorr_edge_nonsyn_scram = np.zeros(nbins)
    acorr_syn_scram = np.zeros(nbins)
    acorr_nonsyn_scram = np.zeros(nbins)
    normalizing_vector_scram = np.zeros(nbins)
    
    corr_charge_reinforced = np.zeros(nbins)
    corr_charge_compensated = np.zeros(nbins)    
    crosscorr = np.zeros(2*nbins)
    crosscorr_edge = np.zeros(2*nbins)
    acorr_edge_syn = np.zeros(nbins)
    acorr_edge_nonsyn = np.zeros(nbins)
    acorr_syn = np.zeros(nbins)
    acorr_nonsyn = np.zeros(nbins)
        
    arr_uni_syn_scram = np.zeros(length_arr)
    arr_uni_nonsyn_scram = np.zeros(length_arr)
    arr_uni_syn = np.zeros(length_arr)
    arr_uni_nonsyn = np.zeros(length_arr)
    for edge in dict_gene_id['edges']:
        #print edge
        arr_delta_charge_scram = np.zeros(length_arr)
        arr_edge_uni_syn_scram = np.zeros(length_arr)
        arr_edge_uni_nonsyn_scram = np.zeros(length_arr)
        arr_delta_charge = np.zeros(length_arr)
        arr_edge_uni_syn = np.zeros(length_arr)
        arr_edge_uni_nonsyn = np.zeros(length_arr)
        for subs in dict_gene_id['edges'][edge]:
            arr_delta_charge[int(subs[4])] = int(subs[5])   
                        
            #print subs
            if str(subs[2])!=str(subs[3]):                    ## Checking for non-synonymous substitutions
                arr_edge_uni_nonsyn[int(subs[4])] = 1
                arr_uni_nonsyn[int(subs[4])] = 1
            else:
                arr_edge_uni_syn[int(subs[4])] = 1            ## Synonymous substitutions
                arr_uni_syn[int(subs[4])] = 1
        
        new_reinforced, new_compensated = charge_corr(arr_delta_charge, nbins)

        corr_charge_reinforced+=new_reinforced
        corr_charge_compensated+=new_compensated      
        
        crosscorr_edge+=crosscorrelation(arr_edge_uni_syn, arr_edge_uni_nonsyn, nbins)
        acorr_edge_syn+=autocorrelation(arr_edge_uni_syn,nbins)
        acorr_edge_nonsyn+=autocorrelation(arr_edge_uni_nonsyn,nbins)
    
    crosscorr+=crosscorrelation(arr_uni_syn, arr_uni_nonsyn, nbins)         
    acorr_syn+=autocorrelation(arr_uni_syn,nbins)
    acorr_nonsyn+=autocorrelation(arr_uni_nonsyn,nbins)
    
    #print np.sum(acorr_edge_syn), np.sum(acorr_edge_nonsyn)
    #print np.sum(acorr_syn), np.sum(acorr_nonsyn)
    normalizing_vector+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))

    return acorr_syn,acorr_nonsyn,acorr_edge_syn,acorr_edge_nonsyn, corr_charge_reinforced, corr_charge_compensated, crosscorr, crosscorr_edge, normalizing_vector


charge_dict = {'R':1,'K':1,'D':-1,'E':-1,'G' : 0,'P' : 0,'A' : 0,'V' : 0,'L' : 0,'I' : 0,'M' : 0,'C' : 0,'F' : 0,'Y' : 0,'W': 0,'H' : 0,'Q' : 0,'N' : 0,'S' : 0,'T' : 0, '*': 0}

nbins = 150

crosscorr_scram = np.zeros(2*nbins)
crosscorr_edge_scram = np.zeros(2*nbins)
corr_charge_reinforced_scram = np.zeros(nbins)
corr_charge_compensated_scram = np.zeros(nbins)
acorr_syn_scram = np.zeros(nbins)
acorr_nonsyn_scram = np.zeros(nbins)
acorr_edge_syn_scram = np.zeros(nbins)
acorr_edge_nonsyn_scram = np.zeros(nbins)

crosscorr = np.zeros(2*nbins)
crosscorr_edge = np.zeros(2*nbins)
corr_charge_reinforced = np.zeros(nbins)
corr_charge_compensated = np.zeros(nbins)
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

with open('temp_mutations.pickle','r') as f:
    mutations = pickle.load(f)
            
for count in range(len(tree_list[:14000])):
    
    print tree_list[count]
    
    if (new_protein_dict['data_contig/'+tree_list[count].split('/')[1].split('-')[0]+'_protein/'+tree_list[count].split('/')[1].split('_')[0]+'.fa']==1): ## filtering elegans  proteins
    
        if tree_list[count].split('/')[-1].split('_')[0]==ancestor_list[count].split('/')[-1].split('_')[0]:
        
            gene_id = tree_list[count].split('/')[-1].split('_')[0]            
            
            new_acorr_syn,new_acorr_nonsyn,new_edge_acorr_syn,new_edge_acorr_nonsyn, new_corr_charge_reinforced, new_corr_charge_compensated, new_crosscorr, new_crosscorr_edge, new_normalizing_vector=make_arr_uni(mutations[gene_id]) ## Making uni_array and calculating the autocorrelation
            
            crosscorr+=new_crosscorr
            crosscorr_edge+=new_crosscorr_edge
            corr_charge_reinforced+=new_corr_charge_reinforced
            corr_charge_compensated+=new_corr_charge_compensated
            acorr_syn+=new_acorr_syn
            acorr_nonsyn+=new_acorr_nonsyn
            acorr_edge_syn+=new_edge_acorr_syn
            acorr_edge_nonsyn+=new_edge_acorr_nonsyn
            normalizing_vector+=new_normalizing_vector
        
        else:
            print 'Well, this was unexpected'
            break
        
crosscorr[:nbins]/=normalizing_vector
crosscorr[nbins:2*nbins]/=normalizing_vector
crosscorr_edge[:nbins]/=normalizing_vector
crosscorr_edge[nbins:2*nbins]/=normalizing_vector

corr_charge_reinforced/=normalizing_vector
corr_charge_compensated/=normalizing_vector

acorr_edge_syn/=normalizing_vector
acorr_edge_nonsyn/=normalizing_vector
acorr_syn/=normalizing_vector
acorr_nonsyn/=normalizing_vector


## plotting ##
plt.figure()
plt.plot(acorr_syn[1:]/np.mean(acorr_syn[-10:]), label = 'syn'+' asym: '+ str(round(np.mean(acorr_syn[-10:]), 4)))
plt.plot(acorr_nonsyn[1:]/np.mean(acorr_nonsyn[-10:]), label = 'nonsyn'+' asym: '+ str(round(np.mean(acorr_nonsyn[-10:]), 4)))
plt.plot(acorr_edge_syn[1:]/np.mean(acorr_edge_syn[-10:]), label = 'syn'+' asym_edge: '+ str(round(np.mean(acorr_edge_syn[-10:]), 4)))
plt.plot(acorr_edge_nonsyn[1:]/np.mean(acorr_edge_nonsyn[-10:]), label = 'nonsyn'+' asym_edge: '+ str(round(np.mean(acorr_edge_nonsyn[-10:]), 4)))
plt.plot(corr_charge_reinforced[1:]/np.mean(corr_charge_reinforced[-10:]), label = 'charge_reinforced'+' asym: '+ str(round(np.mean(corr_charge_reinforced[-10:]), 4)))
plt.plot(corr_charge_compensated[1:]/np.mean(corr_charge_compensated[-10:]), label = 'charge_compensated'+' asym: '+ str(round(np.mean(corr_charge_compensated[-10:]), 4)))
plt.plot(crosscorr[1:]/np.mean(crosscorr[-10+nbins:nbins]), label = 'syn'+' asym_edge: '+ str(round(np.mean(acorr_edge_syn[-10:]), 4)))
plt.legend()
plt.legend()
plt.show()
   