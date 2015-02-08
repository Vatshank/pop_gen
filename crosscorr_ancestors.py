import pickle
import glob
import numpy as np
import matplotlib.pyplot as plt

tree_list = glob.glob('njtrees_transcript/*tree')

ancestor_list = glob.glob('njtrees_transcript/*sequences')

with open('new_protein_dict.pickle','r') as f:
    new_protein_dict = pickle.load(f)

with open('temp_mutations.pickle','r') as f:
    mutations = pickle.load(f)

def crosscorrelation(arr_syn,arr_nonsyn,dmax):
    crosscorr = np.zeros(2*dmax)
    for d in range(1,min(dmax,arr_syn.shape[0])):
        crosscorr[d]+=np.sum(arr_syn[d:] * arr_nonsyn[:-d])
        crosscorr[d+dmax]+=np.sum(arr_nonsyn[d:] * arr_syn[:-d])
     
    return crosscorr

def make_arr_uni(dict_gene_id):                
    length_arr = dict_gene_id['length']
    
    nbins = 150
    
    crosscorr = np.zeros(2*nbins)
    crosscorr_edge = np.zeros(2*nbins)
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
        
        crosscorr_edge+=crosscorrelation(arr_edge_uni_syn, arr_edge_uni_nonsyn, nbins)
    
    crosscorr+=crosscorrelation(arr_uni_syn, arr_uni_nonsyn, nbins)    
    normalizing_vector+=np.maximum(length_arr-np.arange(0,nbins), np.zeros(nbins))

    return crosscorr, crosscorr_edge, normalizing_vector
    
nbins = 150
crosscorr = np.zeros(2*nbins)
crosscorr_edge = np.zeros(2*nbins)
normalizing_vector = np.zeros(nbins)           
                                    
for count in range(len(tree_list[:13000])) :
    
    print tree_list[count]
    
    if (new_protein_dict['data_contig/'+tree_list[count].split('/')[1].split('-')[0]+'_protein/'+tree_list[count].split('/')[1].split('_')[0]+'.fa']==1): ## filtering elegans  proteins
    
        if tree_list[count].split('/')[-1].split('_')[0]==ancestor_list[count].split('/')[-1].split('_')[0]:
        
            gene_id = tree_list[count].split('/')[-1].split('_')[0]
            
            new_crosscorr, new_crosscorr_edge, new_normalizing_vector=make_arr_uni(mutations[gene_id])
            
            crosscorr_edge+=new_crosscorr_edge
            crosscorr+=new_crosscorr
            normalizing_vector+=new_normalizing_vector
            
        else:
            print 'Well, this was unexpected'
            break
    
crosscorr[:nbins]/=normalizing_vector
crosscorr[nbins:2*nbins]/=normalizing_vector
crosscorr_edge[:nbins]/=normalizing_vector
crosscorr_edge[nbins:2*nbins]/=normalizing_vector


plt.figure()
plt.plot(crosscorr[1:nbins]/np.mean(crosscorr[-10+nbins:nbins]), label = 'crosscorr'+' asym_edge: '+ str(round(np.mean(crosscorr[-10:]), 4)))
plt.plot(crosscorr_edge[1:nbins]/np.mean(crosscorr_edge[-10+nbins:nbins]), label = 'crosscorr_edge'+' asym: '+ str(round(np.mean(crosscorr_edge[-10:]), 4)))
plt.axhline(y=1,linestyle='--',color='k')
plt.xlabel('Separation(codons)',fontsize=15)
plt.ylabel('Autocorrelation',fontsize=15)
plt.legend()
plt.show()
plt.savefig('plots/summary/trees_crosscorr.png')


            

            
            