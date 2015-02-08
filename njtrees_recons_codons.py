import glob
import numpy as np
from cogent.phylo.maximum_likelihood import ML
from cogent.evolve.models import F81,H04G
from cogent import LoadTree, LoadSeqs
from cogent.phylo import nj
import pickle
import datetime
from cogent.core.genetic_code import DEFAULT as standard_code

def calc_hamming_distance(aln):
    dists= {}
    aln_dict={sr.Name:np.asarray(list(str(sr))) for sr in aln.iterSeqs()}
    keys = aln_dict.keys()
    for ki,key1 in enumerate(keys):
        for key2 in keys[:ki]:
            dists[(key1,key2)]=np.mean(aln_dict[key1]!= aln_dict[key2])
            dists[(key2,key1)] = dists[(key1,key2)]
    return dists

#def autocorrelation(arr, dmax):                                ## defining the autocorrelation function
#    acorr = np.zeros(dmax)
#    for d in range(1,min(dmax,arr.shape[0])):
#        acorr[d]+=np.sum(arr[d:] * arr[:-d])
#     
#    return acorr 
#    
#def make_arr_uni(dict_gene_id):
#    length_arr = dict_gene_id['length']
#    
#    nbins = 150
#    
#    acorr_syn = np.zeros(nbins)
#    acorr_nonsyn = np.zeros(nbins)
#    normalizing_vector = np.zeros(nbins)
#    
#    for edge in dict_gene_id['edges']:
#        print edge
#        arr_uni_syn = np.zeros(length_arr)
#        arr_uni_nonsyn = np.zeros(length_arr)
#        for subs in dict_gene_id['edges'][edge]:
#            #print subs
#            if str(subs[2])!=str(subs[3]):
#                arr_uni_nonsyn[int(subs[4])] = 1
#            else:
#                arr_uni_syn[int(subs[4])] = 1
#        
#        acorr_syn+=autocorrelation(arr_uni_syn,nbins)
#        acorr_nonsyn+=autocorrelation(arr_uni_nonsyn,nbins)
#        normalizing_vector+=np.maximum(len(arr_uni_nonsyn)-np.arange(0,nbins), np.zeros(nbins))
#        
#    return acorr_syn,acorr_nonsyn,normalizing_vector
#
#
#
#def rec_mutations_tree(T,ancestors):
#    sequence_root = ancestors.getSeq(T.Name)
#    
#    codons_root = []
#    for codon in range(len(sequence_root)/3):
#        codons_root.append(str(sequence_root[(3*codon):(codon+1)*3]))
#    
#    codons_root = np.array(codons_root)
#        
#    length_aln = len(codons_root)
#    mutations[gene_id]['length'] = length_aln    
#    
#    
#    
#    for child in T.Children:
#        if (child.istip()!=True):
#            print child.Name
#            ##convert these into codons!
#            
#            sequence_child = ancestors.getSeq(child.Name)
#            codons_child = []
#            for codon in range(len(sequence_child)/3):
#                codons_child.append(str(sequence_child[(3*codon):(codon+1)*3]))
#            codons_child = np.array(codons_child)
#                
#            #print np.where(sequence_root!=sequence_child)            
#            #mutations = np.where(sequence_root!=sequence_child)
#            mutations[gene_id]['edges'][child.Name]= []
#            
#            for where in np.where(codons_root!=codons_child)[0]:
#                mutations[gene_id]['edges'][child.Name].append([str(codons_root[where]),str(codons_child[where]), standard_code[str(codons_root[where])], standard_code[str(codons_child[where])],str(where)])
#                
#            rec_mutations_tree(child,ancestors)
    


#with open('mutations_ancestors.pickle','r') as f:
#    mutations = pickle.load(f)        
                        
nbins = 150

acorr_syn = np.zeros(nbins)
acorr_nonsyn = np.zeros(nbins)
normalizing_vector = np.zeros(nbins)

mutations = {}

transcript_list = glob.glob('data_contig/*trans*/*fa')

protein_list = glob.glob('data_contig/*protein*/*fa')


for count in range(19598,28666): #range(len(transcript_list)):
    
    print transcript_list[count]
    
    gene_id = transcript_list[count].split('/')[-1].split('.fa')[0]
    
    #mutations[gene_id] = {}
    #
    #mutations[gene_id]['edges'] = {}
    
    aln_trans = LoadSeqs(transcript_list[count])
    aln_trans =  aln_trans.takeSeqs(aln_trans.Names[:105])##[:-3] ##alignment without the outgroup
    
    distances = calc_hamming_distance(aln_trans)
    njtree = nj.nj(distances)
   
    lf = F81().makeLikelihoodFunction(njtree, digits=3, space=2)
    #lf = H04G().makeLikelihoodFunction(njtree, digits=2, space=3)
    lf.setAlignment(aln_trans)
    lf.optimise(show_progress=False, local=True)
    
    ancestors = lf.likelyAncestralSeqs()
    
    njtree.writeToFile('njtrees_transcript/'+gene_id+'_tree')
    with open('njtrees_transcript/'+gene_id+'_ancestral_sequences','w') as f:
        pickle.dump(ancestors,f)
    
    #rec_mutations_tree(njtree,ancestors)
    #
    #new_acorr_syn,new_acorr_nonsyn,new_normalizing_vector=make_arr_uni(mutations[gene_id])
    #
    #acorr_syn+=new_acorr_syn
    #acorr_nonsyn+=new_acorr_nonsyn
    #normalizing_vector+=new_normalizing_vector
    #
#with open(str(datetime.today())+'_mutations_ancestors.pickle','w') as f:
#    pickle.dump(mutations,f)    
