#### makes the tree and ancestral_sequences(pickle) files ####

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


nbins = 150

acorr_syn = np.zeros(nbins)
acorr_nonsyn = np.zeros(nbins)
normalizing_vector = np.zeros(nbins)

mutations = {}

transcript_list = glob.glob('data_contig/*trans*/*fa')

protein_list = glob.glob('data_contig/*protein*/*fa')


for count in range(25019,28666): #range(len(transcript_list)):
    
    print transcript_list[count]
    
    gene_id = transcript_list[count].split('/')[-1].split('.fa')[0]
        
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
    