import glob
import pickle
from collections import defaultdict
import numpy as np


contig_list = glob.glob('annotations/*pickle')

coverage_values = {} #defaultdict(list)

coverage_by_worm = {}##

for contig in contig_list[0:50]:
    contig_name = contig.split('/')[1].split('_')[0]
    print contig_name    
    coverage_by_worm[contig_name] = {}##
    
    coverage_info ={}


    
    with open('coverage/'+contig_name+'.pickle') as f:
        coverage_info = pickle.load(f)
    
    for worm in coverage_info:
        worm_name = worm.split('/')[1].split('_RLA')[0]
        coverage_by_worm[contig_name][worm_name] = {}##
        for gene_id in coverage_info[worm]:
            coverage_values[gene_id] = []
            coverage_by_worm[contig_name][worm_name][gene_id] = []##    
            
                    
    for worm in coverage_info:
        worm_name = worm.split('/')[1].split('_RLA')[0]
        #print worm
        if len(coverage_info[worm])!=0:            
            for gene_id in coverage_info[worm]:
                #print gene_id
                coverage_by_worm[contig_name][worm_name][gene_id].append(np.mean([X[5] for X in coverage_info[worm][gene_id]]))##
                coverage_values[gene_id].extend([X[5] for X in coverage_info[worm][gene_id]])
                #for i in range(len(coverage_info[worm][gene_id])):
                #    coverage_values[gene_id].append(coverage_info[worm][gene_id][i][5]) 
                    
                    
    for gene_id in coverage_values:
        coverage_values[gene_id] = np.mean(coverage_values[gene_id])
        
with open('mean_coverage_values.pickle','w') as f:
    pickle.dump(coverage_values,f)
    
with open('mean_coverage_by_worm.pickle','w') as f:
    pickle.dump(coverage_by_worm,f)

    