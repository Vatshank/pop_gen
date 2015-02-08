import glob
import pickle
from collections import defaultdict
import numpy as np


contig_list = glob.glob('annotations/*pickle')

for contig in contig_list[2:3]:
    
    coverage_info ={}
    coverage_values = {}#defaultdict(list)
    contig_name = contig.split('/')[1].split('_')[0]
    
    with open('coverage/'+contig_name+'.pickle') as f:
        coverage_info = pickle.load(f)
    
    for worm in coverage_info:
        for gene_id in coverage_info[worm]:
            coverage_values[gene_id] = []    
            
                    
    for worm in coverage_info:
        print worm
        if len(coverage_info[worm])!=0:            
            for gene_id in coverage_info[worm]:
                print gene_id
                for i in range(len(coverage_info[worm][gene_id])):
                    coverage_values[gene_id].append(coverage_info[worm][gene_id][i][5]) 
                    
                    
    map(np.mean,coverage_values.itervalues())