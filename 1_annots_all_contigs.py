import pickle

import os

outdir = 'annotations'

try:
    os.mkdir(outdir)                       ##Making directories (one for mRNA and the other for proteins) for the Contig
except:
    print "Cannot make directory"


annotations = {}

# loop over annotations file, record all coding regions
with open('pacificus_complete_phase.gff3','r') as f_anno:
    line  = f_anno.readline()
   
    entries = line.strip().split()
    
    while line!='':
        if entries[1]=='Coding_transcript':  # found coding sequence
            gene_id = entries[-1].split('=')[1]  
            contig = entries[0]
       	    if contig not in annotations:
                annotations[contig] = {}
    
            annotations[contig][gene_id] = {}
            annotations[contig][gene_id]['type']=entries[2]
            annotations[contig][gene_id]['range']=[int(entries[3]), int(entries[4])]
            annotations[contig][gene_id]['strand']=entries[6]
            annotations[contig][gene_id]['exons']=[]

            # Continue over all exons of the coding sequence
            line  = f_anno.readline()
         
            entries = line.strip().split()
            while line!='' and gene_id == entries[-1].split('=')[1]: #while the gene_id does not change
                annotations[contig][gene_id]['exons'].append([int(entries[3]), int(entries[4])])               
                line  = f_anno.readline()              
                entries = line.strip().split()
                            
        else:
            line  = f_anno.readline()
            
            entries = line.strip().split()
            
for contig in annotations:
    annotations_list = annotations[contig].items()
    annotations_list.sort(key = lambda x:x[1]['range'][0])
    with open(outdir+'/'+contig+'_annotation.pickle', 'w') as f:
        pickle.dump(annotations_list, f)     
    
