import pickle

contig = 'Contig1'
annotations = {}
line_count=0
maxlines=10000000
with open('pacificus_complete_phase.gff3','r') as f_anno:
    line  = f_anno.readline()
    line_count+=1
    entries = line.strip().split()
    while line!='' and line_count<maxlines:
        if entries[0]==contig and entries[1]=='Coding_transcript':
            gene_id = entries[-1].split('=')[1]
            annotations[gene_id] = {}
            annotations[gene_id]['type']=entries[2]
            annotations[gene_id]['range']=[int(entries[3]), int(entries[4])]
            annotations[gene_id]['strand']=entries[6]
            annotations[gene_id]['exons']=[]

            line  = f_anno.readline()
            line_count+=1
            entries = line.strip().split()
            while gene_id == entries[-1].split('=')[1]:
                annotations[gene_id]['exons'].append([int(entries[3]), int(entries[4])])               
                line  = f_anno.readline()
                line_count+=1
                entries = line.strip().split()
        else:
            line  = f_anno.readline()
            line_count+=1
            entries = line.strip().split()
                    
annotations_list = annotations.items()
annotations_list.sort(key = lambda x:x[1]['range'][0])

with open(contig+'_annotation.pickle', 'w') as f:
    pickle.dump(annotations_list, f)

