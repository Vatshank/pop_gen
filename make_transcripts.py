# -*- coding: utf-8 -*-
import os
import pickle
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import glob

### splice ###
def splice(mRNA, exons, start_pos):
    spliced_seq = ''
    for exon in exons:
        spliced_seq+=mRNA[(exon[0]-start_pos):(exon[1]-start_pos+1)]

    return spliced_seq


contig='Contig1'
with open(contig+'_annotation.pickle', 'r') as f_anno:
    annotations_list = pickle.load(f_anno)

worm_list = glob.glob('ori_worms/*fa')

outdir = contig+'_mRNA'
try:
    os.mkdir(outdir)
except:
    print "Cannot make directory"

for (gene_id, anno) in annotations_list[:3]:
    all_mRNAs = {}
    all_spliced_transcripts = {}
    all_proteins = {}
    with open(gene_id+'.fa','w') as f1:
        my_sequences = []
        for worm_file_name in worm_list:
            worm = worm_file_name.split('/')[-1][:-11]
            for seq_record in SeqIO.parse(worm_file_name, "fasta"):
                print seq_record.id
                if seq_record.id==contig:
                    mRNA_seq = str(seq_record.seq[anno['range'][0]:anno['range'][1]+1])
                    spliced_transcript = splice(mRNA_seq, anno['exons'], start_pos =anno['range'][0]) 
                    all_mRNAs[worm] = mRNA_seq
                    all_spliced_transcripts[worm] = Seq(spliced_transcript, IUPAC.unambiguous_dna)
                    all_proteins[worm] = all_spliced_transcripts[worm].translate()
                    break

    with open(outdir+'/'+gene_id+'.fa', 'w') as gene_tmp:        
        for worm, seq in all_mRNAs.iteritems():
            gene_tmp.write('>'+worm+'\n')
            gene_tmp.write(str(seq)+'\n')

                
                
### get annotation ###


