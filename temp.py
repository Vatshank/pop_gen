# -*- coding: utf-8 -*-
import os
import pickle
from Bio import SeqIO                        ##Importing the biopython modules to help translate 
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import glob                                  ##Importing glob to read our multiple files in the directory
import time

def splice(mRNA, exons, start_pos):          ##Splice Function
    spliced_seq = ''
    
    for exon in exons:
        spliced_seq+=mRNA[(exon[0]-start_pos):(exon[1]-start_pos+1)]       ##Building the spliced sequence and adjusting the positions

    return spliced_seq

def reverse_complement(seq):                                ##reverse complement function for '-' strands
    r_seq = list(seq[::-1])                                 ##Converting str to list

    for count, n in enumerate(r_seq):
        if n == 'A':
            r_seq[count] = 'T'
        elif n == 'T':
            r_seq[count] = 'A'
        elif n == 'G':
            r_seq[count] = 'C'
        else:
            r_seq[count] = 'G'

        
    rc_seq = ''.join(r_seq)                                 ##Converting list back to string
    return rc_seq

annotations_list = {}

all_mRNAs = {}                   	                          ##Creating dictionaries
all_spliced_transcripts = {}
all_proteins = {}

contig_list = glob.glob('annotations/*pickle')


for contig in contig_list[:10] + contig_list[500:510]:

    worm_list = glob.glob('ori_worms/*fa')                                    ##Listing all the worms in the directory 'ori_worms'

#    with open(contig, 'r') as f_anno:                    ##Opening the pickled file
#        annotations_list[contig] = pickle.load(f_anno)                                ##Read as a list
#
#    outdir_mRNA = 'data_contig1/'+contig.split('/')[1].split('_')[0]+'_mRNA'                                              ##Naming directory for a particular Contig
#    outdir_protein = 'data_contig1/'+contig.split('/')[1].split('_')[0]+'_protein'
#    outdir_trans = 'data_contig1/'+contig.split('/')[1].split('_')[0]+'_transcripts'
#
#    try:
#        os.mkdir(outdir_mRNA), os.mkdir(outdir_protein), os.mkdir(outdir_trans)                       ##Making directories (one for mRNA and the other for proteins) for the Contig
#    except:
#        print "Cannot make directory"
#
#    all_mRNAs[contig] = {}                   	                          ##Creating dictionaries
#    all_spliced_transcripts[contig] = {}
#    all_proteins[contig] = {}
#    
#    for (gene_id, anno) in annotations_list[contig]:
#        all_mRNAs[contig][gene_id] = {}                   	                          ##Creating dictionaries
#        all_spliced_transcripts[contig][gene_id] = {}
#        all_proteins[contig][gene_id] = {}
#    
    
    for worm_file_name in worm_list[:1]:
        worm = worm_file_name.split('/')[-1][:-11]                          ##Modifying the name a bit
        print worm,
        
        t0 = time.time()
        
        for seq_record in SeqIO.parse(worm_file_name, "fasta"):                                           ##SeqIO.parse() returns an iterator giving SeqRecord objects
        
            if seq_record.id==contig.split('/')[1].split('_')[0]:
                
                t1 = time.time()
                
                print seq_record.id, str(t1 - t0)
                
                #for (gene_id, anno) in annotations_list[contig]:                                               
                #    mRNA_seq = str(seq_record.seq[(anno['range'][0]-1):anno['range'][1]])                     ##Getting the mRNA sequence
                #
                #    if anno['strand']=='-':
                #        spliced_transcript = reverse_complement(splice(mRNA_seq, anno['exons'][::-1], start_pos =anno['range'][0]))     ##Calling Splice function for '-'strand
                #    else:                  
                #        spliced_transcript = splice(mRNA_seq, anno['exons'], start_pos =anno['range'][0])     ##Calling Splice function for '+'strand
                #
                #    all_mRNAs[contig][gene_id][worm] = mRNA_seq                                                                ##Storing all the mRNAs 
                #    all_spliced_transcripts[contig][gene_id][worm] = Seq(spliced_transcript, IUPAC.ambiguous_dna)            ##Converting spliced transcripts ro Seq objects and storing them
                #    all_proteins[contig][gene_id][worm] = all_spliced_transcripts[contig][gene_id][worm].translate()                        ##Translating the spliced sequences and storing the proteins
                break

    
for (gene_id, anno) in annotations_list[contig]:    
    with open(outdir_mRNA+'/'+gene_id+'.fa', 'w') as gene_tmp:           ##Creating a file for a particular mRNA in the Contig
    
        for worm, seq in all_mRNAs[contig][gene_id].iteritems():                          ##.iteritems() iterates over the dictionary
            gene_tmp.write('>'+worm+'\n')
            gene_tmp.write(str(seq)+'\n')

    with open(outdir_protein+'/'+gene_id+'.fa', 'w') as pro_tmp:        ##Creating a file for a particular protein in the Contig
        
        for worm, seq in all_proteins[contig][gene_id].iteritems():                      ##.iteritems() iterates over the dictionary
            pro_tmp.write('>'+worm+'\n')
            pro_tmp.write(str(seq)+'\n')                

    with open(outdir_trans+'/'+gene_id+'.fa', 'w') as trans_tmp:        ##Creating a file for a particular protein in the Contig
    
        for worm, seq in all_spliced_transcripts[contig][gene_id].iteritems():                      ##.iteritems() iterates over the dictionary
            trans_tmp.write('>'+worm+'\n')
            trans_tmp.write(str(seq)+'\n')                


