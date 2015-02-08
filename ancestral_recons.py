import subprocess
import glob
import os

list_transcript_alignment = glob.glob ('test_transcript/*fa')
list_transcript_tree = glob.glob('test_transcript_tree/*')

list_protein_alignment =glob.glob('data_contig/*pro*/*fa')
list_protein_tree = glob.glob('tree_contig_protein/*/*')

path_transcript = '/Users/vatshank/test_transcript/'
path_transcript_tree = '/Users/vatshank/test_transcript_tree/'
 


for count in range(len(list_transcript_alignment)):
    
    if (list_transcript_alignment[count].split('/')[1].split('.fa')[0] == list_transcript_tree[count].split('/')[1]) :
    
        print list_transcript_alignment[count]
        
        name = list_transcript_alignment[count].split('/')[1].split('.fa')[0]
        
        try:
            os.mkdir('ancestral_seq_reconstruction/'+name)
        except:
            print 'Cannot make directory'
        
        os.chdir('ancestral_seq_reconstruction/'+name)
        
        subprocess.call(['fastml','-s',path_transcript+list_transcript_alignment[count].split('/')[1]]) #,'-t',path_transcript_tree+list_transcript_tree[count].split('/')[1]])
        
        os.chdir('../../')
        
    else:
        print 'Well, this was unexpected'
        break