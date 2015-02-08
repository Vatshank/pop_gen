import subprocess
import glob
import os


all_transcript_list = glob.glob ('data_contig/*tran*/*fa')

all_pro_list = glob.glob ('data_contig/*pro*/*fa')




entry = all_transcript_list[0].split('/')[1].split('_')[0]  

outdir = entry
    
try:
    os.mkdir('tree_contig_transcript/'+outdir)
except:
    print "Cannot make directory"        
    
for transcript in all_transcript_list:
    
    print transcript
    
    output = subprocess.check_output(['./FastTree','-nt',transcript])

    outdir = transcript.split('/')[1].split('_')[0]
    
    if entry != transcript.split('/')[1].split('_')[0]:
         try:
             os.mkdir('tree_contig_transcript/'+outdir)
         except:
             print "Cannot make directory"       
                         
                        
    with open('tree_contig_transcript/'+ outdir + '/' + transcript.split('/')[2].split('.fa')[0],'w') as f:
             f.write(str(output))
    
    index_no = all_transcript_list.index(transcript) 
         
    entry = transcript.split('/')[1].split('_')[0]





entry = all_pro_list[0].split('/')[1].split('_')[0]

outdir = entry

try:
    os.mkdir('tree_contig_protein/'+outdir)
except:
    print "Cannot make directory"        
    

for protein in all_pro_list:
    
    print protein
    
    output = subprocess.check_output(['./FastTree',protein])

    outdir = protein.split('/')[1].split('_')[0]
    
    if entry != protein.split('/')[1].split('_')[0]:
         try:
             os.mkdir('tree_contig_protein/'+outdir)
         except:
             print "Cannot make directory"       
                         
                        
    with open('tree_contig_protein/'+ outdir + '/' + protein.split('/')[2].split('.fa')[0],'w') as f:
             f.write(str(output))
    
    index_no = all_pro_list.index(protein) 
         
    entry = protein.split('/')[1].split('_')[0]                      