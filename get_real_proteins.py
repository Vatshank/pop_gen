import os
import pickle

real_proteins_list = []


with open('pacificus_complete_protein_vs_c_elegans_blast.txt','r') as f:
    line = f.readline()
    entries = line.strip().split()   
    
    gene = entries[0] 
    
    real_proteins_list.append(entries[0])    
    
    while line!= '':             
        
        if entries[0]!=gene:
            gene = entries[0]
            real_proteins_list.append(entries[0])
            line = f.readline()
            entries = line.strip().split()
            
            
        else:
            line = f.readline()
            entries = line.strip().split()

                  
with open('real_proteins.pickle','w') as f:
    pickle.dump(real_proteins_list,f)
    
                                                