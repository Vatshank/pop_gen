import glob
import pickle

protein_list = glob.glob('data_contig/*pro*/*fa')

with open('real_proteins.pickle','r') as f:
    elegans_list = pickle.load(f)


new_protein_dict = {}

for protein in protein_list:
    print protein
    
    switch= 0
    
    for count in range(len(elegans_list)):
        if (str(protein.split('/')[2].split('.fa')[0]) == str(elegans_list[count])):
            switch = 1
            break
         
    if switch==1:
        new_protein_dict[protein]=1
    else:
        new_protein_dict[protein]=0
         
with open('new_protein_dict.pickle','w') as f:
    pickle.dump(new_protein_dict,f)  