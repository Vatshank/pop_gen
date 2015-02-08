import os
import pickle
import glob

contig_list = glob.glob('annotations/*pickle')

coverage_list = glob.glob('worms/*vcf')



for contig in contig_list[484:]:
    coverage_info ={}    
    contig_name = contig.split('_')[0].split('/')[-1]
    print contig_name
    
    #try:
    #    os.mkdir('coverage/'+contig_name)
    #except:
    #    print 'Cannot make directory'

      
    with open(contig,'r') as f_anno:
       annotations_list = pickle.load(f_anno)
       
    
    for worm in coverage_list:
        print worm
        coverage_info[worm] = {}
        
        
        with open(worm,'r') as f_coverage:
            
            line = f_coverage.readline()

            while line!= '':
                if line.split()[0] == contig_name: 
                    break
                line = f_coverage.readline()                                            #while True:
                                                            #    line = f_coverage.readline()
                                                            #                                            #for line in f_coverage:
                                                            #    
            if line=='':
                continue                                        #    if line.split()[0]==contig_name:
                                                            #        break
                                                            #    if not line: break
                       
            for gene_id,anno in annotations_list:
                tstart, tstop = anno['range']
                #line =  f_coverage.readline()
                entries = line.strip().split()            
                position = int(entries[1])
                coverage_info[worm][gene_id] = []
               
                while contig_name==entries[0] and position<tstart:
                    line =  f_coverage.readline()
                    entries = line.strip().split()            
                    position = int(entries[1])
                    
                if contig_name==entries[0]:
                    while contig_name==entries[0] and position<tstop:
                        #line =  f_coverage.readline()
                        #entries = line.strip().split()            
                        #position = int(entries[1])
                        #### action happens ####
                    #coverage_info[worm][gene_id].append((position, entries[0], entries[3], entries[4], float(entries[5])))
                        #coverage_info[worm][gene_id][position] = (entries[0], entries[3], entries[4], float(entries[5]))
                        #coverage_info[worm][gene_id][position] = {}
                        #coverage_info[worm][gene_id][position]['contig'] = entries[0]
               	        #coverage_info[worm][gene_id][position]['reference'] = entries[3]
                        #coverage_info[worm][gene_id][position]['substitution'] = entries[4]
                        #coverage_info[worm][gene_id][position]['quality'] = float(entries[5])
                        
                        if entries[7].split(';')[0] == 'INDEL':
                            coverage_info[worm][gene_id].append((position, entries[0], entries[3], entries[4], float(entries[5]), float(entries[7].split(';')[1].split('=')[1])))
                            
                            #coverage_info[worm][gene_id][position]['coverage'] = float(entries[7].split(';')[1].split('=')[1])
                            
                        else:
                            
                            coverage_info[worm][gene_id].append((position, entries[0], entries[3], entries[4], float(entries[5]), float(entries[7].split(';')[0].split('=')[1])))
                            
                            #coverage_info[worm][gene_id][position]['coverage'] = float(entries[7].split(';')[0].split('=')[1])
                        
                        line =  f_coverage.readline()
                        entries = line.strip().split()            
                        position = int(entries[1])
                        
                        
                else:
                    break
        
        worm_name=worm.split('/')[1].split('_')[0]
        #with open('coverage/'+contig_name+'/'+worm_name+'.pickle','w') as f:
        #    pickle.dump(coverage_info[worm],f)                                
    with open('coverage/'+contig_name+'.pickle','w') as f:
        pickle.dump(coverage_info,f)
              
                
                
        
        
        
        
        
        
        
        
        
        
       
       