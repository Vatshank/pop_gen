import numpy as np
import pickle
import matplotlib.pyplot as plt

categories = { 'length': [ [0, 113], [113,188], [188,291], [291,431], [431,99999]],
                'GC': [[0,0.48], [0.48,0.51], [0.51,0.54], [0.54,0.58], [0.58,1]],
                'syn': [[0,0.05], [0.05,0.10],[0.10,0.15],[0.15,0.20],[0.20,1]],
                'nonsyn':[[0,0.05], [0.05,0.10],[0.10,0.15],[0.15,0.20],[0.20,1]]}


with open('2013-06-05_polymorphism_correlation.pickle','r') as f:
    acorr_syn, acorr_nonsyn, normalizing_vector, acorr_syn_no_singleton, acorr_nonsyn_no_singleton, normalizing_vector_no_singleton, acorr_syn_elegans, acorr_nonsyn_elegans, normalizing_vector_elegans, acorr_syn_no_elegans, acorr_nonsyn_no_elegans, normalizing_vector_no_elegans = pickle.load(f)
    
    

for key in acorr_syn:                                                        ## plotting   
    #plt.figure()
    #plt.plot(acorr_syn[key][1:]/np.mean(acorr_syn[key][-10:]), label = 'syn'+str(key)+' asym: '+ str(round(np.mean(acorr_syn[key][-10:]), 4)))
    #plt.plot(acorr_nonsyn[key][1:]/np.mean(acorr_nonsyn[key][-10:]), label = 'nonsyn'+str(key)+' asym: '+ str(round(np.mean(acorr_nonsyn[key][-10:]), 4)))    
    #plt.plot(acorr_syn_no_singleton[key][1:]/np.mean(acorr_syn_no_singleton[key][-10:]), label = 'syn_no_singleton'+str(key)+' asym: '+ str(round(np.mean(acorr_syn_no_singleton[key][-10:]), 4)))
    #plt.plot(acorr_nonsyn_no_singleton[key][1:]/np.mean(acorr_nonsyn_no_singleton[key][-10:]), label = 'nonsyn_no_singleton'+str(key)+' asym: '+ str(round(np.mean(acorr_nonsyn_no_singleton[key][-10:]), 4)))
    #plt.legend()
    #plt.savefig('plots/'+'all_'+str(key[0])+'-'+str(key[1])+'-'+str(key[2])+'.png')
    # 
    
    plt.figure()
    plt.plot(acorr_syn_elegans[key][1:]/np.mean(acorr_syn_elegans[key][-10:]), label = 'syn_elegans'+str(key)+' asym: '+ str(round(np.mean(acorr_syn_elegans[key][-10:]), 4)))
    plt.plot(acorr_nonsyn_elegans[key][1:]/np.mean(acorr_nonsyn_elegans[key][-10:]), label = 'nonsyn_elegans'+str(key)+' asym: '+ str(round(np.mean(acorr_nonsyn_elegans[key][-10:]), 4)))    
    #plt.plot(acorr_syn_no_elegans[key][1:]/np.mean(acorr_syn_no_elegans[key][-10:]), label = 'syn_no_elegans'+str(key)+' asym: '+ str(round(np.mean(acorr_syn_no_elegans[key][-10:]), 4)))
    #plt.plot(acorr_nonsyn_no_elegans[key][1:]/np.mean(acorr_nonsyn_no_elegans[key][-10:]), label = 'nonsyn_no_elegans'+str(key)+' asym: '+ str(round(np.mean(acorr_nonsyn_no_elegans[key][-10:]), 4)))
    plt.xlabel('Separation(codons)',fontsize=15)
    plt.ylabel('Autocorrelation',fontsize=15)
    plt.legend()
    plt.savefig('plots/summary/'+'elegans_'+str(key[0])+'-'+str(key[1])+'-'+str(key[2])+'.png')
       
    
    
    #GC_strata = []
    #
    #for k in acorr_syn:
    #    if k[0]=='GC': GC_strata.append(k)
    #
    #GC_strata.sort(key=lambda x:x[1])
    #
    #for x in GC_strata:
    #    plt.plot(acorr_syn[x]/np.mean(acorr_syn[x][100:]), label=x)

    
#acorr={}
#
#acorr['syn'] = acorr_syn
#acorr['syn_no_singleton'] = acorr_syn_no_singleton
#acorr['syn_elegans'] = acorr_syn_elegans 
#acorr['syn_no_elegans'] = acorr_syn_no_elegans
#
#acorr['nonsyn'] = acorr_nonsyn
#acorr['nonsyn_no_singleton'] = acorr_nonsyn_no_singleton
#acorr['nonsyn_elegans'] = acorr_nonsyn_elegans 
#acorr['nonsyn_no_elegans'] = acorr_nonsyn_no_elegans
#        
#                
#strata = {}
#
#for categ in categories:
#        strata[categ]= []
#        
#for key in strata:
#    for k in acorr_syn:
#        if k[0] == key:strata[key].append(k)
#        
#        
#for key in strata:
#    strata[key].sort(key =lambda x:x[1])
#    
#    for acorr_type in acorr:
#        plt.figure()
#        for x in strata[key]:
#            plt.plot(acorr[acorr_type][x][1:]/np.mean(acorr[acorr_type][x][100:]),label = acorr_type+str(x))
#        plt.legend()
#        plt.savefig('plots/strata/'+acorr_type+'-'+key+'.png')           
#            