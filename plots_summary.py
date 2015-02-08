import numpy as np
import matplotlib.pyplot as plt
import pickle

with open('2013-06-21_polymorphism_correlation.pickle','r') as f:
    acorr_syn, acorr_nonsyn, normalizing_vector, acorr_syn_no_singleton, acorr_nonsyn_no_singleton, normalizing_vector_no_singleton, acorr_syn_elegans, acorr_nonsyn_elegans, normalizing_vector_elegans, acorr_syn_no_elegans, acorr_nonsyn_no_elegans, normalizing_vector_no_elegans, acorr_syn_total, acorr_nonsyn_total, normalizing_vector_total, acorr_syn_no_singleton_total, acorr_nonsyn_no_singleton_total, normalizing_vector_no_singleton_total, acorr_syn_elegans_total, acorr_nonsyn_elegans_total, normalizing_vector_elegans_total, acorr_syn_no_elegans_total, acorr_nonsyn_no_elegans_total, normalizing_vector_no_elegans_total = pickle.load(f)


nbins = 150

plt.figure()
plt.plot(acorr_syn_elegans_total[1:]/np.mean(acorr_syn_elegans_total[-10:]), label = 'syn_elegans'+' asym: '+ str(round(np.mean(acorr_syn_elegans_total[-10:]), 4)))
plt.plot(acorr_nonsyn_elegans_total[1:]/np.mean(acorr_nonsyn_elegans_total[-10:]), label = 'nonsyn_elegans'+' asym: '+ str(round(np.mean(acorr_nonsyn_elegans_total[-10:]), 4)))
plt.plot(acorr_syn_no_elegans_total[1:]/np.mean(acorr_syn_no_elegans_total[-10:]), label = 'syn_non_elegans'+' asym: '+ str(round(np.mean(acorr_syn_no_elegans_total[-10:]), 4)))
plt.plot(acorr_nonsyn_no_elegans_total[1:]/np.mean(acorr_nonsyn_no_elegans_total[-10:]), label = 'nonsyn_non_elegans'+' asym: '+ str(round(np.mean(acorr_nonsyn_no_elegans_total[-10:]), 4)))
plt.axhline(y=1,linestyle='--',color='k')
plt.xlabel('Separation(codons)',fontsize=15)
plt.ylabel('Autocorrelation',fontsize=15)
plt.legend()
plt.show()
plt.savefig('plots/summary/comparison_elegans.png')

plt.figure()
plt.plot(acorr_syn_total[1:]/np.mean(acorr_syn_total[-10:]), label = 'syn_all'+' asym: '+ str(round(np.mean(acorr_syn_total[-10:]), 4)))
plt.plot(acorr_nonsyn_total[1:]/np.mean(acorr_nonsyn_total[-10:]), label = 'nonsyn_all'+' asym: '+ str(round(np.mean(acorr_nonsyn_total[-10:]), 4)))
plt.plot(acorr_syn_no_singleton_total[1:]/np.mean(acorr_syn_no_singleton_total[-10:]), label = 'syn_no_singletons'+' asym: '+ str(round(np.mean(acorr_syn_no_singleton_total[-10:]), 4)))
plt.plot(acorr_nonsyn_no_singleton_total[1:]/np.mean(acorr_nonsyn_no_singleton_total[-10:]), label = 'nonsyn_no_singletons'+' asym: '+ str(round(np.mean(acorr_nonsyn_no_singleton_total[-10:]), 4)))
plt.axhline(y=1,linestyle='--',color='k')
plt.xlabel('Separation(codons)',fontsize=15)
plt.ylabel('Autocorrelation',fontsize=15)
plt.legend()
plt.show()
plt.savefig('plots/summary/comparison_singletons.png')