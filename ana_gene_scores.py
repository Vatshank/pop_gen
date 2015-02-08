import pickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats


with open('corr_vectors.pickle', 'r') as infile:
    gene_acorrs = pickle.load(infile)
    
    
genes = np.array(gene_acorrs.keys())

gene_scores = []

first_wd = 20
last_wd = 50
nbins=100

for gene in genes:
    tmp = []
    gs = gene_acorrs[gene]
    for ii in range(3):
        tmp.append(np.sum(gs[ii][:first_wd])/np.sum(gs[3][:first_wd]))
        tmp.append(np.sum(gs[ii][last_wd:nbins])/np.sum(gs[3][last_wd:nbins]))
    tmp.extend(gs[4:])
    gene_scores.append(tmp)

gene_scores =np.array(gene_scores)

for ii in range(3):
    y,x = np.histogram(gene_scores[:,2*ii]-gene_scores[:,2*ii+1], bins=np.linspace(-0.01,0.02,101))
    plt.plot((x[1:]+x[:-1])*0.5, y, label=str(ii))



for perc in [0, 50,80,90,95,99,100]:
    ii=0
    threshold = stats.scoreatpercentile(gene_scores[:,2*ii]-gene_scores[:,2*ii+1],perc)

    ii=0
    high_syn = gene_scores[:,2*ii]-gene_scores[:,2*ii+1]>threshold
    print np.sum(high_syn)
    ii=1
    high_nonsyn = gene_scores[:,2*ii]-gene_scores[:,2*ii+1]>threshold
    
    acorr_nonsyn = np.zeros(nbins)
    acorr_syn = np.zeros(nbins)
    crosscorr = np.zeros(2*nbins)
    normalizing_vector = np.zeros(nbins)
    
    index_set  = np.where(~high_syn)[0]
    for gene in genes[index_set]:
        gs = gene_acorrs[gene]
        acorr_nonsyn+=gs[1]
        acorr_syn+=gs[0]
        crosscorr+=gs[2]
        normalizing_vector+=gs[3]
        
    acorr_nonsyn/=normalizing_vector
    acorr_syn/=normalizing_vector
    crosscorr[:nbins]/=normalizing_vector
    crosscorr[nbins:2*nbins]/=normalizing_vector
    
    plt.figure()
    plt.plot(acorr_nonsyn/np.mean(acorr_nonsyn[nbins-10:]))
    plt.plot(acorr_syn/np.mean(acorr_syn[nbins-10:]))
    plt.plot(crosscorr[:nbins]/np.mean(crosscorr[nbins-10:nbins]))
    
    
