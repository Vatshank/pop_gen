def crosscorrelation(arr_syn.arr_nonsyn,dmax):
    crosscorr = np.zeros(dmax)
    for d in range(1,min(dmax,arr_syn.shape[0])):
        crosscorr[d]+=np.sum(arr_syn[d:] * arr_non[:-d])
     
    return crosscorr   