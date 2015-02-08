import numpy as np
from Bio import AlignIO
import glob

def autocorrelation(arr, dmax):
    acorr = np.zeros(dmax)
    #acorr[0]+=np.mean(arr**2)
    for d in range(1,min(dmax,arr.shape[0])):
        acorr[d]+=np.sum(arr[d:] * arr[:-d])
     
    return acorr               
            

    


indir_trans = 'Contig1_transcripts'
indir_pro = 'Contig1_pro'

transcript_list = glob.glob(indir_trans+'/Contig1*fa')
protein_list = glob.glob(indir_pro+'/Contig1*fa')


codons = {}

acorr_transcript = {}
acorr = np.zeros(150)
mean_poly_density = []
for transcript in transcript_list[:500]:
    print transcript
    alignment = AlignIO.read(transcript,"fasta")

    matrix = np.array(alignment)

    uni= []
    
    for n in range(matrix.shape[1]):
        uni.append(np.unique(matrix[:,n]).shape[0]) 
    
    arr_uni = np.array(uni)>1   # same as->  arr_uni[arr_uni==1] = 0 #arr_uni[arr_uni>1] = 1  
    mean_poly_density.append(np.mean(arr_uni))
    

    acorr+=autocorrelation(arr_uni, len(acorr))
    acorr_transcript[transcript] = autocorrelation(arr_uni, len(acorr))
    

    
    