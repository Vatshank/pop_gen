import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
import glob

def autocorrelation(arr, dmax):                                ## defining the autocorrelation function
    acorr = np.zeros(dmax)
    for d in range(1,min(dmax,arr.shape[0])):
        acorr[d]+=np.sum(arr[d:] * arr[:-d])
     
    return acorr               
            

    


indir_trans = 'Contig1_transcripts'                            ## the directories from the data is read
indir_pro = 'Contig1_pro'

transcript_list = glob.glob(indir_trans+'/Contig1*fa')        ## making a list of all files in the respective directories
protein_list = glob.glob(indir_pro+'/Contig1*fa')


codons = {}

acorr_transcript = {}
acorr = np.zeros(50)
normalizing_vector = np.zeros_like(acorr)
mean_poly_density = []


for transcript in transcript_list:                    ## looping over files in the transcript list
    print transcript
    alignment = AlignIO.read(transcript,"fasta")            ## reading the alignment
    
    codons=[]
    
    for codon in range(alignment.get_alignment_length()/3):                        ## getting codons from the transcript
        
        codons.append(map(str,[X.seq for X in alignment[:,(3*codon):(codon+1)*3]]))
        
    matrix = np.array(codons).T                                ## using transpose of the matrix

    uni= []
    
    for n in range(matrix.shape[1]):                            ## getting uni, which has the number of different types codons at each site
        uni.append(np.unique(matrix[:,n]).shape[0])	        
    
    arr_uni = np.array(uni)>1                                ## same as->  arr_uni[arr_uni==1] = 0 #arr_uni[arr_uni>1] = 1  
                                                        
    mean_poly_density.append(np.mean(arr_uni))
    

    acorr+=autocorrelation(arr_uni, len(acorr))                ## adding the autocorrelation values for the transripts
    normalizing_vector+=len(arr_uni)-np.arange(0,len(acorr))
    acorr_transcript[transcript] = autocorrelation(arr_uni, len(acorr))

acorr/=normalizing_vector
plt.figure()
plt.plot(acorr)
plt.plot(np.arange(len(acorr)), np.mean(acorr[-10:])*np.ones_like(acorr))
