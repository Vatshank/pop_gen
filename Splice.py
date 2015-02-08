exons = [[5,6],[8,10],[11,15]]

def splice(mRNA,exons,mRNA_start) :
    
    with open('first_exons','w') as f_exons:
        
    
        for n in exons:
            data = mRNA[n[0]-mRNA_start:n[1]-mRNA_start+1] # adjusting the position of the mRNA
            f_exons.write(data)
            print data
        

splice('AUCGUAGUGACUAUAUCGCGGCGGGCUACGUAGCUGACAUCGCGUACGUAGC',exons,5)                