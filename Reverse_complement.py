def reverse_complement(seq):
    
    r_seq = list(seq[::-1])
    print seq
    
    print ''.join(r_seq)
    
    count= 0
    for n in r_seq:
        if n == 'A':
            r_seq[count] = 'T'
        elif n == 'T':
            r_seq[count] = 'A'
        elif n == 'G':
            r_seq[count] = 'C'
        else:
            r_seq[count] = 'G'
        
        count  +=1    
    
    rc_seq = ''.join(r_seq)
    
    print rc_seq
    
    return rc_seq

