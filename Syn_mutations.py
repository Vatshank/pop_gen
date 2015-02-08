import numpy as np
import random

a = [random.randint(0,1) for i in range(15)]
b = [random.randint(0,1) for i in range(15)]

pos_syn = []
pos_nonsyn = []

for n in range(len(a)):
    if(a[n]==1) and (b[n]==0):
        pos_syn.append(n)
    elif(a[n]==1) and (b[n]==1):
        pos_nonsyn.append(n)
    print "$$" 
    print pos_syn
    print pos_nonsyn
       
        
        