import numpy as np
import random
a = [random.randint(0,1) for i in range(15)]
arr_a = np.array(a)
print arr_a
acorr = np.zeros(100)
for d in range(1,arr_a.shape[0]):
    if d<=100:
        acorr[d]+=np.mean(arr_a[d:] * arr_a[:-d])
       
    else:
        break
    
