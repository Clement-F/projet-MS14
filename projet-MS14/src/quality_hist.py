import numpy as np
import matplotlib.pyplot as plt

with open('quality.sol') as f:
    lines = f.readlines()

qual=[]
for k in range(9,len(lines)-2):
    qual.append(float(lines[k]))
n = len(qual)
Quality = np.array(qual) 

for k in range(0,n): 
    if(Quality[k]>10e3): 
        print(" element ",k," is of quality :",Quality[k],"\n")


plt.hist(Quality,20, log=True)
plt.xlabel("qualité")
plt.ylabel("quantité")
plt.title("histogramme de la qualité du maillage")
plt.show()