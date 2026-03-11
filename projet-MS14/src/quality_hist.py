import numpy as np
import matplotlib.pyplot as plt

with open('quality.sol') as f:
    lines = f.readlines()

qual=[]
for k in range(9,len(lines)-2):
    qual.append(float(lines[k]))
n = len(qual)
Quality = np.array(qual) 


plt.hist(Quality,20, log=True)
plt.xlabel("qualité")
plt.ylabel("quantité")
plt.title("histogramme de la qualité du maillage")
plt.show()