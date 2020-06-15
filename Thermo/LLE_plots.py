import numpy as np 
import matplotlib.pyplot as plt
from LLESolver import LLESolvers

Length = [2585/54.1,1672/104.1]
Species = ["PB","PS"]
Method = "PCSAFT"
Solver = "GTP"

Temp = np.linspace(270,370,50)
r = LLESolvers(Solver,Method,Species,Length)

x1 = []
x2 = []
for i in range(len(Temp)):
    if i==0:
        X=r.LLE(Temp[i])
    else:
        X=r.LLE(Temp[i],1e5,[x1[i-1],x2[i-1]])
        # X = r.LLE(Temp[i])
    x1.append(X[0])
    x2.append(X[1])
    print(i)

print(x1,x2)
plt.rc('font', family='serif',size=16)
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
fig = plt.figure(figsize=(7, 5.25))
ax = fig.add_subplot(1, 1, 1)

ax.plot(x1,Temp,color='k', ls='solid')
ax.plot(x2,Temp,color='k', ls='solid')
ax.set_xlabel(r'\textit{x}')
ax.set_ylabel(r'\textit{T}')
ax.set_xlim(0,1)
plt.savefig('LLE.png')