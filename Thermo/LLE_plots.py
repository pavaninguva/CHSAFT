import numpy as np 
import matplotlib.pyplot as plt
from LLESolver import LLESolvers

Length = [1100/54.1,1672/104.1]
Species = ["PB","PS"]
Method = "UNIFAC"
Solver = "GTP"

Temp = np.linspace(270,310,50)
rSAFT = LLESolvers(Solver,"PCSAFT",Species,Length)
rUNIFAC = LLESolvers(Solver,"UNIFAC",Species,Length)

x1_SAFT = []
x2_SAFT = []
x1_UNIFAC = []
x2_UNIFAC = []
for i in range(len(Temp)):
    if i==0:
        X_SAFT=rSAFT.LLE(Temp[i])
        X_UNIFAC=rUNIFAC.LLE(Temp[i])
    else:
        X_SAFT=rSAFT.LLE(Temp[i],1e5,[x1_SAFT[i-1],x2_SAFT[i-1]])
        X_UNIFAC=rUNIFAC.LLE(Temp[i],1e5,[x1_UNIFAC[i-1],x2_UNIFAC[i-1]])
        # X = r.LLE(Temp[i])
    x1_UNIFAC.append(X_UNIFAC[0])
    x2_UNIFAC.append(X_UNIFAC[1])
    x1_SAFT.append(X_SAFT[0])
    x2_SAFT.append(X_SAFT[1])
    print(X_UNIFAC,X_SAFT)
    print(i)
plt.rc('font', family='serif',size=16)
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
fig = plt.figure(figsize=(7, 5.25))
ax = fig.add_subplot(1, 1, 1)

ax.plot(x1_SAFT,Temp,color='k', ls='solid')
ax.plot(x2_SAFT,Temp,color='k', ls='solid')
ax.plot(x1_UNIFAC,Temp,color='k', ls='dashed')
ax.plot(x2_UNIFAC,Temp,color='k', ls='dashed')
ax.plot([0.7,0.5,0.3],[294.15,301.05,298.75],'rx')
ax.set_xlabel(r'\textit{x}')
ax.set_ylabel(r'\textit{T}')
ax.set_xlim(0,1)
plt.savefig('LLE.png')