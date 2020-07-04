import numpy as np 
import matplotlib.pyplot as plt
from LLESolver import LLESolvers

Length = [2585/54.1,1672/104.1]
Species = ["PB","PS"]
Solver = "GTP"

Temp = np.linspace(250,370,50)
rSAFT = LLESolvers(Solver,"PCSAFT",Species,Length,k=-0.000173665)
rUNIFAC = LLESolvers(Solver,"UNIFAC",Species,Length)
rFH     = LLESolvers(Solver,"FH",Species,Length)

x1_SAFT = []
x2_SAFT = []
x1_UNIFAC = []
x2_UNIFAC = []
x1_FH = []
x2_FH = []
for i in range(len(Temp)):
    if i==0:
        X_SAFT=rSAFT.LLE(Temp[i])
        X_UNIFAC=rUNIFAC.LLE(Temp[i])
        X_FH    = rFH.LLE(Temp[i])
    else:
        X_SAFT=rSAFT.LLE(Temp[i],1e5,[x1_SAFT[i-1],x2_SAFT[i-1]])
        X_UNIFAC=rUNIFAC.LLE(Temp[i],1e5,[x1_UNIFAC[i-1],x2_UNIFAC[i-1]])
        X_FH    = rFH.LLE(Temp[i],1e5,[x1_FH[i-1],x2_FH[i-1]])
    x1_UNIFAC.append(X_UNIFAC[0])
    x2_UNIFAC.append(X_UNIFAC[1])
    x1_SAFT.append(X_SAFT[0])
    x2_SAFT.append(X_SAFT[1])
    x1_FH.append(X_FH[0])
    x2_FH.append(X_FH[1])
    # print(X_UNIFAC)
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
ax.plot(x1_FH,Temp,color='k', ls='dotted')
ax.plot(x2_FH,Temp,color='k', ls='dotted')
ax.plot([0.7,0.6,0.4,0.2],[341.95,351.15,360.45,355.75],'rx')
ax.set_xlabel(r'\textit{x}')
ax.set_ylabel(r'\textit{T}')
ax.set_xlim(0,1)
plt.savefig('LLE.png')