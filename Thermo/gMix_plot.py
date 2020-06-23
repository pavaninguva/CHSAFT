from Thermo import ThermoMix
from LLESolver import LLESolvers
import matplotlib.pyplot as plt
import numpy as np

xComp = np.linspace(0.000001,0.999999,100)
Length = [1100/54.1,1672/104.1]
Species = ["PB","PS"]
Temp = 270

gMix_FH = []
gMix_UNIFAC = []
gMix_SAFT = []
for i in range(len(xComp)):
    FH  = ThermoMix("FH",Species,Length,[xComp[i],1-xComp[i]])
    gMix_FH.append(FH.GibbsFreeMixing())
    
    UNIFAC  = ThermoMix("UNIFAC",Species,Length,[xComp[i],1-xComp[i]])
    gMix_UNIFAC.append(UNIFAC.GibbsFreeMixing())

    SAFT  = ThermoMix("PCSAFT",Species,Length,[xComp[i]/2585/(xComp[i]*(1/2585-1/1672)+1/1672),1-(xComp[i]/2585/(xComp[i]*(1/2585-1/1672)+1/1672))],[Temp])
    gMix_SAFT.append(SAFT.GibbsFreeMixing())
plt.rc('font', family='serif',size=16)
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
fig = plt.figure(figsize=(7, 5.25))
ax = fig.add_subplot(1, 1, 1)
ax.plot(xComp,gMix_FH,color='k', ls='solid',label=r"Flory-Huggins")
ax.plot(xComp,gMix_UNIFAC,color='0.50', ls='dashed',label=r"UNIFAC")
ax.plot(xComp,gMix_SAFT,color='k', ls='solid', linewidth=0.5,label=r"\textit{s}PC-SAFT")
# ax.legend(loc="lower left")
ax.set_xlabel(r'\textit{x}')
ax.set_ylabel(r'\textit{g}')
ax.set_xlim(0,1)
plt.savefig('gMix.png')