from Thermo import ThermoMix
from LLESolver import LLESolvers
import matplotlib.pyplot as plt
import numpy as np

xComp = np.linspace(0.000001,0.999999,100)
Length = [100,100]
Species = ["PB","PS"]
Temp = 370

gMix_FH = []
gMix_UNIFAC = []
gMix_SAFT = []
gMix_SAFT_CH = []
FH  = ThermoMix("FH",Species,Length,[Temp])
UNIFAC  = ThermoMix("UNIFAC",Species,Length,[Temp])
SAFT  = ThermoMix("PCSAFT",Species,Length,[Temp])
SAFT_CH = ThermoMix("PCSAFT",Species,Length,[Temp],CH="On")
for i in range(len(xComp)):
    
    gMix_FH.append(FH.GibbsFreeMixing(xComp[i]))
    
    
    gMix_UNIFAC.append(UNIFAC.GibbsFreeMixing(xComp[i]))

    gMix_SAFT.append(SAFT.GibbsFreeMixing(xComp[i]))

    gMix_SAFT_CH.append(SAFT_CH.GibbsFreeMixing(xComp[i]))
plt.rc('font', family='serif',size=16)
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
fig = plt.figure(figsize=(7, 5.25))
ax = fig.add_subplot(1, 1, 1)
ax.plot(xComp,gMix_FH,color='k', ls='solid',label=r"Flory-Huggins")
ax.plot(xComp,gMix_UNIFAC,color='0.50', ls='dashed',label=r"UNIFAC")
ax.plot(xComp,gMix_SAFT,color='k', ls='solid', linewidth=0.5,label=r"\textit{s}PC-SAFT")
ax.plot(xComp,gMix_SAFT_CH,color='0.75', ls='solid', linewidth=0.5,label=r"\textit{s}PC-SAFT")

ax.legend(loc="lower left")
ax.set_xlabel(r'\textit{x}')
ax.set_ylabel(r'\textit{g}')
ax.set_xlim(0,1)
plt.savefig('gMix.png')