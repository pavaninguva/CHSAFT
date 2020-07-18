from Thermo import ThermoMix
from LLESolver import LLESolvers
import matplotlib.pyplot as plt
import numpy as np

xComp = np.linspace(0.000001,0.999999,200)
Length = [50,50]
Species = ["PB","PS"]
Temp = 298

gMix_FH = []
gMix_UNIFAC = []
gMix_SAFT = []
gMix_SAFT_CH = []
FH  = ThermoMix("FH",Species,Length,[Temp])
UNIFAC  = ThermoMix("UNIFAC",Species,Length,[Temp])
SAFT  = ThermoMix("PCSAFT",Species,Length,[Temp])
# SAFT_CH = ThermoMix("PCSAFT",Species,Length,[Temp],CH="On",k=0.000720986)
for i in range(len(xComp)):
    
    gMix_FH.append(FH.GibbsFreeMixing(xComp[i]))
    
    
    gMix_UNIFAC.append(UNIFAC.GibbsFreeMixing(xComp[i]))

    gMix_SAFT.append(SAFT.GibbsFreeMixing(xComp[i]))

    # gMix_SAFT_CH.append(SAFT_CH.GibbsFreeMixing(xComp[i]))

x_eq_FH = LLESolvers("GTP","FH",Species,Length).LLE(Temp)
g_eq_FH = [FH.GibbsFreeMixing(x_eq_FH[0]),FH.GibbsFreeMixing(x_eq_FH[1])]
x_eq_UNIFAC = LLESolvers("GTP","UNIFAC",Species,Length).LLE(Temp)
g_eq_UNIFAC = [UNIFAC.GibbsFreeMixing(x_eq_UNIFAC[0]),UNIFAC.GibbsFreeMixing(x_eq_UNIFAC[1])]
x_eq_SAFT = LLESolvers("GTP","PCSAFT",Species,Length).LLE(Temp)
g_eq_SAFT = [SAFT.GibbsFreeMixing(x_eq_SAFT[0]),SAFT.GibbsFreeMixing(x_eq_SAFT[1])]
plt.rc('font', family='serif',size=16)
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
fig = plt.figure(figsize=(7, 5.25))
ax = fig.add_subplot(1, 1, 1)
for direction in ["left","bottom"]:
    ax.spines[direction].set_position('zero')
    # ax.spines[direction].set_smart_bounds(True)
for direction in ["top"]:
    ax.spines[direction].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.plot(xComp,gMix_SAFT,color='k', ls='solid',label=r"PC-SAFT")
ax.plot(xComp,gMix_UNIFAC,color='0.50', ls='dashed',label=r"UNIFAC")
ax.plot(xComp,gMix_FH,color='k', ls='dotted',label=r"Flory-Huggins")
# ax.plot([0,1],[-0.099,-0.099],color='k', ls='solid', linewidth=0.55,label=r"\textit{s}PC-SAFT")
# ax.plot(x_eq_FH,g_eq_FH,'--k',linewidth=0.5,label=r"Tangent Plane")
# ax.plot(x_eq_FH,g_eq_FH,'Dk',fillstyle='none',label=r"Binodal point")
# ax.plot(x_eq_UNIFAC,g_eq_UNIFAC,'--k',linewidth=0.5,label=r"Tangent Plane")
# ax.plot(x_eq_UNIFAC,g_eq_UNIFAC,'Dk',fillstyle='none',label=r"Binodal point")
# ax.plot(x_eq_SAFT,g_eq_SAFT,'--k',linewidth=0.5,label=r"Tangent Plane")
# ax.plot(x_eq_SAFT,g_eq_SAFT,'Dk',fillstyle='none',label=r"Binodal point")

# ax.plot(xComp,gMix_SAFT_CH,color='0.75', ls='solid', linewidth=0.5,label=r"\textit{s}PC-SAFT")

ax.legend(loc="upper left",frameon=False)
ax.set_xlabel(r'$x_1$')
ax.set_ylabel(r'$\Delta G_{\mathrm{mix}}/(Nk_\mathrm{B}T)$')
ax.set_xlim(0,1)
# ax.set_ylim(-1.2,3.7)
plt.tight_layout()
plt.savefig('gMix.png')