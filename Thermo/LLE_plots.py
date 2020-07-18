import numpy as np 
import matplotlib.pyplot as plt
from LLESolver import LLESolvers

Length = [6350/101.1,1390/104.1]
Species = ["PMMA","PS"]
Solver = "GTP"

Temp = np.linspace(360,440,200)
rSAFT = LLESolvers(Solver,"PCSAFT",Species,Length)
rUNIFAC = LLESolvers(Solver,"UNIFAC",Species,Length)
rFH     = LLESolvers(Solver,"FH",Species,Length)

x1_SAFT = []
x2_SAFT = []
x1_UNIFAC = []
x2_UNIFAC = []
x1_FH = []
x2_FH = []
T_SAFT  = []
T_UNIFAC  = []
T_FH  = []
for i in range(len(Temp)):
    if i==0:
        X_SAFT=rSAFT.LLE(Temp[i])
        X_UNIFAC=rUNIFAC.LLE(Temp[i])
        X_FH    = rFH.LLE(Temp[i])
    else:
        X_SAFT=rSAFT.LLE(Temp[i],1e5,[X_SAFT[0],X_SAFT[1]])
        X_UNIFAC=rUNIFAC.LLE(Temp[i],1e5,[X_UNIFAC[0],X_UNIFAC[1]])
        X_FH    = rFH.LLE(Temp[i],1e5,[X_FH[0],X_FH[1]])
    if abs(X_SAFT[0]-X_SAFT[1])>1e-3:
        x1_SAFT.append(X_SAFT[0])
        x2_SAFT.append(X_SAFT[1])
        T_SAFT.append(Temp[i])
    if abs(X_UNIFAC[0]-X_UNIFAC[1])>1e-3:
        x1_UNIFAC.append(X_UNIFAC[0])
        x2_UNIFAC.append(X_UNIFAC[1])
        T_UNIFAC.append(Temp[i])
    if abs(X_FH[0]-X_FH[1])>1e-3:
        x1_FH.append(X_FH[0])
        x2_FH.append(X_FH[1])
        T_FH.append(Temp[i])
    # x1_FH.append(X_FH[0])
    # x2_FH.append(X_FH[1])
    # print(X_UNIFAC)
    print(i)
plt.rc('font', family='serif',size=16)
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
fig = plt.figure(figsize=(7, 5.25))
ax = fig.add_subplot(1, 1, 1)
x2_SAFT = x2_SAFT[::-1]
for x in x2_SAFT:
  x1_SAFT.append(x)

# x2_UNIFAC = x2_UNIFAC[::-1]
# for x in x2_UNIFAC:
#   x1_UNIFAC.append(x)

x2_FH = x2_FH[::-1]
for x in x2_FH:
  x1_FH.append(x)

T_SAFT_2 = T_SAFT[::-1]
T_SAFT.extend(T_SAFT_2)

# T_UNIFAC_2 = T_UNIFAC[::-1]
# T_UNIFAC.extend(T_UNIFAC_2)

T_FH_2 = T_FH[::-1]
T_FH.extend(T_FH_2)

ax.plot(x1_SAFT,T_SAFT,color='k', ls='solid',label=r"PC-SAFT")
ax.plot(x1_UNIFAC,T_UNIFAC,color='k', ls='dashed',label=r"UNIFAC")
ax.plot(x2_UNIFAC,T_UNIFAC,color='k', ls='dashed')
# ax.plot(x1_FH,T_FH,color='k', ls='dotted',label=r"Flory-Huggins")
ax.plot([0.112796209,0.112796209,0.112796209,0.112796209,0.112796209,0.112796209,0.112796209,0.222429907,0.222429907,0.222429907,0.222429907,0.329032258,0.329032258,0.432727273,0.432727273,0.533632287,0.533632287,0.533632287,0.631858407,0.631858407,0.631858407,0.727510917,0.727510917,0.727510917,0.727510917,0.820689655,0.820689655,0.820689655,0.820689655,0.820689655,0.820689655,0.820689655],[423.15,418.15,413.15,403.15,393.15,383.15,373.15,423.15,418.15,413.15,403.15,423.15,418.15,423.15,418.15,423.15,418.15,413.15,423.15,418.15,413.15,423.15,418.15,413.15,403.15,423.15,418.15,413.15,403.15,393.15,383.15,373.15],'sk',fillstyle='none',label=r"one-phase")
ax.plot([0.222429907,0.222429907,0.222429907,0.329032258,0.329032258,0.329032258,0.329032258,0.329032258,0.432727273,0.432727273,0.432727273,0.432727273,0.432727273,0.533632287,0.533632287,0.533632287,0.533632287,0.631858407,0.631858407,0.631858407,0.631858407,0.727510917,0.727510917,0.727510917],[393.15,383.15,373.15,413.15,403.15,393.15,383.15,373.15,413.15,403.15,393.15,383.15,373.15,403.15,393.15,383.15,373.15,403.15,393.15,383.15,373.15,393.15,383.15,373.15],'xk',fillstyle='none',label=r"two-phase")
ax.legend(loc="upper left",frameon=False,ncol=2)
ax.set_xlabel(r'weight fraction of PMMA')
ax.set_ylabel(r'Temperature / K')
ax.set_xlim(0,1)
ax.set_ylim(360,440)
plt.tight_layout()
plt.savefig('LLE.pdf')