import matplotlib.pyplot as plt
import numpy as np

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
plt.rc("axes", titlesize="small")
plt.rc("legend", fontsize="small")

fig = plt.figure(figsize=(4, 3))
ax = fig.add_subplot(1, 1, 1)
ax.set_ylabel("Temperature")
ax.set_xlabel(r"$x_{1}$")
ax.set_ylim([0,2])
ax.axes.yaxis.set_ticklabels([])

#Generate plots

x = np.linspace(0,1,100)

binodal = -10*(x-0.9)*(x-0.10)
spinodal = -25.6*(x-0.75)*(x-0.25)

ax.plot(x,binodal,"k-",label="Binodal")
ax.plot(x,spinodal,"k--",label="Spinodal")
ax.plot(0.5,8/5, marker="x", color="r")

# ax.annotate(r"$T_{C}$", (0.53,1.63))
ax.annotate("", xytext=(0.19,0.6), xy=(0.3, 0.6), arrowprops=dict(arrowstyle="<->"), ha="center")
ax.annotate(r"Metastable" "\n" "     Zone", xy=(0.25,0.6), xytext=(0.0,1.4), arrowprops=dict(arrowstyle="-"))

ax.legend(frameon=False)

plt.tight_layout()
# plt.show()
plt.savefig("spinodal-schematic.png", dpi=300)