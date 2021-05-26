import matplotlib.pyplot as plt  
import numpy as np

plt.rcParams["text.usetex"] = True

plt.rc('font', family='serif')
# plt.rc('xtick', labelsize='small')
# plt.rc('ytick', labelsize='small')
# plt.rc("axes", titlesize="medium")
# plt.rc("legend", fontsize="small")

fig = plt.figure(figsize=(4, 3))
ax = fig.add_subplot(1, 1, 1)
ax.set_ylabel(r"$\frac{\Lambda}{L_{0}  \rho_{m}RT}$")
ax.set_ylabel(r"$\Lambda / \mathrm{N m^{-1}}$")
ax.set_xlabel(r"PS Molecular Weight / $\mathrm{kg \ mol^{-1}}$")
# ax.axes.yaxis.set_ticklabels([])

#Temperature in K
T = 273+199

#Data from plot in dynes / cm
experimental = np.array([0.505, 1.0, 1.113, 1.226, 1.313, 1.31, 1.38])

#Molecular weight of PS
experimental_MW = (104.0/1000.0)*np.array([17.0, 39.0, 85.0, 98.0, 332.0, 396.0, 418.0])

# rhoRT = (1000*8.314*T)*np.reciprocal(experimental_MW)

rhoRT = (1000*8.314*(1/0.104)*T)

#Data in SI units
experimental_SI = 0.001*experimental
# experimental_SI = (71.99e-3)*experimental

# experimental_RG = np.array([6.36e-9,6.36e-9,6.36e-9,6.36e-9,1.158e-8, 1.27e-8, 1.3e-8])

experimental_RG = 6.36e-9

scaled_experimental = np.divide(experimental_SI, rhoRT)
scaled_experimental_final = np.divide(scaled_experimental, experimental_RG)
# print(scaled_experimental)

sim_N_Chain_unifac = np.array([17.0, 39.0, 85.0, 98.0, 110.0, 120.0, 130.0, 140.0])

sim_N_Chain_saft = np.array([17.0, 39.0, 85.0, 98.0, 110.0, 120.0, 200.0, 300, 400.0])

sim_N_Chain_saft_CR = np.array([17.0,39.0,85.0, 98.0, 120.0, 150.0, 200.0])

sim_MW_unifac = 0.104*sim_N_Chain_unifac
sim_MW_saft = 0.104*sim_N_Chain_saft
sim_MW_SAFT_CR = 0.104*sim_N_Chain_saft_CR

sim_RG_UNIFAC = (6.36e-10)*np.sqrt(np.sqrt(100*sim_N_Chain_unifac))
sim_RG_SAFT = (6.36e-10)*np.sqrt(np.sqrt(100*sim_N_Chain_saft))
sim_RG_SAFT_CR = (6.36e-10)*np.sqrt(np.sqrt(100*sim_N_Chain_saft_CR))

sim_scale_UNIFAC = rhoRT*sim_RG_UNIFAC
sim_scale_SAFT = rhoRT*sim_RG_SAFT
sim_scale_SAFT_CR = rhoRT*sim_RG_SAFT_CR

# Simulation results
# sim_UNIFAC = np.array([0.02136*(2.61e-9), 0.03636*(3.95e-9), 0.04318*(5.83e-9), 0.04362*(6.26e-9), 0.04385*(6.63e-9), 0.04393*(6.92e-9), 0.04394*(7.21e-9), 0.04390*(7.48e-9)])
# sim_SAFT = np.array([0.004085*(2.61e-9), 0.003901*(3.95e-9), 0.005226*(5.83e-9), 0.005666*(6.26e-9), 0.005976*(6.63e-9), 0.006176*(6.92e-9), 0.006647*(8.94e-9), 0.006609*(1.09e-8),0.006312*(1.24e-8)])

sim_UNIFAC = np.array([0.02136, 0.03636, 0.04318, 0.04362, 0.04385, 0.04393, 0.04394, 0.04390])
sim_SAFT = np.array([0.004085, 0.003901, 0.005226, 0.005666, 0.005976, 0.006176, 0.006647, 0.006609,0.006312])

sim_SAFT_CR = np.array([0.01629,0.03684, 0.04744, 0.04826,0.04896,0.04915, 0.0482])

sim_UNIFAC_SI = np.multiply(sim_UNIFAC, sim_scale_UNIFAC)
sim_SAFT_SI = np.multiply(sim_SAFT, sim_scale_SAFT)
sim_SAFT_CR_SI = np.multiply(sim_SAFT_CR, sim_scale_SAFT_CR)

ax.semilogy(experimental_MW, experimental_SI, "xr", label="Experimental")
ax.semilogy(sim_MW_unifac, sim_UNIFAC_SI, "^k:",label="UNIFAC")
ax.semilogy(sim_MW_SAFT_CR, sim_SAFT_CR_SI, "+k--", label="PC-SAFT (CR)")
ax.semilogy(sim_MW_saft, sim_SAFT_SI,"sk-",label="PC-SAFT (Fitted)")

plt.tight_layout()
plt.legend(frameon=False)
# plt.show()
plt.savefig("surfacetension.png",dpi=300)

 
