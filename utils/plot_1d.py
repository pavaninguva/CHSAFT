import numpy as np 
import pandas as pd 

import matplotlib.pyplot as plt

#Load CSV File
PS_PB_csv_loc = "../PSPB_1D.csv"
PS_PMMA_csv_loc = "../PSPMMA_1D.csv"

PSPB_df = pd.read_csv(PS_PB_csv_loc)
PSPMMA_df = pd.read_csv(PS_PMMA_csv_loc)


plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
plt.rc("axes", titlesize="small")
plt.rc("legend", fontsize="small")

fig1, ax1 = plt.subplots(figsize=(4, 4))
ax1.set_ylabel(r'$\phi_{1}$')
ax1.set_xlabel(r'$\tilde{x}_{0}$')
# ax1.set_xlim([0,400000])
ax1.xaxis.set_major_locator(plt.MaxNLocator(5))

ax1.plot(PSPMMA_df["Points:0"], PSPMMA_df["f_5-0"] ,color="black")


fig2, ax2 = plt.subplots(figsize=(4, 4))
ax2.set_ylabel(r'$\phi_{1}$')
ax2.set_xlabel(r'$\tilde{x}_{0}$')
# ax2.set_xlim([0,400000])
ax2.xaxis.set_major_locator(plt.MaxNLocator(5))

ax2.plot(PSPB_df["Points:0"], PSPB_df["f_5-0"] , color="black")


# Save figs
fig1.savefig("PSPMMA_1D.png", dpi=300)
fig2.savefig("PSPB_1D.png", dpi=300)



