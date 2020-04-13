import numpy as np 
import pandas as pd 
import glob as glob
import re

import matplotlib.pyplot as plt

# Load Gibbs CSV files from resuts directory
gibbs_csv_list = glob.glob("*.csv")


plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
plt.rc("axes", titlesize="small")
plt.rc("legend", fontsize="small")

fig = plt.figure(figsize=(4, 4))
ax = fig.add_subplot(1, 1, 1)
ax.set_ylabel(r'System Gibbs Energy $G_{system}$')
ax.set_xlabel(r'Time $\tilde{t}$')
ax.set_xlim([0,500000])


# We slice the list to just work with the 1st n items
for run in gibbs_csv_list:
    # Read in the csv
    df = pd.read_csv(run)
    # Extract the run number for label
    run_name_extract = re.search(r"gibbs_(.+?).csv", run)
    run_name = run_name_extract.group(1)
    print (run_name)
    if "Vas94_constant" in run_name:
        label_1023= r"Case 1, Constant mobility "
        ax.plot(df["time"][:1600000], df["gibbs"][:1600000], label = label_1023, color="k")
        # Add the label to the last point of each line. Makes it easier to view
        # ax.text(df["time"].iloc[-1], df["gibbs"].iloc[-1], run_name)
    if "Vas94_var" in run_name:
        label_1029= r"Case 1, Full model"
        ax.plot(df["time"], df["gibbs"], label = label_1029, color="c")
        # Add the label to the last point of each line. Makes it easier to view
        # ax.text(df["time"].iloc[-1], df["gibbs"].iloc[-1], run_name)
    if "He97_constant" in run_name:
        label_0= r"Case 2, Constant mobility"
        ax.plot(df["time"], df["gibbs"], label = label_0, color="b")
        # Add the label to the last point of each line. Makes it easier to view
        # ax.text(df["time"].iloc[-1], df["gibbs"].iloc[-1], run_name)
    if "He97_var" in run_name:
        label_1003= r"Case 2, Full model"
        ax.plot(df["time"], df["gibbs"], label = label_1003, color="r")
        # Add the label to the last point of each line. Makes it easier to view
        # ax.text(df["time"].iloc[-1], df["gibbs"].iloc[-1], run_name)
    
    
    ax.legend()

plt.tight_layout()
# plt.show()
plt.savefig("gibbsplot.png", dpi=600)