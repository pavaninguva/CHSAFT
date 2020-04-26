import numpy as np 
import pandas as pd 
from math import pi
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, least_squares
from Thermo import PCSAFT,GibbsMixingUNIFAC

# r=PCSAFT([26.95,34.235],[1100,1670],[288.84, 348.2],[4.097e-10,4.152e-10],[6.02214086e23],[298],[1e-4],[0.5,0.5],[[0.,0.00497],[0.00497,0.]])
# r=PCSAFT([26.95,34.235],[1100,1670],[288.84, 348.2],[4.097e-10,4.152e-10],[6.02214086e23],[298],[1e-3],[0.5,0.5],[[0.,0.00497],[0.00497,0.]])
r = GibbsMixingUNIFAC(["PMMA","PS"],[0.5,0.5],298)
print(r.ln_gamma())

