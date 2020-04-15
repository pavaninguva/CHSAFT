import numpy as np 
import pandas as pd 
from math import pi
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, least_squares
from Thermo import PCSAFT

# r=PCSAFT([26.95,34.235],[1100,1670],[288.84, 348.2],[4.097e-10,4.152e-10],[6.02214086e23],[298],[1e-4],[0.5,0.5],[[0.,0.00497],[0.00497,0.]])
r=PCSAFT([26.95,34.235],[1100,1670],[288.84, 348.2],[4.097e-10,4.152e-10],[6.02214086e23],[298],[1e-3],[0.5,0.5],[[0.,0.00497],[0.00497,0.]])

print('{:18s}'.format('a_res: '), r.a_res())
print('{:18s}'.format('a_disp: '), r.a_disp())
print('{:18s}'.format('a_hc: '), r.a_hc())
print('{:18s}'.format('a_hs: '), r.a_hs())
print('{:18s}'.format('Z: '), r.Z())
print('{:18s}'.format('Z_disp: '), r.Z_Disp())
print('{:18s}'.format('Z_HC: '), r.Z_HC())
print('{:18s}'.format('Z_HS: '), r.Z_HS())