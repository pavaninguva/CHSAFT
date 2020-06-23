import timeit

import numpy as np 
import pandas as pd 
from math import pi
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, least_squares
from Thermo import ThermoMix,PCSAFT
code = """
from LLESolver import LLESolvers
Temp = 275
Pre = 1e5
Length = [1100/54.1,1340/104.1]
Species = ["PB","PS"]
Method = "PCSAFT"
Solver = "GTP"

r=LLESolvers(Solver,Method,Species,Length)
print(r.LLE(Temp,Pre))
"""
# print(timeit.timeit(code, number=10)/10)
code = """
from LLESolver import LLESolvers
Temp = 275
Pre = 1e5
Length = [1100/54.1,1340/104.1]
Species = ["PB","PS"]
Method = "PCSAFT"
Solver = "vonSolms"

r=LLESolvers(Solver,Method,Species,Length)
r.LLE(Temp,Pre)
"""
# print(timeit.timeit(code, number=10)/10)
# r=LLESolvers("vonSolms",Method,Species,Length)
# print(r.LLE(Temp,Pre),t1-t0)

# r=PCSAFT([26.95,34.235],[1100,1670],[288.84, 348.2],[4.097e-10,4.152e-10],[6.02214086e23],[298],[1e-3],[0.5,0.5],[[0.,0.00497],[0.00497,0.]])
r=ThermoMix("UNIFAC",["PS","PMMA"],[1100/104.1,1670/100],[0.5,0.5],[298])
print(r.dGibbsFreeMixing())
# r = GibbsMixingUNIFAC(["PMMA","PS"],[0.5,0.5],298)
# print('ln gamma: =======')
# print('{:18s}'.format('value: '), r.ln_gamma_comb()+r.ln_gamma_res())
# print('ln gamma comb: =======')
# print('{:18s}'.format('value: '), r.ln_gamma_comb())
# print('{:18s}'.format('phi: '), r.phi())
# print('{:18s}'.format('r: '), r.r_i())
# print('ln gamma FV: =======')
# print('{:18s}'.format('value: '), r.ln_gamma_fv())
# print('{:18s}'.format('V: '), r.volume())
# print('{:18s}'.format('Ci: '), r.C_i_coeff())
# print('{:18s}'.format('Vi: '), r.red_volume())
# print('{:18s}'.format('Vm: '), r.red_volume_mix())
# print('ln gamma res: =======')
# print('{:18s}'.format('value: '), r.ln_gamma_res())
# print('{:18s}'.format('GammaK: '), r.lnGammaK())
# print('{:18s}'.format('X: '), r.Xm())
# print('{:18s}'.format('H: '), r.Hm())
# print('{:18s}'.format('GammaKi: '), r.lnGammaKi())
# print('{:18s}'.format('Xi: '), r.Xi())
# print('{:18s}'.format('Hi: '), r.Hmi())