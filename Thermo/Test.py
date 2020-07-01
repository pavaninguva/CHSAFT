import timeit

import numpy as np 
import pandas as pd 
from math import pi
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, least_squares
from Thermo import ThermoMix,PCSAFT
from LLESolver import LLESolvers
# print(timeit.timeit(code, number=10)/10)
# r=LLESolvers("vonSolms",Method,Species,Length)
# print(r.LLE(Temp,Pre),t1-t0)

# r=PCSAFT([26.95,34.235],[1100,1670],[288.84, 348.2],[4.097e-10,4.152e-10],[6.02214086e23],[298],[1e-3],[0.5,0.5],[[0.,0.00497],[0.00497,0.]])
Length = [10,10]
Species = ["PB","PS"]
Temp = 370
SAFT  = ThermoMix("PCSAFT",Species,Length,[Temp])
SAFT_CH = ThermoMix("PCSAFT",Species,Length,[Temp],CH="On")
print(SAFT.GibbsFreeMixing(0.5),SAFT_CH.GibbsFreeMixing(0.5))
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