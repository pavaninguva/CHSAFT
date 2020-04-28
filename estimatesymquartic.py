"""
This code needs to be run outside the fenics environment since neither;
the singularity container and the fenics conda environment are incompatible;
with Scipy. 
"""
import numpy as np
from scipy.optimize import curve_fit

def fhfull(N_1, N_2, chi, x):
    f = (x*np.log(x) / N_1) + ((1.0-x)*np.log(1.0-x) / N_2) + x*(1.0-x)*chi
    return f

def symquarticdoublewell (x, A):
    g = A * (x**2.0) * (1.0- x)**2.0       
    return g

def symquarticdoublewell_fullFH (N_1, N_2, chi, x):

    # Define a linspace between 0.0 - 1.0: 
    var = np.linspace(1e-10, (1.0-(1e-10)), 2000)

    g_list = []
    for mol in var:
        g = fhfull(N_1, N_2, chi, mol)
        g_list.append(g)

    popt, pcov = curve_fit(symquarticdoublewell, var, g_list)

    f = popt[0] * (x**2.0) * (1.0 -x)**2.0  

    return f, popt[0]

a, b = symquarticdoublewell_fullFH (800, 1400, 0.004, 0.5)


# This gives us the optimal A value 
print (b)