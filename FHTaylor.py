# import matplotlib.pyplot as plt

import numpy as np

from sympy import *

import math

# from scipy.optimize import curve_fit


"""
This section is for doing the full Taylor approximation for the FH equation
"""

x = symbols("x")

N_1, N_2, chi = symbols("N_1 N_2 chi", real=True, constant=True)

f = (x*log(x) / N_1) + ((1.0-x)*log(1.0-x) / N_2) + x*(1.0-x)*chi

dfdx = diff(f, x)

d2fdx2 = diff(dfdx, x)
d3fdx3 = diff(d2fdx2, x)
d4fdx4 = diff(d3fdx3, x)
d5fdx5 = diff(d4fdx4, x)
d6fdx6 = diff(d5fdx5, x)
d7fdx7 = diff(d6fdx6, x)
d8fdx8 = diff(d7fdx7, x)
d9fdx9 = diff(d8fdx8, x)
d10fdx10 = diff(d9fdx9, x)
d11fdx11 = diff(d10fdx10, x)
d12fdx12 = diff(d11fdx11, x)
d13fdx13 = diff(d12fdx12, x)
d14fdx14 = diff(d13fdx13, x)
d15fdx15 = diff(d14fdx14, x)
d16fdx16 = diff(d15fdx15, x)
d17fdx17 = diff(d16fdx16, x)
d18fdx18 = diff(d17fdx17, x)
d19fdx19 = diff(d18fdx18, x)
d20fdx20 = diff(d19fdx19, x)




def fhfull(N_1, N_2, chi, x):
    f = (x*np.log(x) / N_1) + ((1.0-x)*np.log(1.0-x) / N_2) + x*(1.0-x)*chi
    return f


def dfdx_fh(N_1, N_2, chi, x):
    f = -chi*x + chi*(1.0 - x) - np.log(1.0 - x)/N_2 - \
        1/N_2 + np.log(x)/N_1 + 1/N_1
    return f


def d2fdx2_fh(N_1, N_2, chi, x):
    f = -2*chi + 1/(N_2*(1.0 - x)) + 1/(N_1*x)
    return f


def d3fdx3_fh(N_1, N_2, x):
    f = 1/(N_2*(1.0 - x)**2) - 1/(N_1*x**2)
    return f


def d4fdx4_fh(N_1, N_2, x):
    f = 2/(N_2*(1.0 - x)**3) + 2/(N_1*x**3)
    return f


def d5fdx5_fh(N_1, N_2, x):
    f = 6/(N_2*(1.0 - x)**4) - 6/(N_1*x**4)
    return f


def d6fdx6_fh(N_1, N_2, x):
    f = 24/(N_2*(1.0 - x)**5) + 24/(N_1*x**5)
    return f


def d7fdx7_fh(N_1, N_2, x):
    f = 120/(N_2*(1.0 - x)**6) - 120/(N_1*x**6)
    return f


def d8fdx8_fh(N_1, N_2, x):
    f = 720/(N_2*(1.0 - x)**7) + 720/(N_1*x**7)
    return f


def d9fdx9_fh(N_1, N_2, x):
    f = 5040/(N_2*(1.0 - x)**8) - 5040/(N_1*x**8)
    return f

def d10fdx10_fh(N_1, N_2, x):
    f = 40320/(N_2*(1.0 - x)**9) + 40320/(N_1*x**9)
    return f

def d11fdx11_fh(N_1, N_2, x):
    f = 362880/(N_2*(1.0 - x)**10) - 362880/(N_1*x**10)
    return f

def d12fdx12_fh(N_1, N_2, x):
    f = 3628800/(N_2*(1.0 - x)**11) + 3628800/(N_1*x**11)
    return f

def d13fdx13_fh(N_1, N_2, x):
    f = 39916800/(N_2*(1.0 - x)**12) - 39916800/(N_1*x**12)
    return f

def d14fdx14_fh(N_1, N_2, x):
    f = 479001600/(N_2*(1.0 - x)**13) + 479001600/(N_1*x**13)
    return f

def d15fdx15_fh(N_1, N_2, x):
    f = 6227020800/(N_2*(1.0 - x)**14) - 6227020800/(N_1*x**14)
    return f

def d16fdx16_fh(N_1, N_2, x):
    f = 87178291200/(N_2*(1.0 - x)**15) + 87178291200/(N_1*x**15)
    return f

def d17fdx17_fh(N_1, N_2, x):
    f = 1307674368000/(N_2*(1.0 - x)**16) - 1307674368000/(N_1*x**16)
    return f

def d18fdx18_fh(N_1, N_2, x):
    f = 20922789888000/(N_2*(1.0 - x)**17) + 20922789888000/(N_1*x**17)
    return f

def d19fdx19_fh(N_1, N_2, x):
    f = 355687428096000/(N_2*(1.0 - x)**18) - 355687428096000/(N_1*x**18)
    return f

def d20fdx20_fh(N_1, N_2, x):
    f = 6402373705728000/(N_2*(1.0 - x)**19) + 6402373705728000/(N_1*x**19)
    return f

# def d10fdx10_fh(N_1, N_2, x):
#     f = 
#     return f


def taylorapprox_fullFH(N_1, N_2, chi, x, a=0.5):
    f0 = fhfull(N_1, N_2, chi, a)  
    f1 = (dfdx_fh(N_1, N_2, chi, a) / (math.factorial(1))) * (x-a) 
    f2 = (d2fdx2_fh(N_1, N_2, chi, a) / (math.factorial(2)))*(x-a)**2 
    f3 = (d3fdx3_fh(N_1, N_2, a) / (math.factorial(3)))*(x-a)**3 
    f4 = (d4fdx4_fh(N_1, N_2, a) / (math.factorial(4)))*(x-a)**4
    f5 = (d5fdx5_fh(N_1, N_2, a) / (math.factorial(5)))*(x-a)**5
    f6 = (d6fdx6_fh(N_1, N_2, a) / (math.factorial(6)))*(x-a)**6
    f7 = (d7fdx7_fh(N_1, N_2, a) / (math.factorial(7)))*(x-a)**7 
    f8 = (d8fdx8_fh(N_1, N_2, a) / (math.factorial(8)))*(x-a)**8 
    f9 = (d9fdx9_fh(N_1, N_2, a) / (math.factorial(9)))*(x-a)**9
    f10 = (d10fdx10_fh(N_1, N_2, a) / (math.factorial(10)))*(x-a)**10
    f11 = (d11fdx11_fh(N_1, N_2, a) / (math.factorial(11)))*(x-a)**11
    f12 = (d12fdx12_fh(N_1, N_2, a) / (math.factorial(12)))*(x-a)**12
    f13 = (d13fdx13_fh(N_1, N_2, a) / (math.factorial(13)))*(x-a)**13
    f14 = (d14fdx14_fh(N_1, N_2, a) / (math.factorial(14)))*(x-a)**14
    f15 = (d15fdx15_fh(N_1, N_2, a) / (math.factorial(15)))*(x-a)**15
    f16 = (d16fdx16_fh(N_1, N_2, a) / (math.factorial(16)))*(x-a)**16
    f17 = (d17fdx17_fh(N_1, N_2, a) / (math.factorial(17)))*(x-a)**17
    f18 = (d18fdx18_fh(N_1, N_2, a) / (math.factorial(18)))*(x-a)**18
    f19 = (d19fdx19_fh(N_1, N_2, a) / (math.factorial(19)))*(x-a)**19
    f20 = (d20fdx20_fh(N_1, N_2, a) / (math.factorial(20)))*(x-a)**20

    f = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 +f10 + f11 + f12 + f13 + f14 + f15 +f16 +f17 + f18 +f19 +f20

    return f



"""
Now we implement the case where we do a Taylor approximation on just the ln function
"""

# write expressions for the logarithms and evaluating the series expansion. 

g1 = log(x)

g2 = log(1-x)

# We need to find a way to slot a "Rational" into this declaration to enable the output to be nice workable fractions.
series = g2.series(x, Rational(1/2), 11).removeO()

# print (series)

def taylorapprox_logonlyFH (N_1, N_2, chi, x):
    combinatorial_1 = x*(2.0*x - 512.0*(x - 0.5)**10.0/5.0 + 512.0*(x - 1.0/2.0)**9.0/9.0 - 32.0*(x - 0.5)**8.0 + 128*(x - 0.5)**7.0/7.0 - 32.0*(x - 0.5)**6.0/3.0 + 32.0*(x - 0.5)**5.0/5.0 - 4.0*(x - 0.5)**4.0 + 8.0*(x - 0.5)**3.0/3.0 - 2.0*(x - 0.5)**2.0 - 1.0 - np.log(2.0)) / N_1

    combinatorial_2 = (1.0 - x)*(-2.0*x - 512.0*(x - 0.5)**10.0/5.0 - 512.0*(x - 0.5)**9.0/9.0 - 32.0*(x - 0.5)**8.0 - 128.0*(x - 0.5)**7.0/7.0 - 32.0*(x - 0.5)**6.0/3.0 - 32.0*(x - 0.5)**5.0/5.0 - 4.0*(x - 0.5)**4.0 - 8.0*(x - 0.5)**3.0/3.0 - 2.0*(x - 0.5)**2.0 - np.log(2) + 1)/ N_2

    residual = x*(1.0-x)*chi

    f = combinatorial_1 + combinatorial_2 + residual

    return f



def polynomialfit_fullFH (N_1, N_2, chi, x, order):

    # Define a linspace between 0.0 - 1.0: 
    var = np.linspace(1e-10, (1.0-(1e-10)), 2000)

    g_list = []
    for mol in var:
        g = fhfull(N_1, N_2, chi, mol)
        g_list.append(g)
    
    f = np.poly1d(np.polyfit(var, g_list, order))

    f_x = f(x)

    return f_x

def symquarticdoublewell (x, A):
    g = A * (x**2.0) * (1.0- x)**2.0       
    return g

def heatofmixing (chi, x):
    g = x *(1.0 - x)* chi

    return g

# def symquarticdoublewell_fullFH (N_1, N_2, chi, x):

#     # Define a linspace between 0.0 - 1.0: 
#     var = np.linspace(1e-10, (1.0-(1e-10)), 2000)

#     g_list = []
#     for mol in var:
#         g = fhfull(N_1, N_2, chi, mol)
#         g_list.append(g)

#     popt, pcov = curve_fit(symquarticdoublewell, var, g_list)

#     f = popt[0] * (x**2.0) * (1.0 -x)**2.0  

#     return f, popt[0]


# a, b = symquarticdoublewell_fullFH (2960, 3500, 0.0064, 0.5)

# print (b)


    
     
                    

# print(symquarticdoublewell_fullFH(500, 500, 0.006, 0.4))

# fig, ax = plt.subplots()


# fhlist = []
# approxlist = []
# x_a = np.linspace(0, 1.0, 1000)

# for i in x_a:
#     fh = fhfull(500, 5000, 0.0064, i)
#     approx = symquarticdoublewell_fullFH(500, 5000, 0.0064, i)
#     fhlist.append(fh)
#     approxlist.append(approx)

# ax.plot(x_a, fhlist, label="analytical")
# ax.plot(x_a, approxlist, label="approx")
# ax.legend()


# plt.show()
