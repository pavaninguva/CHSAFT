from LLESolver import LLESolvers
from scipy.optimize import minimize,least_squares,brentq
import numpy as np

def Objective(k,*arg):
    Species,Length,Temp,Pre,x = arg
    # with open('Thermo.py', 'r') as file:
    #     data = file.readlines()
    # data[1058] = '"PB": {"PS":'+str(k)+'},\n'
    # data[1059] = '"PS": {"PB":'+str(k)+'}\n'
    # with open('Thermo.py', 'w') as file:
    #     file.writelines( data )
    r = LLESolvers("GTP","PCSAFT",Species,Length)
    err = 0.
    T = np.linspace(270,np.min(Temp),10)
    X = [0.05,0.95]
    for i in range(len(T)):
        X=r.LLE(T[i],Pre,[X[0],X[1]])
        
    for i in range(len(Temp)):
        X = r.LLE(Temp[i],Pre,[X[0],X[1]])
        print(X)
        err += np.min([abs(1-X[0]-x[i])/x[i],abs(1-X[1]-x[i])/x[i]])
    return err

Temp = [341.95,351.15,360.45,355.75]
x = [0.3,0.4,0.6,0.8]
Pre = 101325.
Length = [2585/54.1,1672/104.1]
Species = ["PB","PS"]


# k = minimize(Objective,0.0004,args=(Species,Length,Temp,Pre,x))
print(Objective(0.00045,Species,Length,Temp,Pre,x))