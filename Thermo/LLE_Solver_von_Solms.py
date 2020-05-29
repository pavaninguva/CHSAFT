from Thermo import ThermoMix
from scipy.optimize import minimize

def Spinodal(x,*arg):
    Method,Species,Length,Temp,Pre = arg
    r  = ThermoMix(Method,Species,Length,[x[0],1-x[0]],[Temp],[Pre])
    return r.dGibbsFreeMixing()

Temp = 298
Pre = 1e5
Length = [1100/104.1,1670/100]
Species = ["PS","PMMA"]
Method = "PCSAFT"

bnds = ((0,1),)

xSpn = minimize(Spinodal,(0.5),method='SLSQP',bounds=bnds,args=(Method,Species,Length,Temp,Pre))
print(xSpn)