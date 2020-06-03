from Thermo import ThermoMix
from scipy.optimize import minimize,least_squares,brentq
import numpy as np
import matplotlib.pyplot as plt

class LLESolvers(object):
 def __init__(self,Solver,Method,Species,Length):
     self.Method = Method
     self.Solver = Solver
     self.Species = Species
     self.Length = Length
 def LLE(self,Temp=298,Pre=1e5):
     if self.Solver=="vonSolms":
         r = vonSolmsLLE(self.Method,self.Species,self.Length)
         x = r.LLE(Temp,Pre)
     else:
         r = GibbsTP(self.Method,self.Species,self.Length)
         x = r.LLE(Temp,Pre)
     return x

class vonSolmsLLE(object):
 def __init__(self,Method,Species,Length):
     self.Method = Method
     self.Species = Species
     self.Length = Length
 
 def Spinodal(self,x,*arg):
     Method,Species,Length,Temp,Pre,m = arg
     r  = ThermoMix(Method,Species,Length,[x[0],1-x[0]],[Temp],[Pre])
     return m*r.dGibbsFreeMixing()

 def Equilibrium(self,x,*arg):
     Method,Species,Length,Temp,Pre,xSpn = arg
     r0 = ThermoMix(Method,Species,Length,[xSpn,1-xSpn],[Temp],[Pre])
     g0 = r0.GibbsFreeMixing()

     r  = ThermoMix(Method,Species,Length,[x[0],1-x[0]],[Temp],[Pre])
     g  = r.GibbsFreeMixing()
     dg  = r.dGibbsFreeMixing()
     return (g-g0)/(x[0]-xSpn)-dg
 def LLE(self,Temp=298,Pre=1e5):
     Method = self.Method
     Species = self.Species
     Length = self.Length
     
     # Finding Spinodal points first
     bnds = ((0,1),)
     xSpn1 = minimize(self.Spinodal,(0.5),method='SLSQP',bounds=bnds,args=(Method,Species,Length,Temp,Pre,-1.))
     xSpn2 = minimize(self.Spinodal,(0.5),method='SLSQP',bounds=bnds,args=(Method,Species,Length,Temp,Pre,1.))

     x1 = [xSpn1.x[0],xSpn1.x[0]]
     x2 = [xSpn2.x[0],1.]
     i=0.
     while (abs(x2[0]-x2[1])>1e-6 or abs(x1[0]-x1[1])>1e-6) and i<100:
         if (i % 2)==0:
             x2[0] = x2[1]
             x = least_squares(self.Equilibrium,(x2[0]+1)/2,bounds=(xSpn2.x[0],1),args=(Method,Species,Length,Temp,Pre,x1[1]))
             x2[1] = x.x[0]
         else:
             x1[0] = x1[1]
             x = least_squares(self.Equilibrium,(x1[0])/2,bounds=(0,xSpn1.x[0]),args=(Method,Species,Length,Temp,Pre,x2[1]))
             x1[1] = x.x[0]
         i+=1.
     return [x1[1],x2[1]]
     
# Temp = 350
# Pre = 1e5
# Length = [1100/54.1,1670/104.1]
# Species = ["PMMA","PS"]
# Method = "FH"

# bnds = ((0,1),)

# xSpn2 = minimize(Spinodal,(0.5),method='SLSQP',bounds=bnds,args=(Method,Species,Length,Temp,Pre,1.))
# xSpn1 = minimize(Spinodal,(0.5),method='SLSQP',bounds=bnds,args=(Method,Species,Length,Temp,Pre,-1.))
# print(xSpn1.x,xSpn2.x)

# x = np.linspace(0.00001,0.9999,100)
# f = []
# for i in range(100):
#     f.append(ThermoMix(Method,Species,Length,[x[i],1.-x[i]],[Temp],[Pre]).GibbsFreeMixing())

# plt.plot(x,f,color='0.75', ls='dotted',label=r"\textit{s}PC-SAFT")

# x1 = xSpn1.x[0]
# x2 = xSpn2.x[0]
# fSpn1 = ThermoMix(Method,Species,Length,[x1,1.-x1],[Temp],[Pre]).GibbsFreeMixing()
# fSpn2 = ThermoMix(Method,Species,Length,[x2,1.-x2],[Temp],[Pre]).GibbsFreeMixing()
# plt.plot(x1,fSpn1,color='0.5',  marker="+")
# plt.plot(x2,fSpn2,color='0.5',  marker="+")
# for i in range(20):
#     if (i % 2)==0:
#         x2 = least_squares(Equilibrium,(x2+1)/2,bounds=(xSpn2.x[0],1),args=(Method,Species,Length,Temp,Pre,x1))
#         x2 = x2.x[0]
#         # x2 = brentq(Equilibrium,xSpn2.x[0],0.9999,args=(Method,Species,Length,Temp,Pre,x1))
#         print(x2)
#     else:
#         x1 = least_squares(Equilibrium,(x1)/2,bounds=(0,xSpn1.x[0]),args=(Method,Species,Length,Temp,Pre,x2))
#         x1 = x1.x[0]
#         # x1 = brentq(Equilibrium,0.0001,xSpn1.x[0],args=(Method,Species,Length,Temp,Pre,x2))
#         print(x1)

# fBi1 = ThermoMix(Method,Species,Length,[x1,1.-x1],[Temp],[Pre]).GibbsFreeMixing()
# fBi2 = ThermoMix(Method,Species,Length,[x2,1.-x2],[Temp],[Pre]).GibbsFreeMixing()  
# plt.plot([x1,x2],[fBi1,fBi2],color='k', ls='solid')  
# # print(Equilibrium(0.5,Method,Species,Length,Temp,Pre,xSpn.x))
# plt.savefig('gMix.png')

