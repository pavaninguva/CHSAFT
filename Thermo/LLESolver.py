from Thermo import ThermoMix
from scipy.optimize import minimize,least_squares,root,fsolve
import numpy as np
import matplotlib.pyplot as plt

class LLESolvers(object):
 def __init__(self,Solver,Method,Species,Length):
     self.Method = Method
     self.Solver = Solver
     self.Species = Species
     self.Length = Length
 def LLE(self,Temp=298,Pre=1e5,x0=None):
     if self.Solver=="vonSolms":
         r = vonSolmsLLE(self.Method,self.Species,self.Length)
         x = r.LLE(Temp=Temp,Pre=Pre)
     if self.Solver=="GTP":
         r = TangentPlaneSolver(self.Method,self.Species,self.Length)
         x = r.LLE(Temp=Temp,Pre=Pre,x0=x0)
     if self.Solver=="Crude":
         r = CrudeSolver(self.Method,self.Species,self.Length)
         x = r.LLE(Temp=Temp,Pre=Pre,x0=x0)
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
     g0 = r0.g_res()

     r  = ThermoMix(Method,Species,Length,[x[0],1-x[0]],[Temp],[Pre])
     g  = r.g_res()
    #  dg  = r.dGibbsFreeMixing()
     return (g[0]-g0[0])/(x[0]-xSpn)-g[1]
 def LLE(self,Temp=298,Pre=1e5):
     Method = self.Method
     Species = self.Species
     Length = self.Length
     
     # Finding Spinodal points first
     bnds = ((0,1),)
     try:
        xSpn1 = minimize(self.Spinodal,(0.5),method='SLSQP',bounds=bnds,args=(Method,Species,Length,Temp,Pre,-1.))
        xSpn2 = minimize(self.Spinodal,(0.5),method='SLSQP',bounds=bnds,args=(Method,Species,Length,Temp,Pre,1.))

        x1 = [xSpn1.x[0],xSpn1.x[0]]
        x2 = [xSpn2.x[0],1.]
        i=0.
        while (abs(x2[0]-x2[1])>1e-10 or abs(x1[0]-x1[1])>1e-10) and i<100:
            if (i % 2)==0:
                x2[0] = x2[1]
                x = least_squares(self.Equilibrium,(x2[0]+1)/2,bounds=(xSpn2.x[0],1),args=(Method,Species,Length,Temp,Pre,x1[1]))
                x2[1] = x.x[0]
            else:
                x1[0] = x1[1]
                x = least_squares(self.Equilibrium,(x1[0])/2,bounds=(0,xSpn1.x[0]),args=(Method,Species,Length,Temp,Pre,x2[1]))
                x1[1] = x.x[0]
            i+=1.
        A = [x1[1],x2[1]]
     except:
        A = [0.5,0.5]
     return A
class TangentPlaneSolver(object):
 def __init__(self,Method,Species,Length):
     self.Method = Method
     self.Species = Species
     self.Length = Length
 def Equilibrium(self,x,*arg):
     Method,Species,Length,Temp,Pre = arg
     Mr = [Length[0]*54.1,Length[1]*104.1]

     r1  = ThermoMix(Method,Species,Length,[x[0]/Mr[0]/(x[0]*(1/Mr[0]-1/Mr[1])+1/Mr[1]),1-x[0]/Mr[0]/(x[0]*(1/Mr[0]-1/Mr[1])+1/Mr[1])],[Temp],[Pre])
     g1  = r1.g_res()
    #  g1  = [r1.GibbsFreeMixing(),r1.dGibbsFreeMixing()]

     r2  = ThermoMix(Method,Species,Length,[x[1]/Mr[0]/(x[1]*(1/Mr[0]-1/Mr[1])+1/Mr[1]),1-x[1]/Mr[0]/(x[1]*(1/Mr[0]-1/Mr[1])+1/Mr[1])],[Temp],[Pre])
     g2  = r2.g_res()
    #  g2  = [r2.GibbsFreeMixing(),r2.dGibbsFreeMixing()]
     return [g2[1]-g1[1],(g2[1]*x[1]/Mr[0]/(x[1]*(1/Mr[0]-1/Mr[1])+1/Mr[1])-g1[1]*x[0]/Mr[0]/(x[0]*(1/Mr[0]-1/Mr[1])+1/Mr[1]))+g1[0]-g2[0]]
 def Jac(self,x,*arg):
     Method,Species,Length,Temp,Pre = arg
     dx  = 1e-8
     r1  = ThermoMix(Method,Species,Length,[x[0]-dx,1-x[0]+dx],[Temp],[Pre])
     g1  = r1.g_res()
    #  g1  = [r1.GibbsFreeMixing(),r1.dGibbsFreeMixing()]

     r2  = ThermoMix(Method,Species,Length,[x[1]-dx,1-x[1]+dx],[Temp],[Pre])
     g2  = r2.g_res()
    #  g2  = [r2.GibbsFreeMixing(),r2.dGibbsFreeMixing()]
     
     r3  = ThermoMix(Method,Species,Length,[x[0]+dx,1-x[0]-dx],[Temp],[Pre])
     g3  = r3.g_res()
    #  g1  = [r1.GibbsFreeMixing(),r1.dGibbsFreeMixing()]

     r4  = ThermoMix(Method,Species,Length,[x[1]+dx,1-x[1]-dx],[Temp],[Pre])
     g4  = r4.g_res()

     Jac = [[-(g3[1]-g1[1])/2/dx,(g4[1]-g2[1])/2/dx],[-x[0]*(g3[1]-g1[1])/2/dx,x[1]*(g4[1]-g2[1])/2/dx]]   
     return Jac
 def LLE(self,Temp=298,Pre=1e5,x0=None):
     Method = self.Method
     Species = self.Species
     Length = self.Length
     
     # Finding Spinodal points first
     #  try:
    #  x = least_squares(self.Equilibrium,[0.05,0.95],bounds=((0,0),(1,1)),args=(Method,Species,Length,Temp,Pre))
     if not x0:
        x = least_squares(self.Equilibrium,[0.05,0.95],bounds=((0,0),(1,1)),args=(Method,Species,Length,Temp,Pre))
        # x = root(self.Equilibrium,[0.1,0.90],args=(Method,Species,Length,Temp,Pre))
     else:
        x = least_squares(self.Equilibrium,x0,bounds=((0,0),(1,1)),args=(Method,Species,Length,Temp,Pre))
        # x = root(self.Equilibrium,x0,args=(Method,Species,Length,Temp,Pre))
     A = x.x
    #  except:
    #     A = [0.5,0.5]
     return A
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

class CrudeSolver(object):
 def __init__(self,Method,Species,Length):
    self.Method = Method
    self.Species = Species
    self.Length = Length

 def equations(self,x,*arg):
    Method,Species,Length,Temp,Pre = arg
    Mr = [Length[0]*54.1,Length[1]*104.1]

    r1  = ThermoMix(Method,Species,Length,[x[0]/Mr[0]/(x[0]*(1/Mr[0]-1/Mr[1])+1/Mr[1]),1-x[0]/Mr[0]/(x[0]*(1/Mr[0]-1/Mr[1])+1/Mr[1])],[Temp],[Pre])
    g1  = r1.GibbsFreeMixing()
    mu1 = r1.dGibbsFreeMixing()

    r2  = ThermoMix(Method,Species,Length,[x[1]/Mr[0]/(x[1]*(1/Mr[0]-1/Mr[1])+1/Mr[1]),1-x[1]/Mr[0]/(x[1]*(1/Mr[0]-1/Mr[1])+1/Mr[1])],[Temp],[Pre])
    g2  = r2.GibbsFreeMixing()
    # print(x[0], x[1])
    mu2 = r2.dGibbsFreeMixing()
    return (g2-g1-(x[1]-x[0])*mu1, mu2-mu1)

 def LLE(self,Temp=298,Pre=1e5,x0=None):
    Method = self.Method
    Species = self.Species
    Length = self.Length
    return fsolve(self.equations,[0.05,0.95],args=(Method,Species,Length,Temp,Pre))
    # return least_squares(self.equations,[0.05,0.95],bounds=((0,0),(1,1)),args=(Method,Species,Length,Temp,Pre)).x

