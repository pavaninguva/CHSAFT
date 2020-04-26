import numpy as np 
import pandas as pd 
from math import pi
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, least_squares
from Thermo import GibbsMixingPCSAFT,GibbsMixingFH,PCSAFT,GibbsMixingUNIFAC

def PressureSolver(xComp,Temp):
    V1 = np.polyval([7.35086529028252e-10,-1.03249580765898e-06,0.000775563883778852,-3.06499793800575],Temp)
    V2 = np.polyval([7.82715931536264e-10,-1.24621026088872e-06,0.000900911392438747,-2.98549207494120],Temp)

    RK1 = np.polyval([-1.21165605150984e-10,4.20398328285820e-08,-2.62459781413283e-06],Temp)
    RK2 = np.polyval([-4.56383914456231e-11,2.80296255928667e-08,-4.97738063465419e-06],Temp)

    VE  = xComp[1]*(1-xComp[1])*(RK1+RK2*(1-2*xComp[1]))

    return VE+(10**V1)*xComp[0]+(10**V2)*xComp[1]

def GibbsSolver(xComp,Temp):

    RK1 = np.polyval([-6.51228073979852e-08,0.000111074948331728,-0.0673265285969701,16.3051829504820],Temp)
    RK2 = np.polyval([5.74691987750363e-09,-1.02564201990007e-05,0.00661753267719602,-1.75161871071347],Temp)
    RK3 = np.polyval([-8.27913510636852e-10,1.45035546661912e-06,-0.000921029728818829,0.237298198603018],Temp)
    RK4 = np.polyval([1.13263755156765e-10,-1.97504851618488e-07,0.000125473368989465,-0.0316722753241965],Temp)

    gE  = xComp[0]*(1-xComp[0])*(RK1+RK2*(1-2*xComp[0])+RK3*(1-2*xComp[0])**2+RK4*(1-2*xComp[0])**3)

    return gE+xComp[0]*np.log(xComp[0])+xComp[1]*np.log(xComp[1])

def funPCSAFT(x,*arg):
    p,q=arg
    return [q(x[0])-q(x[1]),p(x[1])-p(x[0])+x[0]*q(x[0])-x[1]*q(x[1])]

def funPCSAFTApx(xComp,Temp,*arg):
    k_B=1.38064852e-23
    nsegment, Mr, epsilon, sigma, NParticles, k = arg
    
    Vol1   = PressureSolver([1,0],Temp)
    Vol2   = PressureSolver([0,1],Temp)
    VolMix = PressureSolver(xComp,Temp)

    print(Vol2[0],VolMix[0])

    r1 = PCSAFT([nsegment[0]],[Mr[0]],[epsilon[0]],[sigma[0]],NParticles,Temp,Vol1[0],[1.],[[0.]])
    r2 = PCSAFT([nsegment[1]],[Mr[1]],[epsilon[1]],[sigma[1]],NParticles,Temp,Vol2[0],[1.],[[0.]])
    r  = PCSAFT(nsegment,Mr,epsilon,sigma,NParticles,Temp,VolMix,xComp,k)

    tideal = xComp[0]*np.log(xComp[0])+xComp[1]*np.log(xComp[1])
    tres   = xComp[0]*(r1.Z()+r1.a_res())+xComp[1]*(r2.Z()+r2.a_res())-1.
    tmix   = r.Z()+r.a_res()-1.

    return (tideal+tmix-tres)*NParticles*k_B*Temp



def funFH(x,*arg):
    A,B,C,Temp,Nmono,Vmono=arg
    r1=GibbsMixingFH(Temp,[x[0],1-x[0]],Nmono,Vmono,A,B,C)
    r2=GibbsMixingFH(Temp,[x[1],1-x[1]],Nmono,Vmono,A,B,C)
    return [r1.dGibbsFreeMixing()-r2.dGibbsFreeMixing(),r2.GibbsFreeMixing()-r1.GibbsFreeMixing()+x[0]*r1.GibbsFreeMixing()-x[1]*r2.GibbsFreeMixing()]

k_B=1.38064852e-23
NParticles = np.array([6.02214086e23])
nsegment = np.array([26.95,34.235])
Mr = np.array([1100,1670])
epsilon = np.array([288.84, 348.2])
sigma = np.array([4.097e-10,4.152e-10])
k = np.array([[0.,0.00497],[0.00497,0.]])

xComp1=np.linspace(0.01,0.99,20)
xComp2=1-xComp1

Temp=np.linspace(400,650,20)

kPC=np.zeros(20)
kFH=np.zeros(20)
kUN=np.zeros(20)

xLLEFH1=np.zeros(len(Temp))
xLLEFH2=np.zeros(len(Temp))

xLLEPC1=np.zeros(len(Temp))
xLLEPC2=np.zeros(len(Temp))

xLLEPCA1=np.zeros(len(Temp))
xLLEPCA2=np.zeros(len(Temp))

x0=[0.01,0.99]
x1=[0.01,0.99]
x2=[0.01,0.99]
for j in range(20):
    for i in range(len(xComp1)):
        # rPC=GibbsMixingPCSAFT(nsegment,Mr,epsilon,sigma,NParticles,[Temp[j]],[1e5],[xComp1[i],xComp2[i]],k)
        rFH=GibbsMixingFH(Temp[j],[xComp1[i],xComp2[i]],[1000,1000],[0.179,0.111],-1.57e-2,18.7,0)
        rUN=GibbsMixingUNIFAC(["PMMA","PS"],[xComp1[i],xComp2[i]],Temp[j])
        # kPC[i]=rPC.GibbsFreeMixing()
        kFH[i]=rFH.GibbsFreeMixing()
        kUN[i]=rUN.GibbsFreeMixing()/8.314/Temp[j]
    # gMixPC=plt.plot(xComp2,kPC)
    # plt.savefig('gMixPC.png')
    # gMixFH=plt.plot(xComp2,kFH)
    # plt.savefig('gMixFH.png')
    gMixUN=plt.plot(xComp2,kUN)
    plt.savefig('gMixUN.png')
    
    # p=np.poly1d(np.polyfit(xComp2,kPC,15))
    # q=np.polyder(p,1)

    # n = np.poly1d(np.polyfit(xComp2,kPCA,15))
    # m = np.polyder(n,1)

    # x0=fsolve(funFH,x0,args=(-1.57e-2,18.7,0,Temp[j],[1000,1000],[0.179,0.111]))
    # xLLEFH1[j]=x0[0]
    # xLLEFH2[j]=x0[1]
    
    # x1=fsolve(funPCSAFT,x1,args=(p,q))
    # xLLEPC1[j]=x1[0]
    # xLLEPC2[j]=x1[1]
    
    # x2=fsolve(funPCSAFT,x2,args=(n,m))
    # xLLEPCA1[j]=x2[0]
    # xLLEPCA2[j]=x2[1]
    print(j)

# LLE=plt.plot(xLLEFH1*1100/(1670+xLLEFH1*(1100-1670)),Temp,'-',xLLEFH2*1100/(1670+xLLEFH2*(1100-1670)),Temp,'-',xLLEPC1*1100/(1670+xLLEPC1*(1100-1670)),Temp,'--',xLLEPC2*1100/(1670+xLLEPC2*(1100-1670)),Temp,'--',xLLEPCA1*1100/(1670+xLLEPCA1*(1100-1670)),Temp,'--',xLLEPCA2*1100/(1670+xLLEPCA2*(1100-1670)),Temp,'--')
# plt.savefig('LLE.png')
