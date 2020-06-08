import numpy as np 
import pandas as pd 
import copy
from math import pi
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, least_squares


class ThermoMix(object):
 def __init__(self,Method,Species,Length,xComp,Temp=[298],Pre=[1e5]):
         self.Method = Method
         self.Species = Species
         self.Length = Length
         self.xComp = xComp
         self.Temp = Temp
         self.Pre = Pre
 def GibbsFreeMixing(self):
       if self.Method == "FH":
              r = GibbsMixingFH(self.Species,self.xComp,self.Length)
              gMix = r.GibbsFreeMixing()
       elif self.Method == "UNIFAC":
              r = GibbsMixingUNIFAC(self.Species,self.Length,self.xComp,self.Temp)
              gMix = r.GibbsFreeMixing()
       elif self.Method == "PCSAFT":
              r = GibbsMixingPCSAFT(self.Species,self.Length,self.Temp,self.Pre,self.xComp)
              A = r.GibbsFreeMixing()
              gMix = A[0]
       return gMix
 def dGibbsFreeMixing(self):
       if self.Method == "FH":
              r = GibbsMixingFH(self.Species,self.xComp,self.Length)
              dgMix = r.dGibbsFreeMixing()
       elif self.Method == "UNIFAC":
              r = GibbsMixingUNIFAC(self.Species,self.Length,self.xComp,self.Temp)
              dgMix = r.dGibbsFreeMixing()
       elif self.Method == "PCSAFT":
              r = GibbsMixingPCSAFT(self.Species,self.Length,self.Temp,self.Pre,self.xComp)
              A = r.dGibbsFreeMixing()
              dgMix = A[0]
       return dgMix
 def g_res(self):
       if self.Method == "FH":
              print('Not supported')
       elif self.Method == "UNIFAC":
              print('Not supported')
       elif self.Method == "PCSAFT":
              r = GibbsMixingPCSAFT(self.Species,self.Length,self.Temp,self.Pre,self.xComp)
              A = r.g_res()
       return A

class PCSAFT(object):
 def __init__(self, nsegment, Mr, epsilon,sigma,NParticles,Temperature,Volume,xComp,k):
        self.epsilon = np.array(epsilon)
        self.sigma = np.array(sigma)
        self.k=np.array(k)
        self.nsegment = np.array(nsegment)
        self.Mr = np.array(Mr)
        self.xComp=np.array(xComp)
        self.NParticles=np.array(NParticles)
        self.Temp=np.array(Temperature)
        self.Vol=np.array(Volume)
    ## Residual Helmholtz Free Energy ##
 def a_res(self):
        a_hc=self.a_hc()
        a_disp=self.a_disp()
        A=a_hc+a_disp
        return A

    ## Hard Chain Term ##
 def a_hc(self):
        result=0.
        xComp=self.xComp
        nsegment=self.nsegment
        g_hs=self.g_hs()
        a_hs=self.a_hs()
        m_mean=self.m_mean()
        
        # PC-SAFT
       #  for i in range(len(xComp)):
       #      B=xComp[i]*(nsegment[i]-1)*np.log(g_hs[i])
       #      result=result+B
       #  A=m_mean*a_hs-result
        
        # sPC-SAFT
        A=m_mean*a_hs-(m_mean-1)*np.log(g_hs)
        return A

 def m_mean(self):
        result=0.
        xComp=self.xComp
        nsegment=self.nsegment
        for i in range(len(xComp)):
            B=xComp[i]*nsegment[i]
            result=result+B
        return result

 def a_hs(self):
        zeta0=self.zetan(0)
        zeta1=self.zetan(1)
        zeta2=self.zetan(2)
        zeta3=self.zetan(3)
        
        # PC-SAFT
       #  A=1/zeta0*((3*zeta1*zeta2)/(1-zeta3)+pow(zeta2,3)/(zeta3*pow((1-zeta3),2))+(pow(zeta2,3)/pow(zeta3,2)-zeta0)*np.log(1-zeta3))
        
        # sPC-SAFT
        A=(4*zeta3-3*pow(zeta3,2))/pow(1-zeta3,2)
        return A

 def g_hs(self):
        zeta2=self.zetan(2)
        zeta3=self.zetan(3)
        HSd=self.HSd()
        
        # PC-SAFT
       #  A=HSd
       #  for i in range(len(HSd)):
       #      A[i]=1/(1-zeta3)+pow(HSd[i],2)*3*zeta2/(2*HSd[i]*pow((1-zeta3),2))+pow(HSd[i],4)*2*pow(zeta2,2)/((pow(2*HSd[i],2))*pow((1-zeta3),3))
        
        # sPC-SAFT
        A=(1-zeta3/2)/pow(1-zeta3,3)    
        return A

 def zetan(self,m):
        result=0.
        HSd=self.HSd()
        Vol=self.Vol
        NParticles=self.NParticles
        nsegment=self.nsegment
        xComp=self.xComp
        for i in range(len(xComp)):
            result=xComp[i]*nsegment[i]*pow(HSd[i],m)+result
        return pi*NParticles/(6*Vol)*result

 def HSd(self):
        sigma=self.sigma
        epsilon=self.epsilon
        Temp=self.Temp
        return sigma*(1-0.12*np.exp(-3*epsilon/(Temp)))
 def sHSd(self):
        HSd=self.HSd()
        nsegment=self.nsegment
        xComp=self.xComp
        A=0
        B=0
        for i in range(len(xComp)):
               A+=xComp[i]*nsegment[i]*pow(HSd[i],3)
               B+=xComp[i]*nsegment[i]
        return pow(A/B,1/3)

    ##Dispersion Term
 def a_disp(self):
        I1=self.I1()
        I2=self.I2()
        m2eo3=self.m2eno3(1)
        m2e2o3=self.m2eno3(2)
        C1=self.C1()
        m_mean=self.m_mean()
        NParticles=self.NParticles
        Vol=self.Vol
        return -2*pi*NParticles*I1*m2eo3/Vol-pi*m_mean*NParticles*C1*I2*m2e2o3/Vol

 def C1(self):
        pfrac=self.zetan(3)
        m_mean=self.m_mean()
        return pow(1+m_mean*((8*pfrac-2*pow(pfrac,2))/pow(1-pfrac,4))+(1-m_mean)*(20*pfrac-27*pow(pfrac,2)+12*pow(pfrac,3)-2*pow(pfrac,4))/pow((1-pfrac)*(2-pfrac),2),-1)

 def m2eno3(self,m):
        xComp=self.xComp
        nsegment=self.nsegment
        epsilon=self.epsilon
        sigma=self.sigma
        Temp=self.Temp
        k=self.k
        result=0.
        for i in range(len(xComp)):
            for j in range(len(xComp)):
                result=xComp[i]*xComp[j]*nsegment[i]*nsegment[j]*pow((np.sqrt(epsilon[i]*epsilon[j])*(1-k[i,j]))/(Temp),m)*pow(0.5*(sigma[i]+sigma[j]),3)+result
        return result

 def I1(self):
        acorr=np.array([[0.9105631445,-0.3084016918,-0.0906148351],
                        [0.6361281449,0.1860531159,0.4527842806],
                        [2.6861347891,-2.5030047259,0.5962700728],
                        [-26.547362491,21.419793629,-1.7241829131],
                        [97.759208784,-65.255885330,-4.1302112531],
                        [-159.59154087,83.318680481,13.776631870],
                        [91.297774084,-33.746922930,-8.6728470368]])
        result=0.
        m_mean=self.m_mean()
        pfrac=self.zetan(3)
        a=np.array([0.,0.,0.,0.,0.,0.,0.])
        for i in range(7):
            a[i]=acorr[i,0]+(m_mean-1)/m_mean*acorr[i,1]+(m_mean-1)*(m_mean-2)/pow(m_mean,2)*acorr[i,2]
            result=a[i]*pow(pfrac,i)+result
        return result

 def I2(self):
        bcorr=np.array([[0.7240946941,-0.5755498075,0.0976883116],
                        [2.2382791861,0.6995095521,-0.2557574982],
                        [-4.0025849485,3.8925673390,-9.1558561530],
                        [-21.003576815,-17.215471648,20.642075974],
                        [26.855641363,192.67226447,-38.804430052],
                        [206.55133841,-161.82646165,93.626774077],
                        [-355.60235612,-165.20769346,-29.666905585]])
        result=0.
        m_mean=self.m_mean()
        pfrac=self.zetan(3)
        b=np.array([0.,0.,0.,0.,0.,0.,0.])
        for i in range(7):
            b[i]=bcorr[i,0]+(m_mean-1)/m_mean*bcorr[i,1]+(m_mean-1)*(m_mean-2)/pow(m_mean,2)*bcorr[i,2]
            result=b[i]*pow(pfrac,i)+result
        return result
 ## Compressibility Factor ##
 def Z(self):
    Z_HC=self.Z_HC()
    Z_Disp=self.Z_Disp()
    return 1+Z_HC+Z_Disp

 def Z_HC(self):
    m_mean=self.m_mean()
    Z_HS=self.Z_HS()
    xComp=self.xComp
    nsegment=self.nsegment
    g_hs=self.g_hs()
    dg_hs=self.dg_hs()
    t1=m_mean*Z_HS

    # PC-SAFT
#     t2=0.
#     for i in range(len(xComp)):
#         t2=xComp[i]*(nsegment[i]-1)*pow(g_hs[i],-1)*dg_hs[i]+t2
#     A = t1-t2
    
    # sPC-SAFT
    A = t1+(1-m_mean)*pow(g_hs,-1)*dg_hs
    return A

 def Z_HS(self):
    zeta0=self.zetan(0)
    zeta1=self.zetan(1)
    zeta2=self.zetan(2)
    zeta3=self.zetan(3)
    
    # sPC-SAFT
    A=2*zeta3*(2-zeta3)/pow(1-zeta3,3)
    
    # PC-SAFT
#     A=zeta3/(1-zeta3)+3*zeta1*zeta2/(zeta0*pow((1-zeta3),2))+(3*pow(zeta2,3)-zeta3*pow(zeta2,3))/(zeta0*pow((1-zeta3),3))
    return A

 def dg_hs(self):
    zeta2=self.zetan(2)
    zeta3=self.zetan(3)
    HSd=self.HSd()
    
    # PC-SAFT
#     A=HSd*0.
#     for i in range(len(HSd)):
#         A[i]=zeta3/pow((1-zeta3),2)+HSd[i]/2*(3*zeta2/pow((1-zeta3),2)+6*zeta2*zeta3/pow((1-zeta3),3))+pow(HSd[i]/2,2)*(4*pow(zeta2,2)/pow((1-zeta3),3)+6*pow(zeta2,2)*zeta3/pow((1-zeta3),4))
    
    # sPC-SAFT
    A=zeta3*(5-2*zeta3)/(2*pow(1-zeta3,4))
    return A 
 def Z_Disp(self):
    m2eo3=self.m2eno3(1)
    m2e2o3=self.m2eno3(2)
    I2=self.I2()
    C1=self.C1()
    dI1=self.dI1()
    dI2=self.dI2()
    m_mean=self.m_mean()
    NParticles=self.NParticles
    Vol=self.Vol
    pfrac=self.zetan(3)
    C2=-pow(C1,2)*(m_mean*(-4*pow(pfrac,2)+20*pfrac+8)/(pow((1-pfrac),5))+(1-m_mean)*(2*pow(pfrac,3)+12*pow(pfrac,2)-48*pfrac+40)/(pow((1-pfrac)*(2-pfrac),3)))
    return -2*pi*NParticles*dI1*m2eo3/Vol-pi*m_mean*NParticles*(C1*dI2+C2*pfrac*I2)*m2e2o3/Vol
 def dI1(self):
    acorr=np.array([[0.9105631445,-0.3084016918,-0.0906148351],
                    [0.6361281449,0.1860531159,0.4527842806],
                    [2.6861347891,-2.5030047259,0.5962700728],
                    [-26.547362491,21.419793629,-1.7241829131],
                    [97.759208784,-65.255885330,-4.1302112531],
                    [-159.59154087,83.318680481,13.776631870],
                    [91.297774084,-33.746922930,-8.6728470368]])
    result=0.
    m_mean=self.m_mean()
    pfrac=self.zetan(3)
    a=np.array([0.,0.,0.,0.,0.,0.,0.])
    for i in range(7):
        a[i]=acorr[i,0]+(m_mean-1)/m_mean*acorr[i,1]+(m_mean-1)*(m_mean-2)/pow(m_mean,2)*acorr[i,2]
        result=a[i]*(i+1)*pow(pfrac,i)+result
    return result

 def dI2(self):
    bcorr=np.array([[0.7240946941,-0.5755498075,0.0976883116],
                    [2.2382791861,0.6995095521,-0.2557574982],
                    [-4.0025849485,3.8925673390,-9.1558561530],
                    [-21.003576815,-17.215471648,20.642075974],
                    [26.855641363,192.67226447,-38.804430052],
                    [206.55133841,-161.82646165,93.626774077],
                    [-355.60235612,-165.20769346,-29.666905585]])
    result=0.
    m_mean=self.m_mean()
    pfrac=self.zetan(3)
    b=np.array([0.,0.,0.,0.,0.,0.,0.])
    for i in range(7):
        b[i]=bcorr[i,0]+(m_mean-1)/m_mean*bcorr[i,1]+(m_mean-1)*(m_mean-2)/pow(m_mean,2)*bcorr[i,2]
        result=b[i]*(i+1)*pow(pfrac,i)+result
    return result
 def Pressure(self):
    Z=self.Z()
    Vol=self.Vol
    Temp=self.Temp
    NParticles=self.NParticles
    k_B=1.38064852e-23
    return Z*NParticles*k_B*Temp/Vol
# Residual chemical potential
 def mu_res(self):
     xComp = self.xComp
     a_res = self.a_res()
     Z_res = self.Z()-1.
     da_dx = self.da_res_dx()
     return a_res+Z_res+da_dx-sum(xComp*da_dx)
 def da_res_dx(self):
     da_HC_dx=self.da_HC_dx()
     da_disp_dx=self.da_disp_dx()
     return da_HC_dx+da_disp_dx
 def da_HC_dx(self):
     nsegment = self.nsegment
     xComp    = self.xComp
     m_mean   = self.m_mean()
     g_hs     = self.g_hs()
     a_hs     = self.a_hs()
     da_hs_dx = self.da_hs_dx()
     dg_hs_dx = self.dg_hs_dx()
     
     # PC-SAFT
#      A = nsegment*a_hs-(nsegment-1.)*np.log(g_hs)+m_mean*da_hs_dx-np.sum(xComp*(nsegment-1.)*pow(g_hs,-1)*dg_hs_dx,axis=1)
     
     # sPC-SAFT
     A = nsegment*a_hs-(nsegment-1.)*np.log(g_hs)+m_mean*da_hs_dx-(m_mean-1)*pow(g_hs,-1)*dg_hs_dx
     return A
 
 def da_hs_dx(self):
     zeta0_dx = self.zetan_dx(0)
     zeta1_dx = self.zetan_dx(1)
     zeta2_dx = self.zetan_dx(2)
     zeta3_dx = self.zetan_dx(3)
     
     zeta0=self.zetan(0)
     zeta1=self.zetan(1)
     zeta2=self.zetan(2)
     zeta3=self.zetan(3)

     a_hs = self.a_hs()
     Z_HS = self.Z_HS()

     t1 = -zeta0_dx/zeta0*a_hs
     
     t20 = 1/zeta0
     t21 = 3*(zeta1_dx*zeta2+zeta1*zeta2_dx)/(1-zeta3)
     t22 = 3*zeta1*zeta2*zeta3_dx*pow((1-zeta3),-2)
     t23 = 3*pow(zeta2,2)*zeta2_dx/(zeta3*pow(1-zeta3,2))
     t24 = pow(zeta2,3)*zeta3_dx*(3*zeta3-1)/(pow(zeta3,2)*pow(1-zeta3,3))
     t251 = (3*pow(zeta2,2)*zeta2_dx*zeta3-2*pow(zeta2,3)*zeta3_dx)*pow(zeta3,-3)-zeta0_dx
     t252 = np.log(1-zeta3)
     t26 = (zeta0-pow(zeta2,3)*pow(zeta3,-2))*zeta3_dx/(1-zeta3)

     # PC-SAFT
#      A = t1+t20*(t21+t22+t23+t24+t251*t252+t26)

     # sPC-SAFT
     A = Z_HS/zeta3*zeta3_dx

     return A

 def dg_hs_dx(self):
     zeta0_dx = self.zetan_dx(0)
     zeta1_dx = self.zetan_dx(1)
     zeta2_dx = self.zetan_dx(2)
     zeta3_dx = self.zetan_dx(3)
     
     zeta0=self.zetan(0)
     zeta1=self.zetan(1)
     zeta2=self.zetan(2)
     zeta3=self.zetan(3)

     dg_eta = self.dg_hs()

     # PC-SAFT  
#      HSd=self.HSd()
#      NComp = len(HSd)
#      dg_hs_dx = np.zeros((NComp,NComp))
#      for i in range(NComp):
#             t1 = zeta3_dx*pow(1-zeta3,-2)
#             t2 = HSd[i]/2*(3*zeta2_dx*pow(1-zeta3,-2)+6*zeta2*zeta3_dx*pow(1-zeta3,-3))
#             t3 = pow(HSd[i]/2,2)*(4*zeta2*zeta2_dx*pow(1-zeta3,-3)+6*pow(zeta2,2)*zeta3_dx*pow(1-zeta3,-4))
#             dg_hs_dx[:,i] = t1+t2+t3

     # sPC-SAFT
     dg_hs_dx =  dg_eta*zeta3_dx/zeta3
     return dg_hs_dx
 
 def da_disp_dx(self):
     nsegment=self.nsegment
     NParticles=self.NParticles
     Vol=self.Vol

     m_mean=self.m_mean()

     I1_dx = self.I1_dx()
     I2_dx = self.I2_dx()
     C1_dx = self.C1_dx()
     
     I1=self.I1()
     I2=self.I2()
     C1=self.C1()

     m2eo3_dx  = self.m2eno3_dx(1)  
     m2e2o3_dx = self.m2eno3_dx(2)  

     m2eo3  = self.m2eno3(1)  
     m2e2o3 = self.m2eno3(2)

     t1  = 2*(I1_dx*m2eo3+I1*m2eo3_dx) 
     t21 = (nsegment*C1*I2+m_mean*C1_dx*I2+m_mean*C1*I2_dx)*m2e2o3
     t22 = m_mean*C1*I2*m2e2o3_dx

     return -pi*NParticles/Vol*(t1+t21+t22)
 
 def m2eno3_dx(self,m):
     nsegment=self.nsegment
     xComp=self.xComp
     epsilon=self.epsilon
     sigma=self.sigma
     Temp = self.Temp
     k = self.k
     result = 0*xComp
     for j in range(len(xComp)):
        result=2*xComp[j]*nsegment*nsegment[j]*pow((np.sqrt(epsilon*epsilon[j])*(1-k[:,j]))/(Temp),m)*pow(0.5*(sigma+sigma[j]),3)+result
     return result
 
 def C1_dx(self):
     nsegment=self.nsegment
     pfrac = self.zetan(3)
     zeta3_dx=self.zetan_dx(3)
     m_mean=self.m_mean()
     C1 =self.C1()

     C2=-pow(C1,2)*(m_mean*(-4*pow(pfrac,2)+20*pfrac+8)/(pow((1-pfrac),5))+(1-m_mean)*(2*pow(pfrac,3)+12*pow(pfrac,2)-48*pfrac+40)/(pow((1-pfrac)*(2-pfrac),3)))
     
     t1 = C2*zeta3_dx
     t21 = nsegment*(8*pfrac-2*pow(pfrac,2))*pow(1-pfrac,-4)
     t22 = nsegment*(20*pfrac-27*pow(pfrac,2)+12*pow(pfrac,3)-2*pow(pfrac,4))*pow((1-pfrac)*(2-pfrac),-2)

     return t1-pow(C1,2)*(t21-t22)

 def I1_dx(self):
     acorr=np.array([[0.9105631445,-0.3084016918,-0.0906148351],
                    [0.6361281449,0.1860531159,0.4527842806],
                    [2.6861347891,-2.5030047259,0.5962700728],
                    [-26.547362491,21.419793629,-1.7241829131],
                    [97.759208784,-65.255885330,-4.1302112531],
                    [-159.59154087,83.318680481,13.776631870],
                    [91.297774084,-33.746922930,-8.6728470368]])
     result=0.
     m_mean=self.m_mean()
     pfrac=self.zetan(3)
     zeta3_dx =self.zetan_dx(3)
     nsegment = self.nsegment
     a=np.zeros((7))
     a_dx = np.zeros((7,len(nsegment)))
     for i in range(7):
        a[i]=acorr[i,0]+(m_mean-1)/m_mean*acorr[i,1]+(m_mean-1)*(m_mean-2)/pow(m_mean,2)*acorr[i,2]
        a_dx[i]=nsegment*pow(m_mean,-2)*acorr[i,1]+nsegment*pow(m_mean,-2)*(3-4*pow(m_mean,-1))*acorr[i,2]
        result=a[i]*i*pow(pfrac,i-1)*zeta3_dx+a_dx[i]*pow(pfrac,i)+result
     return result
 def I2_dx(self):
     bcorr=np.array([[0.7240946941,-0.5755498075,0.0976883116],
                    [2.2382791861,0.6995095521,-0.2557574982],
                    [-4.0025849485,3.8925673390,-9.1558561530],
                    [-21.003576815,-17.215471648,20.642075974],
                    [26.855641363,192.67226447,-38.804430052],
                    [206.55133841,-161.82646165,93.626774077],
                    [-355.60235612,-165.20769346,-29.666905585]])
     result=0.
     m_mean=self.m_mean()
     pfrac=self.zetan(3)
     zeta3_dx =self.zetan_dx(3)
     nsegment = self.nsegment
     b=np.zeros((7))
     b_dx = np.zeros((7,len(nsegment)))
     for i in range(7):
        b[i]=bcorr[i,0]+(m_mean-1)/m_mean*bcorr[i,1]+(m_mean-1)*(m_mean-2)/pow(m_mean,2)*bcorr[i,2]
        b_dx[i]=nsegment*pow(m_mean,-2)*bcorr[i,1]+nsegment*pow(m_mean,-2)*(3-4*pow(m_mean,-1))*bcorr[i,2]
        result=b[i]*i*pow(pfrac,i-1)*zeta3_dx+b_dx[i]*pow(pfrac,i)+result
     return result
 def zetan_dx(self,m):
     Vol=self.Vol
     NParticles=self.NParticles
     nsegment=self.nsegment
     HSd = self.HSd()
     return pi*NParticles*nsegment*pow(HSd,m)/Vol/6

class  GibbsMixingPCSAFT(object):
 def __init__(self, Species, Length,Temperature,Pressure,xComp):
        NSpecies = len(Species)
        epsilon = []
        sigma = []
        nsegment = []
        Mr = []
        k = np.zeros((NSpecies,NSpecies))
        for i in range(len(Species)):
               epsilon.append(EPSILON[Species[i]])
               sigma.append(SIGMA[Species[i]])
               nsegment.append(SEGMENT[Species[i]]*Length[i])
               Mr.append(Segment_mw[Species[i]]*Length[i])
               for j in range(len(Species)):
                      if Species[i]==Species[j]:
                            k[i,j]=0
                      else:
                            k[i,j]=Binary_k[Species[i]][Species[j]]          
        self.epsilon = np.array(epsilon)
        self.sigma = np.array(sigma)
        self.nsegment = np.array(nsegment)
        self.Mr = np.array(Mr)
        self.k  = k
        self.xComp=np.array(xComp)
        self.NParticles=6.02e23
        self.Temp=np.array(Temperature)
        self.Pre=np.array(Pressure)
 def Pressure(self,Vol,*arg):
    k_B=1.38064852e-23
    Pre, nsegment, Mr, epsilon, sigma, NParticles, Temp, xComp, k = arg
    r=PCSAFT(nsegment,Mr,epsilon,sigma,NParticles,Temp,pow(10,Vol),xComp,k)
    Z=r.Z()
    return (Z*NParticles*k_B*Temp/pow(10,Vol)-Pre)

 ## Gibbs Free Energy of Mixing ##
 def GibbsFreeMixing(self):
    xComp=self.xComp
    Temp=self.Temp
    NParticles=self.NParticles
    nsegment=self.nsegment
    epsilon=self.epsilon
    sigma=self.sigma
    Pre=self.Pre
    Mr=self.Mr
    k=self.k
    tideal=0.
    tres=0.
    k_B=1.38064852e-23
    for i in range(len(xComp)):
        tideal=xComp[i]*np.log(xComp[i])+tideal
        Vol1=least_squares(self.Pressure,np.log10(NParticles*nsegment[i]*pow(sigma[i],3)*pow(1-0.12*np.exp(-3*epsilon[i]/(Temp)),3))*1.05,bounds=(np.log10(sum(NParticles*nsegment[i]*pow(sigma[i],3)*pow(1-0.12*np.exp(-3*epsilon[i]/(Temp)),3)))*1.3,np.log10(sum(NParticles*nsegment[i]*pow(sigma[i],3)*pow(1-0.12*np.exp(-3*epsilon[i]/(Temp)),3)))*0.6),args=(Pre,[nsegment[i]],[Mr[i]],[epsilon[i]],[sigma[i]],NParticles,Temp,[1.],[[0.]]))
        r=PCSAFT([nsegment[i]],[Mr[i]],[epsilon[i]],[sigma[i]],NParticles,Temp,pow(10,Vol1.x[0]),[1.],[[0.]])
        Z_res=r.Z()-1.
        a_res=r.a_res()
        tres=(a_res+Z_res+np.log(Z_res+1))*xComp[i]+tres
    Vol=least_squares(self.Pressure,np.log10(sum(NParticles*xComp*nsegment*pow(sigma,3)*pow(1-0.12*np.exp(-3*epsilon/(Temp)),3)))*1.05,bounds=(np.log10(sum(NParticles*xComp*nsegment*pow(sigma,3)*pow(1-0.12*np.exp(-3*epsilon/(Temp)),3)))*1.3,np.log10(sum(NParticles*xComp*nsegment*pow(sigma,3)*pow(1-0.12*np.exp(-3*epsilon/(Temp)),3)))*0.6),args=(Pre,nsegment,Mr,epsilon,sigma,NParticles,Temp,xComp,k))
    # Vol=least_squares(self.Pressure,-2.7,bounds=(-3.2,-2.5),args=(Pre,nsegment,Mr,epsilon,sigma,NParticles,Temp,xComp,k))
    r=PCSAFT(nsegment,Mr,epsilon,sigma,NParticles,Temp,pow(10,Vol.x[0]),xComp,k)
    g_res=r.a_res()+r.Z()-1+np.log(r.Z())
    

    return g_res-tres+tideal


 def dGibbsFreeMixing(self):
    xComp=self.xComp
    Temp=self.Temp
    NParticles=self.NParticles
    nsegment=self.nsegment
    epsilon=self.epsilon
    sigma=self.sigma
    Pre=self.Pre
    Mr=self.Mr
    k=self.k
    tideal=0.
    tres=0.
    k_B=1.38064852e-23
    
    i=0
    Vol1=least_squares(self.Pressure,np.log10(NParticles*nsegment[i]*pow(sigma[i],3)*pow(1-0.12*np.exp(-3*epsilon[i]/(Temp)),3))*1.05,bounds=(np.log10(sum(NParticles*nsegment[i]*pow(sigma[i],3)*pow(1-0.12*np.exp(-3*epsilon[i]/(Temp)),3)))*1.3,np.log10(sum(NParticles*nsegment[i]*pow(sigma[i],3)*pow(1-0.12*np.exp(-3*epsilon[i]/(Temp)),3)))*0.6),args=(Pre,[nsegment[i]],[Mr[i]],[epsilon[i]],[sigma[i]],NParticles,Temp,[1.],[[0.]]))
    r1 = PCSAFT([nsegment[i]],[Mr[i]],[epsilon[i]],[sigma[i]],NParticles,Temp,pow(10,Vol1.x[0]),[1.],[[0.]])
    mu_pure1 = r1.mu_res()
    
    i=1
    Vol2=least_squares(self.Pressure,np.log10(NParticles*nsegment[i]*pow(sigma[i],3)*pow(1-0.12*np.exp(-3*epsilon[i]/(Temp)),3))*1.05,bounds=(np.log10(sum(NParticles*nsegment[i]*pow(sigma[i],3)*pow(1-0.12*np.exp(-3*epsilon[i]/(Temp)),3)))*1.3,np.log10(sum(NParticles*nsegment[i]*pow(sigma[i],3)*pow(1-0.12*np.exp(-3*epsilon[i]/(Temp)),3)))*0.6),args=(Pre,[nsegment[i]],[Mr[i]],[epsilon[i]],[sigma[i]],NParticles,Temp,[1.],[[0.]]))
    r2 = PCSAFT([nsegment[i]],[Mr[i]],[epsilon[i]],[sigma[i]],NParticles,Temp,pow(10,Vol2.x[0]),[1.],[[0.]])
    mu_pure2 = r2.mu_res()

    VolMix=least_squares(self.Pressure,np.log10(sum(NParticles*xComp*nsegment*pow(sigma,3)*pow(1-0.12*np.exp(-3*epsilon/(Temp)),3)))*1.05,bounds=(np.log10(sum(NParticles*xComp*nsegment*pow(sigma,3)*pow(1-0.12*np.exp(-3*epsilon/(Temp)),3)))*1.3,np.log10(sum(NParticles*xComp*nsegment*pow(sigma,3)*pow(1-0.12*np.exp(-3*epsilon/(Temp)),3)))*0.6),args=(Pre,nsegment,Mr,epsilon,sigma,NParticles,Temp,xComp,k))
    
    r_mix=PCSAFT(nsegment,Mr,epsilon,sigma,NParticles,Temp,pow(10,VolMix.x[0]),xComp,k)
    mu_mix = r_mix.mu_res()
    
    return np.log(xComp[0])-np.log(xComp[1])+((mu_mix[0]-mu_mix[1]))-((mu_pure1-mu_pure2))
 def g_res(self):
     xComp=self.xComp
     Temp=self.Temp
     NParticles=self.NParticles
     nsegment=self.nsegment
     epsilon=self.epsilon
     sigma=self.sigma
     Pre=self.Pre
     Mr=self.Mr
     k=self.k
     tideal=0.
     tres=0.
     k_B=1.38064852e-23
     
     VolMix=least_squares(self.Pressure,np.log10(sum(NParticles*xComp*nsegment*pow(sigma,3)*pow(1-0.12*np.exp(-3*epsilon/(Temp)),3)))*1.05,bounds=(np.log10(sum(NParticles*xComp*nsegment*pow(sigma,3)*pow(1-0.12*np.exp(-3*epsilon/(Temp)),3)))*1.3,np.log10(sum(NParticles*xComp*nsegment*pow(sigma,3)*pow(1-0.12*np.exp(-3*epsilon/(Temp)),3)))*0.6),args=(Pre,nsegment,Mr,epsilon,sigma,NParticles,Temp,xComp,k))
    
     r_mix=PCSAFT(nsegment,Mr,epsilon,sigma,NParticles,Temp,pow(10,VolMix.x[0]),xComp,k)
     mu_mix = r_mix.mu_res()
     g_mix  = r_mix.a_res()+r_mix.Z()-1+np.log(r_mix.Z())
     dg_res = np.log(xComp[0])-np.log(xComp[1])+(mu_mix[0]-mu_mix[1])
     g_res  =sum(xComp*(np.log(xComp)+g_mix))
     return [g_res,dg_res]
 
class GibbsMixingFH(object):
 def __init__(self,Species,xComp,Length):
        self.xComp=np.array(xComp)
        self.Nmono=np.array(Length)
        self.chi = CHI[Species[0]][Species[1]]
 def GibbsFreeMixing(self):
       chi = self.chi
       Nmono=self.Nmono
       xComp=self.xComp
       return xComp[0]/Nmono[0]*np.log(xComp[0])+xComp[1]/Nmono[1]*np.log(xComp[1])+chi*xComp[0]*xComp[1]
 def dGibbsFreeMixing(self):
       chi = self.chi
       Nmono=self.Nmono
       xComp=self.xComp
       return np.log(xComp[0])/Nmono[0]-np.log(xComp[1])/Nmono[1]+chi*(1-2*xComp[0])+1/Nmono[0]-1/Nmono[1]

class GibbsMixingUNIFAC(object):
 def __init__(self,Species,Length,xComp,Temp):
       self.xComp = xComp
       self.Temp = Temp[0]
       self.Species = Species
       self.Length = Length
 def GibbsFreeMixing(self):
       ln_gamma      = self.ln_gamma()
       xComp         = self.xComp
       Temp          = self.Temp
       A = 0.
       for i in range(len(xComp)):
              A+= xComp[i]*(np.log(xComp[i])+ln_gamma[i])
       return A
 def ln_gamma(self):
       Species = self.Species
       ln_gamma_comb = self.ln_gamma_comb()
       ln_gamma_fv   = self.ln_gamma_fv()
       ln_gamma_res  = self.ln_gamma_res()
       ln_gamma      = np.zeros(len(ln_gamma_comb))
       for i in range(len(Species)):
              ln_gamma[i] = ln_gamma_comb[i]+ln_gamma_fv[i]+ln_gamma_res[i]
       return ln_gamma
# Combinatorial contribution
 def ln_gamma_comb (self):
       xComp = self.xComp

       phi = self.phi()

       ln_gamma_comb = np.log(phi/xComp) + 1 - (phi/xComp)
       return ln_gamma_comb

 def  phi(self):
       xComp = self.xComp
       ri    = self.r_i()
       phi   = np.zeros((len(ri)))
       A=0.
       for i in range(len(xComp)):
              A += ((xComp[i]*(ri[i]**(0.75))))
       phi = xComp*ri**(0.75)/A
       return phi

 def r_i(self):
       Species = self.Species
       Length = self.Length
       r  = np.zeros((len(Species)))
       for i in range(len(Species)):
              for group in list(Polymer_groups[Species[i]].keys()):
                     r[i] += Polymer_groups[Species[i]][group]*R_k[group]*Length[i]
       return r

 # Free volume contribution
 def volume(self):
       Species = self.Species
       Temp    = self.Temp
       volume  = np.zeros((len(Species)))
       for i in range(len(Species)):
              if Species[i] in list(T_G.keys()):
                     if Temp <= T_G[Species[i]]:
                            coeff = expansion_coefficients[Species[i]]["LEQ_TG"]
                     else:
                            coeff = expansion_coefficients[Species[i]]["GEQ_TG"]
              else:
                     coeff = expansion_coefficients[Species[i]]
              rho = rho_reference[Species[i]] / (1 + coeff*(Temp - rho_reference_temp[Species[i]]))
              volume[i] = 1 / rho
       return volume

 def C_i_coeff (self):
       Species = self.Species
       C_1 = np.zeros((len(Species)))
       Length = self.Length
       for i in range(len(Species)):
              if Species[i] == "PMMA":
                     C_1[i] = c_i_coeff[Species[i]]*Length[i]*Segment_mw[Species[i]]
              elif Species[i] == "PS":
                     C_1[i] = c_i_coeff[Species[i]]*Length[i]*Segment_mw[Species[i]]
              elif Species[i] == "PB":
                     C_1[i] = -0.640 + ((0.6744*0.146*2.0) + (1.1167*0.304))*Length[i]
       return C_1

 def red_volume(self):
       vol = self.volume()
       ri  = self.r_i()
       red_volume = vol/(15.17*1.28*ri)
       return red_volume

 def red_volume_mix (self):
       xComp = self.xComp
       vol = self.volume()
       ri = self.r_i()
       Species = self.Species
       Length = self.Length
       A = np.zeros((len(Species)))
       wComp = np.zeros((len(Species)))
       # Calculate weight fractions from mole fractions
       for i in range(len(Species)):
              A[i] = xComp[i]*Length[i]*Segment_mw[Species[i]]
       wComp  = A / sum(A)
       red_vol_mix = (sum(vol*wComp))/(15.17*1.28*sum(ri*wComp))
       return red_vol_mix

 def ln_gamma_fv(self):
       red_vol     = self.red_volume()
       red_vol_mix = self.red_volume_mix()
       Ci          = self.C_i_coeff()
       ln_gamma_fv = np.zeros((len(Ci)))
       for i in range(len(red_vol)):
              ln_gamma_fv[i] = 3.0*Ci[i]*np.log((red_vol[i]**(1/3)- 1.0)/(red_vol_mix**(1/3) - 1.0)) - Ci[i]*((red_vol[i]/red_vol_mix)-1.0)*((1.0 - (1.0/red_vol[i]**(1/3)))**(-1.0))
       return ln_gamma_fv
 
 # Residual Contribution
 def ln_gamma_res(self):
       Species = self.Species
       lnGammaK = self.lnGammaK()
       lnGammaKi = self.lnGammaKi()
       ln_gamma_res = np.zeros((len(Species)))
       Length = self.Length
       for i in range(len(Species)):
              for group in list(Polymer_groups[Species[i]].keys()):
                     ln_gamma_res[i]+=Length[i]*Polymer_groups[Species[i]][group]*(lnGammaK[group]-lnGammaKi[Species[i]][group])
       return ln_gamma_res
 def lnGammaK(self):
       Temp =self.Temp
       H = self.Hm()
       lnGammaK = copy.deepcopy(R_k)
       for k in list(R_k):
              B = 0.
              C = 0.
              for m in list(R_k):
                     A = 0.
                     for n in list(R_k):
                            A+= H[n]*np.exp(-a_nm[n][m]/Temp)
                     B += H[m]*np.exp(-a_nm[m][k]/Temp)
                     C += H[m]*np.exp(-a_nm[k][m]/Temp)/A
              lnGammaK[k] = QK[k]*(1-np.log(B)-C)
       return lnGammaK
 def Hm(self):
       B = 0.
       A = copy.deepcopy(R_k)
       X = self.Xm()
       for m in list(R_k):
              A[m] = QK[m]*X[m]
              B   += A[m]
       H = copy.deepcopy(A)
       for m in list(A):
              H[m] = A[m]/B
       return H

 def Xm(self):
       Species = self.Species
       xComp   = self.xComp
       Length = self.Length
       C = 0.
       B = copy.deepcopy(R_k)
       for m in list(R_k):
              A = 0.
              for j in range(len(Species)):
                     if m in list(Polymer_groups[Species[j]].keys()):
                            A += xComp[j]*Polymer_groups[Species[j]][m]*Length[j]
                     else:
                            A += 0
              B[m] = A
              C   += A
       X = copy.deepcopy(B)
       for m in list(B):
              X[m] = B[m]/C
       return X

 def lnGammaKi(self):
       Temp = self.Temp
       Hi = self.Hmi()
       Species = self.Species
       lnGammaKi = copy.deepcopy(Polymer_groups)
       for i in range(len(Species)):
              for k in list(Polymer_groups[Species[i]].keys()):
                     B = 0.
                     C = 0.
                     for m in list(Polymer_groups[Species[i]].keys()):
                            A = 0.
                            for n in list(Polymer_groups[Species[i]].keys()):
                                   A+= Hi[Species[i]][n]*np.exp(-a_nm[n][m]/Temp)
                            B += Hi[Species[i]][m]*np.exp(-a_nm[m][k]/Temp)
                            C += Hi[Species[i]][m]*np.exp(-a_nm[k][m]/Temp)/A
                     lnGammaKi[Species[i]][k] = QK[k]*(1-np.log(B)-C)
       return lnGammaKi

 def Hmi(self):
       Species = self.Species
       Xi      = self.Xi()
       Hi = copy.deepcopy(Polymer_groups)
       for i in range(len(Species)):
              A = copy.deepcopy(R_k)
              B = 0.
              for m in list(Polymer_groups[Species[i]].keys()):
                     A[m] = QK[m]*Xi[Species[i]][m]
                     B    += QK[m]*Xi[Species[i]][m]
              for m in list(Polymer_groups[Species[i]].keys()):
                     Hi[Species[i]][m] = A[m]/B
       return Hi

 def Xi(self):
       Xi = copy.deepcopy(Polymer_groups)
       Species = self.Species
       Length = self.Length
       for i in range(len(Species)):
              B = 0.
              A = copy.deepcopy(R_k)
              for m in list(Polymer_groups[Species[i]].keys()):
                     A[m]  = Length[i]*Polymer_groups[Species[i]][m]
                     B    += Length[i]*Polymer_groups[Species[i]][m]
              for m in list(Polymer_groups[Species[i]].keys()):
                     Xi[Species[i]][m] = A[m]/B
       return Xi

       
Polymer_length = {
    "PMMA": 930,
    "PS":961,
    "PB":490,
    "EtOH": 1,
    "Benzene": 1
}

# Molecular weight of each segment: 
Segment_mw = {
    "PMMA": 86,
#     "PMA": 100,
    "PS": 104.1,
    "PB": 54.1,
    "EtOH": 46.07,
    "Benzene": 78.11
}

# C_{i} per unit M_{w} 
c_i_coeff = {
    "PMMA": 0.0092,
    "PS": 0.0066,
}

# Molar group volumes of all relevant groups
R_k = {
    "CH=CH":1.1167,
    "CH2":0.6744,
    "C":0.2195,
    "CH3COO":1.9031,
    "CH3":0.9011,
    "ACH":0.5313,
    "ACCH":0.8121,
    "OH": 1.0
}

# A_mn
a_nm = {
    "CH=CH":{"CH2": 2520,
             "CH=CH":0,
             "C": 2520,
             "CH3": 2520,
             "CH3COO": 71.23,
             "ACH": 340.7,
             "ACCH": 4102,
             "OH": 636.1},
    "CH2":{"CH2": 0,
             "CH=CH":-200,
             "C": 0,
             "CH3": 0,
             "CH3COO": 232.1,
             "ACH": 61.13,
             "ACCH": 76.50,
             "OH":986.5},
    "C":{"CH2": 0,
             "CH=CH":-200,
             "C": 0,
             "CH3": 0,
             "CH3COO": 232.1,
             "ACH": 61.13,
             "ACCH": 76.50,
             "OH": 636.1},
    "CH3COO":{"CH2": 114.8,
             "CH=CH":269.3,
             "C": 114.8,
             "CH3": 114.8,
             "CH3COO": 0,
             "ACH": 85.84,
             "ACCH": -170,
             "OH": 636.1},
    "CH3":{"CH2": 0,
             "CH=CH":-200,
             "C": 0,
             "CH3": 0,
             "CH3COO": 232.1,
             "ACH": 61.13,
             "ACCH": 76.50,
             "OH":986.5},
    "ACH":{"CH2": -11.12,
             "CH=CH":-94.78,
             "C": -11.12,
             "CH3": -11.12,
             "CH3COO": 5.994,
             "ACH": 0,
             "ACCH": 167.0,
             "OH": 636.1},
    "ACCH":{"CH2": -69.7,
             "CH=CH":-269.7,
             "C": -69.7,
             "CH3": -69.7,
             "CH3COO": 5688,
             "ACH": -146.8,
             "ACCH": 0,
             "OH": 636.1},
    "OH":{"CH2": 156.4,
             "CH=CH":0,
             "C": 156.4,
             "CH3": 156.4,
             "CH3COO": 0,
             "ACH": 89.6,
             "ACCH": 0,
             "OH": 0},
}

# QK
QK = {
    "CH=CH":0.867,
    "CH2":0.540,
    "C":0.0,
    "CH3COO":1.728,
    "CH3":0.848,
    "ACH":0.400,
    "ACCH":0.384,
    "OH": 1.2
}

# Number of each group in the polymer repeating unit
Polymer_groups = {
    "PB":{
        "CH=CH":1,
        "CH2":2
    },
    "PMMA":{
        "CH3":1,
        "CH2":1,
        "C":1,
        "CH3COO":1
    },
    "PS":{
        "ACH": 5,
        "ACCH":1,
        "CH2":1
    },
    "EtOH":{
        "OH": 1,
        "CH3": 1,
        "CH2": 1
    },
    "Benzene":{
        "ACH": 6,
    },
}

# Setting up volume expansion dicts
# rho_reference has units of g/ cm^3
rho_reference = {
    "PMMA":1.190,
    "PS":1.05,
    "PB":0.915
}

#units of Celcius
rho_reference_temp = {
    "PMMA":20.0,
    "PS":20.0,
    "PB":25.0
}

# Glass transition temperatures
T_G = {
    "PMMA": 105,
    "PS": 80
}

# volume thermal expansion coefficients: 
expansion_coefficients = {
    "PMMA":{
        "LEQ_TG":2.60e-4,
        "GEQ_TG":5.60e-4
    },
    "PS":{
        "LEQ_TG":1.7e-4,
        "GEQ_TG":5.1e-4
    },
    "PB":6.99e-4
}

EPSILON = { 
       "PMMA": 264.6,
       "PS": 348.2,
       "PB": 288.84
}

SIGMA = { 
       "PMMA": 3.553e-10,
       "PS": 4.152e-10,
       "PB": 4.097e-10
}

SEGMENT = { 
       "PMMA": 0.0270*100,
       "PS": 0.0205*104.1,
       "PB": 0.0245*54.1
}

Binary_k = {
"PB": {"PS":0.00539},
"PS": {"PB":0.00539}
}

CHI = {
       "PMMA": {"PS": 0.315},
       "PS": {"PMMA": 0.315}
}