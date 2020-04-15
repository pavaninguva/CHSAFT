import numpy as np 
import pandas as pd 
from math import pi
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, least_squares

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
        # for i in range(len(xComp)):
        #     B=xComp[i]*(nsegment[i]-1)*np.log(g_hs[i])
        #     result=result+B
        # A=m_mean*a_hs-result
        
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
        # A=1/zeta0*((3*zeta1*zeta2)/(1-zeta3)+pow(zeta2,3)/(zeta3*pow((1-zeta3),2))+(pow(zeta2,3)/pow(zeta3,2)-zeta0)*np.log(1-zeta3))
        
        # sPC-SAFT
        A=(4*zeta3-3*pow(zeta3,2))/pow(1-zeta3,2)
        return A

 def g_hs(self):
        zeta2=self.zetan(2)
        zeta3=self.zetan(3)
        HSd=self.HSd()
        
        # PC-SAFT
        # A=HSd
        # for i in range(len(HSd)):
        #     A[i]=1/(1-zeta3)+pow(HSd[i],2)*3*zeta2/(2*HSd[i]*pow((1-zeta3),2))+pow(HSd[i],4)*2*pow(zeta2,2)/((pow(2*HSd[i],2))*pow((1-zeta3),3))
        
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
    # t2=0.
    # for i in range(len(xComp)):
    #     t2=xComp[i]*(nsegment[i]-1)*pow(g_hs[i],-1)*dg_hs[i]+t2
    # A = t1-t2
    
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
    # A=zeta3/(1-zeta3)+3*zeta1*zeta2/(zeta0*pow((1-zeta3),2))+(3*pow(zeta2,3)-zeta3*pow(zeta2,3))/(zeta0*pow((1-zeta3),3))
    return A

 def dg_hs(self):
    zeta2=self.zetan(2)
    zeta3=self.zetan(3)
    HSd=self.HSd()
    
    # PC-SAFT
    # A=HSd*0.
    # for i in range(len(HSd)):
    #     A[i]=zeta3/pow((1-zeta3),2)+HSd[i]/2*(3*zeta2/pow((1-zeta3),2)+6*zeta2*zeta3/pow((1-zeta3),3))+pow(HSd[i]/2,2)*(4*pow(zeta2,2)/pow((1-zeta3),3)+6*pow(zeta2,2)*zeta3/pow((1-zeta3),4))
    
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

class  GibbsMixingPCSAFT(object):
 def __init__(self, nsegment, Mr, epsilon,sigma,NParticles,Temperature,Pressure,xComp,k):
        self.epsilon = np.array(epsilon)
        self.sigma = np.array(sigma)
        self.k=np.array(k)
        self.nsegment = np.array(nsegment)
        self.Mr = np.array(Mr)
        self.xComp=np.array(xComp)
        self.NParticles=np.array(NParticles)
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
        tres=(a_res+Z_res)*xComp[i]+tres
    Vol=least_squares(self.Pressure,np.log10(sum(NParticles*xComp*nsegment*pow(sigma,3)*pow(1-0.12*np.exp(-3*epsilon/(Temp)),3)))*1.05,bounds=(np.log10(sum(NParticles*xComp*nsegment*pow(sigma,3)*pow(1-0.12*np.exp(-3*epsilon/(Temp)),3)))*1.3,np.log10(sum(NParticles*xComp*nsegment*pow(sigma,3)*pow(1-0.12*np.exp(-3*epsilon/(Temp)),3)))*0.6),args=(Pre,nsegment,Mr,epsilon,sigma,NParticles,Temp,xComp,k))
    # Vol=least_squares(self.Pressure,-2.7,bounds=(-3.2,-2.5),args=(Pre,nsegment,Mr,epsilon,sigma,NParticles,Temp,xComp,k))
    r=PCSAFT(nsegment,Mr,epsilon,sigma,NParticles,Temp,pow(10,Vol.x[0]),xComp,k)
    a_res=r.a_res()
    Z_res=r.Z()-1.
    print(10**Vol1.x[0],10**Vol.x[0])
    
    gMix=(a_res+Z_res)

    return (tideal+gMix-tres)*NParticles*Temp*k_B

class GibbsMixingFH(object):
 def __init__(self,Temp,xComp,Nmono,Vmono,A,B,C):
        self.xComp=np.array(xComp)
        self.Nmono=np.array(Nmono)
        self.Vmono=np.array(Vmono)
        self.A=np.array(A)
        self.B=np.array(B)
        self.C=np.array(C)
        self.Temp=np.array(Temp)
 def GibbsFreeMixing(self):
       A=self.A
       B=self.B
       C=self.C
       Temp=self.Temp
       Nmono=self.Nmono
       Vmono=self.Vmono
       xComp=self.xComp
       N_A=6.02e23
       chi=(A+B/Temp+C/Temp**2)
       vComp=np.zeros(len(xComp))
       vComp[0]=Nmono[0]*Vmono[0]*xComp[0]/(sum(Nmono*Vmono*xComp))
       vComp[1]=Nmono[1]*Vmono[1]*xComp[1]/(sum(Nmono*Vmono*xComp))
       Vol=Nmono*Vmono
       t1=vComp[0]/Vol[0]*np.log(vComp[0])
       t2=(1-vComp[0])/Vol[1]*np.log(1-vComp[0])
       t3=chi*vComp[0]*(1-vComp[0])
       return t1+t2+t3

 def dGibbsFreeMixing(self):
       A=self.A
       B=self.B
       C=self.C
       Temp=self.Temp
       Nmono=self.Nmono
       Vmono=self.Vmono
       xComp=self.xComp
       N_A=6.02e23
       chi=(A+B/Temp+C/Temp**2)
       vComp=np.zeros(len(xComp))
       vComp[0]=Nmono[0]*Vmono[0]*xComp[0]/(sum(Nmono*Vmono*xComp))
       vComp[1]=Nmono[1]*Vmono[1]*xComp[1]/(sum(Nmono*Vmono*xComp))
       Vol=Nmono*Vmono
       return np.log(vComp[0])/Vol[0]-np.log(1-vComp[0])/Vol[1]-2*chi*vComp[0]+chi-pow(Vol[1],-1)+pow(Vol[0],-1)