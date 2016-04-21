# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 21:24:44 2016

@author: KILLER
"""

import scipy as sc
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import pylab

Mh=1
Mc=2
Thin=373.15 #K
Tcin=303.15#K
U=300
A=100
n=10
Thguess=n*[Thin]
Tcguess=n*[Tcin]
Tguess=sc.array(Thguess+Tcguess)

### heat balance

Ht1=238 ## enthalpy at t1 , Btu/ lb 
Ht2=252 ## enthalpy at t2 , Btu/ lb 
Ht3=378 ## enthalpy of vapour at t2


def Cph(T):
    Cp=4184+10**(-4)*T+10**(-6)*T**2+10**(-9)*T**3
    return Cp
def Cpc(T):
    Cp=4184+10**(-4)*T+10**(-6)*T**2+10**(-9)*T**3
    return Cp    
#the array whose elements you want to modify is the first argument to the function
def residuals(T,U,A,Thin,Tcin,Mh,Mc):
    n=len(T)
    Th=T[:n/2]
    Tc=T[n/2:]
    dA=A/((n/2)-1)
    errHL=U*(Thin-Tc[0])/(Mh*Cph(Thin))+(Th[1]-Thin)/dA
    errCL=U*(Thin-Tc[0])/(Mc*Cpc(Tc[0]))+(Tc[1]-Tc[0])/dA
    errHR=U*(Th[-1]-Tcin)/(Mh*Cph(Th[-1]))+(Th[-1]-Th[-2])/dA
    errCR=U*(Tc[-1]-Tcin)/(Mh*Cpc(Tcin))+(Tcin-Tc[-2])/dA
    errH=sc.zeros(n/2)
    errC=sc.zeros(n/2)
    errH[0]=errHL;errH[-1]=errHR
    errC[0]=errCL;errC[-1]=errCR
    errH[1:-1]=U*(Th[1:-1]-Tc[1:-1])/(Mh*Cph(Th[1:-1]))+(Th[2:]-Th[1:-1])/dA
    #errH[1:-1]=U*(Th[1:-1]-Tc[1:-1])/(Mh*Cph(Th[1:-1]))+(Th[2:]-Th[0:-2])/dA for central difference
    errC[1:-1]=U*(Th[1:-1]-Tc[1:-1])/(Mc*Cpc(Tc[1:-1]))+(Tc[2:]-Tc[1:-1])/dA
    return sc.concatenate((errH,errC))

soln=opt.leastsq(residuals,Tguess,args=(U,A,Thin,Tcin,Mh,Mc))
Tsoln=soln[0]
Thsoln=Tsoln[:20/2];Thsoln[0]=Thin
Tcsoln=Tsoln[20/2:];Tcsoln[-1]=Tcin
print Thsoln    
print Tcsoln
pylab.plot(range(0,100,10),Thsoln)
pylab.plot(range(0,100,10),Tcsoln)
pylab.show()

## for preheat
qv=( Mh *( Ht3 -Ht2))
qs=Mc*( Ht2 - Ht1)
Q=qs+qv ## total heat required for naphtha

Nt=float(raw_input("no.of tube on hot fluid side="))
n=float(raw_input("no.of passes on hot fluid side="))
L=float(raw_input("length of tube="))   ##ft
at1=0.546## flow area,in^2
at=((Nt*at1)/(144.0*n))  ##  total area,ft^2
print ('total area=',at)
Gt=(Mc/(at))## mass velocity,lb/(hr)*(ft^2)
mu1 =1.09## at 451F,lb/(f t)*(hr)
D=0.0695## f t
Ret=(( D)*( Gt)/ mu1 )## reynolds number
print ('Reynolds no. =',Ret)
jH=168
Z =0.142## // Z=k ((c)(mu1)/k )^(1/3 )
Hi=(( jH) *(1/ D)*(Z))## individual heat transfer coefficient 
Hio =(( Hi) *(0.834/1))##  Correct Hio to the surface at the OD is
ho1 =300  ## assumption
tw =(Tcsoln[-1])+(((Hio)/(Hio+ho1))*(Thsoln[1] -Tcsoln[-1]))
deltw =(tw -Tcsoln[-1])
Av=(qv/300)
As=qs/60
A1=As+Av
ho=(Q/A1)
Uc=(( Hio)*(ho)/(Hio +ho))## clean overall coefficient
A2=0.2618## actual surface supplied for each tube
A=( Nt*L*A2)##ft^2/// total surface area
UD=((Q)/((A)*(tw)))## actual design  overall coeffcient
print (' actual design  overall coeffcient',UD)
