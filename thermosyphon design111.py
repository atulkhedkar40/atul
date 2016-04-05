# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 21:24:44 2016

@author: KILLER
"""

import scipy as sc
import matplotlib.pyplot as plt
import scipy.optimize as opt
import pylab
Mh=38500
Mc=51000
Thin=525 #inlet cold fluid(K)
Tcin=315#inlet gasoil(K)
U=300
A=100
n=20
Thguess=n*[Thin]
Tcguess=n*[Tcin]
Tguess=sc.array(Thguess+Tcguess)

### heat balance

Ht1=238 ## enthalpy at t1 , Btu/ lb 
Ht2=252  ## enthalpy at t2 , Btu/ lb 
Ht3=378 ## enthalpy of vapour at t2

def Cph(T):
    Cp=5100+10**(-4)*T+10**(-6)*T**2+10**(-12)*T**6
    return Cp
def Cpc(T):
    Cp=3600+10**(-4)*T+10**(-26)*T**2+10**(-12)*T**7
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
##pylab.plot(range(0,100,20),Thsoln)
##pylab.plot(range(0,100,20),Tcsoln)
pylab.show() 

## for hot side fluid
Nt=raw_input("no.of tube on hot fluid side=")
n=raw_input("no.of passes on hot fluid side=")
L=12##ft
at1=0.546## flow area,in^2
at=((Nt*at1)/(144*n))##  total area,ft^2
print at
Gt=(Mc/(at))## mass velocity,lb/(hr)*(ft^2)
mu1 =1.09## at 451F,lb/(f t)*(hr)
D=0.0695## f t
Ret=(( D)*( Gt)/ mu1 )## reynolds number
print (Ret)