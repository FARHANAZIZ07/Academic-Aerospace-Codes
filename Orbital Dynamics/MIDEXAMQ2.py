# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 12:31:07 2022

@author: dell
"""

import numpy as np
a=6
e=0.75
nu=np.pi/2
print(nu)
mu=1
P=2*np.pi*(((a**3)/mu)**0.5)
print('P=',P)
T=((1-e)/(1+e))**0.5
x=np.tan(nu/2)
E=2*np.arctan(T*x)
print('E=',E)
M=E-e*np.sin(E)
print('M=',M)
n=M/2*np.pi
print('n=',n)
TOF=n*P
print('TOF=',TOF)