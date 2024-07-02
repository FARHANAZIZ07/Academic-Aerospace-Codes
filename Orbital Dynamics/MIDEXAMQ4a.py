# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 17:53:32 2022

@author: dell
"""
import numpy as np
Energy=(((0.7156)**2)/2)-(1/2)
print('Energy=',Energy)
a=-1/(2*(Energy))
print('a=',a)
r=np.array([0,2,0])
v=np.array([0.6422,0,0.3708])
k=np.array([0,0,1])
h=np.cross(r,v)
#print('h=',h)
z=np.cross(v,h)
ev=(z-1/2*(r))/1
e=np.linalg.norm(ev)
print('e=',e)
y=np.linalg.norm(np.dot(k,h))
print(y)
i=np.arccos((np.dot(k,h))/np.linalg.norm(h))
print('i=',i)
n=np.cross(k,h)/np.linalg.norm(np.cross(k,h))
print(n)
Bomg=np.arctan2(n[1],n[0])
print('Bomg=',Bomg)
Somg=np.arccos((np.linalg.norm(np.dot(n,ev)))/(e))
print('Somg=',Somg)
mu=np.arccos(np.linalg.norm(np.dot(r,ev))/(np.linalg.norm(r)*e))
print('mu=',mu)

T=2*np.pi*((a)**3)**0.5
print('Operiod=',T)

#at half the period

MU=180
A=a
P=A*(1-e**2)
R=P/(1+(e*np.cos(MU)))
print('P=',P)
print('R=',R)
pp=2

print(np.cos(MU))
#position and velocity
RR=R*np.cos(MU)+R*np.sin(MU)
V=((1/pp)**0.5)*(-np.sin(MU)+(e+np.cos(MU)))
print('RR=',RR)
print('V=',V)
