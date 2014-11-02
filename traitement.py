# -*- coding: utf-8 -*-
"""
Created on Sat Nov  1 11:02:40 2014

@author: manu
"""

from numpy import *



res=load("Estirrer_MC10.npz")
Et=res["Et"]
freq=res["freq"]

step=10000
rangef=range(0,len(freq)-step,step)
Nx=zeros((len(Et[:,0,0,0]),len(rangef)))
Ny=zeros((len(Et[:,0,0,0]),len(rangef)))
Nz=zeros((len(Et[:,0,0,0]),len(rangef)))



for i in range(len(Et[:,0,0,0])):
    for j in range(len(rangef)):
      Nx[i,j]=360**2/sum(abs(corrcoef(real(Et[i,:,rangef[j]:rangef[j]+step,0])))**2)
      Ny[i,j]=360**2/sum(abs(corrcoef(real(Et[i,:,rangef[j]:rangef[j]+step,1])))**2)
      Nz[i,j]=360**2/sum(abs(corrcoef(real(Et[i,:,rangef[j]:rangef[j]+step,2])))**2)
    
    