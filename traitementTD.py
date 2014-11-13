# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 14:05:20 2014

@author: e68972
"""
from __future__ import division
from numpy import *
from pylab import *


Et=load("Estirrerfast_MC5.npz")["Et"]
freq=1e9
Lt=.1e-6
N=int(3*Lt*freq)
t=linspace(Lt/N,Lt,N)

plot(t,Et[0,0,:,2])

def nextpow2(v):
    v -= 1
    v |= v >> 1
    v |= v >> 2
    v |= v >> 4
    v |= v >> 8
    v |= v >> 16
    return v + 1

Fs = N/Lt #sampling frequency
periode = 1/Fs                     # sample length
L = N                  #Number of points
tt = array(arange(0,L-1,1))*periode                # time...
NFFT = 2^nextpow2(N) # Next power of 2 from length of y

mc=len(Et[:,0,0,0])
nphi=len(Et[0,:,0,0])
nf=len(f)

FFTEt=zeros((mc,nphi,NFFT/2-1,3),'complex')

for i in range(mc):
    for j in range(nphi):  
        Yx = fft(Et[i,j,:,0],NFFT)
        Yy = fft(Et[i,j,:,1],NFFT)
        Yz = fft(Et[i,j,:,2],NFFT)
        f = Fs/2*linspace(0,1,NFFT/2)
        FFTEt[i,j,:,0] = Yx[1:NFFT/2]
        FFTEt[i,j,:,1] = Yy[1:NFFT/2]
        FFTEt[i,j,:,2] = Yz[1:NFFT/2]


df=50
nx=zeros((mc,nf-df))
ny=zeros((mc,nf-df))
nz=zeros((mc,nf-df))
for i in range(mc):
    for l in range(nf-df):
        ax=corrcoef(abs(FFTEt[i,:,l:l+df,0]))
        ay=corrcoef(abs(FFTEt[i,:,l:l+df,1]))
        az=corrcoef(abs(FFTEt[i,:,l:l+df,2]))
        nx[i,l]=nphi**2/sum(abs(ax)**2)
        ny[i,l]=nphi**2/sum(abs(ay)**2)
        nz[i,l]=nphi**2/sum(abs(az)**2)

plot(f[:-df]/1e6,nx.mean(axis=0))
plot(f[:-df]/1e6,ny.mean(axis=0))
plot(f[:-df]/1e6,nz.mean(axis=0))
grid()



df=100
nfx=zeros((mc,nf-df))
nfy=zeros((mc,nf-df))
nfz=zeros((mc,nf-df))
for i in range(mc):
    for l in range(nf-df):
        ax=corrcoef(abs(FFTEt[i,:,l:l+df,0].T))
        ay=corrcoef(abs(FFTEt[i,:,l:l+df,1].T))
        az=corrcoef(abs(FFTEt[i,:,l:l+df,2].T))
        nfx[i,l]=df**2/sum(abs(ax)**2)
        nfy[i,l]=df**2/sum(abs(ay)**2)
        nfz[i,l]=df**2/sum(abs(az)**2)

plot(f[:-df]/1e6,nfx.mean(axis=0))
plot(f[:-df]/1e6,nfy.mean(axis=0))
plot(f[:-df]/1e6,nfz.mean(axis=0))
grid()



df=100
npx=zeros((nf))
npy=zeros((nf))
npz=zeros((nf))

for l in range(nf):
    ax=corrcoef(abs(FFTEt[:,:,l,0]))
    ay=corrcoef(abs(FFTEt[:,:,l,1]))
    az=corrcoef(abs(FFTEt[:,:,l,2]))
    npx[l]=mc**2/sum(abs(ax)**2)
    npy[l]=mc**2/sum(abs(ay)**2)
    npz[l]=mc**2/sum(abs(az)**2)

plot(f/1e6,npx)
plot(f/1e6,npy)
plot(f/1e6,npz)
grid()

