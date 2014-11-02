include("ImageCreator.jl")
include("dipole.jl")

using NPZ
using PyPlot
figure()
const c = 299792458.
const mu0 = 4*pi*1e-7
const eps0 = 1/(mu0*c^2)

f0=100e6
f1=3e9
nf=30
freq=linspace(f0,f1,nf)#10.^(linspace(log10(f0),log10(f1),nf))
nphi=360
phi=linspace(2*pi/nphi,2*pi,nphi)

n_paddles=10

L=8.7
l=3.7
h=2.9
order=20
losses=0.99*ones(nf)

MC=1
#x=rand(MC)*L-1+1
#y=rand(MC)*l-1+1
#z=rand(MC)*h-1+1

x=npzread("posxyz.npz")["x"][1:MC]
y=npzread("posxyz.npz")["y"][1:MC]
z=npzread("posxyz.npz")["z"][1:MC]

Et=zeros(Complex64,MC,nphi,nf,3)
#Bt=zeros(Complex64,MC,nphi,nf,3)

#paddle positions
#axis:
X=7*ones(n_paddles)
Y=1.5*ones(n_paddles)
#height
Z=rand(n_paddles)*3
az=rand(n_paddles)*2*pi

#dipole at phi=0
radius = 1
tilt=acos(2*rand(n_paddles)-1)#[pi/2-acos(sqrt(2/3)),pi/6,pi/3]
azimut=rand(n_paddles)*2*pi
phase=rand(n_paddles)*2*pi
#dipole moment
#total time averaged radiated power P= 1 W dipole moment => |p|=sqrt(12πcP/µOω⁴)
Pow=1
amplitude=sqrt(12*pi*c*Pow./(mu0*(2*pi*freq).^4))/n_paddles
for i=1:nphi
  POS=[]
  for k=1:n_paddles
    #println("$i/$nphi")
    xx=X[k]+radius*cos(phi[i]+az[k])
    yy=Y[k]+radius*sin(phi[i]+az[k])
    zz=Z[k]
    POSp=IC(L,l,h,xx,yy,zz,tilt[k],azimut[k]+phi[i],phase[k],1,order)
    Lm=(order+2)*L
    lm=(order+2)*l
    hm=(order+2)*h
    dmax=maximum([Lm,lm,hm])
    dist=sqrt(POSp[:,4].^2+POSp[:,5].^2+POSp[:,6].^2)
    U=find(dist.<dmax)
    if k==1
      POS=POSp[U,:]
    else
      POS=vcat(POS,POSp[U,:])
    end
  end
  numberofimages=length(POS[:,1])
  #field computation
  for mc=1:MC
    for f=1:nf
      r=[x[mc],y[mc],z[mc]]
      E=zeros(Complex128,3)
      #B=zeros(Complex128,3)
      for m=1:numberofimages
        p=vec(POS[m,1:3])*amplitude[f]*losses[f]^POS[m,7]
        R=vec(POS[m,4:6])
        ord=POS[m,7]
        Ed,Bd=Hertz_dipole_ff (r, p, R, phase[k], freq[f])
        E+=Ed
        #B+=Bd
      end
      Et[mc,i,f,:]=E
      #Bt[mc,i,f,:]=B
    end
    println("$mc/$MC")
  end
  println("$i/$nphi")
  clf()
  suptitle("$i/$nphi")
  subplot(311)
  plot(phi,real(squeeze(Et[1,:,:,1],1)))
  ylabel("V/m")
  subplot(312)
  plot(phi,real(squeeze(Et[1,:,:,2],1)))
  ylabel("V/m")
  subplot(313)
  plot(phi,real(squeeze(Et[1,:,:,3],1)))
  ylabel("V/m")
  xlabel("\$\\phi\$")
end
npzwrite("Estirrer_MC$n_paddles.npz", ["Et" => Et, "freq" => freq])


cor(real(squeeze(Et[1,:,:,3],1)))




