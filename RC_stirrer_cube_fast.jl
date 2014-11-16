include("ImageCreator.jl")
#include("dipole.jl")

using NPZ
#using PyPlot
#figure()
const c = 299792458.
const mu0 = 4*pi*1e-7
const eps0 = 1/(mu0*c^2)
epsr=1
#simulation parameters
const Lt=.5e-6
const freq=2e9

L=8.7
l=3.7
h=2.9
losses=0.99
dmax=Lt*c
order=int(round(dmax/minimum([L;l;h]))+1)
const w=2*pi*freq  # omega
const k=w/c     # wave number
const N=int(5*Lt*freq)
const t=linspace(Lt/N,Lt,N)


nphi=360
phi=linspace(2*pi/nphi,2*pi,nphi)

n_paddles=5
nx=3
ny=3
nz=3
x0=2
y0=1
z0=1.5

x=linspace(0,1,nx)+x0
y=linspace(0,1,ny)+y0
z=linspace(0,1,nz)+z0

Et=zeros(nx,ny,nz,nphi,N,3)
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
  for kp=1:n_paddles
    #println("$i/$nphi")
    xx=X[kp]+radius*cos(phi[i]+az[kp])
    yy=Y[kp]+radius*sin(phi[i]+az[kp])
    zz=Z[kp]
    POSp=IC(L,l,h,xx,yy,zz,tilt[kp],azimut[kp]+phi[i],phase[kp],1,order)
    Lm=(order+2)*L
    lm=(order+2)*l
    hm=(order+2)*h
    dmax=maximum([Lm,lm,hm])
    dist=sqrt(POSp[:,4].^2+POSp[:,5].^2+POSp[:,6].^2)
    U=find(dist.<dmax)
    if kp==1
      POS=POSp[U,:]
    else
      POS=vcat(POS,POSp[U,:])
    end
  end
  numberofimages=length(POS[:,1])
  #field computation
  for xxx=1:nx
    for yyy=1:ny
      for zzz=1:nz
        r=[x[xxx],y[yyy],z[zzz]] #position of the receiver
        E=zeros(N,3)
        #B=zeros(N,3)
        for m=1:numberofimages
          ord=POS[m,7]#order of the dipole
          p=vec(POS[m,1:3])*losses^ord #image dipole moment
          R=vec(POS[m,4:6]) #image dipole position
          rprime=r-R  # r'=r-R
          magrprime=sqrt(sum(rprime.^2)) # |r-R|
          #krp=k*magrprime  # k*|r-R|
          rprime_cross_p = cross(rprime, p) # (r-R) x p
          rp_c_p_c_rp = cross(rprime_cross_p, rprime) # ((r-R) x p) x (r-R)
          ta=int(magrprime/c/Lt*N)
          if ta == 0
            ta=1
          end
          if ta<N+1
            #expfac=exp(1im*(-w*t+krp-phi))
            Ep=1/(4*pi*eps0*epsr)*(1/(c^2*magrprime^3)*rp_c_p_c_rp) #+(1/magrprime^3-w*im/(c*magrprime^2))*(3*rprime*dot(rprime,p)/magrprime^2-p))
            E[ta,:]+=Ep'
            #Bp=1/(4*pi*eps0*epsr)*1/(magrprime^2*c^3)*(w^2*rprime_cross_p)
            #B[:,ta]+=Bp
          end
        end
        Et[xxx,yyy,zzz,i,:,:]=E
        println("$xxx $yyy $zzz")
      end
    end
  end
end
npzwrite("Estirrerfast_cube$n_paddles.npz", ["Et" => Et, "t" => t])






