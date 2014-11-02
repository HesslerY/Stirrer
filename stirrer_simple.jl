#include("ImageCreator.jl")
include("dipole.jl")

const c = 299792458.
const mu0 = 4*pi*1e-7
const eps0 = 1/(mu0*c^2)

f0=100e6
f1=1e9
nf=20
freq=10.^(linspace(log10(f0),log10(f1),nf))

n_paddles=100
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

nphi=360
phi=linspace(2*pi/nphi,2*pi,nphi)

Et=zeros(Complex64,nphi,nf,3)
Bt=zeros(Complex64,nphi,nf,3)
using NPZ

for i=1:nphi
  x=X+radius*cos(phi[i]+az)
  y=Y+radius*sin(phi[i]+az)
  z=Z.*ones(length(phi[i]))
  i0=[sin(tilt).*cos(azimut+phi[i]) sin(tilt).*sin(azimut+phi[i]) cos(tilt)]
  POS=[i0 x y z zeros(length(X)) phase]
  numberofimages=length(POS[:,1])

  #field computation
  x=1
  y=1
  z=1
  for f=1:nf
    r=[x,y,z]
    E=zeros(Complex128,3)
    B=zeros(Complex128,3)
    for m=1:numberofimages
      p=vec(POS[m,1:3])*amplitude[f]#*losses[f]^POS[m,7]
      R=vec(POS[m,4:6])
      ord=POS[m,7]
      Ed,Bd=Hertz_dipole_ff (r, p, R, POS[m,8], freq[f])
      E+=Ed
      B+=Bd
    end
    Et[i,f,:]=E
    Bt[i,f,:]=B
  end
end
#Back up
npzwrite("Estirrer_MC.npz", ["Et" => Et,"Bt" => Bt, "freq" => freq])

using PyPlot
plot(phi,real(Et[:,:,2]))



figure()
plt.hist(abs(Et[:,nf,2]),30)




#faire une stat de MC en déduire l'entropie et associer une valeur d'ordre de Tchebyschef
