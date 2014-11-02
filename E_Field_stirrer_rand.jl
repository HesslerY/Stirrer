include("ImageCreator.jl")
include("dipole.jl")

const c = 299792458.
const mu0 = 4*pi*1e-7
const eps0 = 1/(mu0*c^2)

f0=100e6
f1=3e9
nf=50
freq=10.^(linspace(log10(f0),log10(f1),nf))

#cavity
const L=8.7
const l=3.7
const h=2.9
const losses=0.998*ones(nf)  #losses can be frequency dependent...
#crit=0.01 #convergence criterion
#order=int(log(crit)/log(maximum(losses)))
const order=20

##dipole
#dipole position
const X=1
const Y=1
const Z=1
const tilt=pi/2-acos(sqrt(2/3));
const azimut=pi/4
const phase=0
#dipole moment
#total time averaged radiated power P= 1 W dipole moment => |p|=sqrt(12πcP/µOω⁴)
const Pow=1
const amplitude=sqrt(12*pi*c*Pow./(mu0*(2*pi*freq).^4))

#Images
POS=IC(L,l,h,X,Y,Z,tilt,azimut,phase,1,order)
Lm=(order+2)*L
lm=(order+2)*l
hm=(order+2)*h

dmax=maximum([Lm,lm,hm])
dist=sqrt(POS[:,4].^2+POS[:,5].^2+POS[:,6].^2)
U=find(dist.<dmax)
POS=POS[U,:]



const numberofimages=length(POS[:,1])
const M=30
const na=180
const phi=linspace(2*pi/na,2*pi,na)

Et=zeros(Complex64,M,nf,na,3)
Bt=zeros(Complex64,M,nf,na,3)

using NPZ


for i=1:M
  #random observation circle
  radius=1
  th=acos(2*rand()-1)
  ph=rand()*2*pi
  Cx=rand()*(L-3*radius)+1.5*radius
  Cy=rand()*(l-3*radius)+1.5*radius
  Cz=rand()*(h-3*radius)+1.5*radius
  x=radius*cos(phi)*(-sin(ph))+radius*sin(phi)*cos(th)*cos(ph)+Cx
  y=radius*cos(phi)*cos(ph)+radius*sin(phi)*cos(th)*sin(ph)+Cy
  z=radius*sin(phi)*(-sin(th))
  for f=1:nf
    println("$i/$M, $f/$nf")
    for n=1:na
      perc=round(n/na*1000)/10
      r=[x[n],y[n],z[n]]
      E=zeros(Complex128,3)
      B=zeros(Complex128,3)
      for m=1:numberofimages
        p=vec(POS[m,1:3])*amplitude[f]*losses[f]^POS[m,7]
        R=vec(POS[m,4:6])
        ord=POS[m,7]
        Ed,Bd=Hertz_dipole_ff (r, p, R, phase, freq[f])
        E+=Ed
        B+=Bd
      end
      Et[i,f,n,:]=E
      Bt[i,f,n,:]=B
    end
  end
  npzwrite("Estirrer_part.npz", ["Et" => Et[i,:,:,:],"Bt" => Bt[i,:,:,:], "freq" => freq, "x" => x, "y" => y, "z" => z])
end
#Back up
npzwrite("Estirrer_MC.npz", ["Et" => Et,"Bt" => Bt, "freq" => freq, "x" => x, "y" => y, "z" => z])

