const c = 299792458.
const mu0 = 4*pi*1e-7
const eps0 = 1/(mu0*c^2)

function Hertz_dipole (r, p, R, phi, f, t=0, epsr=1.)
    #=
    Calculate E field and B field strength of hertzian dipole(s)
    p: array of dipole moments
    R: array of dipole positions
    r: observation point
    f: frequency
    t: time
    phi: array with dipole phase angles (0..2pi)
    return: field strength at observation point r at time t (3-tuple: Ex, Ey, Ez)
    =#
    rprime=r-R  # r'=r-R
    magrprime=sqrt(sum(rprime.^2)) # |r-R|
    w=2*pi*f  # omega
    k=w/c     # wave number
    krp=k*magrprime  # k*|r-R|
    rprime_cross_p = cross(rprime, p) # (r-R) x p
    rp_c_p_c_rp = cross(rprime_cross_p, rprime) # ((r-R) x p) x (r-R)
    expfac=exp(1im*(-w*t+krp-phi))/(4*pi*eps0*epsr)
    E=expfac*(w^2/(c^2*magrprime^3)*rp_c_p_c_rp+(1/magrprime^3-w*im/(c*magrprime^2))*(3*rprime*dot(rprime,p)/magrprime^2-p))
	  B=expfac/(magrprime*c^3)*(w^2*rprime_cross_p)/magrprime*(1-c/(im*w*magrprime))
    return E,B
  end


function Hertz_dipole_ff (r, p, R, phi, f, t=0, epsr=1.)
  #=
  Calculate E field and B field strength in far field of hertzian dipole(s)
  p: array of dipole moments
  R: array of dipole positions
  r: observation point
  f: frequency
  t: time
  phi: array with dipole phase angles (0..2pi)
  return: field strength at observation point r at time t (3-tuple: Ex, Ey, Ez)
  =#
  rprime=r-R  # r'=r-R
  magrprime=sqrt(sum(rprime.^2)) # |r-R|
  w=2*pi*f  # omega
  k=w/c     # wave number
  krp=k*magrprime  # k*|r-R|
  rprime_cross_p = cross(rprime, p) # (r-R) x p
  rp_c_p_c_rp = cross(rprime_cross_p, rprime) # ((r-R) x p) x (r-R)
  expfac=exp(1im*(-w*t+krp-phi))/(4*pi*eps0*epsr)
  E=expfac*(w^2/(c^2*magrprime^3)*rp_c_p_c_rp)
  B=expfac/(magrprime^2*c^3)*(w^2*rprime_cross_p)
  return E,B
end


function Hertz_dipole_nf (r, p, R, phi, f, t=0, epsr=1.)
    #=
  Calculate E field and B field strength in near field of hertzian dipole(s)
  p: array of dipole moments
  R: array of dipole positions
  r: observation point
  f: frequency
  t: time
  phi: array with dipole phase angles (0..2pi)
  return: field strength at observation point r at time t (3-tuple: Ex, Ey, Ez)
  =#
    rprime=r-R  # r'=r-R
    magrprime=sqrt(sum(rprime.^2)) # |r-R|
    w=2*pi*f  # omega
    k=w/c     # wave number
    krp=k*magrprime  # k*|r-R|
	  rprime_cross_p = cross(rprime, p) # (r-R) x p
    expfac=exp(1im*(-w*t+krp-phi))/(4*pi*eps0*epsr)
    E=expfac*((1/magrprime^3-w*im/(c*magrprime^2))*(3*rprime*dot(rprime,p)/magrprime^2-p))
	  B=expfac/(magrprime*c^2)*(w^2*rprime_cross_p)/magrprime*(-1/(im*w*magrprime))
    return E,B
	end


