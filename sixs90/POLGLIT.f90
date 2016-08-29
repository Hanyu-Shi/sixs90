subroutine POLGLIT(xts,xtv,phi,wspd,azw,ropq,ropu)
    implicit none
    real(8) :: xts,xtv,phi,azw,wspd,ropq,ropu
!      compute the polarized components of the surface
!      to an agitated surface model

!   xts is the sun zenith angle in degrees
!   xtv is the view zenith angle in degrees
!   phi is the relative azimuth between sun and view in degrees
!   wspd is the wind speed for use in ponderating
!       azw=azim. of the sun - azim. of the wind (in deg.)

!   the polarization by the actual glint reflectance
!   ropq and ropu and the stokes parameter (Q,U) of the polarized
!   surface reflectance
    real(8) :: pi,dtr,csca,sca,alpha,m,rl,rr,factor
!   csca is the cosine of the scattering angle (sca)
!   alpha is the incidence angle used for specular reflection
!   alphap is the refraction angle after specular reflection
!   m is the refractive index of water (fixed to 1.33)
!   rr and rl and parameters used in the reflection mattrix
!   computation.
!   factor is a multiplicative factor to account for agitated surface

! following is a set of intermediate variable used in the computation
!  of factor
    real(8) :: phw
    real(8) :: cs,cv,ss,sv,zx,zy,tantilt,tilt,proba,xe,xn,xe2,xn2
    real(8) :: coef
    real(8) :: sigmaC,sigmaU,C21,C03,C40,C04,C22
    real(8) :: mus,muv,sinv,cksi,ksi
    real(8) :: nr,ni
    real(8) :: r1,r2,r3

    m=1.33
    nr=1.33
    ni=0.0
    pi=acos(0.d0)*2.0
    dtr=pi/180.0
    csca=-cos(xts*dtr)*cos(xtv*dtr)-sin(xts*dtr)*sin(xtv*dtr)*cos(phi*dtr)
    sca=acos(csca)
    alpha=(pi-sca)/2.0

    rl=(sqrt(m*m-sin(alpha)*sin(alpha))-m*m*cos(alpha))/   &
          (sqrt(m*m-sin(alpha)*sin(alpha))+m*m*cos(alpha))

    rr=(cos(alpha)-sqrt(m*m-sin(alpha)*sin(alpha)))/    &
         (cos(alpha)+sqrt(m*m-sin(alpha)*sin(alpha)))

    r1=(rl*rl+rr*rr)/2.
    r2=(rl*rl-rr*rr)/2.
    r3=rl*rr
    r3=0.0

    phw=azw*dtr
    cs=cos(xts*dtr)
    cv=cos(xtv*dtr)
    ss=sin(xts*dtr)
    sv=sin(xtv*dtr)
    Zx=-sv*sin(phi*dtr)/(cs+cv)
    Zy=(ss+sv*cos(phi*dtr))/(cs+cv)
    tantilt=sqrt(zx*zx+zy*zy)
    tilt=atan(tantilt)
    sigmaC=0.003+0.00192*wspd
    sigmaU=0.00316*wspd
    C21=0.01-0.0086*wspd
    C03=0.04-0.033*wspd
    C40=0.40
    C22=0.12
    C04=0.23
    xe=(cos(phw)*Zx+sin(phw)*Zy)/sqrt(SigmaC)
    xn=(-sin(phw)*Zx+cos(phw)*Zy)/sqrt(SigmaU)
    xe2=xe*xe
    xn2=xn*xn
    coef=1-C21/2.*(xe2-1)*xn-C03/6.*(xn2-3)*xn
    coef=coef+c40/24.*(xe2*xe2-6*xe2+3)
    coef=coef+C04/24.*(xn2*xn2-6*xn2+3)
    coef=coef+C22/4.*(xe2-1)*(xn2-1)
    proba=coef/2./pi/sqrt(sigmaU)/sqrt(sigmaC)*exp(-(xe2+xn2)/2.)
!      write(6,*) "probapol:",proba
!      write(6,*) "coefpol:",coef
!      write(6,*) "tiltpol:",tilt
!      write(6,*) "phiw pol:",phw
! Compute Fresnel's coefficient R1
!      cos2chi=cv*cs+sv*ss*cos(phi*dtr)
!      if (cos2chi.gt.1.0)cos2chi=0.99999999999
!      if (cos2chi.lt.-1.0)cos2chi=-0.99999999999
!      coschi=sqrt(0.5*(1+cos2chi))
!      sinchi=sqrt(0.5*(1-cos2chi))
!      Call Fresnel(nr,ni,coschi,sinchi,R1f)
!      Call pfresnel(nr,ni,coschi,sinchi,r1f,r2f)
!      write(6,*) "R1 fresnel:",R1f," r1 actual:",r1
!      write(6,*) "R2 fresnel:",R2f," r2 actual:",r2
! Compute Reflectance of the sun glint
!      Rog=pi*R1*proba/4./cs/cv/(cos(tilt)**4)
        factor=pi*proba/4./cs/cv/(cos(tilt)**4)
! compute rotation angle for Q and U
    muv=cos(xtv*dtr)
    mus=cos(xts*dtr)
    sinv=sin(xtv*dtr)
    if (xtv.gt.0.5) then
    if (sin(phi*dtr).lt.0) then
    cksi=(muv*csca+mus)/sqrt(1.-csca*csca)/sinv
    else
    cksi=-(muv*csca+mus)/sqrt(1.-csca*csca)/sinv
    endif
    else
    cksi=1.0
    endif
    if (cksi.gt.1.) cksi=1.
    if (cksi.lt.-1.) cksi=-1.
    ksi=acos(cksi)/dtr

    ropq=r2*(2.*cksi*cksi-1.)*factor
    ropu=-r2*2.*cksi*sqrt(1.-cksi*cksi)*factor
    return
end

subroutine pfresnel(nr,ni,coschi,sinchi,r1,r2)
!
! to compute the Fresnel's coefficient of reflection (see for
! example M. Born and E. Wolf, Principles of Optics, Pergamon Press, fifth
! edition, 1975, pp 628
! input parameters: nr=index of refraction of the sea water
!                   ni=extinction coefficient of the sea water
!                   coschi & sinchi=cosine and sine of the incident radiation
!                                   with respect of the wave facet normal.
! output parameter: R1=Fresnel's coefficient for reflection
!
    implicit none
    real(8) :: nr,ni,a1,a2,u,v,Rr2,Rl2,b1,b2,R1,coschi,sinchi,R2
! absolute value for a1 to get v=0 when ni=0
    a1=abs(nr*nr-ni*ni-sinchi*sinchi)
    a2=sqrt((nr*nr-ni*ni-sinchi*sinchi)**2.+4*nr*nr*ni*ni)
    u=sqrt(0.5*(a1+a2))
    v=sqrt(0.5*(-a1+a2))
    Rr2=((coschi-u)**2+v*v)/((coschi+u)**2+v*v)
    b1=(nr*nr-ni*ni)*coschi
    b2=2*nr*ni*coschi
    Rl2=((b1-u)**2+(b2+v)**2)/((b1+u)**2+(b2-v)**2)
    R1=(Rr2+Rl2)/2.
    R2=(Rl2-Rr2)/2.
    return
end
