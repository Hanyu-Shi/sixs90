subroutine oceaalbe(pws,paw,xsal,pcl,pwl,brdfalbe)
! INPUT:  pws=speed of wind (in m/s)
!         paw=azim. of sun - azim. of wind (in deg.)
!     xsal=salinity (in ppt)
!     pcl=pigment concentration (in mg.m-3)
!         pwl=wavelength of the computation (in micrometer)
! OUTPUT: brdfalbe=the spherical albedo of the ocean
!
    implicit none
    real(8) :: Ref(39)
    real(8) :: pwl,azw,pcl,wl,wspd,C,pws,brdfalbe,w,wlp,paw
    real(8) :: ref_i,rwc,rw,rogalbe,a,rwb,xsal
    real(8) :: nr,ni
    integer :: iwl

! effective reflectance of the whitecaps (Koepke, 1984)
    data Ref/                                                         &
     0.220,0.220,0.220,0.220,0.220,0.220,0.215,0.210,0.200,0.190,     &
     0.175,0.155,0.130,0.080,0.100,0.105,0.100,0.080,0.045,0.055,     &
     0.065,0.060,0.055,0.040,0.000,0.000,0.000,0.000,0.000,0.000,     &
     0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/
! conversion of parameter
    C=pcl
    wspd=pws
    azw=paw
    wl=pwl

! COMPUTE WHITECAPS REFLECTANCE (LAMBERTIAN)
    W=2.95e-06*(wspd**3.52)
    iwl=1+int((wl-0.2)/0.1)
    wlp=0.5+(iwl-1)*0.1
    Ref_i=ref(iwl+1)+(wl-wlp)/0.1*(ref(iwl)-ref(iwl+1))
    Rwc=W*Ref_i
! COMPUTE WATER REFRACTION INDEX
    call indwat(wl,xsal,nr,ni)
! COMPUTE BACKSCATTERED REFLECTANCE FROM THE SEA WATER (LAMBERTIAN)
!  water reflectance below the sea surface
    call morcasiwat(wl,C,Rw)
! SUNGLINT spherical albedo
    call glitalbe(wspd,nr,ni,azw,rogalbe)
!  water reflectance above the sea surface, (albedo re=0.485)
! albedo is a=re is taken from table 2 of Austin,1974,The remote sensing
! of spectral radiance from below the ocean surface, in Optical
! Aspects of Oceanography (N.G. Jerlov and E. Steeman Nielsen,Eds),
! Academic,London,pp. 317-344
    a=0.485
    Rwb=(1.-Rogalbe)*(1.-a)*Rw/(1-a*Rw)
! SPHERICAL ALBEDO OF SEA WATER
    brdfalbe=Rwc+(1-W)*Rogalbe+(1-Rwc)*Rwb
    return
end
