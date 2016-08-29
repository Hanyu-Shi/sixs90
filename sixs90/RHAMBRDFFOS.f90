subroutine RAHMBRDFFOS(rho0, af, xk, mu, rm, rosur, wfisur, fisur)
!***********************************************************************
!
!           a semi-empirical bidirectional reflectance model
!
!  purpose:
!
!  to generate a single bidirectional reflectance factor value for a
!  semi-infinite medium, given the illumination and viewing geometry,
!  the optical and the structural properties of the scatterers.
!
!  definitions:
!   geometrical conditions
!     mu1          : illumination zenith angle, in radians
!     mu2          : observation zenith angle, in radians
!     fi           : relative azimuth angle, in radians
!   optical characteristics of the scatterers:
!     rho0         : intensity of the reflectance of the surface cover.
!                    n/d value greater or equal to 0.0 (rho_0)
!     af           : phase function parameter:
!                    asymmetry factor, n/d value between -1.0 and 1.0
!     xk           : structural parameter of the medium
!
!  references:
!
!        [1] r. rahman, verstraete m., and pinty b., 1993, `a coupled
!            surface-atmosphere reflectance (csar) model. part 1: model
!            description and inversion on synthetic data', submitted to
!            jgr.
!        [2] r. rahman, pinty b., and verstraete m., 1993, `a coupled
!            surface-atmosphere reflectance (csar) model. part 2: a
!            semi-empirical model usable with noaa/avhrr data',
!            submitted to jgr.
!
!***********************************************************************
!
  use paramdef
  implicit none
  integer :: mu, j, k, i
  real(8) :: xk, af, rho0
  real(8) :: rm(-mu:mu)
  real(8) :: coef1, coef2, cospha, geofac
  real(8) :: fi, mu1, mu2, phafun, tante1, tante2
  real(8) :: rosur(0:mu_p, mu_p, 83), fisur(83), wfisur(83), pisp !psip !! HyS
!
  pisp = acos(-1.)
!      mu=mu_p
  call gauss(0.d0, pisp, fisur, wfisur, 83)
  do i = 0, mu
    do j = 1, mu
      do k = 1, 83
        if (i==0) then
          mu1 = rm(0)
        else
          mu1 = rm(i)
        end if
        mu2 = rm(j)
        fi = fisur(k) + pisp
! compute various trigonometric expressions:
        cospha = mu1*mu2 + sqrt(1.-mu1*mu1)*sqrt(1.-mu2*mu2)*cos(fi)
        if (cospha<=-1.0) cospha = -1.0
        if (cospha>=1.0) cospha = 1.0
        tante1 = sqrt(1.-mu1*mu1)/mu1
        tante2 = sqrt(1.-mu2*mu2)/mu2
        geofac = sqrt(tante1*tante1+tante2*tante2-2.0*tante1*tante2*cos(fi))
! compute the first term
        coef1 = (mu1**(xk-1.))*(mu2**(xk-1.))/((mu1+mu2)**(1.-xk))
! compute the phase function:
        phafun = (1.0-af*af)/((1.0+af*af-2.0*af*(-cospha))**1.5)

! compute the opposition (hot spot) function:
        coef2 = 1. + (1.-rho0)/(1.+geofac)
! compute the bidirectional reflectance factor:
        rosur(i, j, k) = rho0*coef1*phafun*coef2
!      write(6,*) "rosur", mu1,mu2,fi,rosur(i,j,k)
      end do
    end do
  end do
!
  return
end subroutine rahmbrdffos
