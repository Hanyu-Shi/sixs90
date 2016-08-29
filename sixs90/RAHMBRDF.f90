subroutine rahmbrdf(rho0,af,xk,mu,np,rm,rp,brdfint)
!***********************************************************************
!
!           A semi-empirical bidirectional reflectance model
!
!  Purpose:
!
!  To generate a single bidirectional reflectance factor value for a
!  semi-infinite medium, given the illumination and viewing geometry,
!  the optical and the structural properties of the scatterers.
!
!  Definitions:
!   geometrical conditions
!     mu1          : Illumination zenith angle, in radians
!     mu2          : Observation zenith angle, in radians
!     fi           : Relative azimuth angle, in radians
!   optical characteristics of the scatterers:
!     Rho0         : Intensity of the reflectance of the surface cover.
!                    N/D value greater or equal to 0.0 (Rho_0)
!     af           : Phase function parameter:
!                    Asymmetry factor, N/D value between -1.0 and 1.0
!     xk           : Structural parameter of the medium
!
!  References:
!
!        [1] R. Rahman, Verstraete M., and Pinty B., 1993, `A coupled
!            surface-atmosphere reflectance (CSAR) model. Part 1: Model
!            Description and Inversion on Synthetic Data', submitted to
!            JGR.
!        [2] R. Rahman, Pinty B., and Verstraete M., 1993, `A coupled
!            surface-atmosphere reflectance (CSAR) model. Part 2: a
!            semi-empirical model usable with NOAA/AVHRR data',
!            submitted to JGR.
!
!***********************************************************************
!
    implicit none
    integer :: mu,np,j,k
    real(8) :: xk,af,rho0
    real(8) :: rm(-mu:mu),rp(np),brdfint(-mu:mu,np)
    real(8) :: coef1,coef2,cospha,geofac
    real(8) :: fi,mu1,mu2,phafun,pi,tante1,tante2
!
    pi =acos(-1.d0)
    mu1=rm(0)
    do k=1,np
        do j=1,mu
            mu2=rm(j)

!    if (abs(mu1).lt.0.03) then ! Svetlana - Alex' correction
!     mu1=0.03
!    endif
!    if (abs(mu2).lt.0.03) then ! Svetlana - Alex' correction
!     mu2=0.03
!    endif

            if (j.eq.mu)then
                fi=rm(-mu)
            else
                fi=rp(k)+rm(-mu)
            endif
! Compute various trigonometric expressions:
            cospha=mu1*mu2+sqrt(1.-mu1*mu1)*sqrt(1.-mu2*mu2)*cos(fi)
            tante1=sqrt(1.-mu1*mu1)/mu1
            tante2=sqrt(1.-mu2*mu2)/mu2
            geofac=sqrt(tante1*tante1+tante2*tante2-2.0*tante1*tante2*cos(fi))
! Compute the first term
            coef1=(mu1**(xk-1.))*(mu2**(xk-1.))/((mu1+mu2)**(1.-xk))
! Compute the phase function:
            phafun=(1.0-af*af)/((1.0+af*af-2.0*af*(-cospha))**1.5)
! Compute the opposition (hot spot) function:
            coef2=1.+(1.-rho0)/(1.+geofac)
! Compute the bidirectional reflectance factor:
            brdfint(j,k)=rho0*coef1*phafun*coef2
        enddo
    enddo
    return
end
