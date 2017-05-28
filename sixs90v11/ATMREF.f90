subroutine atmref(iaer,iaer_prof,tamoy,taer,trmoy,pizmoy,piza,                      &
                  tamoyp,taerp,trmoyp,palt,phi,xmus,xmuv,phirad,nt,mu,np,rm,gb,rp,  &
                  rorayl,roaero,romix,rqrayl,rqaero,rqmix,rurayl,ruaero,rumix,      &
                  ipol,xlm1,xlm2,rorayl_fi,romix_fi,nfi,                            &
                  nfilut,filut,rolut,rolutq,rolutu)
    implicit none
    integer :: mu,np,nfi
    integer :: nfilut(mu)
    integer :: iaer,nt,ipol,iaer_prof
    integer :: i,ifi,j
	
    real(8) :: rolut(mu,41),rolutq(mu,41),rolutu(mu,41),rolutd(mu,41)
    real(8) :: filut(mu,41)
    real(8) :: rm(-mu:mu),rp(np),gb(-mu:mu)
    real(8) :: tamoy,trmoy,pizmoy,tamoyp,trmoyp,palt
    real(8) :: xlm2(-mu:mu,np)
    real(8) :: xlphim(nfi)
    real(8) :: rorayl_fi(nfi),romix_fi(nfi)
    real(8) :: phi,xmus,xmuv,phirad
    real(8) :: delta,sigma,tamol,tamolp
    real(8) :: rorayl,roaero,romix,rqrayl,rqaero,rqmix
    real(8) :: rurayl,ruaero,rumix
    real(8) :: xlm1(-mu:mu,np),xqm1(-mu:mu,np),xum1(-mu:mu,np)
    real(8) :: taerp,piza,taer

    common /sixs_del/ delta,sigma

!   atmospheric reflectances
    rorayl=0.
    roaero=0.
    romix=0.
    rqrayl=0.
    rqaero=0.
    rqmix=0.
    rurayl=999.
    ruaero=999.
    rumix=999.

! 3 possible cases (satellite,plane,ground)
    if(palt.gt.0.0)then
        rm(-mu)=-xmuv
        rm(mu)=xmuv
        rm(0)=-xmus

!  -----rayleigh reflectance = rorayl,rprayl
        tamol=0.
        tamolp=0.
        call ospol(iaer_prof,tamol,trmoy,piza,tamolp,trmoyp,palt,      &
                   phirad,nt,mu,np,rm,gb,rp,xlm1,xqm1,xum1,xlphim,nfi, &
                   nfilut,filut,rolutd,rolutd,rolutd)
        if (ipol.ne.1)then
            rorayl=xlm1(-mu,1)/xmus
            romix=rorayl
            do ifi=1,nfi
                rorayl_fi(ifi)=(xlphim(ifi)/xmus)
                romix_fi(ifi)=(xlphim(ifi)/xmus)
            enddo
        endif
        if (ipol.ne.0)then
!              -> here we define 2 reflectances from Stockes' parameters
!                 but they don't have any physical interpretations, this is
!                 just to be coherent with the reflectance rorayl
!              -> parameters rorayl2,roaero2,romix2 have been introduced
!                 to compute the degrees of polarization.
            rorayl=xlm1(-mu,1)/xmus
            rqrayl=xqm1(-mu,1)/xmus
            rqmix=rqrayl
            rurayl=xum1(-mu,1)/xmus
            rumix=rurayl
            do ifi=1,nfi
                rorayl_fi(ifi)=(xlphim(ifi)/xmus)
            enddo
        endif

        if(iaer.eq.0) then
            romix=rorayl
            rqmix=rqrayl
            rumix=rurayl
            roaero=0.0
            rqaero=0.0
            ruaero=0.0
            return
        endif

!  -----aerosol reflectance = roaero,rpaero
        tamol=0.
        tamolp=0.
        if (ipol.ne.1)then
            call os(iaer_prof,tamoy,tamol,pizmoy,tamoyp,tamolp,palt, &
                    phirad,nt,mu,np,rm,gb,rp,xlm1,xlphim,nfi,rolutd)
            roaero=(xlm1(-mu,1)/xmus)
        endif
        if (ipol.ne.0)then
        call ospol(iaer_prof,taer,tamol,piza,taerp,tamolp,palt,         &
                   phirad,nt,mu,np,rm,gb,rp,xlm1,xqm1,xum1,xlphim,nfi,  &
                   nfilut,filut,rolutd,rolutd,rolutd)
        rqaero=xqm1(-mu,1)/xmus
        ruaero=xum1(-mu,1)/xmus
        if (ipol.eq.1)roaero=xlm1(-mu,1)/xmus
        endif

!  -----rayleigh+aerosol reflectance = romix,rpmix
        if (ipol.ne.1)then
            call os(iaer_prof,tamoy,trmoy,pizmoy,tamoyp,trmoyp,palt,  &
                    phirad,nt,mu,np,rm,gb,rp,xlm1,xlphim,nfi,rolutd)
            romix=(xlm1(-mu,1)/xmus)
            do ifi=1,nfi
                romix_fi(ifi)=(xlphim(ifi)/xmus)
            enddo
        endif
        if (ipol.ne.0)then
            call ospol(iaer_prof,taer,trmoy,piza,taerp,trmoyp,palt,        &
                       phirad,nt,mu,np,rm,gb,rp,xlm1,xqm1,xum1,xlphim,nfi, &
                       nfilut,filut,rolut,rolutq,rolutu)
            rqmix=xqm1(-mu,1)/xmus
            rumix=xum1(-mu,1)/xmus
            if (ipol.eq.1)then
                romix=xlm1(-mu,1)/xmus
                    do ifi=1,nfi
                        romix_fi(ifi)=(xlphim(ifi)/xmus)
                    enddo
            endif
        endif

        do i=1,mu
            do j=1,41
                rolut(i,j)=rolut(i,j)/xmus
                rolutq(i,j)=rolutq(i,j)/xmus
                rolutu(i,j)=rolutu(i,j)/xmus
            enddo
        enddo

    endif
    return
end
