subroutine iapibrdf(pild,pxlt,prl,ptl,prs,pihs,pc,mu,np,rm,rp,brdfint)
!
! interface between the computer code of the model of Iaquinta and Pinty
! the computer code is courtesy of Jean Ianquinta
! see module IAPITOOLS.f for a complete description
!
!
    implicit none
    integer :: np,mu
    real(8) :: rm(-mu:mu),rp(np),brdfint(-mu:mu,np)
    real(8) :: ro_1
    integer :: iwr,k,j
    integer :: pild,pihs
    real(8) :: pxlt,prl,ptl,prs,pc
    logical:: ier
    common/sixs_ier/iwr,ier

    real(8) :: mu1,mu2,fi
    real(8) :: pi
! begin of Iaquinta and Pinty model parameter and declaration
    parameter (Pi=3.141592653589793)
    common /gauss_m/xgm (20),wgm (20),n
    real(8) :: xgm,wgm
    integer :: n
    common /p/xLt,Rl,Tl,Rs,c,ild
    real(8) :: xLt,Rl,Tl,Rs,c
    integer :: ild
    common /ld/a_ld,b_ld,c_ld,d_ld
    real(8) :: a_ld,b_ld,c_ld,d_ld
    common /Ro/Ro_1_c,Ro_1_s,Ro_mult
    real(8) :: Ro_1_c,Ro_1_s,Ro_mult
    real(8) :: Theta_i,Phi_i
    real(8) :: Theta_v,Phi_v
    integer :: ihs
!       xLt =  leaf area index
!       Rl = leaf reflection coefficient (Bi-Lambertian)
!       Tl = leaf transmission coefficient (Bi-Lambertian)
!       ild = leaf angle distribution :
!                                       1 = planophile
!                                       2 = erectophile
!                                       3 = plagiophile
!                                       4 = extremophile
!                                       5 = uniform
!       Rs = soil albedo
!       c = 2*r*Lambda
!
!       Ro_1_c  = single scattering by the canopy term
!       Ro_1_s  = uncollided by the leaves (or singly scattered by
!                 the soil) radiation
!                (Ro_1 = Ro_1_c + Ro_1_s)
!       Ro_mult = multiple scattering
! transfer paramater to common / / parameter struture
    ild=pild
    Xlt=pXlt
    Rl=pRl
    Tl=pTl
    Rs=pRs
    ihs=pihs
    c=pc

! Check parameter validity
    if ((ild.ne.1).and.(ild.ne.2).and.(ild.ne.3).and.(ild.ne.4).and.(ild.ne.5)) then
        print*,'Leaf angle distribution !'
        stop
    endif
    if (xlt.le.0.) then
        print*,'Leaf area index < 0. !'
        stop
    endif
    if (xlt.lt.1.) then
        print*,'Leaf area index < 1. !'
    endif
    if (xlt.gt.15.) then
        print*,'Leaf area index > 15. !'
    endif
    if (Rl.lt.0.) then
        print*,'Leaf reflectance < 0. !'
        stop
    endif
    if (Rl.gt..99) then
        print*,'Leaf reflectance > .99 !'
        stop
    endif
    if (Tl.lt.0.) then
        print*,'Leaf transmittance < 0. !'
        stop
    endif
    if (Tl.gt..99) then
        print*,'Leaf transmittance > .99 !'
        stop
    endif
    if (Rl+Tl.gt..99) then
        print*,'Single scattering albedo > .99 !'
        stop
    endif
    if (Rs.lt.0.) then
        print*,'Soil albedo < 0. !'
        stop
    endif
    if (Rs.gt..99) then
        print*,'Soil albedo > .99 !'
        stop
    endif
    if (c.lt.0.) then
        print*,'Hot-spot parameter < 0. !'
        stop
    endif
    if (c.gt.2.) then
        print*,'Hot-spot parameter > 2. !'
        stop
    endif
! compute leaf area angle distribution
    call lad
!
! - Hot-spot parameter
!
    if (ihs.eq.0) c=0.
!
    mu1=rm(0)
    Theta_i=acos(mu1)
    Theta_i=Pi-Theta_i
!
!
! - Gauss's quadrature (n points)
!
    n=10
    call gauleg (-1.d0,1.d0,xgm,wgm,n)
!
! - Computation of the multiple scattering (Ro_mult)
!
    call solve (Theta_i)

    do k=1,np
        do j=1,mu
            mu2=rm(j)
            if (j.eq.mu) then
                fi=rm(-mu)
                else
                fi=rp(k)+rm(-mu)
            endif
            Theta_v=acos(mu2)
            if (fi.lt.0.) fi=fi+2.*pi
            if (fi.gt.(2.*pi)) fi=fi-2.*pi
            Phi_i=fi
            Phi_v=0.
            brdfint(j,k)=Ro_1(Theta_i,Phi_i,Theta_v,Phi_v)+Ro_mult
        enddo
    enddo
    return
end
