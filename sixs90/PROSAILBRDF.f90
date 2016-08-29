!***********************************************************************
subroutine prosailbrdf(TypeLidf,LiDFa,LIDFb,Cab,Car,Cbrown, &
            Cw,Cm,N,lai,hspot,psoil, &
            wlmoy,mu,np,rm,rp,brdfint)
    implicit none
    integer :: mu,np
    real(8) :: wlmoy,rm(-mu:mu),rp(np),brdfint(-mu:mu,np)
    real(8) :: pir,dr,rtv,rts,rfi,r_lamda,bq
    integer :: TypeLidf
    real(8) :: LIDFa,LIDFb,Cab,Car,Cbrown,Cw,Cm,N,lai,hspot,psoil

    integer :: k,j

    r_lamda = wlmoy*1000.d0
    pir = acos(-1.d0)
    dr = pir/180.0
    rts = acos(rm(0))
    do k=1,np
        do j=1,mu
            rtv = acos(rm(j))
            if(j.eq.mu) then
                rfi = rm(-mu)
            else
                rfi = rp(k) + rm(-mu)
            endif
            if (rfi .lt. 0.) rfi = rfi + 2.*pir
            if (rfi .gt. (2.*pir)) rfi = rfi - 2.*pir
            if (rfi .gt. pir) rfi = 2.*pir - rfi

            call pro_sail(TypeLidf,LiDFa,LIDFb,Cab,Car,Cbrown, &
                    Cw,Cm,N,lai,hspot,psoil, &
                    rts/dr,rtv/dr,rfi/dr,r_lamda,bq) ! trans to degree
            brdfint(j,k) = bq
        enddo
    enddo
end subroutine prosailbrdf
