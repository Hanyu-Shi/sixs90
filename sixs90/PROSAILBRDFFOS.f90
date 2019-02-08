subroutine prosailbrdffos(TypeLidf,LiDFa,LIDFb,Cab,Car,Anth,Cbrown, &
            Cw,Cm,N,lai,hspot,psoil, &
            wlmoy,mu,rm,rosur,wfisur,fisur)
    use paramdef
    implicit none
    integer :: mu
    real(8) :: wlmoy,rm(-mu:mu),rosur(0:mu_p,mu_p,83),fisur(83),wfisur(83),pisp
    real(8) :: mu1,mu2
    real(8) :: dr,rtv,rts,rfi,r_lamda,bq
    integer :: TypeLidf
    real(8) :: LIDFa,LIDFb,Cab,Car,Anth,Cbrown,Cw,Cm,N,lai,hspot,psoil

    integer :: k,j,i

    r_lamda = wlmoy*1000.d0
    pisp = acos(-1.d0)
    dr = pisp/180.0

    call gauss(0.d0,pisp,fisur,wfisur,83)

    do i = 0, mu
        do j = 1, mu
            do k = 1, 83
                if (i == 0) then
                    mu1 = rm(0)
                else
                    mu1 = rm(i)
                end if
                mu2 = rm(j)
                rfi = fisur(k) + pisp
                rts = acos(mu1)
                rtv = acos(mu2)

                if (rfi < 0.) rfi = rfi + 2.*pisp
                if (rfi > (2.*pisp)) rfi = rfi - 2.*pisp
                if (rfi > pisp) rfi = 2.*pisp - rfi

                call pro_sail(TypeLidf,LiDFa,LIDFb,Cab,Car,Anth,Cbrown, &
                    Cw,Cm,N,lai,hspot,psoil, &
                    rts/dr,rtv/dr,rfi/dr,r_lamda,bq) ! trans to degree
                rosur(i,j,k) = bq
            end do
        end do
    end do
end subroutine prosailbrdffos
