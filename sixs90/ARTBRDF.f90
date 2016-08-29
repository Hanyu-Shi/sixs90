subroutine artbrdf(dia,M,wlmoy,mu,np,rm,rp,brdfint)
    implicit none
    integer :: mu,np
    real(8) :: wlmoy,rm(-mu:mu),rp(np),brdfint(-mu:mu,np)
    real(8) :: dia, M
    real(8) :: pisp,dr,rtv,rts,rfi,r_lamda,bq
    integer :: k,j

    r_lamda = wlmoy * 1000.d0
    pisp = acos(-1.d0)
    dr = pisp / 180.d0
    rts = acos(rm(0))
    do k = 1, np
        do j = 1, mu
            rtv = acos(rm(j))
            if (j == mu) then
                rfi = rm(-mu)
            else
                rfi = rp(k) + rm(-mu)
            end if

            if (rfi < 0.) rfi = rfi + 2.*pisp
            if (rfi > (2.*pisp)) rfi = rfi - 2.*pisp
            if (rfi > pisp) rfi = 2.*pisp - rfi

            call ART(rts,rtv,rfi,dia,M,r_lamda,bq)
            brdfint(j,k) = bq
        end do
    end do
end subroutine
