subroutine brdfgrid(mu,np,rm,rp,brdfdat,angmu,angphi,brdfint)
    implicit none
    integer :: mu,np
    real(8) :: rp(np),brdfint(-mu:mu,np),rm(-mu:mu)  &
        ,angmu(10),angphi(13),brdfdat(10,13)
    real(8) :: brdftemp(10,13)
    real(8) :: gaussmu,gaussphi,y
    integer :: j,k

    do j=1,np
        do k=1,mu
            brdfint(k,j)=0.
        enddo
    enddo
        call splie2(angphi,brdfdat,10,13,brdftemp)
    do j=1,np
        do k=1,mu
            gaussmu=rm(k)
            gaussphi=rp(j)
            call splin2(angmu,angphi,brdfdat,brdftemp,10,13,gaussmu,gaussphi,y)
            brdfint(k,j)=y
        enddo
    enddo
    return
end
