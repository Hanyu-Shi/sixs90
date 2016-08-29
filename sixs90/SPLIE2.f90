subroutine splie2(x2a,ya,m,n,y2a)
    implicit none
    integer,parameter::nn=100
    integer :: m,n,j,k
    real(8) :: x2a(n),ya(m,n),y2a(m,n),ytmp(nn),y2tmp(nn)
    do j=1,m
        do k=1,n
            ytmp(k)=ya(j,k)
        enddo
        call spline(x2a,ytmp,n,1.d30,1.d30,y2tmp)
        do k=1,n
            y2a(j,k)=y2tmp(k)
        enddo
    enddo
    return
end
