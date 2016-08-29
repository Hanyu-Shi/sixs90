subroutine splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
    implicit none
    integer,parameter::nn=100
    integer :: m,n,j,k
    real(8) :: x1,x2,y
    real(8) :: x1a(m),x2a(n),ya(m,n),y2a(m,n),ytmp(nn),y2tmp(nn)
    real(8) :: yytmp(nn)
    do j=1,m
        do  k=1,n
            ytmp(k)=ya(j,k)
            y2tmp(k)=y2a(j,k)
        enddo
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
    enddo
    call spline(x1a,yytmp,m,1.d30,1.d30,y2tmp)
    call splint(x1a,yytmp,y2tmp,m,x1,y)
    return
end
