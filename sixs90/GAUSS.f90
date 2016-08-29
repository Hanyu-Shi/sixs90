subroutine gauss(x1,x2,x,w,n)
    implicit none
    integer :: n
    real(8) :: x1,x2,x(n),w(n)
    real(8) :: xm,xl,z,p1,p2,p3,pp,z1
    integer :: m,i,j
    real(8),parameter::eps = 3.0e-14
    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    do i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
        p1=1.d0
        p2=0.d0
        do j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
        enddo
            pp=n*(z*p1-p2)/(z*z-1.d0)
            z1=z
            z=z1-p1/pp
        if(abs(z-z1).gt.eps)go to 1
        if(abs(z).lt.eps) z=0.
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
    enddo
    return
end
