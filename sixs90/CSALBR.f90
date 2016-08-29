subroutine csalbr(xtau,xalb)
    implicit none
    real(8) :: xtau,xalb,fintexp3
    xalb=(3*xtau-fintexp3(xtau)*(4+2*xtau)+2*exp(-xtau))
    xalb=xalb/(4.+3*xtau)
    return
end


real(8) function fintexp3(xtau)
    implicit none
    real(8) :: xx,xtau,fintexp1
    xx=(exp(-xtau)*(1.-xtau)+xtau*xtau*fintexp1(xtau))/2.
    fintexp3=xx
    return
end

real(8) function fintexp1(xtau)
! accuracy 2e-07... for 0<xtau<1
    implicit none
    real(8) :: xx,a(0:5),xtau,xftau
    integer :: i
    data (a(i),i=0,5) /-.57721566,0.99999193,-0.24991055,&
                       0.05519968,-0.00976004,0.00107857/
    xx=a(0)
    xftau=1.
    do i=1,5
        xftau=xftau*xtau
        xx=xx+a(i)*xftau
    enddo
    fintexp1=xx-log(xtau)
    return
end
