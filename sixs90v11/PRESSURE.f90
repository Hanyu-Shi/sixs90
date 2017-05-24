subroutine pressure(uw,uo3,xps)
    implicit none
    common /sixs_atm/z(34),p(34),t(34),wh(34),wo(34)
    real(8) :: z,p,t,wh,wo,xa,xb,xalt,xtemp,xwo,xwh,g
    real(8) :: air,ro3,roair,ds
    integer :: i,isup,iinf,l,k
    real(8) :: rmo3(34),rmwh(34)
    real(8) :: ps,xps,uo3,uw
    real(8) :: xps2
!   log linear interpolation
        xps2=-xps
        if (xps2.ge.100.) xps2=99.99
        i=0
10      i=i+1
        if (z(i).le.xps2) goto 10
        isup=i
        iinf=i-1
        xa=(z(isup)-z(iinf))/log(p(isup)/p(iinf))
        xb=z(isup)-xa*log(p(isup))
        ps=exp((xps2-xb)/xa)
! interpolating temperature wator vapor and ozone profile versus altitud
        xalt=xps2

        xtemp=(t(isup)-t(iinf))/(z(isup)-z(iinf))
        xtemp=xtemp*(xalt-z(iinf))+t(iinf)
        xwo=(wo(isup)-wo(iinf))/(z(isup)-z(iinf))
        xwo=xwo*(xalt-z(iinf))+wo(iinf)
        xwh=(wh(isup)-wh(iinf))/(z(isup)-z(iinf))
        xwh=xwh*(xalt-z(iinf))+wh(iinf)
! uptading atmospheric profile
!  1rst level: target     , complete to 34
!  with interpolated layers
    z(1)=xalt
    p(1)=ps
    t(1)=xtemp
    wh(1)=xwh
    wo(1)=xwo
    do i=2,33-iinf+1
        z(i)=z(i+iinf-1)
        p(i)=p(i+iinf-1)
        t(i)=t(i+iinf-1)
        wh(i)=wh(i+iinf-1)
        wo(i)=wo(i+iinf-1)
    enddo
    l=33-iinf+1
    do i=l+1,34
        z(i)=(z(34)-z(l))*(i-l)/(34-l)+z(l)
        p(i)=(p(34)-p(l))*(i-l)/(34-l)+p(l)
        t(i)=(t(34)-t(l))*(i-l)/(34-l)+t(l)
        wh(i)=(wh(34)-wh(l))*(i-l)/(34-l)+wh(l)
        wo(i)=(wo(34)-wo(l))*(i-l)/(34-l)+wo(l)
    enddo

! compute modified h2o and o3 integrated content
    uw=0.
    uo3=0.
    g=98.1
    air=0.028964/0.0224
    ro3=0.048/0.0224
    do k=1,33
        roair=air*273.16*p(k)/(1013.25*t(k))
        rmwh(k)=wh(k)/(roair*1000.)
        rmo3(k)=wo(k)/(roair*1000.)
    enddo
    do k=2,33
        ds=(p(k-1)-p(k))/p(1)
        uw=uw+((rmwh(k)+rmwh(k-1))/2.)*ds
        uo3=uo3+((rmo3(k)+rmo3(k-1))/2.)*ds
    enddo
    uw=uw*p(1)*100./g
    uo3=uo3*p(1)*100./g
    uo3=1000.*uo3/ro3
    return
end
