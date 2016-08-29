subroutine artalbe(dia,M,wlmoy,brdfalb)
    implicit none
    integer, parameter :: nta = 24, nfa = 48
    real(8) :: teta1,teta2,phi1,phi2,ta(nta),fa(nfa),wta(nta),wfa(nfa)
    real(8) :: dia,M,pi,wlmoy,brdfalb,summ,wl
    real(8) :: mu1,mu2,si1,si2,rts,rtv,phi,pond,y
    integer :: k,j,l

    wl = wlmoy * 1000.d0

    pi = atan(1.d0)*4.d0
    teta1 = 0.d0
    teta2 = pi / 2.d0
    call gauss(teta1,teta2,ta,wta,nta)
    phi1 = 0.d0
    phi2 = 2.d0*pi
    call gauss(phi1,phi2,fa,wfa,nfa)
    brdfalb = 0.d0
    summ = 0.d0

    do k = 1, nfa
        do j = 1, nta
            do l = 1, nta
                si2=sin(ta(j))
                si1=sin(ta(l))
                mu2=cos(ta(j))
                mu1=cos(ta(l))
                rts=ta(j)
                rtv=ta(l)
                phi=fa(k)
                call ART(rts,rtv,phi,dia,M,wl,y)
                pond=mu1*mu2*si1*si2*wfa(k)*wta(j)*wta(l)
                brdfalb = brdfalb + y*pond
                summ = summ + pond
            end do
        end do
    end do
    brdfalb = brdfalb / summ


!    integer :: i,j,k,npoints
!    real(8) :: wlmoy,brdfalb,wl
!    real(8) :: psi,thetaS,thetaV,y,y1
!    real(8) :: x(48),weight(48),tmp
!    real(8) :: dia,M
!    real(8) :: pi,dr
!
!    pi = acos(-1.d0)
!    dr = pi / 180.0
!
!    npoints = 6
!
!    call gauss(-1.d0,1.d0,x,weight,npoints)
!
!    wl = wlmoy*1000.d0
!    y1 = 0.d0
!    do k=1,npoints
!        psi = (x(k)+1)*pi/2.0
!        do j=1,npoints
!            thetaS = (x(j)+1)*pi/4.0
!            tmp = weight(k)*weight(j)*sin(2*thetaS)
!            do i=1,npoints
!                thetaV = (x(i)+1)*pi/4.0
!                call ART(thetaS,thetaV,psi,dia,M,y)
!                y1 = y1 + tmp*weight(i)*sin(2*thetaV)*y
!            enddo
!        enddo
!    enddo
!    brdfalb = y1*pi*pi/32.0
end subroutine
