subroutine rahmalbe (rho0,af,xk,brdfalb)
    implicit none
    integer,parameter::nta=24,nfa=48
!
! see RAHMBRDF.f for description
!
    integer :: j,k,l
    real(8) :: teta1,teta2,phi1,phi2,ta(nta),fa(nfa),wta(nta),wfa(nfa)
    real(8) :: xk,af,rho0
    real(8) :: brdfalb,summ,si1,si2,pond
    real(8) :: coef1,coef2,cospha,geofac
    real(8) :: fi,mu1,mu2,phafun,pi,tante1,tante2
!
    pi =acos(-1.d0)
    teta1=0.
    teta2=pi/2.
    call gauss(teta1,teta2,ta,wta,nta)
    phi1=0.
    phi2=2.*pi
    call gauss(phi1,phi2,fa,wfa,nfa)
    brdfalb=0.
    summ=0.
    do k=1,nfa
        do j=1,nta
            do l=1,nta
                mu2=cos(ta(j))
                mu1=cos(ta(l))
                si2=sin(ta(j))
                si1=sin(ta(l))
                fi=fa(k)
! Compute various trigonometric expressions:
                cospha=mu1*mu2+sqrt(1.-mu1*mu1)*sqrt(1.-mu2*mu2)*cos(fi)
                tante1=sqrt(1.-mu1*mu1)/mu1
                tante2=sqrt(1.-mu2*mu2)/mu2
                geofac=sqrt(tante1*tante1+tante2*tante2-2.0*tante1*tante2*cos(fi))
! Compute the first term
                coef1=(mu1**(xk-1.))*(mu2**(xk-1.))/((mu1+mu2)**(1.-xk))
! Compute the phase function:
                phafun=(1.0-af*af)/((1.0+af*af-2.0*af*(-cospha))**1.5)
! Compute the opposition (hot spot) function:
                coef2=1.+(1.-rho0)/(1.+geofac)
! Compute the bidirectional reflectance factor:
                pond=mu1*mu2*si1*si2*wfa(k)*wta(j)*wta(l)
                brdfalb=brdfalb+rho0*coef1*phafun*coef2*pond
                summ=summ+pond
            enddo
        enddo
    enddo
    brdfalb=brdfalb/summ
    return
end
