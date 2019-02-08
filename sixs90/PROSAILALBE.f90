subroutine prosailalbe(TypeLidf,LiDFa,LIDFb,Cab,Car,Anth,Cbrown, &
            Cw,Cm,N,lai,hspot,psoil,wlmoy,brdfalb)
    implicit none
    integer, parameter :: nta = 24, nfa = 48
    real(8) :: teta1,teta2,phi1,phi2,ta(nta),fa(nfa),wta(nta),wfa(nfa)
    real(8) :: pi,wlmoy,brdfalb,summ,wl,dr
    real(8) :: mu1,mu2,si1,si2,rts,rtv,phi,pond,y
    integer :: k,j,l
    integer :: TypeLidf
    real(8) :: LIDFa,LIDFb,Cab,Car,Anth,Cbrown,Cw,Cm,N,lai,hspot,psoil

    wl = wlmoy * 1000.d0

    pi = atan(1.d0)*4.d0
    dr = pi / 180.d0
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

                call pro_sail(TypeLidf,LiDFa,LIDFb,Cab,Car,Anth,Cbrown, &
                    Cw,Cm,N,lai,hspot,psoil, &
                    rts/dr,rtv/dr,phi/dr,wl,y)
                pond=mu1*mu2*si1*si2*wfa(k)*wta(j)*wta(l)
                brdfalb = brdfalb + y*pond
                summ = summ + pond
            end do
        end do
    end do
    brdfalb = brdfalb / summ
    return
end subroutine prosailalbe
