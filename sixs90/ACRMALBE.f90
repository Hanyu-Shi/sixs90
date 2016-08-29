subroutine acrmalbe(lai2,sl2,clmp2,szz,eln2,thm2,nratio2,slw2,lmod2,&
            ncomp2,ccomp2,N2,dcell2,asp2, &
            lai1,sl1,clmp1,eln1,thm1,nratio1,slw1,lmod1, &
            ncomp1,ccomp1,N1,dcell1,asp1,&
            s1,s2,s3,s4, &
            wlmoy,brdfalb)
    implicit none
    integer, parameter :: nta = 24, nfa = 48
    real(8) :: teta1,teta2,phi1,phi2,ta(nta),fa(nfa),wta(nta),wfa(nfa)
    real(8) :: pi,wlmoy,brdfalb,summ,wl,dr
    real(8) :: mu1,mu2,si1,si2,rts,rtv,phi,pond,y
    integer :: k,j,l
    character(30)::lmod1,lmod2
    integer :: ncomp2,ncomp1
    real(8) :: LAI2,sl2,clmp2,szz,eln2,thm2,nratio2,slw2,N2,dcell2,asp2, &
         LAI1,sl1,clmp1,eln1,thm1,nratio1,slw1,N1,dcell1,asp1, &
         s1,s2,s3,s4,ccomp1(10), ccomp2(10)

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

                call ACRM(lai2,sl2,clmp2,szz,eln2,thm2,nratio2,slw2,lmod2,&
                    ncomp2,ccomp2,N2,dcell2,asp2, &
                    lai1,sl1,clmp1,eln1,thm1,nratio1,slw1,lmod1, &
                    ncomp1,ccomp1,N1,dcell1,asp1,&
                    s1,s2,s3,s4,&
                    rts,rtv,phi,wl,y)
                pond=mu1*mu2*si1*si2*wfa(k)*wta(j)*wta(l)
                brdfalb = brdfalb + y*pond
                summ = summ + pond
            end do
        end do
    end do
    brdfalb = brdfalb / summ
    return
end
