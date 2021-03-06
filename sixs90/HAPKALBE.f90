subroutine hapkalbe(om,af,s0,h,brdfalb)
    implicit none
    integer,parameter::nta=24,nfa=48
    real(8) :: teta1,teta2,phi1,phi2,ta(nta),fa(nfa),wta(nta),wfa(nfa)
    real(8) :: om,af,s0,h,mu1,mu2
    real(8) :: brdfalb,summ,si2,si1
    real(8) :: fi,f,cg,h1,h2,h1h2,pg,p0,g,bg,pond
    real(8) :: pi
    integer :: k,j,l
    pi=atan(1.d0)*4.
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
                f=om/4./(mu2+mu1)
                cg=mu1*mu2+sqrt(1.-mu1*mu1)*sqrt(1.-mu2*mu2)*cos(fi)
                h1=(1.+2*mu1)/(1.+2.*sqrt(1.-om)*mu1)
                h2=(1.+2*mu2)/(1.+2.*sqrt(1.-om)*mu2)
                h1h2=h1*h2
                pg=(1-af*af)/((1+af*af+2*af*cg)**1.5)
                p0=(1-af*af)/((1+af*af+2*af)**1.5)
                g=acos(cg)
                bg=(s0/(om*p0))/(1.+tan(g/2.)/h)
                pond=mu1*mu2*si1*si2*wfa(k)*wta(j)*wta(l)
                brdfalb=brdfalb+f*((1.+bg)*pg+h1h2-1.)*pond
                summ=summ+pond
            enddo
        enddo
    enddo
    brdfalb=brdfalb/summ
    return
end
