subroutine iso(iaer_prof,tamoy,trmoy,pizmoy,tamoyp,trmoyp,palt,nt,mu,rm,gb,xf)

    use paramdef
    implicit none
    real(8) :: psl(-1:nqmax_p,-mu:mu)
!  dimension for gauss integration
    integer :: mu
    real(8) :: rm(-mu:mu),gb(-mu:mu)
!  dimension for os computation
    real(8) ::  xf(-1:1)
! array for sos computation
    real(8) :: xpl(-mu:mu),bp(0:mu,-mu:mu),ch(0:nt)
    real(8) :: xdel(0:nt),ydel(0:nt),h(0:nt),altc(0:nt)
    real(8) :: i1(0:nt,-mu:mu),i2(0:nt,-mu:mu),i3(-mu:mu)
    real(8) :: in(-mu:mu),inm1(-mu:mu),inm2(-mu:mu)
    real(8) :: acu,acu2,ta,piz
    real(8) :: tr,trp,tap,hr,ha,zx,yy,dd,ppp2,ppp1,ca
    real(8) :: cr,ratio,taup,th,xt1,xt2,aaaa,ron,beta0,beta2
    real(8) :: tavion0,tavion1,tavion2,tavion,zi1,xpk,ypk,x,y
    real(8) :: z,xi1,xi2,bpjk,bpjmk,f,a,b,c,d,xx,a1,d1,g1
    real(8) :: y1,xpj,xxx,ii1,ii2
    real(8) :: tamoy,trmoy,pizmoy
    real(8) :: tamoyp,trmoyp,palt
    real(8) :: delta,sigma
    integer :: snt,nt,iplane,ntp,j,it,itp,i,ig,k,index,iwr,m
    integer :: jj,l
    logical:: ier
    integer :: igmax,iaer_prof

    ! HyS
    real(8) :: xmus

    common/sixs_del/delta,sigma
    common/sixs_ier/iwr,ier
    common /multorder/ igmax

    ! HyS WuTeng
    xmus = rm(0)

    snt=nt
    iplane=0
    acu=1.e-20
    acu2=1.e-3
    ta=tamoy
    piz=pizmoy
    tr=trmoy
    do m=-1,1
        xf(m)=0.
    enddo
!
!     molecular ratio within the layer
!     computations are performed assuming a scale of 8km for
!     molecules and 2km for aerosols
!
! the optical thickness above plane are recomputed to give o.t above pla
    trp=trmoy-trmoyp
    tap=tamoy-tamoyp
!     print *, 'tamoy,trmoy,pizmoy,tap,trp,palt,nt'
!     print *,tamoy,trmoy,pizmoy,tap,trp,palt,nt
    acu=1.e-20
! if plane observations recompute scale height for aerosol knowing:
! the aerosol optical depth as measure from the plane   = tamoyp
! the rayleigh scale   height =             = hr (8km)
! the rayleigh optical depth  at plane level        = trmoyp
! the altitude of the plane                 = palt
! the rayleigh optical depth for total atmos        = trmoy
! the aerosol  optical depth for total atmos        = tamoy
! if not plane observations then ha is equal to 2.0km
! ntp local variable: if ntp=nt     no plane observation selected
!                        ntp=nt-1   plane observation selected

    hr=8.0
!   it's a mixing rayleigh+aerosol
    if(palt.le.900..and.palt.gt.0.0)then
        if (tap.gt.1.e-03) then
            ha=-palt/log(tap/ta)
        else
            ha=2.
        endif
        ntp=nt-1
    else
        ha=2.0
        ntp=nt
    endif

    ta=tamoy
    tr=trmoy
    piz=pizmoy
!
! compute mixing rayleigh, aerosol
! case 1: pure rayleigh
! case 2: pure aerosol
! case 3: mixing rayleigh-aerosol
!
    if((ta.le.acu2).and.(tr.gt.ta)) then
        do j=0,ntp
            h(j)=j*tr/ntp
            ydel(j)=1.0
            xdel(j)=0.0
        enddo
    endif
    if((tr.le.acu2).and.(ta.gt.tr)) then
        do j=0,ntp
            h(j)=j*ta/ntp
            ydel(j)=0.0
            xdel(j)=piz
        enddo
    endif

    if(tr.gt.acu2.and.ta.gt.acu2.and.iaer_prof.eq.0)then
        ydel(0)=1.0
        xdel(0)=0.0
        h(0)=0.
        altc(0)=300.
        zx=300.
        iplane=0
        do it=0,ntp
            if (it.eq.0) then
                yy=0.
                dd=0.
                goto 111
            endif
            yy=h(it-1)
            dd=ydel(it-1)
 111        ppp2=300.0
            ppp1=0.0
            itp=it
            call discre(ta,ha,tr,hr,itp,ntp,yy,dd,ppp2,ppp1,zx)
            if(ier)return
            xxx=-zx/ha
            if (xxx.lt.-18) then
                ca=0.
            else
                ca=ta*dexp(xxx)
            endif
            xxx=-zx/hr
            cr=tr*dexp(xxx)
            h(it)=cr+ca
            altc(it)=zx
            cr=cr/hr
            ca=ca/ha
            ratio=cr/(cr+ca)
            xdel(it)=(1.e+00-ratio)*piz
            ydel(it)=ratio
        enddo
    endif

    if(tr.gt.acu2.and.ta.gt.acu2.and.iaer_prof.eq.1)then
        call aero_prof(ta,piz,tr,hr,ntp,xmus,h,ch,ydel,xdel,altc)
    endif

! update plane layer if necessary
    if (ntp.eq.(nt-1)) then
! compute position of the plane layer
        taup=tap+trp
        iplane=-1
        do i=0,ntp
            if (taup.ge.h(i)) iplane=i
        enddo

! update the layer from the end to the position to update if necessary
        th=0.005
        xt1=abs(h(iplane)-taup)
        xt2=abs(h(iplane+1)-taup)
        if ((xt1.gt.th).and.(xt2.gt.th)) then
            do i=nt,iplane+1,-1
                xdel(i)=xdel(i-1)
                ydel(i)=ydel(i-1)
                h(i)=h(i-1)
                altc(i)=altc(i-1)
            enddo
        else
            nt=ntp
            if (xt2.lt.xt1) iplane=iplane+1
        endif
        h(iplane)=taup
        if ( tr.gt.acu2.and.ta.gt.acu2) then
            ca=ta*exp(-palt/ha)
            cr=tr*exp(-palt/hr)
            cr=cr/hr
            ca=ca/ha
            ratio=cr/(cr+ca)
            xdel(iplane)=(1.e+00-ratio)*piz
            ydel(iplane)=ratio
            altc(iplane)=palt
        endif
        if ( tr.gt.acu2.and.ta.le.acu2) then
            ydel(iplane)=1.
            xdel(iplane)=0.
            altc(iplane)=palt
        endif
        if ( tr.le.acu2.and.ta.gt.acu2) then
            ydel(iplane)=0.
            xdel(iplane)=1.*piz
            altc(iplane)=palt
        endif
    endif

    aaaa=delta/(2-delta)
    ron=(1-aaaa)/(1+2*aaaa)
!
!     rayleigh phase function
!
    beta0=1.
    beta2=0.5*ron
!
!    primary scattering
!
    ig=1
    tavion0=0.
    tavion1=0.
    tavion2=0.
    tavion=0.
    do j=-mu,mu
        i3(j)=0.
    enddo
!
!   kernel computations
!
    call kernel(0,mu,rm,xpl,psl,bp)
    do j=-mu,mu
        do k=0,nt
            i2(k,j)=0.0000
        enddo
    enddo
!
!   vertical integration, primary upward radiation
!
    do k=1,mu
        i1(nt,k)=1.0
        zi1=i1(nt,k)
        yy=rm(k)
        do i=nt-1,0,-1
            i1(i,k)=exp(-(ta+tr-h(i))/yy)
        enddo
    enddo


    do k=-mu,-1
        do i=0,nt
            i1(i,k)=0.00
        enddo
    enddo

    do 20 k=-mu,mu
    if(k) 21,20,23
21  index=nt
    go to 25
23  index=0
25  continue
    inm1(k)=i1(index,k)
    inm2(k)=i1(index,k)
    i3(k)=i1(index,k)
20  continue
    tavion=i1(iplane,mu)
    tavion2=i1(iplane,mu)

403 ig=ig+1

    do k=1,mu
        xpk=xpl(k)
        ypk=xpl(-k)
        do i=0,nt
            ii1=0.
            ii2=0.
            x=xdel(i)
            y=ydel(i)
            do j=1,mu
                xpj=xpl(j)
                z=gb(j)
                xi1=i1(i,j)
                xi2=i1(i,-j)
                bpjk=bp(j,k)*x+y*(beta0+beta2*xpj*xpk)
                bpjmk=bp(j,-k)*x+y*(beta0+beta2*xpj*ypk)
                ii2=ii2+z*(xi1*bpjk+xi2*bpjmk)
                ii1=ii1+z*(xi1*bpjmk+xi2*bpjk)
            enddo
            i2(i,k)=ii2
            i2(i,-k)=ii1
        enddo
    enddo

    do k=1,mu
        i1(nt,k)=0.0
        zi1=i1(nt,k)
        yy=rm(k)
        do i=nt-1,0,-1
            jj=i+1
            f=h(jj)-h(i)
            a=(i2(jj,k)-i2(i,k))/f
            b=i2(i,k)-a*h(i)
            c=exp(-f/yy)
            d=1.e+00-c
            xx=h(i)-h(jj)*c
            zi1=c*zi1+(d*(b+a*yy)+a*xx)*0.5e+00
            i1(i,k)=zi1
        enddo
    enddo

    do k=-mu,-1
        i1(0,k)=0.
        zi1=i1(0,k)
        yy=rm(k)
        do i=1,nt
            jj=i-1
            f=h(i)-h(jj)
            c=exp(f/yy)
            d=1.e+00-c
            a=(i2(i,k)-i2(jj,k))/f
            b=i2(i,k)-a*h(i)
            xx=h(i)-h(jj)*c
            zi1=c*zi1+(d*(b+a*yy)+a*xx)*0.5e+00
            i1(i,k)=zi1
        enddo
    enddo

    do 30 k=-mu,mu
        if(k) 31,30,33
31      index=nt
        go to 34
33      index=0
34      continue
        in(k)=i1(index,k)
30  continue
    tavion0=i1(iplane,mu)

    if(ig.gt.2) then
        z=0.
        a1=tavion2
        d1=tavion1
        g1=tavion0
        if(a1.ge.acu.and.d1.ge.acu.and.tavion.ge.acu) then
            y=((g1/d1-d1/a1)/((1.-g1/d1)**2)*(g1/tavion))
            y=abs(y)
            z=max(y,z)
        endif
        do 89 l=-mu,mu
            if (l.eq.0) goto 89
            a1=inm2(l)
            d1=inm1(l)
            g1=in(l)
            if(a1.eq.0.) go to 89
            if(d1.eq.0.) go to 89
            if(i3(l).eq.0.) go to 89
            y=((g1/d1-d1/a1)/((1-g1/d1)**2)*(g1/i3(l)))
            y=abs(y)
            z=max(y,z)
89      continue
        if(z.lt.0.0001) then
            do 506 l=-mu,mu
                if (l.eq.0) goto 506
                y1=1.
                d1=inm1(l)
                g1=in(l)
                if(d1.eq.0.0) go to 506
                y1=1-g1/d1
                g1=g1/y1
                i3(l)=i3(l)+g1
506         continue
            d1=tavion1
            g1=tavion0
            y1=1.
            if (d1.ge.acu) then
                if (abs(g1-d1).ge.acu) then
                    y1=1.-g1/d1
                    g1=g1/y1
                endif
                tavion=tavion+g1
            endif
            go to 405
        endif
        do k=-mu,mu
            inm2(k)=inm1(k)
        enddo
        tavion2=tavion1
    endif

    do k=-mu,mu
        inm1(k)=in(k)
    enddo
    tavion1=tavion0

    do l=-mu,mu
      i3(l)=i3(l)+in(l)
    enddo
    tavion=tavion+tavion0

    z=0.
    do l=-mu,mu
        if(i3(l).ne.0)then
            y=abs(in(l)/i3(l))
            z=max(z,y)
        endif
    enddo
    if(z.lt.0.00001) go to 405

    if(ig-igmax) 403,403,405
405 continue

    xf(1)=xf(1)+i3(mu)
    xf(-1)=tavion
    do k=1,mu
        xf(0)=xf(0)+rm(k)*gb(k)*i3(-k)
    enddo
    nt=snt
    return
end
