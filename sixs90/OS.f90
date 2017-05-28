subroutine os (iaer_prof,tamoy,trmoy,pizmoy,tamoyp,trmoyp,palt, &
               phirad,nt,mu,np,rm,gb,rp,xl,xlphim,nfi,rolut)

! - to vary the number of quadratures
    use paramdef
    implicit none
    integer :: nquad,ifi
    common /num_quad/ nquad
    real(8) :: pha,qha,uha,alphal,betal,gammal,zetal
    common /sixs_polar/ pha(nqmax_p),qha(nqmax_p),uha(nqmax_p),&
            alphal(0:nqmax_p),betal(0:nqmax_p),gammal(0:nqmax_p),zetal(0:nqmax_p)
    integer :: nbmu
! - to vary the number of quadratures

!  dimension for gauss integration
    integer :: mu,np,nfi
    integer :: snt
    integer :: nt,iwr,iplane,mum1,ntp,j,it,itp,i,l,m,iborm
    integer :: is,isp,ig,k,jj,index
    integer :: igmax,iaer_prof
    logical :: ier
    integer :: nfilut(mu),nbisca
    real(8) :: rm(-mu:mu),gb(-mu:mu),rp(np),xlphim(nfi)
!  dimension for os computation
    real(8) ::  xl(-mu:mu,np)
! array for sos computation
    real(8) :: xpl(-mu:mu),bp(0:mu,-mu:mu),                    &
        xdel(0:nt),ydel(0:nt),ch(0:nt),h(0:nt),altc(0:nt)
    real(8) :: i1(0:nt,-mu:mu),i2(0:nt,-mu:mu),i3(-mu:mu),     &
        i4(-mu:mu),in(-mu:mu),inm1(-mu:mu),inm2(-mu:mu)

!CCC Begin Variable for Look up table generation
! azimuth or scattering angle variable for LUT computation (rolut)
! azimuth table for look up table computation (filut), nb of fi in each case (nfilut)
    real(8) :: luttv,lutmuv,iscama,iscami,its,scaa,cscaa
    real(8) :: rolut(mu,61),filut(mu,61)
    real(8) :: psl(-1:nqmax_p,-mu:mu)
!CCC End Variable for Look up table generation
    real(8) :: tamoy,trmoy,pizmoy
    real(8) :: tamoyp,trmoyp,palt,phirad
    real(8) :: delta,sigma
    real(8) :: hr,ta,tr,trp
    real(8) :: tap,piz,accu,accu2,ha,xmus,zx,yy,dd
    real(8) :: ppp2,ppp1,ca,cr,ratio
    real(8) :: taup,th,xt1,xt2,pi,phi,aaaa,ron
    real(8) :: roavion1,roavion2,roavion,spl,sa1
    real(8) :: beta0,beta2,roavion0
    real(8) :: sa2,c,zi1,f,d,xpk,y
    real(8) :: a1,d1,g1,y1,delta0s
    real(8) :: xx,xdb,bpjk,bpjmk,z,xi1,xi2,x,xpj,ypk,a,b,ii1,ii2
    real(8) :: phimul,cfi

    common/sixs_del/delta,sigma
    common/sixs_ier/iwr,ier
    common /multorder/ igmax

    nbmu=nquad
    snt=nt
    hr=8.d0
    ta=tamoy
    tr=trmoy
    trp=trmoy-trmoyp
    tap=tamoy-tamoyp
    piz=pizmoy
    iplane=0
    accu=1.d-20
    accu2=1.d-3
    mum1=mu-1
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
!     it's a mixing rayleigh+aerosol
    if(palt.le.900..and.palt.gt.0.0) then
        if (tap.gt.1.d-03) then
            ha=-palt/log(tap/ta)
        else
            ha=2.d0
        endif
        ntp=nt-1
    else
        ha=2.d0
        ntp=nt
    endif

    xmus=-rm(0)
!
! compute mixing rayleigh, aerosol
! case 1: pure rayleigh
! case 2: pure aerosol
! case 3: mixing rayleigh-aerosol
!
    if((ta.le.accu2).and.(tr.gt.ta)) then
        do j=0,ntp
            h(j)=j*tr/ntp
            ch(j)=exp(-h(j)/xmus)/2.d0
            ydel(j)=1.d0
            xdel(j)=0.d0
            if (j.eq.0) then
                altc(j)=300.d0
            else
                altc(j)=-log(h(j)/tr)*hr
            endif
        enddo
    endif
    if((tr.le.accu2).and.(ta.gt.tr)) then
        do j=0,ntp
            h(j)=j*ta/ntp
            ch(j)=exp(-h(j)/xmus)/2.d0
            ydel(j)=0.d0
            xdel(j)=piz
            if (j.eq.0) then
                altc(j)=300.d0
            else
                altc(j)=-log(h(j)/ta)*ha
            endif
        enddo
    endif

    if(tr.gt.accu2.and.ta.gt.accu2.and.iaer_prof.eq.0)then
        ydel(0)=1.d0
        xdel(0)=0.d0
        h(0)=0.d0
        ch(0)=0.5d0
        altc(0)=300.d0
        zx=300.d0
        iplane=0
        do it=0,ntp
            if(it.ne.0) then
                yy=h(it-1)
                dd=ydel(it-1)
            else
                yy=0.d0
                dd=0.d0
            endif
            ppp2=300.d0
            ppp1=0.d0
            itp=it
            call discre(ta,ha,tr,hr,itp,ntp,yy,dd,ppp2,ppp1,zx)
            if(ier) return
            xx=-zx/ha
            if (xx.le.-20) then
                ca=0.d0
            else
                ca=ta*dexp(xx)
            endif
            xx=-zx/hr
            cr=tr*dexp(xx)
            h(it)=cr+ca
            altc(it)=zx
            ch(it)=exp(-h(it)/xmus)/2.d0
            cr=cr/hr
            ca=ca/ha
            ratio=cr/(cr+ca)
            xdel(it)=(1.d+00-ratio)*piz
            ydel(it)=ratio
        enddo
    endif

    if(tr.gt.accu2.and.ta.gt.accu2.and.iaer_prof.eq.1)then
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
        th=0.0005d0
        xt1=abs(h(iplane)-taup)
        xt2=abs(h(iplane+1)-taup)
        if ((xt1.gt.th).and.(xt2.gt.th)) then
            do i=nt,iplane+1,-1
                xdel(i)=xdel(i-1)
                ydel(i)=ydel(i-1)
                h(i)=h(i-1)
                altc(i)=altc(i-1)
                ch(i)=ch(i-1)
            enddo
        else
            nt=ntp
            if (xt2.lt.xt1) iplane=iplane+1
        endif
        h(iplane)=taup
        if ( tr.gt.accu2.and.ta.gt.accu2) then
            ca=ta*exp(-palt/ha)
            cr=tr*exp(-palt/hr)
            h(iplane)=ca+cr
            cr=cr/hr
            ca=ca/ha
            ratio=cr/(cr+ca)
            xdel(iplane)=(1.d+00-ratio)*piz
            ydel(iplane)=ratio
            altc(iplane)=palt
            ch(iplane)=exp(-h(iplane)/xmus)/2.d0
        endif
        if ( tr.gt.accu2.and.ta.le.accu2) then
            ydel(iplane)=1.d0
            xdel(iplane)=0.d0
            altc(iplane)=palt
        endif
        if ( tr.le.accu2.and.ta.gt.accu2) then
            ydel(iplane)=0.d0
            xdel(iplane)=1.d0*piz
            altc(iplane)=palt
        endif
    endif

    pi=acos(-1.d0)
    phi=phirad
    do l=1,np
        do m=-mu,mu
            xl(m,l)=0.d0
        enddo
    enddo
    do ifi=1,nfi
        xlphim(ifi)=0.d0
    enddo

!CC initialization of look up table variable
    do i=1,mu
        nfilut(i)=0
        do j=1,61
            rolut(i,j)=0.d0
            filut(i,j)=0.d0
        enddo
    enddo
    its=acos(xmus)*180.d0/pi
    do i=1,mu
        lutmuv=rm(i)
        luttv=acos(lutmuv)*180.d0/pi
        iscama=int(180.d0-abs(luttv-its))
        iscami=int(180.d0-(luttv+its))
        nbisca=int((iscama-iscami)/4)+1
        nfilut(i)=nbisca
        filut(i,1)=0.d0
        filut(i,nbisca)=180.d0
        scaa=iscama
        do j=2,nfilut(i)-1
            scaa=scaa-4.d0
            cscaa=cos(scaa*pi/180.d0)
            cfi=-(cscaa+xmus*lutmuv)/(sqrt(1-xmus*xmus)*sqrt(1.-lutmuv*lutmuv))
            filut(i,j)=acos(cfi)*180.d0/pi
        enddo
    enddo

!
!     ************ incident angle mus *******
!
!
    aaaa=delta/(2-delta)
    ron=(1-aaaa)/(1+2*aaaa)
!     rayleigh phase function
!
    beta0=1.d0
    beta2=0.5d0*ron
!
!     fourier decomposition
!
    do j=-mu,mu
        i4(j)=0
    enddo
    iborm=nbmu-3
    if(abs (xmus-1.d0) .lt.1.d-06)iborm=0
    do 16 is=0,iborm
!
!    primary scattering
!
        ig=1
        roavion0=0.d0
        roavion1=0.d0
        roavion2=0.d0
        roavion=0.d0
        do j=-mu,mu
            i3(j)=0.d0
        enddo
!
!     kernel computations
!
        isp=is
        call kernel(isp,mu,rm,xpl,psl,bp)
        if(is.gt.0)beta0=0.d0
        do j=-mu,mu
            if(is.le.2) then
                spl=xpl(0)
                sa1=beta0+beta2*xpl(j)*spl
                sa2=bp(0,j)
            else
                sa2=bp(0,j)
                sa1=0.d0
            endif

!     primary scattering source function at every level within the layer
!
            do k=0,nt
                c=ch(k)
                a=ydel(k)
                b=xdel(k)
                i2(k,j)=c*(sa2*b+sa1*a)
            enddo
        enddo
!
!     vertical integration, primary upward radiation
!
      do k=1,mu
          i1(nt,k)=0.d0
          zi1=i1(nt,k)
          yy=rm(k)
          do i=nt-1,0,-1
              jj=i+1
              f=h(jj)-h(i)
              a=(i2(jj,k)-i2(i,k))/f
              b=i2(i,k)-a*h(i)
              c=exp(-f/yy)
              d=1.0d+00-c
              xx=h(i)-h(jj)*c
              zi1=c*zi1+(d*(b+a*yy)+a*xx)*0.5d+00
              i1(i,k)=zi1
          enddo
      enddo
!
!     vertical integration, primary downward radiation
!
      do k=-mu,-1
          i1(0,k)=0.d0
          zi1=i1(0,k)
          yy=rm(k)
          do i=1,nt
              jj=i-1
              f=h(i)-h(jj)
              c=exp(f/yy)
              d=1.0d+00-c
              a=(i2(i,k)-i2(jj,k))/f
              b=i2(i,k)-a*h(i)
              xx=h(i)-h(jj)*c
              zi1=c*zi1+(d*(b+a*yy)+a*xx)*0.5d+00
              i1(i,k)=zi1
          enddo
      enddo
!
!   inm2 is inialized with scattering computed at n-2
!   i3 is inialized with primary scattering
!
    do k=-mu,mu
        if(k.ne.0) then
            if(k.lt.0) then
                index=nt
            else
                index=0
            endif
            inm1(k)=i1(index,k)
            inm2(k)=i1(index,k)
            i3(k)=i1(index,k)
        endif
    enddo

    roavion2=i1(iplane,mu)
    roavion=i1(iplane,mu)
!
!   loop on successive order
!
17  ig=ig+1
!
!     successive orders
!
!     multiple scattering source function at every level within the laye
!
!     if is < ou = 2 kernels are a mixing of aerosols and molecules kern
!     if is >2 aerosols kernels only
!
    if(is.le.2) then
        do k=1,mu
            xpk=xpl(k)
            ypk=xpl(-k)
            do i=0,nt
                ii1=0.d0
                ii2=0.d0
                x=xdel(i)
                y=ydel(i)
                do j=1,mu
                    xpj=xpl(j)
                    z=gb(j)
                    xi1=i1(i,j)
                    xi2=i1(i,-j)
                    bpjk=bp(j,k)*x+y*(beta0+beta2*xpj*xpk)
                    bpjmk=bp(j,-k)*x+y*(beta0+beta2*xpj*ypk)
                    xdb=z*(xi1*bpjk+xi2*bpjmk)
                    ii2=ii2+xdb
                    xdb=z*(xi1*bpjmk+xi2*bpjk)
                    ii1=ii1+xdb
                enddo
                if (abs(ii2).lt.1.d-30) ii2=0.d0
                if (abs(ii1).lt.1.d-30) ii1=0.d0
                i2(i,k)=ii2
                i2(i,-k)=ii1
            enddo
        enddo
    else
        do k=1,mu
            do i=0,nt
                ii1=0.d0
                ii2=0.d0
                x=xdel(i)
                do j=1,mu
                    z=gb(j)
                    xi1=i1(i,j)
                    xi2=i1(i,-j)
                    bpjk=bp(j,k)*x
                    bpjmk=bp(j,-k)*x
                    xdb=z*(xi1*bpjk+xi2*bpjmk)
                    ii2=ii2+xdb
                    xdb=z*(xi1*bpjmk+xi2*bpjk)
                    ii1=ii1+xdb
                enddo
                if (abs(ii2).lt.1.d-30) ii2=0.d0
                if (abs(ii1).lt.1.d-30) ii1=0.d0
                i2(i,k)=ii2
                i2(i,-k)=ii1
            enddo
        enddo
    endif

!   vertical integration, upward radiation
!
    do k=1,mu
        i1(nt,k)=0.d0
        zi1=i1(nt,k)
        yy=rm(k)
        do i=nt-1,0,-1
            jj=i+1
            f=h(jj)-h(i)
            a=(i2(jj,k)-i2(i,k))/f
            b=i2(i,k)-a*h(i)
            c=exp(-f/yy)
            d=1.d+00-c
            xx=h(i)-h(jj)*c
            zi1=c*zi1+(d*(b+a*yy)+a*xx)*0.5d+00
            if (abs(zi1).le.1.d-20) zi1=0.d0
            i1(i,k)=zi1
        enddo
    enddo
!
!   vertical integration, downward radiation
!
    do k=-mu,-1
        i1(0,k)=0.d0
        zi1=i1(0,k)
        yy=rm(k)
        do i=1,nt
            jj=i-1
            f=h(i)-h(jj)
            c=exp(f/yy)
            d=1.d+00-c
            a=(i2(i,k)-i2(jj,k))/f
            b=i2(i,k)-a*h(i)
            xx=h(i)-h(jj)*c
            zi1=c*zi1+(d*(b+a*yy)+a*xx)*0.5d+00
            if (abs(zi1).le.1.d-20) zi1=0.d0
            i1(i,k)=zi1
        enddo
    enddo
!
!   in is the nieme scattering order
!
    do k=-mu,mu
        if(k.ne.0) then
            if(k.lt.0) then
                index=nt
            else
                index=0
            endif
            in(k)=i1(index,k)
        endif
    enddo
    roavion0=i1(iplane,mu)
!
!   convergence test (geometrical serie)
!
    if(ig.gt.2) then
        z=0.d0
        a1=roavion2
        d1=roavion1
        g1=roavion0
        if(a1.ge.accu.and.d1.ge.accu.and.roavion.ge.accu) then
            y=((g1/d1-d1/a1)/((1-g1/d1)**2)*(g1/roavion))
            y=abs(y)
            z=dmax1(dble(y),z)
            endif
            do l=-mu,mu
                if(l.ne.0) then
                    a1=inm2(l)
                    d1=inm1(l)
                    g1=in(l)
                    if((a1.gt.accu).and.(d1.gt.accu).and.(i3(l).gt.accu))then
                        y=((g1/d1-d1/a1)/((1-g1/d1)**2)*(g1/i3(l)))
                        y=abs(y)
                        z=dmax1(dble(y),z)
                    endif
                endif
            enddo

            if(z.lt.0.0001d0) then
!
!     successful test (geometrical serie)
!
            do l=-mu,mu
                y1=1.d0
                d1=inm1(l)
                g1=in(l)
                if(d1.gt.accu) then
                    y1=1-g1/d1
                    if(abs(g1-d1).gt.accu) then
                        g1=g1/y1
                        i3(l)=i3(l)+g1
                    endif
                endif
            enddo
            d1=roavion1
            g1=roavion0
            y1=1.d0
            if(d1.ge.accu) then
                if(abs(g1-d1).ge.accu) then
                    y1=1-g1/d1
                    g1=g1/y1
                endif
                roavion=roavion+g1
            endif
            go to 18
        endif
        do k=-mu,mu
            inm2(k)=inm1(k)
        enddo
        roavion2=roavion1
    endif

    do k=-mu,mu
        inm1(k)=in(k)
    enddo
    roavion1=roavion0

    do l=-mu,mu
        i3(l)=i3(l)+in(l)
    enddo

    roavion=roavion+roavion0

    z=0.d0
    do l=-mu,mu
        if (abs(i3(l)).ge.accu) then
            y=abs(in(l)/i3(l))
            z=dmax1(z,dble(y))
        endif
    enddo


    if(z.lt.0.00001d0) go to 18
    if(ig-igmax) 17,17,18
18  continue

    delta0s=1
    if(is.ne.0) delta0s=2
    do l=-mu,mu
        i4(l)=i4(l)+delta0s*i3(l)
    enddo

    do l=1,np
        phi=rp(l)
        do m=-mum1,mum1
            if(m.gt.0) then
                xl(m,l)=xl(m,l)+delta0s*i3(m)*cos(is*(phi+pi))
            else
                xl(m,l)=xl(m,l)+delta0s*i3(m)*cos(is*phi)
            endif
        enddo
    enddo

! Look up table generation
    do m=1,mu
        do l=1,nfilut(m)
            phimul=filut(m,l)*pi/180.d0
            rolut(m,l)=rolut(m,l)+delta0s*i3(m)*cos(is*(phimul+pi))
            enddo
    enddo

    if(is.eq.0) then
        do k=1,mum1
            xl(0,1)=xl(0,1)+rm(k)*gb(k)*i3(-k)
        enddo
    endif
    xl(mu,1)=xl(mu,1)+delta0s*i3(mu)*cos(is*(phirad+pi))
    do ifi=1,nfi
        phimul=(ifi-1)*pi/(nfi-1)
        xlphim(ifi)=xlphim(ifi)+delta0s*roavion*cos(is*(phimul+pi))
    enddo
    xl(-mu,1)=xl(-mu,1)+delta0s*roavion*cos(is*(phirad+pi))
    z=0.d0
    do l=-mu,mu
        if(abs(i4(l)).ge.accu) then
            x=abs(i3(l)/i4(l))
            z=dmax1(z,x)
        endif
    enddo

    if(z.gt.0.001d0) go to 16
    goto 19

16  continue
19  continue
    nt=snt
!      write(6,*) "in oS",tamoy,trmoy,xl(-mu,1)/rm(0)
!      write(6,*) 'reflectance ', xl(mu,1)/xmus
    return
end
