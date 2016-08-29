subroutine ossurf(iaer_prof, tamoy, trmoy, pizmoy, tamoyp, trmoyp, &
    palt, phirad, nt, mu, np, rm, gb, rp, rosur, wfisur, fisur, xl, &
    xlphim, nfi, rolut)

! - to vary the number of quadratures
  use paramdef
  implicit none
  integer :: nquad
  common /num_quad/nquad
  real(8) :: pha, qha, uha, alphal, betal, gammal, zetal
  common /sixs_polar/pha(nqmax_p), qha(nqmax_p), uha(nqmax_p), &
    alphal(0:nqmax_p), betal(0:nqmax_p), gammal(0:nqmax_p), zetal(0:nqmax_p)
  integer :: nbmu
! - to vary the number of quadratures


!  dimension for gauss integration
  integer :: mu, np, nfi
  real(8) :: rm(-mu:mu), gb(-mu:mu), rp(np), xlphim(nfi)
!  dimension for os computation
  real(8) :: xl(-mu:mu, np)
! array for sos computation
  real(8) :: xpl(-mu:mu), bp(0:mu, -mu:mu), xdel(0:nt), &
    ydel(0:nt), ch(0:nt), h(0:nt), altc(0:nt)
  real(8) :: i1(0:nt, -mu:mu), i2(0:nt, -mu:mu), i3(-mu:mu), &
    i4(-mu:mu), in(-mu:mu), inm1(-mu:mu), inm2(-mu:mu)

!ccc begin variable for look up table generation
! azimuth or scattering angle variable for lut computation (rolut)
! azimuth table for look up table computation (filut), nb of fi in each case (nfilut)
  real(8) :: luttv, lutmuv, iscama, iscami, scaa, its
  integer :: nbisca!, its  !! HyS
  real(8) :: rolut(mu, 61), filut(mu, 61)
  real(8) :: psl(-1:nqmax_p, -mu:mu)
  integer ::  nfilut(mu)
!ccc end variable for look up table generation


  real(8) :: tamoy, trmoy, pizmoy
  real(8) :: tamoyp, trmoyp, palt, phirad
  real(8) :: delta, sigma
  real(8) :: hr, ta, tr, trp
  real(8) :: tap, piz, accu, accu2, ha, xmus, zx, yy, dd
  real(8) :: ppp2, ppp1, ca, cr, ratio
  real(8) :: taup, th, xt1, xt2, pi, phi, aaaa, ron
  real(8) :: roavion1, roavion2, roavion, spl, sa1
  real(8) :: beta0, beta2, roavion0
  real(8) :: sa2, c, zi1, f, d, xpk, y
  real(8) :: a1, d1, g1, y1, delta0s
  integer :: snt
  integer :: nt, iwr, iplane, mum1, ntp, j, it, itp, i, l, m, iborm
  integer :: is, isp, ig, k, jj, index
  logical :: ier
  real(8) :: xx, xdb, bpjk, bpjmk, z, xi1, xi2, x, xpj, ypk, a, b, ii1, ii2
  integer :: igmax, iaer_prof

!  variable for ground boundary conditions
  real(8) :: rosur(0:mu_p, mu_p, 83), fisur(83), wfisur(83)
  real(8) :: srosur(0:mu_p, mu_p), pisp

  ! hys add
  real(8) :: acu2,cfi,phimul,cscaa
  integer :: ifi


  common /sixs_del/delta, sigma
  common /sixs_ier/iwr, ier
  common /multorder/igmax

!      igmax=1
  nbmu = nquad
! the optical thickness above plane are recomputed to give o.t above pla
!     write(6,*) 'tamoy,trmoy,tamoyp,trmoyp,palt,pizmoy'
!      write(6,*) tamoy,trmoy,tamoyp,trmoyp,palt,pizmoy
!       write(6,*) 'betal 0:80'
!       do i=0,80
!         write(6,*) i,betal(i)
!       enddo
!       write(6,*) 'phase function 83 terms'
!       do i=1,83
!         write(6,*) pha(i)
!       enddo
  snt = nt
  hr = 8.0
  ta = tamoy
  tr = trmoy
  trp = trmoy - trmoyp
  tap = tamoy - tamoyp
  piz = pizmoy
!      print *, 'ta,tr,piz,tap,trp,palt,nt,piz'
!      print *,ta,tr,piz,tap,trp,palt,nt,piz
  pisp = acos(-1.)
  iplane = 0
  accu = 1.e-20
  accu2 = 1.e-3
  mum1 = mu - 1

  ! HyS
  acu2 = accu2

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
  if (palt<=900. .and. palt>0.0) then
    if (tap>1.e-03) then
      ha = -palt/log(tap/ta)
    else
      ha = 2.
    end if
    ntp = nt - 1
  else
    ha = 2.0
    ntp = nt
  end if
!
  xmus = -rm(0)
!
! compute mixing rayleigh, aerosol
! case 1: pure rayleighnedit
! case 2: pure aerosol
! case 3: mixing rayleigh-aerosol
!
!      write(6,*) "rosur ",rosur
!      stop
  if ((ta<=accu2) .and. (tr>ta)) then
    do j = 0, ntp
      h(j) = j*tr/ntp
      ch(j) = exp(-h(j)/xmus)/2.
      ydel(j) = 1.0
      xdel(j) = 0.0
      if (j==0) then
        altc(j) = 300.
      else
        altc(j) = -log(h(j)/tr)*hr
      end if
    end do
  end if
  if ((tr<=accu2) .and. (ta>tr)) then
    do j = 0, ntp
      h(j) = j*ta/ntp
      ch(j) = exp(-h(j)/xmus)/2.
      ydel(j) = 0.0
      xdel(j) = piz
      if (j==0) then
        altc(j) = 300.
      else
        altc(j) = -log(h(j)/ta)*ha
      end if
    end do
  end if
!
  if (tr>accu2 .and. ta>accu2 .and. iaer_prof==0) then
    ydel(0) = 1.0
    xdel(0) = 0.0
    h(0) = 0.
    ch(0) = 0.5
    altc(0) = 300.
    zx = 300.
    iplane = 0
    do it = 0, ntp
      if (it==0) then
        yy = 0.
        dd = 0.
        goto 111
      end if
      yy = h(it-1)
      dd = ydel(it-1)
      111 ppp2 = 300.0
      ppp1 = 0.0
      itp = it
      call discre(ta, ha, tr, hr, itp, ntp, yy, dd, ppp2, ppp1, zx)
      if (ier) return
      xx = -zx/ha
      if (xx<=-20) then
        ca = 0.
      else
        ca = ta*dexp(xx)
      end if
      xx = -zx/hr
      cr = tr*dexp(xx)
      h(it) = cr + ca
      altc(it) = zx
      ch(it) = exp(-h(it)/xmus)/2.
      cr = cr/hr
      ca = ca/ha
      ratio = cr/(cr+ca)
      xdel(it) = (1.e+00-ratio)*piz
      ydel(it) = ratio
!     print *,'discre ',it,cr,ca,xdel(it),ydel(it),zx
    end do
  end if
! comment the aero_prof vermote 11/20/2010

  if (tr>acu2 .and. ta>acu2 .and. iaer_prof==1) then
    call aero_prof(ta, piz, tr, hr, ntp, xmus, h, ch, ydel, xdel, altc)
  end if

! update plane layer if necessary
  if (ntp==(nt-1)) then
! compute position of the plane layer
    taup = tap + trp
    iplane = -1
    do i = 0, ntp
      if (taup>=h(i)) iplane = i
    end do
! update the layer from the end to the position to update if necessary
    th = 0.0005
    xt1 = abs(h(iplane)-taup)
    xt2 = abs(h(iplane+1)-taup)
    if ((xt1>th) .and. (xt2>th)) then
      do i = nt, iplane + 1, -1
        xdel(i) = xdel(i-1)
        ydel(i) = ydel(i-1)
        h(i) = h(i-1)
        altc(i) = altc(i-1)
        ch(i) = ch(i-1)
      end do
    else
      nt = ntp
      if (xt2<xt1) iplane = iplane + 1
    end if
    h(iplane) = taup
    if (tr>accu2 .and. ta>accu2) then
      ca = ta*exp(-palt/ha)
      cr = tr*exp(-palt/hr)
      h(iplane) = ca + cr
      cr = cr/hr
      ca = ca/ha
      ratio = cr/(cr+ca)
      xdel(iplane) = (1.e+00-ratio)*piz
      ydel(iplane) = ratio
      altc(iplane) = palt
      ch(iplane) = exp(-h(iplane)/xmus)/2.
    end if
    if (tr>accu2 .and. ta<=accu2) then
      ydel(iplane) = 1.
      xdel(iplane) = 0.
      altc(iplane) = palt
    end if
    if (tr<=accu2 .and. ta>accu2) then
      ydel(iplane) = 0.
      xdel(iplane) = 1.*piz
      altc(iplane) = palt
    end if
  end if
!
!
!      print *,ha,hr,palt,ta,tr,tap,trp
!      do i=0,nt
!      print *,i,h(i),ch(i),xdel(i),ydel(i),altc(i)
!      enddo
!
  pi = acos(-1.)
  phi = phirad
  do l = 1, np
    do m = -mu, mu
      xl(m, l) = 0.
    end do
  end do
  do ifi = 1, nfi
    xlphim(ifi) = 0.
  end do

!cc initialization of look up table variable
  do i = 1, mu
    nfilut(i) = 0
    do j = 1, (nbmu-1)/2
      rolut(i, j) = 0.
      filut(i, j) = 0.
    end do
  end do
  its = acos(xmus)*180.0/pi
  do i = 1, mu
    lutmuv = rm(i)
    luttv = acos(lutmuv)*180./pi
    iscama = (180-abs(luttv-its))
    iscami = (180-(luttv+its))
    nbisca = int((iscama-iscami)/4) + 1
    nfilut(i) = nbisca
    filut(i, 1) = 0.0
    filut(i, nbisca) = 180.0
    scaa = iscama
    do j = 2, nfilut(i) - 1
      scaa = scaa - 4.0
      cscaa = cos(scaa*pi/180.)
      cfi = -(cscaa+xmus*lutmuv)/(sqrt(1-xmus*xmus)*sqrt(1.-lutmuv*lutmuv))
      filut(i, j) = acos(cfi)*180.0/pi
    end do
  end do
!ccc check initialization  (debug)
!      do i=1,mu
!        lutmuv=rm(i)
!        luttv=acos(lutmuv)*180./pi
!       do j=1,nfilut(i)
!      cscaa=-xmus*lutmuv-cos(filut(i,j)*pi/180.)*sqrt(1.-xmus*xmus)
!    s  *sqrt(1.-lutmuv*lutmuv)
!      scaa=acos(cscaa)*180./pi
!      write(6,*) its,luttv,filut(i,j),scaa
!      enddo
!      enddo
!ccc check initialization  (debug)
!cc end initialization of look up table variable



!
!     ************ incident angle mus *******
!
!
  aaaa = delta/(2-delta)
  ron = (1-aaaa)/(1+2*aaaa)
!     write(6,*) 'ron ',ron
!
!     rayleigh phase function
!
  beta0 = 1.
  beta2 = 0.5*ron
!
!     fourier decomposition
!      write(6,*) "rosur(0,mu,1)",rosur(0,mu,1),rm(0),rm(mu),fisur(1)
!      write(6,*) "rosur",rosur
!      do i=0,mu
!      do j=1,mu
!      do k=1,83
!      rosur(i,j,k)=0.4
!      enddo
!      enddo
!      enddo
!
  do j = -mu, mu
    i4(j) = 0.
  end do
  iborm = nbmu - 3
  if (abs(xmus-1.000000)<1.e-06) iborm = 0
  do is = 0, iborm
!
! compute fourier component of the surface term for is=0
    do i = 0, mu
      do j = 1, mu
        srosur(i, j) = 0.0
        do k = 1, 83
          srosur(i, j) = srosur(i, j) + 2.*rosur(i, j, k)*wfisur(k)*cos(is*fisur(k))
        end do
        srosur(i, j) = srosur(i, j)/pisp
      end do
    end do

!    primary scattering
!
    ig = 1
    roavion0 = 0.
    roavion1 = 0.
    roavion2 = 0.
    roavion = 0.
    do j = -mu, mu
      i3(j) = 0.
    end do
!
!     kernel computations
!
    isp = is
    call kernel(isp, mu, rm, xpl, psl, bp)
    if (is>0) beta0 = 0.
    do j = -mu, mu
      if (is-2) 200, 200, 201
      200 spl = xpl(0)
      sa1 = beta0 + beta2*xpl(j)*spl
      sa2 = bp(0, j)
      goto 202
      201 sa2 = bp(0, j)
      sa1 = 0.
!
!     primary scattering source function at every level within the layer
!
      202 do k = 0, nt
        c = ch(k)
        a = ydel(k)
        b = xdel(k)
        i2(k, j) = c*(sa2*b+sa1*a)
      end do
    end do
!
!     vertical integration, primary upward radiation
!

    do k = 1, mu
      i1(nt, k) = srosur(0, k)*xmus*exp(-h(nt)/xmus)/2.0
      zi1 = i1(nt, k)
      yy = rm(k)
      do i = nt - 1, 0, -1
        jj = i + 1
        f = h(jj) - h(i)
        a = (i2(jj,k)-i2(i,k))/f
        b = i2(i, k) - a*h(i)
        c = exp(-f/yy)
        d = 1.0e+00 - c
        xx = h(i) - h(jj)*c
        zi1 = c*zi1 + (d*(b+a*yy)+a*xx)*0.5e+00
        i1(i, k) = zi1
      end do
    end do
!
!     vertical integration, primary downward radiation
!
    do k = -mu, -1
      i1(0, k) = 0.
      zi1 = i1(0, k)
      yy = rm(k)
      do i = 1, nt
        jj = i - 1
        f = h(i) - h(jj)
        c = exp(f/yy)
        d = 1.0e+00 - c
        a = (i2(i,k)-i2(jj,k))/f
        b = i2(i, k) - a*h(i)
        xx = h(i) - h(jj)*c
        zi1 = c*zi1 + (d*(b+a*yy)+a*xx)*0.5e+00
        i1(i, k) = zi1
      end do
    end do
!
!     inm2 is inialized with scattering computed at n-2
!     i3 is inialized with primary scattering
!
    do k = -mu, mu
      if (k) 21, 20, 23
      21 index = nt
      goto 25
      23 index = 0
      25 continue
      inm1(k) = i1(index, k)
      inm2(k) = i1(index, k)
      i3(k) = i1(index, k)
    20 end do
    roavion2 = i1(iplane, mu)
    roavion = i1(iplane, mu)
!
!     loop on successive order
!
    503 ig = ig + 1
!      write(6,*) 'ig ',ig
!
!     successive orders
!
!     multiple scattering source function at every level within the laye
!
!     if is < ou = 2 kernels are a mixing of aerosols and molecules kern
!     if is >2 aerosols kernels only
!
    if (is-2) 210, 210, 211
    210 do k = 1, mu
      xpk = xpl(k)
      ypk = xpl(-k)
      do i = 0, nt
        ii1 = 0.
        ii2 = 0.
        x = xdel(i)
        y = ydel(i)
        do j = 1, mu
          xpj = xpl(j)
          z = gb(j)
          xi1 = i1(i, j)
          xi2 = i1(i, -j)
          bpjk = bp(j, k)*x + y*(beta0+beta2*xpj*xpk)
          bpjmk = bp(j, -k)*x + y*(beta0+beta2*xpj*ypk)
          xdb = z*(xi1*bpjk+xi2*bpjmk)
          ii2 = ii2 + xdb
          xdb = z*(xi1*bpjmk+xi2*bpjk)
          ii1 = ii1 + xdb
        end do
        if (abs(ii2)<1.e-30) ii2 = 0.
        if (abs(ii1)<1.e-30) ii1 = 0.
        i2(i, k) = ii2
        i2(i, -k) = ii1
      end do
    end do
    goto 213
    211 do k = 1, mu
      do i = 0, nt
        ii1 = 0.
        ii2 = 0.
        x = xdel(i)
        do j = 1, mu
          z = gb(j)
          xi1 = i1(i, j)
          xi2 = i1(i, -j)
          bpjk = bp(j, k)*x
          bpjmk = bp(j, -k)*x
          xdb = z*(xi1*bpjk+xi2*bpjmk)
          ii2 = ii2 + xdb
          xdb = z*(xi1*bpjmk+xi2*bpjk)
          ii1 = ii1 + xdb
        end do
        if (abs(ii2)<1.e-30) ii2 = 0.
        if (abs(ii1)<1.e-30) ii1 = 0.
        i2(i, k) = ii2
        i2(i, -k) = ii1
      end do
    end do
!
!     vertical integration, upward radiation
!
    213 do k = 1, mu
      i1(nt, k) = 0.
!  the surface contribution at the boundary layer
      do l = 1, mu
        i1(nt, k) = i1(nt, k) + 2.*gb(l)*i1(nt, -l)*rm(l)*srosur(l, k)/2.0
      end do
      zi1 = i1(nt, k)
      yy = rm(k)
      do i = nt - 1, 0, -1
        jj = i + 1
        f = h(jj) - h(i)
        a = (i2(jj,k)-i2(i,k))/f
        b = i2(i, k) - a*h(i)
        c = exp(-f/yy)
        d = 1.e+00 - c
        xx = h(i) - h(jj)*c
        zi1 = c*zi1 + (d*(b+a*yy)+a*xx)*0.5e+00
        if (abs(zi1)<=1.e-20) zi1 = 0.
        i1(i, k) = zi1
      end do
    end do
!
!     vertical integration, downward radiation
!
    do k = -mu, -1
      i1(0, k) = 0.
      zi1 = i1(0, k)
      yy = rm(k)
      do i = 1, nt
        jj = i - 1
        f = h(i) - h(jj)
        c = exp(f/yy)
        d = 1.e+00 - c
        a = (i2(i,k)-i2(jj,k))/f
        b = i2(i, k) - a*h(i)
        xx = h(i) - h(jj)*c
        zi1 = c*zi1 + (d*(b+a*yy)+a*xx)*0.5e+00
        if (abs(zi1)<=1.e-20) zi1 = 0.
        i1(i, k) = zi1
      end do
    end do
!
!     in is the nieme scattering order
!
    do k = -mu, mu
      if (k) 31, 30, 33
      31 index = nt
      goto 34
      33 index = 0
      34 continue
      in(k) = i1(index, k)
    30 end do
    roavion0 = i1(iplane, mu)
!
!   convergence test (geometrical serie)
!
    if (ig>2) then
      z = 0.
      a1 = roavion2
      d1 = roavion1
      g1 = roavion0
      if (a1>=accu .and. d1>=accu .and. roavion>=accu) then
        y = ((g1/d1-d1/a1)/((1-g1/d1)**2)*(g1/roavion))
        y = abs(y)
        z = dmax1(dble(y), z)
      end if
      do l = -mu, mu
        if (l==0) goto 99
        a1 = inm2(l)
        d1 = inm1(l)
        g1 = in(l)
        if (a1<=accu) goto 99
        if (d1<=accu) goto 99
        if (i3(l)<=accu) goto 99
        y = ((g1/d1-d1/a1)/((1-g1/d1)**2)*(g1/i3(l)))
        y = abs(y)
        z = dmax1(dble(y), z)
      99 end do
      if (z<0.0001) then
!
!     successful test (geometrical serie)
!
        do l = -mu, mu
          y1 = 1.
          d1 = inm1(l)
          g1 = in(l)
          if (d1<=accu) goto 606
          y1 = 1 - g1/d1
          if (abs(g1-d1)<=accu) then
            goto 606
          end if
          g1 = g1/y1
          i3(l) = i3(l) + g1
        606 end do
        d1 = roavion1
        g1 = roavion0
        y1 = 1.
        if (d1>=accu) then
          if (abs(g1-d1)>=accu) then
            y1 = 1 - g1/d1
            g1 = g1/y1
          end if
          roavion = roavion + g1
        end if
        goto 505
      end if
!
!     inm2 is the (n-2)ieme scattering order
!
      do k = -mu, mu
        inm2(k) = inm1(k)
      end do
      roavion2 = roavion1
    end if
!
!     inm1 is the (n-1)ieme scattering order
!
    do k = -mu, mu
      inm1(k) = in(k)
    end do
    roavion1 = roavion0
!
!     sum of the n-1 orders
!
    do l = -mu, mu
      i3(l) = i3(l) + in(l)
    end do
    roavion = roavion + roavion0
!
!     stop if order n is less than 1% of the sum
!
    z = 0.
    do l = -mu, mu
      if (abs(i3(l))>=accu) then
        y = abs(in(l)/i3(l))
        z = dmax1(z, dble(y))
      end if
    end do

!     if(z.lt.0.00001) go to 505    # 6sv4.0 choice
    if (z<0.00001) goto 505
!
!      stop if order n is greater than 20 in any case
!
    if (ig-igmax) 503, 503, 505
    505 continue
!
!     sum of the fourier component s
!
    delta0s = 1
    if (is/=0) delta0s = 2
    do l = -mu, mu
      i4(l) = i4(l) + delta0s*i3(l)
    end do
!
!     stop of the fourier decomposition
!
    do l = 1, np
      phi = rp(l)
      do m = -mum1, mum1
        if (m>0) then
          xl(m, l) = xl(m, l) + delta0s*i3(m)*cos(is*(phi+pi))
        else
          xl(m, l) = xl(m, l) + delta0s*i3(m)*cos(is*phi)
        end if
      end do
    end do

! look up table generation
    do m = 1, mu
      do l = 1, nfilut(m)
        phimul = filut(m, l)*pi/180.
        rolut(m, l) = rolut(m, l) + delta0s*i3(m)*cos(is*(phimul+pi))
      end do
    end do
! end of look up table generation



    if (is==0) then
      do k = 1, mum1
        xl(0, 1) = xl(0, 1) + rm(k)*gb(k)*i3(-k)
      end do
    end if
    xl(mu, 1) = xl(mu, 1) + delta0s*i3(mu)*cos(is*(phirad+pi))
    do ifi = 1, nfi
      phimul = (ifi-1)*pi/(nfi-1)
      xlphim(ifi) = xlphim(ifi) + delta0s*roavion*cos(is*(phimul+pi))
    end do
    xl(-mu, 1) = xl(-mu, 1) + delta0s*roavion*cos(is*(phirad+pi))
    z = 0.
    do l = -mu, mu
      if (abs(i4(l))<accu) goto 613
      x = abs(i3(l)/i4(l))
      z = dmax1(z, x)
    613 end do

!     if(z.gt.0.001) go to 24     #6sv4.0 choice
    if (z>0.001) goto 24
    goto 243

  24 end do
  243 continue
  nt = snt
!     write(6,*) "in os",tamoy,trmoy,xl(-mu,1)/rm(0)
!      write(6,*) 'reflectance ', xl(mu,1)/xmus
  return
end subroutine ossurf
