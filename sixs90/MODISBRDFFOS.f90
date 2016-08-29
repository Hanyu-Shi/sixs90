subroutine modisbrdffos(p1, p2, p3, mu, rm, rosur, wfisur, fisur)
  use paramdef
  implicit none
  real(8) ::  p1, p2, p3!, xmu, view
  !real(8) ::  dts, dtv, dfs, dfv, dfi
  real(8) ::  rts, rtv, rfi, rpha
  !real(8) ::  rts, rtv, rfs, rfv, rfi, rpha
  real(8) ::  cts, ctv, cfi, cpha!, ct0
  real(8) ::  sts, stv, sfi
  real(8) ::  tanti, tantv
  real(8) ::  cost, sint, tvar
  real(8) ::  rossthick, rosselt, lispars
  real(8) ::  angdist, angtemp, angover
  integer ::  mu, k, j !,np
  real(8) ::  rm(-mu:mu)
  real(8) ::  rosur(0:mu_p, mu_p, 83), fisur(83), wfisur(83), pisp!, psip
  real(8) ::  mu1, mu2

  real(8) :: fac,pi
  integer :: i

  pisp = acos(-1.)
  pi = atan(1.)*4.
  fac = 180./pi
!      mu=mu_p
  call gauss(0.d0, pisp, fisur, wfisur, 83)
  do i = 0, mu
    do j = 1, mu
      do k = 1, 83
        if (i==0) then
          mu1 = rm(0)
        else
          mu1 = rm(i)
        end if
        mu2 = rm(j)
        rfi = fisur(k) + pisp
        rts = acos(mu1)
        rtv = acos(mu2)
        if ((rts*180./pi)>75.) rts = 75.*pi/180.
        if ((rtv*180./pi)>65.) rtv = 65.*pi/180.
        rfi = abs(rfi)
        cts = cos(rts)
        ctv = cos(rtv)
        sts = sin(rts)
        stv = sin(rtv)
        cfi = cos(rfi)
        sfi = sin(rfi)
        cpha = cts*ctv + sts*stv*cfi
        if (cpha>1.0) cpha = 1.0
        if (cpha<-1.0) cpha = -1.0
        rpha = acos(cpha)

        rosselt = (pi/2-rpha)*cpha + sin(rpha)
        rossthick = (rosselt/(cts+ctv)) - pi/4.

        tanti = tan(rts)
        tantv = tan(rtv)

        angdist = tanti*tanti + tantv*tantv - 2.*tanti*tantv*cfi
        angdist = sqrt(angdist)

        angtemp = 1./cts + 1./ctv
        cost = 2.*sqrt(angdist*angdist+tanti*tanti*tantv*tantv*sfi*sfi)
        cost = cost/angtemp
        if (cost>=1.) cost = 1.
        if (cost<=-1.) cost = -1.
        tvar = acos(cost)
        sint = sqrt(1.-cost*cost)
        angover = (tvar-sint*cost)*angtemp/pi
        lispars = angover - angtemp + 0.5*(1.+cpha)/cts/ctv

        rosur(i, j, k) = p1 + p2*rossthick + p3*lispars
!      write(6,*) "modisbrffos ",acos(mu1)*fac,acos(mu2)*fac,rfi*fac,
!     &  rosur(i,j,k)
        if (rosur(i,j,k)<0.) then
          rosur(i, j, k) = 0.
          write (6, *)
        end if
      end do
    end do
  end do
  return
end subroutine modisbrdffos
