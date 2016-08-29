subroutine rlmaignanbrdf(p1, p2, p3, mu, np, rm, rp, brdfint)
  implicit none
  real(8) :: p1, p2, p3!, xmu, view
  !real(8) :: dts, dtv, dfs, dfv, dfi
  real(8) :: rts, rtv, rfi, rpha
  !real(8) :: rts, rtv, rfs, rfv, rfi, rpha
  real(8) :: cts, ctv, cfi, cpha, ct0
  real(8) :: sts, stv, sfi
  real(8) :: tanti, tantv
  real(8) :: cost, sint, tvar
  real(8) :: rossthick, rosselt, lispars
  real(8) :: angdist, angtemp, angover
  integer :: mu, np, k, j
  real(8) :: rm(-mu:mu), rp(np), brdfint(-mu:mu, np),pi
  rts = acos(rm(0))
  pi = atan(1.)*4.
  if ((rts*180./pi)>75.) rts = 75.*pi/180.
  do k = 1, np
    do j = 1, mu
      rtv = acos(rm(j))
      if ((rtv*180./pi)>65.) rtv = 65.*pi/180.
      if (j==mu) then
        rfi = rm(-mu)
      else
        rfi = rp(k) + rm(-mu)
      end if
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

      ct0 = 2./(3.*pi)
      rosselt = ct0*(((pi-2.*rpha)*cpha+2.*sin(rpha))/(cts+ctv))
      rossthick = rosselt*(1.+1./(1.+rpha/(1.5*pi/180.))) - (1./3.)

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
!      write(6,*) rossthick,rosselt,rpha,cpha
      brdfint(j, k) = p1 + p2*rossthick + p3*lispars
      if (brdfint(j,k)<0.) then
        brdfint(j, k) = 0.
        write (6, *)
      end if
    end do
  end do
  return
end subroutine rlmaignanbrdf
