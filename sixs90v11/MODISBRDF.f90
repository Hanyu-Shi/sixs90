subroutine modisbrdf(p1,p2,p3,mu,np,rm,rp,brdfint)
    implicit none
    real(8) :: p1,p2,p3
    real(8) :: rts,rtv,rfi,rpha
    real(8) :: cts,ctv,cfi,cpha
    real(8) :: sts,stv,sfi
    real(8) :: tanti,tantv
    real(8) :: cost,sint,tvar
    real(8) :: rossthick,rosselt,lispars
    real(8) :: angdist,angtemp,angover
    integer :: mu,np,k,j
    real(8) :: rm(-mu:mu),rp(np),brdfint(-mu:mu,np),pi
    rts=acos(rm(0))
    pi=atan(1.d0)*4.d0

    do k=1,np
        do j=1,mu
            rtv=acos(rm(j))
            if (j.eq.mu) then
                rfi=rm(-mu)
            else
                rfi=rp(k)+rm(-mu)
            endif
            rfi=abs(rfi)
            cts=cos(rts)
            ctv=cos(rtv)
            sts=sin(rts)
            stv=sin(rtv)
            cfi=cos(rfi)
            sfi=sin(rfi)
            cpha=cts*ctv+sts*stv*cfi
            rpha=acos(cpha)

            rosselt=(pi/2-rpha)*cpha+sin(rpha)
            rossthick=(rosselt/(cts+ctv))-pi/4.d0

            tanti=tan(rts)
            tantv=tan(rtv)

            angdist=tanti*tanti+tantv*tantv-2.*tanti*tantv*cfi
            angdist=sqrt(angdist)

            angtemp=1.d0/cts+1.d0/ctv
            cost=2.d0*sqrt(angdist*angdist+tanti*tanti*tantv*tantv*sfi*sfi)
            cost=cost/angtemp
            if (cost.ge.1.d0) cost=1.d0
            if (cost.le.-1.d0) cost=-1.d0
            tvar=acos(cost)
            sint=sqrt(1.-cost*cost)
            angover=(tvar-sint*cost)*angtemp/pi
            lispars=angover-angtemp+0.5*(1.+cpha)/cts/ctv
            brdfint(j,k)=p1+p2*rossthick+p3*lispars
        enddo
    enddo
    return
end
