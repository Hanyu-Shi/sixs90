subroutine POLNAD(xts,xtv,phi,pveg,ropq,ropu)
    implicit none
    real(8) :: xts,xtv,phi,pveg,ropq,ropu
    real(8) :: pi,dtr,csca,sca,alpha,alphap,N,mui,mut
    real(8) :: mus,muv,sinv,cksi,ksi,xf1,xf2
    real(8) :: fpalpha,rpveg,rpsoil,rpsur
    N=1.5
    pi=acos(0.d0)*2.0
    dtr=pi/180.0
    csca=-cos(xts*dtr)*cos(xtv*dtr)-sin(xts*dtr)*sin(xtv*dtr)*cos(phi*dtr)
    sca=acos(csca)
    alpha=(pi-sca)/2.0
    alphap=asin(sin(alpha)/N)
    mui=cos(alpha)
    mut=cos(alphap)

    xf1=(N*mut-mui)/(N*mut+mui)
    xf2=(N*mui-mut)/(N*mui+mut)
    fpalpha=0.5*(xf1*xf1-xf2*xf2)
    rpveg=fpalpha/4./(cos(xts*dtr)+cos(xtv*dtr))
    rpsoil=fpalpha/4./cos(xts*dtr)/cos(xtv*dtr)
    rpsur=rpveg*pveg+rpsoil*(1.-pveg)
    muv=cos(xtv*dtr)
    mus=cos(xts*dtr)
    sinv=sin(xtv*dtr)
    if (xtv.gt.0.5) then
        if (sin(phi*dtr).lt.0) then
            cksi=(muv*csca+mus)/sqrt(1.-csca*csca)/sinv
        else
            cksi=-(muv*csca+mus)/sqrt(1.-csca*csca)/sinv
        endif
    else
        cksi=0.0
    endif

    if (cksi.gt.1.) cksi=1.
    if (cksi.lt.-1.) cksi=-1.
    ksi=acos(cksi)/dtr

    ropq=rpsur*(2.*cksi*cksi-1.)
    ropu=-rpsur*2.*cksi*sqrt(1.-cksi*cksi)
    return
end
