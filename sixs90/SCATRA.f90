subroutine scatra(iaer_prof,taer,taerp,tray,trayp,    &
                  piza,palt,nt,mu,rm,gb,xmus,xmuv,    &
                  ddirtt,ddiftt,udirtt,udiftt,sphalbt,&
                  ddirtr,ddiftr,udirtr,udiftr,sphalbr,&
                  ddirta,ddifta,udirta,udifta,sphalba)
    implicit none
    integer :: mu
    real(8) :: rm(-mu:mu),gb(-mu:mu)
!   computations of the direct and diffuse transmittances
!   for downward and upward paths , and spherical albedo
    real(8) :: xtrans(-1:1)
    real(8) :: taer,taerp,tray,trayp,piza,palt,xmus,xmuv
    real(8) :: udiftt,sphalbt,ddirtr,ddiftr,udirtr,udiftr,sphalbr
    real(8) :: ddirtt,ddiftt,udirtt,ddirta,ddifta,udirta,udifta
    real(8) :: sphalba,tamol,tamolp
    integer :: nt,it,iaer_prof

    ddirtt=1.
    ddiftt=0.
    udirtt=1.
    udiftt=0.
    sphalbt=0.
    ddirtr=1.
    ddiftr=0.
    udirtr=1.
    udiftr=0.
    sphalbr=0.
    ddirta=1.
    ddifta=0.
    udirta=1.
    udifta=0.
    sphalba=0.

    do 1 it=1,3
! it=1 rayleigh only, it=2 aerosol only, it=3 rayleigh+aerosol
        if (it.eq.2.and.taer.le.0.) goto 1
! compute upward,downward diffuse transmittance for rayleigh,aerosol
        if (it.eq.1) then
            if (palt.gt.900) then
                udiftt=(2./3.+xmuv)+(2./3.-xmuv)*exp(-tray/xmuv)
                udiftt=udiftt/((4./3.)+tray)-exp(-tray/xmuv)
                ddiftt=(2./3.+xmus)+(2./3.-xmus)*exp(-tray/xmus)
                ddiftt=ddiftt/((4./3.)+tray)-exp(-tray/xmus)
                ddirtt=exp(-tray/xmus)
                udirtt=exp(-tray/xmuv)
                call csalbr(tray,sphalbt)
            endif
            if (palt.lt.900) then
                tamol=0.
                tamolp=0.
                rm(-mu)=-xmuv
                rm(mu)=xmuv
                rm(0)=xmus
                call iso(iaer_prof,tamol,tray,piza,tamolp,trayp,&
                            palt,nt,mu,rm,gb,xtrans)
                udiftt=xtrans(-1)-exp(-trayp/xmuv)
                udirtt=exp(-trayp/xmuv)
                rm(-mu)=-xmus
                rm(mu)=xmus
                rm(0)=xmus
                ddiftt=(2./3.+xmus)+(2./3.-xmus)*exp(-tray/xmus)
                ddiftt=ddiftt/((4./3.)+tray)-exp(-tray/xmus)
                ddirtt=exp(-tray/xmus)
                udirtt=exp(-tray/xmuv)
                call csalbr(tray,sphalbt)
            endif
            if (palt.le.0.) then
                udiftt=0.
                udirtt=1.
            endif
        endif
        if (it.eq.2) then
            tamol=0.
            tamolp=0.
            rm(-mu)=-xmuv
            rm(mu)=xmuv
            rm(0)=xmus
            call iso(iaer_prof,taer,tamol,piza,taerp,tamolp,&
                        palt,nt,mu,rm,gb,xtrans)
            udiftt=xtrans(-1)-exp(-taerp/xmuv)
            udirtt=exp(-taerp/xmuv)
            rm(-mu)=-xmus
            rm(mu)=xmus
            rm(0)=xmus
            call iso(iaer_prof,taer,tamol,piza,taerp,tamolp,&
                    999.d0,nt,mu,rm,gb,xtrans)
            ddirtt=exp(-taer/xmus)
            ddiftt=xtrans(1)-exp(-taer/xmus)
            sphalbt=xtrans(0)*2.
            if (palt.le.0.) then
                udiftt=0.
                udirtt=1.
            endif
        endif
        if (it.eq.3) then
            rm(-mu)=-xmuv
            rm(mu)=xmuv
            rm(0)=xmus
            call iso(iaer_prof,taer,tray,piza,taerp,trayp,&
                    palt,nt,mu,rm,gb,xtrans)
            udirtt=exp(-(taerp+trayp)/xmuv)
            udiftt=xtrans(-1)-exp(-(taerp+trayp)/xmuv)
            rm(-mu)=-xmus
            rm(mu)=xmus
            rm(0)=xmus
            call iso(iaer_prof,taer,tray,piza,taerp,trayp,&
                    999.d0,nt,mu,rm,gb,xtrans)
            ddiftt=xtrans(1)-exp(-(taer+tray)/xmus)
            ddirtt=exp(-(taer+tray)/xmus)
            sphalbt=xtrans(0)*2.
            if (palt.le.0.) then
                udiftt=0.
                udirtt=1.
            endif
        endif

        if (it.eq.2) goto 2
        if (it.eq.3) goto 1
        ddirtr=ddirtt
        ddiftr=ddiftt
        udirtr=udirtt
        udiftr=udiftt
        sphalbr=sphalbt
        goto 1
2       ddirta=ddirtt
        ddifta=ddiftt
        udirta=udirtt
        udifta=udiftt
        sphalba=sphalbt
1       continue
    return
end
