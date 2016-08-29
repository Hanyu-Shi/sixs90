subroutine equivwl(iinf,isup,step,wlmoy)
    implicit none
    common /sixs_ffu/s(1501),wlinf,wlsup
    real(8) :: step,wlmoy,s,wlinf,wlsup,seb,wlwave,sbor,wl,swl,coef
    integer :: iinf,isup,l

    seb=0.
    wlwave=0.

    do l=iinf,isup
        sbor=s(l)
        if(l.eq.iinf.or.l.eq.isup) sbor=sbor*0.5
        wl=.25+(l-1)*step
        call solirr(wl,swl)
        coef=sbor*step*swl
        seb=seb+coef
        wlwave=wlwave+wl*coef
    enddo
    wlmoy=wlwave/seb
    return
end
