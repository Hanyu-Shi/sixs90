subroutine roujbrdf(k0,k1,k2,mu,np,rm,rp,brdfint)
! model can be found in JGR 1992 paper, Vol. 97,No D18, Page20,445-20,468
    implicit none
    integer :: mu,np
    dimension rm(-mu:mu),rp(np),brdfint(-mu:mu,np)
    real(8) :: k0,k1,k2,pi,psi,cpsi
    real(8) :: rm,rp,brdfint,xmus,xmuv,fi,fr,tts,ttv,f2,f1
    integer :: k,j
    xmus=rm(0)
    pi=atan(1.d0)*4.
    do k=1,np
        do j=1,mu
            xmuv=rm(j)
            if (j.eq.mu) then
                fi=rm(-mu)
            else
                fi=rp(k)+rm(-mu)
            endif
            fr=acos(cos(fi))
            tts=tan(acos(xmus))
            ttv=tan(acos(xmuv))
            cpsi=xmus*xmuv+sin(acos(xmus))*sin(acos(xmuv))*cos(fi)
            if (cpsi.lt.1.) then
                psi=acos(cpsi)
            else
                psi=0.
            endif
            f2=4./(3.*pi*(xmus+xmuv))*((pi/2-psi)*cpsi+sin(psi))-1./3.
            f1=0.5*((pi-fr)*cos(fr)+sin(fr))*tts*ttv-tts-ttv-sqrt(   &
            tts*tts+ttv*ttv-2*tts*ttv*cos(fr))
            f1=f1/pi
            brdfint(j,k)=k0+k1*f1+k2*f2
        enddo
    enddo
    return
end
