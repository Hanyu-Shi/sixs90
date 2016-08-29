subroutine specinterp(wl,taer55,taer55p,tamoy,tamoyp,pizmoy,pizmoyp,ipol)

! - to vary the number of quadratures
    use paramdef
    implicit none
    integer :: nquad
    common /num_quad/ nquad
    real(8) :: phasel,qhasel,uhasel
    common /sixs_phase/ phasel(20,nqmax_p),qhasel(20,nqmax_p),uhasel(20,nqmax_p)
    real(8) :: pha,qha,uha,alphal,betal,gammal,zetal
    common /sixs_polar/ pha(nqmax_p),qha(nqmax_p),uha(nqmax_p), &
            alphal(0:nqmax_p),betal(0:nqmax_p),gammal(0:nqmax_p),  &
            zetal(0:nqmax_p)
    integer :: nbmu
! - to vary the number of quadratures

    real(8) :: wl,taer55,taer55p,tamoy,tamoyp,pizmoy,pizmoyp,roatm
    real(8) :: dtdir,dtdif,utdir,utdif,sphal,wldis,trayl,traypl
    real(8) :: ext,ome,gasym,phase,coef,rpatm,dpatm
    real(8) :: qhase,uhase
    real(8) :: wlinf,alphaa,betaa,tsca,coeff
    integer :: linf,ll,lsup,k,ipol,ipol0

    common /sixs_disc/ roatm(3,20),dtdir(3,20),dtdif(3,20), &
     utdir(3,20),utdif(3,20),sphal(3,20),wldis(20),trayl(20), &
     traypl(20),rpatm(3,20),dpatm(3,20)
    common /sixs_aer/ext(20),ome(20),gasym(20),phase(20),qhase(20),uhase(20)

    real(8) :: test1,test2,test3,coefl

    nbmu=nquad
    linf=1
    do ll=1,19
      if(wl.ge.wldis(ll).and.wl.le.wldis(ll+1)) linf=ll
    enddo
    if(wl.gt.wldis(20)) linf=19
    lsup=linf+1
    coef=log(wldis(lsup)/wldis(linf))
    coefl=(wl-wldis(linf))/(wldis(lsup)-wldis(linf))
    wlinf=wldis(linf)
    alphaa=log(ext(lsup)*ome(lsup)/(ext(linf)*ome(linf)))/coef
    betaa=ext(linf)*ome(linf)/(wlinf**(alphaa))
    tsca=taer55*betaa*(wl**alphaa)/ext(8)
    alphaa=log(ext(lsup)/(ext(linf)))/coef
    betaa=ext(linf)/(wlinf**(alphaa))
    tamoy=taer55*betaa*(wl**alphaa)/ext(8)
    tamoyp=taer55p*betaa*(wl**alphaa)/ext(8)
    pizmoy=tsca/tamoy
    pizmoyp=pizmoy
    do k=1,nbmu
        alphaa=log(phasel(lsup,k)/phasel(linf,k))/coef
        betaa=phasel(linf,k)/(wlinf**(alphaa))
        pha(k)=betaa*(wl**alphaa)
        if (ipol.ne.0)then
            test1=qhasel(linf,k)
            test2=qhasel(lsup,k)
            test3=qhasel(lsup,k)*qhasel(linf,k)
            if((test1.lt.0.001).or.(test2.lt.0.001).or.(test3.lt.0.0)) then
                qha(k)=qhasel(linf,k)+(qhasel(lsup,k)-qhasel(linf,k))*coefl
            else
                alphaa=log(qhasel(lsup,k)/qhasel(linf,k))/coef
                betaa=qhasel(linf,k)/(wlinf**(alphaa))
                qha(k)=betaa*(wl**alphaa)
            endif

            test1=uhasel(linf,k)
            test2=uhasel(lsup,k)
            test3=uhasel(lsup,k)*qhasel(linf,k)
            if((test1.lt.0.001).or.(test2.lt.0.001).or.(test3.lt.0.0)) then
                uha(k)=uhasel(linf,k)+(uhasel(lsup,k)-uhasel(linf,k))*coefl
            else
                alphaa=log(uhasel(lsup,k)/uhasel(linf,k))/coef
                betaa=uhasel(linf,k)/(wlinf**(alphaa))
                uha(k)=betaa*(wl**alphaa)
            endif
        endif
    enddo
! here we don't need coefficients for computation of the polarization
    ipol0=0
!      write(6,*) "tamoy ",tamoy
    call trunca(coeff,ipol0)
    tamoy=tamoy*(1.-pizmoy*coeff)
    tamoyp=tamoyp*(1.-pizmoyp*coeff)
    pizmoy=pizmoy*(1.-coeff)/(1.-pizmoy*coeff)
    return
end
