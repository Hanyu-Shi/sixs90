subroutine interp (iaer,idatmp,wl,taer55,taer55p,xmud,              &
       romix,rorayl,roaero,phaa,phar,rqmix,rqrayl,rqaero,qhaa,qhar, &
       rumix,rurayl,ruaero,uhaa,uhar,                               &
       tsca,tray,trayp,taer,taerp,dtott,utott,astot,asray,asaer,    &
       utotr,utota,dtotr,dtota,ipol,                                &
            roatm_fi,romix_fi,rorayl_fi,nfi,                        &
        roluts,rolut,rolutsq,rolutq,rolutsu,rolutu,nfilut)

    use paramdef
    implicit none
    integer :: nfi,i,ifi,j
    real(8) :: test1,test2,test3
    real(8) :: wl,taer55,taer55p
    real(8) :: xmud,romix,rorayl,roaero,phaa,phar,tsca,tray
    real(8) :: rqmix,rqrayl,rqaero,qhaa,qhar,rqatm,qhase
    real(8) :: rumix,rurayl,ruaero,uhaa,uhar,ruatm,uhase
    real(8) :: trayp,taer,taerp,dtott,utott,astot,asray,asaer,utotr
    real(8) :: utota,dtotr,dtota,ext,ome,gasym,phase,roatm,dtdir
    real(8) :: dtdif,utdir,utdif,sphal,wldis,trayl,traypl,delta,sigma
    real(8) :: alphaa,betaa,alphar,betar,alphac,betac,coef,coefl,wlinf
    real(8) :: drinf,drsup,dtinf,dtsup,dtotc,dainf,dasup,urinf,ursup
    real(8) :: utinf,utsup,utotc,uainf,uasup,arinf,arsup,atinf,atsup
    real(8) :: aainf,aasup,depolar1,depolar2
    real(8) :: romix_fi(nfi),rorayl_fi(nfi),roatm_fi(3,20,nfi)
    real(8) :: rolut(mu_p,61),roluts(20,mu_p,61)
    real(8) :: rolutq(mu_p,61),rolutsq(20,mu_p,61)
    real(8) :: rolutu(mu_p,61),rolutsu(20,mu_p,61)
    integer :: nfilut(mu_p),mu
    integer :: iaer,idatmp,linf,ll,lsup,ipol
    real(8) :: ruaer0

    common /sixs_aer/ext(20),ome(20),gasym(20),phase(20),qhase(20),uhase(20)
    common /sixs_disc/ roatm(3,20),dtdir(3,20),dtdif(3,20),  &
    utdir(3,20),utdif(3,20),sphal(3,20),wldis(20),trayl(20),&
    traypl(20),rqatm(3,20),ruatm(3,20)
    common /sixs_del/ delta,sigma


    mu=mu_p
    linf=1
    do ll=1,19
        if(wl.gt.wldis(ll).and.wl.le.wldis(ll+1)) linf=ll
    enddo
    if(wl.gt.wldis(20)) linf=19
    lsup=linf+1

    alphaa=0.
    betaa=0.
    alphar=0.
    betar=0.
    alphac=0.
    betac=0.
    phaa=0.
    qhaa=0.
    uhaa=0.
    phar=0.
    qhar=0.
    uhar=0.
    roaero=0.
    rqaero=0.
    ruaero=0.
    rorayl=0.
    rqrayl=0.
    rurayl=0.
    romix=0.
    rqmix=0.
    rumix=0.
    dtota=1.
    utota=1.
    asaer=0.
    taer=0.
    taerp=0.
    coef=log(wldis(lsup)/wldis(linf))
    wlinf=wldis(linf)

    depolar1=2.*(1.-delta)/(2.+delta)
    depolar2=3.*delta/(2.+delta)

    if(iaer.ne.0) then
        alphaa=log(phase(lsup)/phase(linf))/coef
        betaa=phase(linf)/(wlinf**(alphaa))
        phaa=betaa*(wl**alphaa)
    endif

    phar=depolar1*.75*(1.+xmud*xmud)+depolar2
    if (idatmp.eq.0) goto 2234

    alphar=log(roatm(1,lsup)/roatm(1,linf))/ coef
    betar=roatm(1,linf)/(wlinf**(alphar))
    rorayl=betar*(wl**alphar)
    do ifi=1,nfi
        alphar=log(roatm_fi(1,lsup,ifi)/roatm_fi(1,linf,ifi))/coef
        betar=roatm_fi(1,linf,ifi)/(wlinf**(alphar))
        rorayl_fi(ifi)=betar*(wl**alphar)
    enddo

    alphac=log(roatm(2,lsup)/roatm(2,linf))/coef
    betac=roatm(2,linf)/(wlinf**(alphac))
    romix=betac*(wl**alphac)
    do ifi=1,nfi
        alphac=log(roatm_fi(2,lsup,ifi)/roatm_fi(2,linf,ifi))/coef
        betac=roatm_fi(2,linf,ifi)/(wlinf**(alphac))
        romix_fi(ifi)=betac*(wl**alphac)
    enddo

    if(iaer.eq.0) goto 2234
    alphaa=log(roatm(3,lsup)/roatm(3,linf))/coef
    betaa=roatm(3,linf)/(wlinf**(alphaa))
    roaero=betaa*(wl**alphaa)

    coefl=(wl-wldis(linf))/(wldis(lsup)-wldis(linf))
    do i=1,mu
        do j=1,nfilut(i)
            alphac=log(roluts(lsup,i,j)/roluts(linf,i,j))/coef
            betac=roluts(linf,i,j)/(wlinf**(alphac))
            rolut(i,j)=betac*(wl**alphac)
            if ((rolutsq(lsup,i,j).gt.0.001).and.(rolutsq(linf,i,j).gt.0.001)) then
                alphac=log(rolutsq(lsup,i,j)/rolutsq(linf,i,j))/coef
                betac=rolutsq(linf,i,j)/(wlinf**(alphac))
                rolutq(i,j)=betac*(wl**alphac)
            else
                rolutq(i,j)=rolutsq(linf,i,j)+(rolutsq(lsup,i,j)-rolutsq(linf,i,j))*coefl
            endif

            if ((rolutsu(lsup,i,j).gt.0.001).and.(rolutsu(linf,i,j).gt.0.001)) then
                alphac=log(rolutsu(lsup,i,j)/rolutsu(linf,i,j))/coef
                betac=rolutsu(linf,i,j)/(wlinf**(alphac))
                rolutu(i,j)=betac*(wl**alphac)
            else
                rolutu(i,j)=rolutsu(linf,i,j)+(rolutsu(lsup,i,j)-rolutsu(linf,i,j))*coefl
            endif

        enddo
    enddo
2234 continue

    if(iaer.eq.0) goto 3240
    if ((qhase(lsup).gt.0.001).and.(qhase(linf).gt.0.001)) then
        alphaa=log(qhase(lsup)/qhase(linf))/coef
        betaa=qhase(linf)/(wlinf**(alphaa))
        qhaa=betaa*(wl**alphaa)
    else
        qhaa=qhase(linf)+(qhase(lsup)-qhase(linf))*coefl
    endif
3240  continue
    qhar=depolar1*.75*(xmud*xmud-1.)
    if (idatmp.eq.0) goto 3234

    test1=abs(rqatm(1,linf))
    test2=abs(rqatm(1,lsup))
    test3=rqatm(1,lsup)*rqatm(1,linf)
    if((test1.lt.0.001).or.(test2.lt.0.001).or.(test3.lt.0.0)) then
        rqrayl=rqatm(1,linf)+(rqatm(1,lsup)-rqatm(1,linf))  &
                *(wl-wldis(linf))/(wldis(lsup)-wldis(linf))
    else
        alphar=log(rqatm(1,lsup)/rqatm(1,linf))/ coef
        betar=rqatm(1,linf)/(wlinf**(alphar))
        rqrayl=betar*(wl**alphar)
    endif

    test1=abs(rqatm(2,linf))
    test2=abs(rqatm(2,lsup))
    test3=rqatm(2,lsup)*rqatm(2,linf)
    if((test1.lt.0.001).or.(test2.lt.0.001).or.(test3.lt.0.0)) then
        rqmix=rqatm(2,linf)+(rqatm(2,lsup)-rqatm(2,linf))  &
                *(wl-wldis(linf))/(wldis(lsup)-wldis(linf))
    else
        alphac=log(rqatm(2,lsup)/rqatm(2,linf))/coef
        betac=rqatm(2,linf)/(wlinf**(alphac))
        rqmix=betac*(wl**alphac)
    endif
    if(iaer.eq.0) goto 3234

    test1=abs(rqatm(3,linf))
    test2=abs(rqatm(3,lsup))
    test3=rqatm(3,lsup)*rqatm(3,linf)
    if((test1.lt.0.001).or.(test2.lt.0.001).or.(test3.lt.0.0)) then
        rqaero=rqatm(3,linf)+(rqatm(3,lsup)-rqatm(3,linf))  &
            *(wl-wldis(linf))/(wldis(lsup)-wldis(linf))
    else
        alphaa=log(rqatm(3,lsup)/rqatm(3,linf))/coef
        betaa=rqatm(3,linf)/(wlinf**(alphaa))
        rqaero=betaa*(wl**alphaa)
    endif

3234  continue

    if(iaer.eq.0) goto 4242
    if ((uhase(lsup).gt.0.001).and.(uhase(linf).gt.0.001)) then
        alphaa=log(uhase(lsup)/uhase(linf))/coef
        betaa=uhase(linf)/(wlinf**(alphaa))
        uhaa=betaa*(wl**alphaa)
    else
        uhaa=uhase(linf)+(uhase(lsup)-uhase(linf))*coefl
    endif
4242   continue
    uhar=depolar1*3./2.*xmud
    if (idatmp.eq.0) goto 4234

    test1=abs(ruatm(1,linf))
    test2=abs(ruatm(1,lsup))
    test3=ruatm(1,lsup)*ruatm(1,linf)

    if((test1.lt.0.001).or.(test2.lt.0.001).or.(test3.lt.0.0)) then
        rurayl=ruatm(1,linf)+(ruatm(1,lsup)-ruatm(1,linf))  &
            *(wl-wldis(linf))/(wldis(lsup)-wldis(linf))
    else
        alphar=log(ruatm(1,lsup)/ruatm(1,linf))/ coef
        betar=ruatm(1,linf)/(wlinf**(alphar))
        rurayl=betar*(wl**alphar)
    endif

    test1=abs(ruatm(2,linf))
    test2=abs(ruatm(2,lsup))
    test3=ruatm(2,lsup)*ruatm(2,linf)
    if((test1.lt.0.001).or.(test2.lt.0.001).or.(test3.lt.0.0)) then
        rumix=ruatm(2,linf)+(ruatm(2,lsup)-ruatm(2,linf)) &
            *(wl-wldis(linf))/(wldis(lsup)-wldis(linf))
    else
        alphac=log(ruatm(2,lsup)/ruatm(2,linf))/coef
        betac=ruatm(2,linf)/(wlinf**(alphac))
        rumix=betac*(wl**alphac)
    endif
    if(iaer.eq.0) goto 4234
    test1=abs(ruatm(3,linf))
    test2=abs(ruatm(3,lsup))
    test3=ruatm(3,lsup)*ruatm(3,linf)
    if((test1.lt.0.001).or.(test2.lt.0.001).or.(test3.lt.0.0)) then
        ruaer0=ruatm(3,linf)+(ruatm(3,lsup)-ruatm(3,linf))  &
            *(wl-wldis(linf))/(wldis(lsup)-wldis(linf))
    else
        alphaa=log(ruatm(3,lsup)/ruatm(3,linf))/coef
        betaa=ruatm(3,linf)/(wlinf**(alphaa))
        ruaero=betaa*(wl**alphaa)
    endif
4234    continue

    alphar=log(trayl(lsup)/trayl(linf))/coef
    betar=trayl(linf)/(wlinf**(alphar))
    tray=betar*(wl**alphar)
    if (idatmp.ne.0.) then
        alphar=log(traypl(lsup)/traypl(linf))/coef
        betar=traypl(linf)/(wlinf**(alphar))
        trayp=betar*(wl**alphar)
    else
        trayp=0.
    endif

    if(iaer.eq.0) goto 1235
    alphaa=log(ext(lsup)*ome(lsup)/(ext(linf)*ome(linf)))/coef
    betaa=ext(linf)*ome(linf)/(wlinf**(alphaa))
    tsca=taer55*betaa*(wl**alphaa)/ext(8)
    alphaa=log(ext(lsup)/ext(linf))/coef
    betaa=ext(linf)/(wlinf**(alphaa))
    taerp=taer55p*betaa*(wl**alphaa)/ext(8)
    taer=taer55*betaa*(wl**alphaa)/ext(8)

1235 drinf=dtdif(1,linf)+dtdir(1,linf)
    drsup=dtdif(1,lsup)+dtdir(1,lsup)
    alphar=log(drsup/drinf)/coef
    betar=drinf/(wlinf**(alphar))
    dtotr=betar*(wl**alphar)
    dtinf=dtdif(2,linf)+dtdir(2,linf)
    dtsup=dtdif(2,lsup)+dtdir(2,lsup)
    alphac=log((dtsup*drinf)/(dtinf*drsup))/coef
    betac=(dtinf/drinf)/(wlinf**(alphac))
    dtotc=betac*(wl**alphac)
    dainf=dtdif(3,linf)+dtdir(3,linf)
    dasup=dtdif(3,lsup)+dtdir(3,lsup)
    if(iaer.eq.0) goto 1236
    alphaa=log(dasup/dainf)/coef
    betaa=dainf/(wlinf**(alphaa))
    dtota=betaa*(wl**alphaa)

1236 dtott=dtotc*dtotr
    urinf=utdif(1,linf)+utdir(1,linf)
    ursup=utdif(1,lsup)+utdir(1,lsup)
    alphar=log(ursup/urinf)/ coef
    betar=urinf/(wlinf**(alphar))
    utotr=betar*(wl**alphar)
    utinf=utdif(2,linf)+utdir(2,linf)
    utsup=utdif(2,lsup)+utdir(2,lsup)
    alphac=log((utsup*urinf)/(utinf*ursup))/ coef
    betac=(utinf/urinf)/(wlinf**(alphac))
    utotc=betac*(wl**alphac)
    uainf=utdif(3,linf)+utdir(3,linf)
    uasup=utdif(3,lsup)+utdir(3,lsup)
    if(iaer.eq.0) goto 1237
    alphaa=log(uasup/uainf)/ coef
    betaa=uainf/(wlinf**(alphaa))
    utota=betaa*(wl**alphaa)

1237 utott=utotc*utotr
    arinf=sphal(1,linf)
    arsup=sphal(1,lsup)
    alphar=log(arsup/arinf)/ coef
    betar=arinf/(wlinf**(alphar))
    asray=betar*(wl**alphar)
    atinf=sphal(2,linf)
    atsup=sphal(2,lsup)
    alphac=log(atsup/atinf)/coef
    betac=atinf/(wlinf**(alphac))
    astot=betac*(wl**alphac)
    aainf=sphal(3,linf)
    aasup=sphal(3,lsup)
    if(iaer.eq.0) goto 1239
    alphaa=log(aasup/aainf)/coef
    betaa=aainf/(wlinf**(alphaa))
    asaer=betaa*(wl**alphaa)
1239 return
end
