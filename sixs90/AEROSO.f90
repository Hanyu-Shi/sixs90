subroutine aeroso(iaer,co,xmud,wldis,FILE,ipol)

! - to vary the number of quadratures
    use paramdef
    implicit none
    integer :: nquad
    common /num_quad/ nquad
    real(8) :: ph,qh,uh
    common /sixs_aerbas/ ph(20,nqmax_p),qh(20,nqmax_p),uh(20,nqmax_p)
    real(8) :: cgaus_S(nqmax_p),pdgs_S(nqmax_p)
    real(8) :: phasel,qhasel,uhasel
    common /sixs_phase/ phasel(20,nqmax_p),qhasel(20,nqmax_p),uhasel(20,nqmax_p)
    integer :: nbmu, nbmu_2
    real(8) :: cosang(nqmax_p),weight(nqmax_p)
! - to vary the number of quadratures

    real(8) :: cij(4),nis,sumni  !,ni(4)
    real(8) :: co(4),dd(4,20),qq(4,20),ci(4)
    real(8) :: pha(5,20,nqmax_p),qha(5,20,nqmax_p),uha(5,20,nqmax_p)
    real(8) :: sca(20),wldis(20),uu(4,20)
    real(8) :: xmud,ext,ome,gasym,phase,qhase,uhase
    real(8) :: coef,sigm,pi
    integer :: i,j,k,l,j1,j2,iaer,icp,ipol
    !character FILE*80
    character(80) :: FILE

    common /sixs_aer/ ext(20),ome(20),gasym(20),phase(20),qhase(20),uhase(20)

    real(8) :: ex,sc,asy,vi
    common /sixs_coef/ ex(4,20),sc(4,20),asy(4,20),vi(4)

    real(8) :: wldisc(20)

    data wldisc /0.350,0.400,0.412,0.443,0.470,0.488,0.515,0.550,   &
                0.590,0.633,0.670,0.694,0.760,0.860,1.240,1.536,    &
                1.650,1.950,2.250,3.750/

    pi=4.*atan(1.d0)

    if (iaer.eq.12) then
        open(10,file=FILE)
        read (10,*) nbmu
        read(10,*)
        do l=1,20
            read(10,'(10x,4(3x,f8.4,3x))')ext(l),sca(l),ome(l),gasym(l)
        enddo
        read(10,'(///)')
        do k=1,nbmu
            read(10,'(8x,20(1x,e10.4))')(phasel(l,k),l=1,20)
        enddo
        if (ipol.ne.0)then
            do k=1,nbmu
                read(10,'(8x,20(1x,e10.4))')(qhasel(l,k),l=1,20)
            enddo
            do k=1,nbmu
                read(10,'(8x,20(1x,e10.4))')(uhasel(l,k),l=1,20)
            enddo
        endif
        nquad=nbmu
        close(10)
    endif

! - calculation of gauss points
    nbmu=nquad
    nbmu_2=(nbmu-3)/2
    call gauss(-1.d0,1.d0,cosang,weight,nbmu-3)
    cgaus_S(1)=-1.0
    pdgs_S(1)=0.0
    do j=1,nbmu_2
        cgaus_S(j+1)=cosang(j)
        pdgs_S(j+1)=weight(j)
    enddo
    cgaus_S(nbmu_2+2)=0.
    pdgs_S(nbmu_2+2)=0.
    do j=nbmu_2+1,nbmu-3
        cgaus_S(j+2)=cosang(j)
        pdgs_S(j+2)=weight(j)
    enddo
    cgaus_S(nbmu)=1.0
    pdgs_S(nbmu)=0.
! - calculation of gauss points

    do k=1,nbmu-1
        if((xmud.ge.cgaus_S(k)).and.(xmud.lt.cgaus_S(k+1))) go to 8
    enddo
    return
8   j1=k
    j2=j1+1
    coef=-(xmud-cgaus_S(j1))/(cgaus_S(j2)-cgaus_S(j1))

    if(iaer.eq.12) then
        do l=1,20
            phase(l)=phasel(l,j1)+coef*(phasel(l,j1)-phasel(l,j2))
        enddo
    if(ipol.ne.0)then
        do l=1,20
            qhase(l)=qhasel(l,j1)+coef*(qhasel(l,j1)-qhasel(l,j2))
            uhase(l)=uhasel(l,j1)+coef*(uhasel(l,j1)-uhasel(l,j2))
        enddo
    endif
    return
    endif

    do l=1,20
        ext(l)=0.
        sca(l)=0.
        if(l.eq.4.and.iaer.eq.0) ext(l)=1.
        ome(l)=0.
        gasym(l)=0.
        phase(l)=0.
        qhase(l)=0.
        uhase(l)=0.
        do k=1,nbmu
            phasel(l,k)=0.
            qhasel(l,k)=0.
            uhasel(l,k)=0.
        enddo
    enddo

    do j=1,4
        ci(j)=co(j)
    enddo

    if(iaer.eq.0) goto 777

    if (iaer.ge.5.and.iaer.le.11) then
!   calling a special aerosol model
!   (background desert model...)
        if (iaer.eq.5) call bdm
!   (biomass burning model...)
        if (iaer.eq.6) call bbm
!   (stratospherique aerosol model...)
        if (iaer.eq.7) call stm
!   (user defined model from size distribution)
        if (iaer.ge.8.and.iaer.le.11) then
            call mie(iaer,wldis,ex,sc,asy,ipol)
        endif

        do l=1,20
            dd(1,l)=ph(l,j1)+coef*(ph(l,j1)-ph(l,j2))
            do k=1,nbmu
                pha(1,l,k)=ph(l,k)
            enddo
        enddo
        if(ipol.ne.0)then
            do l=1,20
                qq(1,l)=qh(l,j1)+coef*(qh(l,j1)-qh(l,j2))
                uu(1,l)=uh(l,j1)+coef*(uh(l,j1)-uh(l,j2))
                do k=1,nbmu
                    qha(1,l,k)=qh(l,k)
                    uha(1,l,k)=uh(l,k)
                enddo
            enddo
        endif
        icp=1
        cij(1)=1.00
!   for normalization of the extinction coefficient
        nis=1.d+00/ex(1,8)
    else
! calling each sra components
        icp=4
! extrapolate each component for wavelength
        do l=1,20
            do j=1,4
                ex(j,l)=0
                sc(j,l)=0.
                asy(j,l)=0.
            enddo
        enddo

!  phase function of 4 components
        do j=1,4
            if (j.eq.1) call dust
            if (j.eq.2) call wate
            if (j.eq.3) call ocea
            if (j.eq.4) call soot
        do l=1,20
            dd(j,l)=ph(l,j1)+coef*(ph(l,j1)-ph(l,j2))
            do k=1,nbmu
                pha(j,l,k)=ph(l,k)
            enddo
            if (ipol.ne.0)then
                qq(j,l)=qh(l,j1)+coef*(qh(l,j1)-qh(l,j2))
                uu(j,l)=uh(l,j1)+coef*(uh(l,j1)-uh(l,j2))
                do k=1,nbmu
                    qha(j,l,k)=qh(l,k)
                    uha(j,l,k)=uh(l,k)
                enddo
            endif
            enddo
        enddo
!   summ of the ci/vi calculation

        sigm=0.
        sumni=0.0
        do i=1,4
            sigm=sigm+ci(i)/vi(i)
        enddo
!   cij coefficients calculation
        do j=1,4
            cij(j)=(ci(j)/vi(j)/sigm)
            sumni=sumni+cij(j)*ex(j,8)
        enddo
        nis=1.d+00/sumni
    endif

!     mixing parameters calculation
    do l=1,20
        do j=1,icp
            ext(l)=ex(j,l)*cij(j)+ext(l)
            sca(l)=sc(j,l)*cij(j)+sca(l)
            gasym(l)=sc(j,l)*cij(j)*asy(j,l)+gasym(l)
            phase(l)=sc(j,l)*cij(j)*dd(j,l)+phase(l)
            do k=1,nbmu
                phasel(l,k)=sc(j,l)*cij(j)*pha(j,l,k)+phasel(l,k)
            enddo
            if (ipol.ne.0)then
                qhase(l)=sc(j,l)*cij(j)*qq(j,l)+qhase(l)
                uhase(l)=sc(j,l)*cij(j)*uu(j,l)+uhase(l)
                do k=1,nbmu
                    qhasel(l,k)=sc(j,l)*cij(j)*qha(j,l,k)+qhasel(l,k)
                    uhasel(l,k)=sc(j,l)*cij(j)*uha(j,l,k)+uhasel(l,k)
                enddo
            endif
        enddo
        ome(l)=sca(l)/ext(l)
        gasym(l)=gasym(l)/sca(l)
        phase(l)=phase(l)/sca(l)
        do k=1,nbmu
            phasel(l,k)=phasel(l,k)/sca(l)
        enddo
        if (ipol.ne.0)then
            qhase(l)=qhase(l)/sca(l)
            uhase(l)=uhase(l)/sca(l)
            do k=1,nbmu
                qhasel(l,k)=qhasel(l,k)/sca(l)
                uhasel(l,k)=uhasel(l,k)/sca(l)
            enddo
        endif
        ext(l)=ext(l)*nis
        sca(l)=sca(l)*nis
    enddo

    if(iaer.ge.8.and.iaer.le.11) then
        open(10,file=FILE)
        write(10,*) nbmu
        write(10,'(3x,A5,1x,5(1x,A10,1x),1x,A10)')'Wlgth','Nor_Ext_Co',    &
                'Nor_Sca_Co','Sg_Sca_Alb','Asymm_Para','Extinct_Co','Scatter_Co'
        do l=1,20
            write(10,'(2x,f8.4,4(3x,f8.4,3x),2(2x,e10.4))')            &
                    wldis(l),ext(l),sca(l),ome(l),gasym(l),ext(l)/nis,sca(l)/nis
        enddo
        write(10,'(//,T20,A16,/,3x,A4,1x,20(3x,f6.4,2x))')   &
                ' Phase Function ','TETA',(wldis(l),l=1,20)
        do k=1,nbmu
            write(10,761)180.*acos(cgaus_S(k))/pi,(phasel(l,k),l=1,20)
        enddo
761     format(2x,f6.2,20(1x,e10.4))
        if(ipol.ne.0)then
            do k=1,nbmu
                write(10,761)180.*acos(cgaus_S(k))/pi,(qhasel(l,k),l=1,20)
            enddo
            do k=1,nbmu
                write(10,761)180.*acos(cgaus_S(k))/pi,(uhasel(l,k),l=1,20)
            enddo
        endif
        close(10)
    endif

777 continue
    return
end
