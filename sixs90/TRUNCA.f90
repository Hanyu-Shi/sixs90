subroutine trunca(coeff,ipol)

! - to vary the number of quadratures
    use paramdef
    implicit none
    real(8) :: coeff
    integer :: nquad,ipol,nbmu_2
    common /num_quad/ nquad
    real(8) :: cgaus_S(nqmax_p), pdgs_S(nqmax_p)
    real(8) :: pl(-1:nqmax_p),pol(0:nqmax_p),deltal(0:nqmax_p)
    real(8) :: pha,qha,uha,alphal,betal,gammal,zetal
    common /sixs_polar/ pha(nqmax_p),qha(nqmax_p),uha(nqmax_p),&
            alphal(0:nqmax_p),betal(0:nqmax_p),gammal(0:nqmax_p),&
            zetal(0:nqmax_p)
! - to vary the number of quadratures

    real(8) :: x,rm,z1p,e,d,co1,co2,co3,xx,c2
    real(8) :: som1,som2,som3,som4
    integer :: nbmu,k,j,i,nn,mm
    real(8) :: cosang(nqmax_p),weight(nqmax_p)

    nbmu=nquad

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


! Computations of Legendre coefficients
    do k=0,nbmu-3
        alphal(k)=0.
        betal(k)=0.
        gammal(k)=0.
        deltal(k)=0.
        zetal(k)=0.
    enddo
    do j=1,nbmu
        x=pha(j)*pdgs_S(j)
        rm=cgaus_S(j)
        pl(-1)=0.
        pl(0)=1.
        do k=0,nbmu-3
            pl(k+1)=((2*k+1.)*rm*pl(k)-k*pl(k-1))/(k+1.)
            betal(k)=betal(k)+x*pl(k)
        enddo
    enddo
    do k=0,nbmu-3
        betal(k)=(2*k+1.)*0.5*betal(k)
! - to put negative coefficients to 0
        if (betal(k).lt.0) then
            do j=k,nbmu-3
                betal(j)=0.0
            enddo
        goto 133
        endif
! - to put negative coefficients to 0
    enddo

133 continue

!   cases of polarization
    if (ipol.ne.0)then
        do j=1,nbmu
            x=qha(j)*pdgs_S(j)
            xx=uha(j)*pdgs_S(j)
            rm=cgaus_S(j)
            pol(0)=0.
            pol(1)=0.
            pol(2)=3.*(1.-rm**2)/2./sqrt(6.0)
            pl(-1)=0.
            pl(0)=1.
            do k=2,nbmu-3
                d=(2.*k+1.)/sqrt((k+3)*(k-1.))
                e=sqrt((k+2.)*(k-2.))/(2.*k+1.)
                pol(k+1)=d*(rm*pol(k)-e*pol(k-1))
                gammal(k)=gammal(k)+x*pol(k)
            enddo
            do k=0,nbmu-3
                pl(k+1)=((2.*k+1.)*rm*pl(k)-k*pl(k-1))/(k+1.)
                deltal(k)=deltal(k)+xx*pl(k)
            enddo
        enddo
        do k=0,nbmu-3
            deltal(k)=deltal(k)*(2.*k+1.)/2.
            gammal(k)=gammal(k)*(2.*k+1.)/2.
        enddo

        do i=2,nbmu-3
            co1=4.*(2.*i+1.)/(i*(i-1.)*(i+1.)*(i+2.))
            co2=i*(i-1.)/((i+1.)*(i+2.))
            co3=co2*deltal(i)
            co2=co2*betal(i)
            nn=i/2
            mm=(i-1)/2
            som1=0.
            som2=0.
            som3=0.
            som4=0.
            do j=1,nn
                c2=(i-1.)*(i-1.)-3.*(2*j-1.)*(i-j)
                som1=som1+c2*betal(i-2*j)
                som2=som2+c2*deltal(i-2*j)
            enddo
            do j=0,mm
                c2=(i-1.)*(i-1.)-3.*j*(2*i-2*j-1.)
                som3=som3+c2*betal(i-2*j-1)
                som4=som4+c2*deltal(i-2*j-1)
            enddo
            zetal(i)=co3-co1*(som2-som3)
            alphal(i)=co2-co1*(som1-som4)
        enddo

        z1p=betal(0)

        do k=0,nbmu-3
            alphal(k)=alphal(k)/z1p
            betal(k)=betal(k)/z1p
            gammal(k)=gammal(k)/z1p
            zetal(k)=zetal(k)/z1p
        enddo
    endif
    coeff=0.0
    return
end
