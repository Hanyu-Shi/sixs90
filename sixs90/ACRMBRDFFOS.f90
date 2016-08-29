!************************************************************************
subroutine acrmbrdffos(lai2,sl2,clmp2,szz,eln2,thm2,nratio2,slw2,lmod2,&
            ncomp2,ccomp2,N2,dcell2,asp2, &
            lai1,sl1,clmp1,eln1,thm1,nratio1,slw1,lmod1, &
            ncomp1,ccomp1,N1,dcell1,asp1,&
            s1,s2,s3,s4, &
            wlmoy,mu,rm,rosur,wfisur,fisur)
    use paramdef
    implicit none
    integer :: mu
    real(8) :: wlmoy,rm(-mu:mu),rosur(0:mu_p,mu_p,83),fisur(83),wfisur(83)
    real(8) :: pisp,dr,rtv,rts,rfi,r_lamda,bq
    character(30)::lmod1,lmod2
    integer :: ncomp2,ncomp1
    real(8) :: LAI2,sl2,clmp2,szz,eln2,thm2,nratio2,slw2,N2,dcell2,asp2, &
         LAI1,sl1,clmp1,eln1,thm1,nratio1,slw1,N1,dcell1,asp1, &
         s1,s2,s3,s4,ccomp1(10), ccomp2(10)
    real(8) :: mu1,mu2
    integer :: k,j,i

    r_lamda = wlmoy*1000.d0
    pisp = acos(-1.d0)
    dr = pisp/180.0

    call gauss(0.d0,pisp,fisur,wfisur,83)
    do i = 0, mu
        do j = 1, mu
            do k = 1, 83
                if (i == 0) then
                    mu1 = rm(0)
                else
                    mu1 = rm(i)
                end if
                mu2 = rm(j)
                rfi = fisur(k) + pisp
                rts = acos(mu1)
                rtv = acos(mu2)

                if (rfi < 0.) rfi = rfi + 2.*pisp
                if (rfi > (2.*pisp)) rfi = rfi - 2.*pisp
                if (rfi > pisp) rfi = 2.*pisp - rfi
                if (abs(rtv) < dr) rfi = 0.d0

                call ACRM(lai2,sl2,clmp2,szz,eln2,thm2,nratio2,slw2,lmod2,&
                    ncomp2,ccomp2,N2,dcell2,asp2, &
                    lai1,sl1,clmp1,eln1,thm1,nratio1,slw1,lmod1, &
                    ncomp1,ccomp1,N1,dcell1,asp1,&
                    s1,s2,s3,s4,&
                    rts,rtv,rfi,r_lamda,bq)
                rosur(i,j,k) = bq
            end do
        end do
    end do
end subroutine acrmbrdffos
