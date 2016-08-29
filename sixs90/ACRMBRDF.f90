!************************************************************************
subroutine acrmbrdf(lai2,sl2,clmp2,szz,eln2,thm2,nratio2,slw2,lmod2,&
            ncomp2,ccomp2,N2,dcell2,asp2, &
            lai1,sl1,clmp1,eln1,thm1,nratio1,slw1,lmod1, &
            ncomp1,ccomp1,N1,dcell1,asp1,&
            s1,s2,s3,s4, &
            wlmoy,mu,np,rm,rp,brdfint)
    implicit none
    integer :: mu,np
    real(8) :: wlmoy,rm(-mu:mu),rp(np),brdfint(-mu:mu,np)
    real(8) :: pisp,dr,rtv,rts,rfi,r_lamda,bq
    character(30)::lmod1,lmod2
    integer :: ncomp2,ncomp1
    real(8) :: LAI2,sl2,clmp2,szz,eln2,thm2,nratio2,slw2,N2,dcell2,asp2, &
         LAI1,sl1,clmp1,eln1,thm1,nratio1,slw1,N1,dcell1,asp1, &
         s1,s2,s3,s4,ccomp1(10), ccomp2(10)

    integer :: k,j

    r_lamda = wlmoy*1000.d0
    pisp = acos(-1.d0)
    dr = pisp/180.0
    rts = acos(rm(0))
    do k=1,np
        do j=1,mu
            rtv = acos(rm(j))
            if(j.eq.mu) then
                rfi = rm(-mu)
            else
                rfi = rp(k) + rm(-mu)
            endif
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
            brdfint(j,k) = bq
        enddo
    enddo
end subroutine acrmbrdf
