subroutine kernel(is,mu,rm,xpl,psl,bp)
    use paramdef
    implicit none
    integer :: mu
    real(8) :: rm(-mu:mu)
    real(8) :: xpl(-mu:mu),bp(0:mu,-mu:mu)
    integer :: is,ip1,j,i,k,ip,ig,l,lp,lm,ij
    real(8) :: xdb,a,b,c,xx,rac3,x,bt,sbp

! - to vary the number of quadratures
    integer :: nquad
    common /num_quad/ nquad
    real(8) :: pha,qha,uha,alphal,betal,gammal,zetal
    common /sixs_polar/ pha(nqmax_p),qha(nqmax_p),uha(nqmax_p),&
        alphal(0:nqmax_p),betal(0:nqmax_p),gammal(0:nqmax_p),&
        zetal(0:nqmax_p)
    real(8) :: psl(-1:nqmax_p,-mu:mu)
! - to vary the number of quadratures
    ip1=nquad-3
    rac3=dsqrt(3.D+00)
    if(is.ne.0) then
        if(is.ne.1) then
            a=1
            do i=1,is
                x=i
                a=a*sqrt((i+is)/x)*0.5
            enddo
            b=a*sqrt(is/(is+1.))*sqrt((is-1.)/(is+2.))
            do j=0,mu
                c=dble(rm(j))
                xx=1.-c*c
                psl(is-1,j)=0.
                xdb=a*xx**(is*0.5)
                if (abs(xdb).lt.1.E-30) xdb=0.0
                psl(is,-j)=xdb
                psl(is,j)=xdb
            enddo
        else
            do j=0,mu
                c=dble(rm(j))
                x=1.-c*c
                psl(0,j)=0.
                psl(0,-j)=0.
                psl(1,-j)=sqrt(x*0.5)
                psl(1,j)=sqrt(x*0.5)
                psl(2,j)=c*psl(1,j)*rac3
                psl(2,-j)=-psl(2,j)
            enddo
            psl(2,0)=-psl(2,0)
        endif
    else
        do j=0,mu
            c=dble(rm(j))
            psl(0,-j)=1.
            psl(0,j)=1.
            psl(1,j)=c
            psl(1,-j)=-c
            xdb=(3.*c*c-1.)*0.5
            if (abs(xdb).lt.1.E-30) xdb=0.0
            psl(2,-j)=xdb
            psl(2,j)=xdb
        enddo
        psl(1,0)=rm(0)
    endif

    k=2
    ip=ip1
    if(is.gt.2)k=is
    if(k.ne.ip) then
        ig=-1
        if(is.eq.1)ig=1
        do l=k,ip-1
            lp=l+1
            lm=l-1
            a=(2*l+1.)/sqrt((l+is+1.)*(l-is+1.))
            b=sqrt(float((l+is)*(l-is)))/(2.*l+1.)
            do j=0,mu
                c=dble(rm(j))
                xdb=a*(c*psl(l,j)-b*psl(lm,j))
                if (abs(xdb).lt.1.E-30) xdb=0.
                psl(lp,j)=xdb
                if(j.ne.0) then
                    psl(lp,-j)=ig*psl(lp,j)
                endif
            enddo
            ig=-ig
        enddo
    endif

    do j=-mu,mu
        xpl(j)=psl(2,j)
    enddo
    ij=ip1
    do j=0,mu
        do k=-mu,mu
            sbp=0.
            if(is.le.ij) then
                do l=is,ij
                    bt=betal(l)
                    sbp=sbp+dble(psl(l,j))*psl(l,k)*bt
                enddo
            endif
            if (abs(sbp).lt.1.E-30) sbp=0.
          bp(j,k)=sbp
        enddo
    enddo
    return
end
