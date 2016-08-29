real(8) function Ro_1 (Theta_i,Phi_i,Theta_e,Phi_e)
    implicit none
    common /gauss_m/xgm (20),wgm (20),n
    real(8) :: xgm,wgm
    integer :: n
    real(8) :: G_f,Geo,h,gamma_f
    common /p/xLt,Rl,Tl,Rs,c,ild
    real(8) :: xLt,Rl,Tl,Rs,c
    integer :: ild
    common /Ro/Ro_1_c,Ro_1_s,Ro_mult
    real(8) :: Ro_1_c,Ro_1_s,Ro_mult
    real(8) :: Theta_i,Phi_i,Theta_e,Phi_e,xmui,xmu,xtmu
    real(8) :: Gi,Ge,Ki,Ke,xLi
    real(8) :: xmm,xrm
    real(8) :: xL
    real(8) :: xdb
    integer :: j

    xtmu=1.E-05
    xmui=abs(cos(Theta_i))
    xmu=cos(Theta_e)
    if (abs(xmu).lt.xtmu) xmu=xtmu

    Gi=G_f (Theta_i)
    Ge=G_f (Theta_e)

    Ki=Gi/xmui
    Ke=Ge/xmu

    xmm=0.5*(xLt+0.)
    xrm=0.5*(xLt-0.)
    Ro_1_c=0.
    xLi=c/Geo (Theta_i,Phi_i,Theta_e,Phi_e)
    do j=1,n
        xL=xmm+xrm*xgm(j)
        xdb=(Ki+Ke*h(xL,xLi))*dble(xL)
        if (abs(xdb).lt.1.E-30) xdb=0.
        if (xdb.le.20) Ro_1_c=Ro_1_c+wgm(j)*xrm*dexp(-xdb)
    enddo

    Ro_1_c=Ro_1_c*Gamma_f (Theta_i,Phi_i,Theta_e,Phi_e)/xmui/xmu

    xdb=(Ki+Ke*h (xLt,xLi))*dble(xLt)
    if (abs(xdb).lt.1.E-30) xdb=0.
    if (xdb.le.20) Ro_1_s=Rs*dexp(-xdb)

    Ro_1=Ro_1_c+Ro_1_s

    return
end

real(8) function Gamma_f (Theta_p,Phi_p,Theta,Phi)
    implicit none
    real(8),parameter::Pi=3.141592653589793
    common /p/xLt,Rl,Tl,Rs,c,ild
    real(8) :: xLt,Rl,Tl,Rs,c,gl
    integer :: ild
    common /gauss_m/xgm (20),wgm (20),n
    real(8) :: xgm,wgm
    integer :: n
    real(8) :: Theta_p,Phi_p,Theta,Phi
    real(8) :: xmm,xrm,xt
    real(8) :: ymm,yrm,yt
    real(8) :: dp,dpp,f
    real(8) :: sum
    integer :: i,j

    xmm=0.5*(Pi/2.+0.)
    xrm=0.5*(Pi/2.-0.)
    ymm=0.5*(2.*Pi+0.)
    yrm=0.5*(2.*Pi-0.)
    Gamma_f = 0.
   do j=1,n
        xt=xmm+xrm*xgm (j)
        sum=0.
        do i=1,n
            yt=ymm+yrm*xgm (i)
            dpp = cos (Theta_p)*cos (xt)+sin (Theta_p)*sin (xt)*cos (Phi_p-yt)
            dp = cos (Theta)*cos (xt)+sin (Theta)*sin (xt)*cos (Phi-yt)
! correction when porting code to HP730
            if (dp*dpp.lt.0.) then
                f=Rl*abs (dp)/Pi
            else
                f=Tl*abs (dp)/Pi
            endif
! end of correction
            sum=sum+wgm (i)*xrm*gl (xt)*f*abs (dpp)
        enddo
        Gamma_f=Gamma_f+wgm (j)*yrm*sum
    enddo
    Gamma_f=Gamma_f/2.
    return
end

real(8) function Geo (Theta_i,Phi_i,Theta_e,Phi_e)
    implicit none
    real(8) :: Theta_i,Phi_i,Theta_e,Phi_e
    Geo=sqrt (abs(tan(Theta_i)**2+tan(Theta_e)**2-2.*tan(Theta_i)*tan(Theta_e)*cos(Phi_i-Phi_e)))
    if (Geo.lt.1.e-35) Geo=1.e-35
    return
end

real(8) function h (xL,xLi)
    implicit none
    real(8),parameter::Pi=3.141592653589793
    real(8) :: xL,xLi
    if (xL.lt.xLi) then
        h=(1.-4./3./Pi)/xLi*xL
    else
        h=1.-4./3./Pi*xLi/xL
    endif
    return
end

real(8) function G_f (Theta)
    implicit none
    real(8),parameter::Pi=3.141592653589793
    common /p/xLt,Rl,Tl,Rs,c,ild
    real(8) :: xLt,Rl,Tl,Rs,c,psi,gl
    integer :: ild
    common /gauss_m/xgm (20),wgm (20),n
    real(8) :: xgm,wgm
    integer :: n
    real(8) :: Theta
    real(8) :: xmm,xrm,xt
    integer :: j

    xmm=0.5*(Pi/2.+0.)
    xrm=0.5*(Pi/2.-0.)
    G_f = 0.
    do j=1,n
        xt=xmm+xrm*xgm (j)
        G_f=G_f+wgm (j)*xrm*Psi (Theta,xt)*gl (xt)
    enddo
    return
end

real(8) function Psi (Theta,xt)
    implicit none
    real(8),parameter::Pi=3.141592653589793
    common /p/xLt,Rl,Tl,Rs,c,ild
    real(8) :: xLt,Rl,Tl,Rs,c
    integer :: ild
    real(8) :: Theta,xt
    real(8) :: cpt,pt
    real(8) :: xmu,smu

    xmu=cos (xt)
    smu=sin (xt)
    if (xmu.eq.1.) then
        Psi=cos (Theta)
    else
        if (sin (Theta).eq.0.) then
            Psi=xmu
        else
            if (smu.eq.0.) then
                cpt=0.
            else
                cpt=1.*xmu/smu*cos (Theta)/sin (Theta)
            endif
            if (abs (cpt).gt.1.) then
                Psi=xmu*cos (Theta)
            else
                pt=acos (-cpt)
                Psi=xmu*cos (Theta)*(2./Pi*pt-1.)+2./Pi*smu*sin (Theta)*sin(pt)
            endif
        endif
    endif
    Psi=abs (Psi)
    return
end

real(8) function gl (Theta)
    implicit none
    real(8),parameter::Pi=3.141592653589793
    common /ld/a_ld,b_ld,c_ld,d_ld
    real(8) :: a_ld,b_ld,c_ld,d_ld
    real(8) :: Theta

    gl=a_ld+b_ld*cos (2.*Theta)+c_ld*cos (4.*Theta)+d_ld*sin (Theta)
    return
end

subroutine gauleg(x1,x2,x,w,n)
    implicit none
    integer :: n
    real(8) :: x1,x2,x(n),w(n)
    real(8) :: eps
    parameter (eps=3.d-14)
    integer :: i,j,m
    real(8) :: p1,p2,p3,pp,xl,xm,z,z1
    m=(n+1)/2
    xm=0.5d00*(x2+x1)
    xl=0.5d00*(x2-x1)
    do i=1,m
        z=cos(3.141592654d00*(i-.25d00)/(n+.5d00))
1       continue
        p1=1.d00
        p2=0.d00
        do j=1,n
            p3=p2
            p2=p1
            p1=((2.d00*j-1.d00)*z*p2-(j-1.d00)*p3)/j
        enddo
        pp=n*(z*p1-p2)/(z*z-1.d00)
        z1=z
        z=z1-p1/pp
        if(abs(z-z1).gt.eps)go to 1
        x(i)=real(xm-xl*z)
        x(n+1-i)=real(xm+xl*z)
        w(i)=real(2.d00*xl/((1.d00-z*z)*pp*pp))
        w(n+1-i)=w(i)
    enddo
    return
end

subroutine solve (Theta_i)
    implicit none
    real(8),parameter::Pi=3.141592653589793
    integer,parameter::m=20
    common /gauss_m/xgm (20),wgm (20),n
    real(8) :: xgm,wgm,g_f
    integer :: n
    common /p/xLt,Rl,Tl,Rs,c,ild
    real(8) :: xLt,Rl,Tl,Rs,c
    integer :: ild
    common /Ro/Ro_1_c,Ro_1_s,Ro_mult
    real(8) :: Ro_1_c,Ro_1_s,Ro_mult
    real(8) :: Theta_i,xmui,Gi
    real(8) :: xdb
    common /l/dL,xL
    real(8) :: dL,xL
    real(8) :: xI0t,xI1t,xImt
    real(8) :: xI (m+1,20)
    real(8) :: Q0d (m),Q0u (m),Q1 (m),S (m),xIf (m+1,20)
    real(8) :: G (20)
    integer :: j,k,l
    real(8) :: xmm,xrm,xmu
    real(8) :: sum
    integer :: nc
    real(8) :: Epsilon

    Epsilon=1.e-4

    xmui=abs (cos (Theta_i))

    dL=xLt/float (m)
!
! - Computation of G-functions
!
    Gi=G_f (Theta_i)
    xmm=0.5*(1.+(-1.))
    xrm=0.5*(1.-(-1.))
    do j=1,n
        xmu=xmm+xrm*xgm (j)
        G (j)=G_f (acos (xmu))
    enddo
!
! - Initialisation of S (k) & xIf (k,j)
!
    do k=1,m
        S (k)=0.
        do j=1,n
            xIf (k,j)=0.
        enddo
    enddo
!
! - Computation of Q0d (k) & Q0u (k) <- first collision source
!
! - (down)
!
    do k=1,m
        xL=(k-.5)*dL
        xdb=Gi/xmui*dble(xL)
        if (abs(xdb).lt.1.E-30) xdb=0.
        if (xdb.lt.20) Q0d(k)=(Rl+Tl)/2.*Gi*dexp(-xdb)
    enddo
!
! - (up)
!
    xdb=Gi/xmui*dble(xLt)
    if (abs(xdb).lt.1.E-30) xdb=0.
    if (xdb.lt.20) xI0t=2.*Rs*xmui*dexp(-xdb)
    do k=m,1,-1
        xL=(k-.5)*dL
        sum=0.
        do j=n/2+1,n
            xmu=xmm+xrm*xgm (j)
            xdb=dble(G(j))/dble(xmu)*(xLt-xL)
            if (abs(xdb).lt.1.E-30) xdb=0.
            if (xdb.lt.20) sum=sum+wgm(j)*xrm*xI0t*(Rl+Tl)/2.*G(j)*dexp(-xdb)
        enddo
        Q0u (k)=sum
    enddo
!
! - Computation of xI (k,j) <- single scattering
!
! - Initialisation of xI (k,j)
!
    do k=1,m+1
        do j=1,n/2
            xI (k,j)=0.
        enddo
    enddo
!
! - (down)
!
    do k=1,m
        do j=1,n/2
            xmu=xmm+xrm*xgm (j)
            xI (k+1,j)=(Q0d (k)-xI (k,j)*(G (j)/2.+xmu/dL))/(G (j)/2.-xmu/dL)
        enddo
    enddo
!
! - (boundary condition)
!
    xI1t=0.
    do j=1,n/2
        xmu=xmm+xrm*xgm (j)
        xI1t=xI1t+wgm (j)*xrm*2.*Rs*abs (xmu)*xI (m+1,j)
    enddo

    do j=n/2+1,n
        xI (m+1,j)=0.
    enddo
!
! - (up)
!
    do k=m,1,-1
        do j=n/2+1,n
            xmu=xmm+xrm*xgm (j)
            xI (k,j)=(Q0d (k)-xI (k+1,j)*(G (j)/2.-xmu/dL))/(G (j)/2.+xmu/dL)
        enddo
    enddo
!
! - Computation of Q1 (k) <- second collision source
!
    do k=1,m
        sum=0.
        do j=1,n
            sum=sum+wgm (j)*xrm*(Rl+Tl)/2.*G (j)*(xI (k+1,j)+xI (k,j))/2.
        enddo
        Q1 (k)=sum
    enddo
!
! - Computation of xI (k,j) <- multiple scattering
!
! - Initialisation of xI (k,j)
!
    do k=1,m+1
        do j=1,n/2
            xI (k,j)=0.
        enddo
    enddo
    l=0
1   l=l+1
!
! - (down)
!
    do k=1,m
        do j=1,n/2
            xmu=xmm+xrm*xgm (j)
            xI (k+1,j)=(S (k)+Q0u (k)+Q1 (k)-xI (k,j)*(G (j)/2.+xmu/dL))/(G (j)/2.-xmu/dL)
        enddo
    enddo
!
!- (boundary condition)
!
    xImt=0.
    do j=1,n/2
        xmu=xmm+xrm*xgm (j)
        xImt=xImt+wgm (j)*xrm*2.*Rs*abs (xmu)*xI (m+1,j)
    enddo
    do j=n/2+1,n
        xI (m+1,j)=xImt+xI1t
    enddo
!
! - (up)
!
    do k=m,1,-1
        do j=n/2+1,n
            xmu=xmm+xrm*xgm (j)
            xI (k,j)=(S (k)+Q0u (k)+Q1 (k)-xI (k+1,j)*(G (j)/2.-xmu/dL))/(G (j)/2.+xmu/dL)
        enddo
    enddo
!
! - End test
!
    nc=0
    do k=1,m+1
        do j=1,n
            if (abs (xIf (k,j)-xI (k,j)).lt.Epsilon) nc=nc+1
            xIf (k,j)=xI (k,j)
        enddo
    enddo
    if ((l.lt.50).and.(nc.ne.(m+1)*n)) then
!
! - Computation of S (k) <- distributed source
!
        do k=1,m
            sum=0.
            do j=1,n
            sum=sum+wgm (j)*xrm*(Rl+Tl)/2.*G (j)*(xI (k+1,j)+xI (k,j))/2.
            enddo
            S (k)=sum
        enddo
        goto 1
    endif
!
! - Computation of Ro_mult
!
    sum=0.
    do j=n/2+1,n
      xmu=xmm+xrm*xgm (j)
      sum=sum+wgm (j)*xrm*xI (1,j)*xmu/xmui
    enddo
    Ro_mult=sum
    return
end


subroutine lad
    implicit none
    real(8),parameter::Pi=3.141592653589793
    common /p/xLt,Rl,Tl,Rs,c,ild
    real(8) :: xLt,Rl,Tl,Rs,c
    integer :: ild
    common /ld/a_ld,b_ld,c_ld,d_ld
    real(8) :: a_ld,b_ld,c_ld,d_ld

    if (ild.eq.1) then
        a_ld=2./Pi
        b_ld=2./Pi
        c_ld=0.
        d_ld=0.
    else
        if (ild.eq.2) then
            a_ld=2./Pi
            b_ld=-2./Pi
            c_ld=0.
            d_ld=0.
        else
            if (ild.eq.3) then
                a_ld=2./Pi
                b_ld=0.
                c_ld=-2./Pi
                d_ld=0.
            else
                if (ild.eq.4) then
                    a_ld=2./Pi
                    b_ld=0.
                    c_ld=2./Pi
                    d_ld=0.
                else
                    a_ld=0.
                    b_ld=0.
                    c_ld=0.
                    d_ld=1.
                endif
            endif
        endif
    endif
    return
end
