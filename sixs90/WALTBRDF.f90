subroutine waltbrdf(a,ap,b,c,mu,np,rm,rp,brdfint)
! this model can be found in applied optics vol 24 / no 3/ pp 383-387
! but it has to be modified (slightly) to match reciprocity principle
    implicit none
    integer :: mu,np,k,j
    real(8) :: rm(-mu:mu),rp(np),brdfint(-mu:mu,np)
    real(8) :: a,ap,b,c
    real(8) :: xmu,ts,view,tv,fi,phi

    xmu=rm(0)
    ts=acos(xmu)
    do k=1,np
        do j=1,mu
            view=rm(j)
            tv=acos(view)
            if (j.eq.mu) then
                fi=rm(-mu)
            else
                fi=rp(k)+rm(-mu)
            endif
            phi=fi
            brdfint(j,k)=a*(ts*ts*tv*tv)+ap*(ts*ts+tv*tv)+b*ts*tv*cos(phi)+c
        enddo
    enddo
    return
end
