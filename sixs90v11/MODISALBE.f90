subroutine modisalbe(p1,p2,p3,brdfalb)
    implicit none
    real(8) :: p1,p2,p3,brdfalb
    brdfalb=p1+p2*0.189184-p3*1.377622
    return
end
