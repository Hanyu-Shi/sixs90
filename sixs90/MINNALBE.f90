subroutine minnalbe(par1,par2,brdfalb)
    implicit none
    real(8) :: par1,par2,brdfalb
    brdfalb=2.*par2/(par1+1.)
    return
end
