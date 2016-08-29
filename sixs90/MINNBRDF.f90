subroutine minnbrdf(par1,par2,mu,np,rm,brdfint)
    implicit none
    real(8) :: par1,par2,xmu,view
    integer :: mu,np,k,j
    real(8) :: rm(-mu:mu),brdfint(-mu:mu,np)
    xmu=rm(0)
    do k=1,np
        do j=1,mu
            view=rm(j)
            brdfint(j,k)=0.5*par2*(par1+1.)*((xmu*view)**(par1-1))
        enddo
    enddo
    return
end
