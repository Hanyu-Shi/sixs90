subroutine poslan(month,jday,tu,xlon,xlat,asol,phi0,avis,phiv)
    implicit none
    logical:: ier
    real(8) :: tu,xlon,xlat,asol,phi0,avis,phiv
    integer :: month,jday,iwr
    common/sixs_ier/iwr,ier

!   landsat5 definition
!   warning !!!
!   xlon and xlat are the coordinates of the scene center.

    avis=0.
    phiv=0.

    call possol(month,jday,tu,xlon,xlat,asol,phi0)
    if(ier)return
    return
end
