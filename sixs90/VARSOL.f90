subroutine varsol(jday,month,dsol)
    implicit none
    real(8) :: dsol,pi,om
    integer :: jday,month,j

!   calculation of the variability of the solar constant during the year.
!   jday is the number of the day in the month
!   dsol is a multiplicative factor to apply to the mean value of
!   solar constant

    if (month.le.2) goto 1
    if (month.gt.8) goto 2
    j=31*(month-1)-((month-1)/2)-2+jday
    goto 3
1   j=31*(month-1)+jday
    goto 3
2   j=31*(month-1)-((month-2)/2)-2+jday

3   pi=2.*acos (0.d0)
    om=(.9856*float(j-4))*pi/180.
    dsol=1./((1.-.01673*cos(om))**2)
    return
end
