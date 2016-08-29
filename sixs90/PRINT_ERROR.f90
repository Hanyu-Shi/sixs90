subroutine print_error(tex)
    implicit none
    !character *(*) tex
    character(*) :: tex
    logical:: ier
    integer :: iwr
    common/sixs_ier/iwr,ier
    ier = .TRUE.
    write(iwr,'(a)')tex
    return
end
