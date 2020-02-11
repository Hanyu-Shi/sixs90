module paramdef
    implicit none
    ! do not change
    integer,parameter::nt_p=30,mu_p=25,mu2_p=48,np_p=49,nfi_p=181,nquad_p=83
    integer,parameter::nt_p_max=100,nqmax_p=1001,nqdef_p=83
    ! Attention
    ! mu2_p has to be equal to (mu_p-1)*2
    ! nquad_p must always be odd
end module paramdef

module paramdef_highaccuracy
    implicit none
    ! do not change
    integer,parameter::nt_p=50,mu_p=75,mu2_p=148,np_p=149,nfi_p=181,nquad_p=200
    integer,parameter::nt_p_max=100,nqmax_p=1000,nqdef_p=83
    ! Attention
    ! mu2_p has to be equal to (mu_p-1)*2
    ! nquad_p must always be odd
end module paramdef_highaccuracy

module paramdef_normalaccuracy
    implicit none
    ! do not change
    integer,parameter::nt_p=30,mu_p=25,mu2_p=48,np_p=49,nfi_p=181,nquad_p=83
    integer,parameter::nt_p_max=100,nqmax_p=1001,nqdef_p=83
    ! Attention
    ! mu2_p has to be equal to (mu_p-1)*2
    ! nquad_p must always be odd
end module paramdef_normalaccuracy
