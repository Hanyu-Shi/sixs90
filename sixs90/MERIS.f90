subroutine meris(iwa)
    implicit none
    real(8) :: s,wlinf,wlsup
    common /sixs_ffu/ s(1501),wlinf,wlsup
    real(8) :: sr(15,1501),wli(15),wls(15)
    integer :: iwa,l,i

    DATA (SR(1,L),L=1,1501)/ 61*0.,                                &
     .0000, .0002, .6080, 1.0000, 1.0000, 1.0000, .3829, .0000,    &
     1432*0./

    DATA (SR(2,L),L=1,1501)/ 73*0.,                                &
     .0000, .0002, .6396, 1.0000, 1.0000, .9990, .3484, .0000,     &
     1420*0./

    DATA (SR(3,L),L=1,1501)/ 92*0.,                                &
     .0000, .0003, .7031, 1.0000, 1.0000, .9990, .2829, .0000,     &
     1401*0./

    DATA (SR(4,L),L=1,1501)/ 100*0.,                               &
     .0000, .0004, .7315, 1.0000, 1.0000, .9990, .2545, .0000,     &
     1393*0./

    DATA (SR(5,L),L=1,1501)/ 120*0.,                               &
     .0000, .0007, .7968, 1.0000, 1.0000, .9990, .1922, .0000,     &
     1373*0./

    DATA (SR(6,L),L=1,1501)/ 144*0.,                               &
     .0000, .0013, .8530, 1.0000, 1.0000, .9980, .1400, .0000,     &
     1349*0./

    DATA (SR(7,L),L=1,1501)/ 162*0.,                               &
     .0000, .0016, .8760, 1.0000, .9990, .9970, .1210, .0000,      &
     1331*0./

    DATA (SR(8,L),L=1,1501)/ 169*0.,                               &
     .0000, .0016, .8800, 1.0000, .9980, .1190, .0000,             &
     1325*0./

    DATA (SR(9,L),L=1,1501)/ 180*0.,                               &
     .0000, .1690, .9990, 1.0000, .9990, .8330, .0008, .0000,      &
     1313*0./

    DATA (SR(10,L),L=1,1501)/ 198*0.,                              &
     .0000, .0011, .8680, 1.0000, .9990, .1360, .0000,             &
     1296*0./

    DATA (SR(11,L),L=1,1501)/ 202*0.,                              &
     .0000, .1430, 1.0000, .1390, .0000,                           &
     1294*0./

    DATA (SR(12,L),L=1,1501)/ 207*0.,                              &
     .0000, .1340, .9990, 1.0000, .9990, .9990, .9990, .8750,      &
     .0012, .0000,                                                 &
     1284*0./

    DATA (SR(13,L),L=1,1501)/ 240*0.,                              &
     .0000, .0002, .7683, 1.0000, .9990, .9980, .9980, .9970,      &
     .9970, .9960, .2618, .0000,                                   &
     1249*0./

    DATA (SR(14,L),L=1,1501)/ 250*0.,                              &
     .0000, .0002, .7294, 1.0000, 1.0000, .9990, .2877, .0000,     &
     1243*0./

    DATA (SR(15,L),L=1,1501)/ 256*0.,                              &
     .0000, .0001, .7029, 1.0000, .9990, .9980, .3152, .0000,      &
     1237*0./

! channel 1 lower and upper wavelength
    wli(1)=0.4025
    wls(1)=0.42
! channel 2 lower and upper wavelength
    wli(2)=0.4325
    wls(2)=0.45
! channel 3 lower and upper wavelength
    wli(3)=0.48
    wls(3)=0.4975
! channel 4 lower and upper wavelength
    wli(4)=0.5
    wls(4)=0.5175
! channel 5 lower and upper wavelength
    wli(5)=0.55
    wls(5)=0.5675
! channel 6 lower and upper wavelength
    wli(6)=0.61
    wls(6)=0.6275
! channel 7 lower and upper wavelength
    wli(7)=0.655
    wls(7)=0.6725
! channel 8 lower and upper wavelength
    wli(8)=0.6725
    wls(8)=0.6875
! channel 9 lower and upper wavelength
    wli(9)=0.7
    wls(9)=0.7175
! channel 10 lower and upper wavelength
    wli(10)=0.7450
    wls(10)=0.76
! channel 11 lower and upper wavelength
    wli(11)=0.755
    wls(11)=0.765
! channel 12 lower and upper wavelength
    wli(12)=0.7675
    wls(12)=0.79
! channel 13 lower and upper wavelength
    wli(13)=0.85
    wls(13)=0.8775
! channel 14 lower and upper wavelength
    wli(14)=0.875
    wls(14)=0.8925
! channel 15 lower and upper wavelength
    wli(15)=0.89
    wls(15)=0.9075

    do i=1,1501
        s(i)=sr(iwa,i)
    enddo

    wlinf=wli(iwa)
    wlsup=wls(iwa)
    return
end
