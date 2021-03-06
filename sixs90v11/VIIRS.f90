subroutine viirs(iwa)
    implicit none
    common /sixs_ffu/ s(1501),wlinf,wlsup
    real(8) :: sr(16,1501),wli(16),wls(16)
    real(8) :: wlinf,wlsup,s
    integer :: iwa,l,i

!
!    VIIRS band M1
!
    data (sr(1,l),l=1,1501)/ 61*0.,  &
      0.7,1.0,1.0,1.0,1.0, &
      1.0,1.0,1.0,0.3, &
     1431*0./
!
!    VIIRS band M2
!
      data (sr(2,l),l=1,1501) / 74*0., &
      0.1,1.0,1.0,1.0,1.0,             &
      1.0,1.0,1.0,0.1,                 &
     1418*0./
!
!    VIIRS band M3
!
      data (sr(3,l),l=1,1501)/ 91*0.,  &
      0.3,1.0,1.0,1.0,1.0,             &
      1.0,1.0,1.0,0.7,                 &
     1401*0./
!
!    VIIRS band M4
!
      data (sr(4,l),l=1,1501)/ 118*0., &
      0.5,1.0,1.0,1.0,1.0,             &
      1.0,1.0,1.0,0.5,                 &
     1374*0./
!
!    VIIRS band M5
!
      data (sr(5,l),l=1,1501)/ 165*0., &
      0.7,1.0,1.0,1.0,1.0,             &
      1.0,1.0,1.0,0.3,                 &
     1327*0./
!
!    VIIRS band M6
!
      data (sr(6,l),l=1,1501)/ 195*0., &
      0.1,1.0,1.0,1.0,1.0,             &
      1.0,0.9,                         &
     1299*0./
!
!    VIIRS band M7
!
      data (sr(7,l),l=1,1501)/ 238*0., &
      0.3,1.0,1.0,1.0,1.0,             &
      1.0,1.0,1.0,1.0,1.0,             &
      1.0,1.0,1.0,1.0,1.0,             &
      1.0,0.3,                         &
     1246*0./
!
!    VIIRS band M8
!
     data (sr(8,l),l=1,1501)/ 392*0.,  &
      0.5,1.0,1.0,1.0,1.0,             &
      1.0,1.0,1.0,0.5,                 &
     1100*0./
!
!    VIIRS band M9
!
     data (sr(9,l),l=1,1501)/ 448*0.,  &
      0.3,1.0,1.0,1.0,1.0,             &
      1.0,0.7,                         &
     1046*0./
!
!    VIIRS band M10
!
     data (sr(10,l),l=1,1501)/ 532*0., &
      0.5,1.0,1.0,1.0,1.0,             &
      1.0,1.0,1.0,1.0,1.0,             &
      1.0,1.0,1.0,1.0,1.0,             &
      1.0,1.0,1.0,1.0,1.0,             &
      1.0,1.0,1.0,1.0,0.5,             &
     944*0./
!
!    VIIRS band M11
!
     data (sr(11,l),l=1,1501)/ 790*0., &
      0.5,1.0,1.0,1.0,1.0,             &
      1.0,1.0,1.0,1.0,1.0,             &
      1.0,1.0,1.0,1.0,1.0,             &
      1.0,1.0,1.0,1.0,1.0,             &
      0.5,                             &
     690*0./
!
!    VIIRS band M12
!
     data (sr(12,l),l=1,1501)/ 1344*0., &
      0.5,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,0.5,                      &
     84*0./
!
!    VIIRS band I1
!
     data (sr(13,l),l=1,1501)/ 140*0.,  &
      0.5,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,0.5,                      &
     1328*0./
!
!    VIIRS band I2
!
     data (sr(14,l),l=1,1501)/ 238*0.,  &
      0.3,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,0.3,                          &
     1246*0./
!
!    VIIRS band I3
!
     data (sr(15,l),l=1,1501)/ 532*0.,  &
      0.5,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,0.5,              &
     944*0./
!
!    VIIRS band I4
!
     data (sr(16,l),l=1,1501)/ 1320*0., &
      0.5,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,1.0,1.0,1.0,              &
      1.0,1.0,0.5,                      &
     28*0./

!
      wli(1)=0.4025
      wls(1)=0.4225
      wli(2)=0.4350
      wls(2)=0.4550
      wli(3)=0.4775
      wls(3)=0.4975
      wli(4)=0.5450
      wls(4)=0.5650
      wli(5)=0.6625
      wls(5)=0.6825
      wli(6)=0.7375
      wls(6)=0.7525
      wli(7)=0.8450
      wls(7)=0.8850
      wli(8)=1.2300
      wls(8)=1.2500
      wli(9)=1.3700
      wls(9)=1.3850
      wli(10)=1.5800
      wls(10)=1.6400
      wli(11)=2.2250
      wls(11)=2.2750
      wli(12)=3.6100
      wls(12)=3.7900
      wli(13)=0.6000
      wls(13)=0.6800
      wli(14)=0.8450
      wls(14)=0.8850
      wli(15)=1.5800
      wls(15)=1.6400
      wli(16)=3.5500
      wls(16)=3.9300
      do i=1,1501
        s(i)=sr(iwa,i)
      enddo
      wlinf=wli(iwa)
      wlsup=wls(iwa)
     return
end
