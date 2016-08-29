subroutine tropic
    implicit none
    integer :: i
    real(8) :: z1(34),p1(34),t1(34),wh1(34),wo1(34)
    real(8) :: z,p,t,wh,wo
    common /sixs_atm/z(34),p(34),t(34),wh(34),wo(34)
!
!   model: tropical mc clatchey
!
    data(z1(i),i=1, 34)/                                              &
        0.,    1.,    2.,    3.,    4.,    5.,    6.,    7.,    8.,   &
        9.,   10.,   11.,   12.,   13.,   14.,   15.,   16.,   17.,   &
       18.,   19.,   20.,   21.,   22.,   23.,   24.,   25.,   30.,   &
       35.,   40.,   45.,   50.,   70.,  100.,99999./
    data (p1(i),i=1,34)/                                              &
     1.013e+03,9.040e+02,8.050e+02,7.150e+02,6.330e+02,5.590e+02,     &
     4.920e+02,4.320e+02,3.780e+02,3.290e+02,2.860e+02,2.470e+02,     &
     2.130e+02,1.820e+02,1.560e+02,1.320e+02,1.110e+02,9.370e+01,     &
     7.890e+01,6.660e+01,5.650e+01,4.800e+01,4.090e+01,3.500e+01,     &
     3.000e+01,2.570e+01,1.220e+01,6.000e+00,3.050e+00,1.590e+00,     &
     8.540e-01,5.790e-02,3.000e-04,0.000e+00/
    data (t1(i),i=1,34)/                                              &
     3.000e+02,2.940e+02,2.880e+02,2.840e+02,2.770e+02,2.700e+02,     &
     2.640e+02,2.570e+02,2.500e+02,2.440e+02,2.370e+02,2.300e+02,     &
     2.240e+02,2.170e+02,2.100e+02,2.040e+02,1.970e+02,1.950e+02,     &
     1.990e+02,2.030e+02,2.070e+02,2.110e+02,2.150e+02,2.170e+02,     &
     2.190e+02,2.210e+02,2.320e+02,2.430e+02,2.540e+02,2.650e+02,     &
     2.700e+02,2.190e+02,2.100e+02,2.100e+02/
    data (wh1(i),i=1,34)/                                             &
     1.900e+01,1.300e+01,9.300e+00,4.700e+00,2.200e+00,1.500e+00,     &
     8.500e-01,4.700e-01,2.500e-01,1.200e-01,5.000e-02,1.700e-02,     &
     6.000e-03,1.800e-03,1.000e-03,7.600e-04,6.400e-04,5.600e-04,     &
     5.000e-04,4.900e-04,4.500e-04,5.100e-04,5.100e-04,5.400e-04,     &
     6.000e-04,6.700e-04,3.600e-04,1.100e-04,4.300e-05,1.900e-05,     &
     6.300e-06,1.400e-07,1.000e-09,0.000e+00/
    data (wo1(i),i=1,34)/                                             &
     5.600e-05,5.600e-05,5.400e-05,5.100e-05,4.700e-05,4.500e-05,     &
     4.300e-05,4.100e-05,3.900e-05,3.900e-05,3.900e-05,4.100e-05,     &
     4.300e-05,4.500e-05,4.500e-05,4.700e-05,4.700e-05,6.900e-05,     &
     9.000e-05,1.400e-04,1.900e-04,2.400e-04,2.800e-04,3.200e-04,     &
     3.400e-04,3.400e-04,2.400e-04,9.200e-05,4.100e-05,1.300e-05,     &
     4.300e-06,8.600e-08,4.300e-11,0.000e+00/

    do i=1,34
        z(i)=z1(i)
        p(i)=p1(i)
        t(i)=t1(i)
        wh(i)=wh1(i)
        wo(i)=wo1(i)
    enddo
    return
end
