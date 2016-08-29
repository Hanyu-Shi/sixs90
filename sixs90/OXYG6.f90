subroutine oxyg6(a,inu)
    implicit none
    real(8) :: a(8)
    real(8) :: acr(8,256)
    integer :: inu,j,k,i

     data ((acr(k,j),k=1,8),j=  1,  8) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15300e+05, 0.15310e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15310e+05, 0.15320e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15320e+05, 0.15330e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15330e+05, 0.15340e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15340e+05, 0.15350e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15350e+05, 0.15360e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15360e+05, 0.15370e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15370e+05, 0.15380e+05/
     data ((acr(k,j),k=1,8),j=  9, 16) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15380e+05, 0.15390e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15390e+05, 0.15400e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15400e+05, 0.15410e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15410e+05, 0.15420e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15420e+05, 0.15430e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15430e+05, 0.15440e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15440e+05, 0.15450e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15450e+05, 0.15460e+05/
     data ((acr(k,j),k=1,8),j= 17, 24) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15460e+05, 0.15470e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15470e+05, 0.15480e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15480e+05, 0.15490e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15490e+05, 0.15500e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15500e+05, 0.15510e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15510e+05, 0.15520e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15520e+05, 0.15530e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15530e+05, 0.15540e+05/
     data ((acr(k,j),k=1,8),j= 25, 32) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15540e+05, 0.15550e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15550e+05, 0.15560e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15560e+05, 0.15570e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15570e+05, 0.15580e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15580e+05, 0.15590e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15590e+05, 0.15600e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15600e+05, 0.15610e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15610e+05, 0.15620e+05/
     data ((acr(k,j),k=1,8),j= 33, 40) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15620e+05, 0.15630e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15630e+05, 0.15640e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15640e+05, 0.15650e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15650e+05, 0.15660e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15660e+05, 0.15670e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15670e+05, 0.15680e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15680e+05, 0.15690e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15690e+05, 0.15700e+05/
     data ((acr(k,j),k=1,8),j= 41, 48) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15700e+05, 0.15710e+05,                              &
     0.15960e-07, 0.18194e-01, 0.38561e-01,-0.15424e-03, 0.36533e-01,    &
    -0.15016e-03, 0.15710e+05, 0.15720e+05,                              &
     0.15504e-07, 0.18194e-01, 0.38520e-01,-0.15408e-03, 0.36492e-01,    &
    -0.15000e-03, 0.15720e+05, 0.15730e+05,                              &
     0.86149e-07, 0.36385e-01, 0.34090e-01,-0.13636e-03, 0.32063e-01,    &
    -0.13228e-03, 0.15730e+05, 0.15740e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15740e+05, 0.15750e+05,                              &
     0.22012e-06, 0.36385e-01, 0.29911e-01,-0.11964e-03, 0.27883e-01,    &
    -0.11556e-03, 0.15750e+05, 0.15760e+05,                              &
     0.52469e-06, 0.36384e-01, 0.26002e-01,-0.10401e-03, 0.23974e-01,    &
    -0.99924e-04, 0.15760e+05, 0.15770e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15770e+05, 0.15780e+05/
     data ((acr(k,j),k=1,8),j= 49, 56) /                                 &
     0.11652e-05, 0.37236e-01, 0.22365e-01,-0.89458e-04, 0.20337e-01,    &
    -0.85374e-04, 0.15780e+05, 0.15790e+05,                              &
     0.24086e-05, 0.38115e-01, 0.19000e-01,-0.75998e-04, 0.16972e-01,    &
    -0.71915e-04, 0.15790e+05, 0.15800e+05,                              &
     0.46300e-05, 0.39407e-01, 0.15907e-01,-0.63629e-04, 0.13880e-01,    &
    -0.59545e-04, 0.15800e+05, 0.15810e+05,                              &
     0.42478e-05, 0.19927e-01, 0.13110e-01,-0.52439e-04, 0.11082e-01,    &
    -0.48357e-04, 0.15810e+05, 0.15820e+05,                              &
     0.40121e-05, 0.19927e-01, 0.13066e-01,-0.52262e-04, 0.11038e-01,    &
    -0.48180e-04, 0.15820e+05, 0.15830e+05,                              &
     0.13653e-04, 0.41148e-01, 0.10543e-01,-0.42171e-04, 0.85153e-02,    &
    -0.38088e-04, 0.15830e+05, 0.15840e+05,                              &
     0.20849e-04, 0.43340e-01, 0.82714e-02,-0.33085e-04, 0.62412e-02,    &
    -0.28988e-04, 0.15840e+05, 0.15850e+05,                              &
     0.29295e-04, 0.44164e-01, 0.62748e-02,-0.25099e-04, 0.42470e-02,    &
    -0.21015e-04, 0.15850e+05, 0.15860e+05/
     data ((acr(k,j),k=1,8),j= 57, 64) /                                 &
     0.37680e-04, 0.45464e-01, 0.45528e-02,-0.18211e-04, 0.25250e-02,    &
    -0.14127e-04, 0.15860e+05, 0.15870e+05,                              &
     0.68735e-04, 0.70552e-01, 0.26882e-02,-0.10606e-04, 0.67182e-03,    &
    -0.66419e-05, 0.15870e+05, 0.15880e+05,                              &
     0.62284e-04, 0.71735e-01, 0.13301e-02,-0.52392e-05,-0.70196e-03,    &
    -0.11788e-05, 0.15880e+05, 0.15890e+05,                              &
     0.39741e-04, 0.77900e-01, 0.35779e-03,-0.11367e-05,-0.13542e-02,    &
     0.45537e-05, 0.15890e+05, 0.15900e+05,                              &
     0.16287e-04, 0.58682e-01, 0.31178e-03, 0.24589e-05, 0.53130e-03,    &
     0.13293e-04, 0.15900e+05, 0.15910e+05,                              &
     0.99855e-04, 0.14228e+00, 0.13456e-02,-0.16726e-05, 0.10356e-02,    &
     0.34749e-05, 0.15910e+05, 0.15920e+05,                              &
     0.19246e-03, 0.31716e+00, 0.57540e-02,-0.16318e-04, 0.52332e-02,    &
    -0.19378e-04, 0.15920e+05, 0.15930e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15930e+05, 0.15940e+05/
     data ((acr(k,j),k=1,8),j= 65, 72) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15940e+05, 0.15950e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15950e+05, 0.15960e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15960e+05, 0.15970e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15970e+05, 0.15980e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15980e+05, 0.15990e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.15990e+05, 0.16000e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16000e+05, 0.16010e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16010e+05, 0.16020e+05/
     data ((acr(k,j),k=1,8),j= 73, 80) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16020e+05, 0.16030e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16030e+05, 0.16040e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16040e+05, 0.16050e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16050e+05, 0.16060e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16060e+05, 0.16070e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16070e+05, 0.16080e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16080e+05, 0.16090e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16090e+05, 0.16100e+05/
     data ((acr(k,j),k=1,8),j= 81, 88) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16100e+05, 0.16110e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16110e+05, 0.16120e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16120e+05, 0.16130e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16130e+05, 0.16140e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16140e+05, 0.16150e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16150e+05, 0.16160e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16160e+05, 0.16170e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16170e+05, 0.16180e+05/
     data ((acr(k,j),k=1,8),j= 89, 96) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16180e+05, 0.16190e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16190e+05, 0.16200e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16200e+05, 0.16210e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16210e+05, 0.16220e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16220e+05, 0.16230e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16230e+05, 0.16240e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16240e+05, 0.16250e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16250e+05, 0.16260e+05/
     data ((acr(k,j),k=1,8),j= 97,104) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16260e+05, 0.16270e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16270e+05, 0.16280e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16280e+05, 0.16290e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16290e+05, 0.16300e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16300e+05, 0.16310e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16310e+05, 0.16320e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16320e+05, 0.16330e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16330e+05, 0.16340e+05/
     data ((acr(k,j),k=1,8),j=105,112) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16340e+05, 0.16350e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16350e+05, 0.16360e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16360e+05, 0.16370e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16370e+05, 0.16380e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16380e+05, 0.16390e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16390e+05, 0.16400e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16400e+05, 0.16410e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16410e+05, 0.16420e+05/
     data ((acr(k,j),k=1,8),j=113,120) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16420e+05, 0.16430e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16430e+05, 0.16440e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16440e+05, 0.16450e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16450e+05, 0.16460e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16460e+05, 0.16470e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16470e+05, 0.16480e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16480e+05, 0.16490e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16490e+05, 0.16500e+05/
     data ((acr(k,j),k=1,8),j=121,128) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16500e+05, 0.16510e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16510e+05, 0.16520e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16520e+05, 0.16530e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16530e+05, 0.16540e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16540e+05, 0.16550e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16550e+05, 0.16560e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16560e+05, 0.16570e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16570e+05, 0.16580e+05/
     data ((acr(k,j),k=1,8),j=129,136) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16580e+05, 0.16590e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16590e+05, 0.16600e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16600e+05, 0.16610e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16610e+05, 0.16620e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16620e+05, 0.16630e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16630e+05, 0.16640e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16640e+05, 0.16650e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16650e+05, 0.16660e+05/
     data ((acr(k,j),k=1,8),j=137,144) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16660e+05, 0.16670e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16670e+05, 0.16680e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16680e+05, 0.16690e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16690e+05, 0.16700e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16700e+05, 0.16710e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16710e+05, 0.16720e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16720e+05, 0.16730e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16730e+05, 0.16740e+05/
     data ((acr(k,j),k=1,8),j=145,152) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16740e+05, 0.16750e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16750e+05, 0.16760e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16760e+05, 0.16770e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16770e+05, 0.16780e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16780e+05, 0.16790e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16790e+05, 0.16800e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16800e+05, 0.16810e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16810e+05, 0.16820e+05/
     data ((acr(k,j),k=1,8),j=153,160) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16820e+05, 0.16830e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16830e+05, 0.16840e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16840e+05, 0.16850e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16850e+05, 0.16860e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16860e+05, 0.16870e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16870e+05, 0.16880e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16880e+05, 0.16890e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16890e+05, 0.16900e+05/
     data ((acr(k,j),k=1,8),j=161,168) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16900e+05, 0.16910e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16910e+05, 0.16920e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16920e+05, 0.16930e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16930e+05, 0.16940e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16940e+05, 0.16950e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16950e+05, 0.16960e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16960e+05, 0.16970e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16970e+05, 0.16980e+05/
     data ((acr(k,j),k=1,8),j=169,176) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16980e+05, 0.16990e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.16990e+05, 0.17000e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17000e+05, 0.17010e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17010e+05, 0.17020e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17020e+05, 0.17030e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17030e+05, 0.17040e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17040e+05, 0.17050e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17050e+05, 0.17060e+05/
     data ((acr(k,j),k=1,8),j=177,184) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17060e+05, 0.17070e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17070e+05, 0.17080e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17080e+05, 0.17090e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17090e+05, 0.17100e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17100e+05, 0.17110e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17110e+05, 0.17120e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17120e+05, 0.17130e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17130e+05, 0.17140e+05/
     data ((acr(k,j),k=1,8),j=185,192) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17140e+05, 0.17150e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17150e+05, 0.17160e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17160e+05, 0.17170e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17170e+05, 0.17180e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17180e+05, 0.17190e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17190e+05, 0.17200e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17200e+05, 0.17210e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17210e+05, 0.17220e+05/
     data ((acr(k,j),k=1,8),j=193,200) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17220e+05, 0.17230e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17230e+05, 0.17240e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17240e+05, 0.17250e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17250e+05, 0.17260e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17260e+05, 0.17270e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17270e+05, 0.17280e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17280e+05, 0.17290e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17290e+05, 0.17300e+05/
     data ((acr(k,j),k=1,8),j=201,208) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17300e+05, 0.17310e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17310e+05, 0.17320e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17320e+05, 0.17330e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17330e+05, 0.17340e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17340e+05, 0.17350e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17350e+05, 0.17360e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17360e+05, 0.17370e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17370e+05, 0.17380e+05/
     data ((acr(k,j),k=1,8),j=209,216) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17380e+05, 0.17390e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17390e+05, 0.17400e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17400e+05, 0.17410e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17410e+05, 0.17420e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17420e+05, 0.17430e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17430e+05, 0.17440e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17440e+05, 0.17450e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17450e+05, 0.17460e+05/
     data ((acr(k,j),k=1,8),j=217,224) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17460e+05, 0.17470e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17470e+05, 0.17480e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17480e+05, 0.17490e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17490e+05, 0.17500e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17500e+05, 0.17510e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17510e+05, 0.17520e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17520e+05, 0.17530e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17530e+05, 0.17540e+05/
     data ((acr(k,j),k=1,8),j=225,232) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17540e+05, 0.17550e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17550e+05, 0.17560e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17560e+05, 0.17570e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17570e+05, 0.17580e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17580e+05, 0.17590e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17590e+05, 0.17600e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17600e+05, 0.17610e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17610e+05, 0.17620e+05/
     data ((acr(k,j),k=1,8),j=233,240) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17620e+05, 0.17630e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17630e+05, 0.17640e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17640e+05, 0.17650e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17650e+05, 0.17660e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17660e+05, 0.17670e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17670e+05, 0.17680e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17680e+05, 0.17690e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17690e+05, 0.17700e+05/
     data ((acr(k,j),k=1,8),j=241,248) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17700e+05, 0.17710e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17710e+05, 0.17720e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17720e+05, 0.17730e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17730e+05, 0.17740e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17740e+05, 0.17750e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17750e+05, 0.17760e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17760e+05, 0.17770e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17770e+05, 0.17780e+05/
     data ((acr(k,j),k=1,8),j=249,256) /                                 &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17780e+05, 0.17790e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17790e+05, 0.17800e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17800e+05, 0.17810e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17810e+05, 0.17820e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17820e+05, 0.17830e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17830e+05, 0.17840e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17840e+05, 0.17850e+05,                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00,    &
     0.00000e+00, 0.17850e+05, 0.17860e+05/
    do i=1,8
        a(i)=acr(i,inu)
    enddo
    return
end