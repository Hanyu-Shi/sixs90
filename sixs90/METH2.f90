subroutine meth2(a,inu)
    implicit none
    real(8) :: a(8)
    real(8) :: acr(8,256)
    integer :: inu,j,k,i
!   methane (5060 - 7610 cm-1)
!
     data ((acr(k,j),k=1,8),j=  1,  8) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.50600e+04, 0.50700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.50700e+04, 0.50800e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.50800e+04, 0.50900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.50900e+04, 0.51000e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.51000e+04, 0.51100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.51100e+04, 0.51200e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.51200e+04, 0.51300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.51300e+04, 0.51400e+04/
     data ((acr(k,j),k=1,8),j=  9, 16) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.51400e+04, 0.51500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.51500e+04, 0.51600e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.51600e+04, 0.51700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.51700e+04, 0.51800e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.51800e+04, 0.51900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.51900e+04, 0.52000e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.52000e+04, 0.52100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.52100e+04, 0.52200e+04/
     data ((acr(k,j),k=1,8),j= 17, 24) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.52200e+04, 0.52300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.52300e+04, 0.52400e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.52400e+04, 0.52500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.52500e+04, 0.52600e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.52600e+04, 0.52700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.52700e+04, 0.52800e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.52800e+04, 0.52900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.52900e+04, 0.53000e+04/
     data ((acr(k,j),k=1,8),j= 25, 32) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.53000e+04, 0.53100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.53100e+04, 0.53200e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.53200e+04, 0.53300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.53300e+04, 0.53400e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.53400e+04, 0.53500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.53500e+04, 0.53600e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.53600e+04, 0.53700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.53700e+04, 0.53800e+04/
     data ((acr(k,j),k=1,8),j= 33, 40) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.53800e+04, 0.53900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.53900e+04, 0.54000e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.54000e+04, 0.54100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.54100e+04, 0.54200e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.54200e+04, 0.54300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.54300e+04, 0.54400e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.54400e+04, 0.54500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.54500e+04, 0.54600e+04/
     data ((acr(k,j),k=1,8),j= 41, 48) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.54600e+04, 0.54700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.54700e+04, 0.54800e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.54800e+04, 0.54900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.54900e+04, 0.55000e+04,                           &
     0.50197e+00, 0.25265e+00, 0.33519e-02,-0.24607e-04, 0.68955e-03, &
    -0.20482e-04, 0.55000e+04, 0.55100e+04,                           &
     0.23899e+01, 0.60596e+00, 0.27671e-04,-0.12307e-04,-0.33058e-02, &
    -0.52945e-05, 0.55100e+04, 0.55200e+04,                           &
     0.24379e+01, 0.30699e+00,-0.60867e-03,-0.90704e-05,-0.32892e-02, &
    -0.50115e-05, 0.55200e+04, 0.55300e+04,                           &
     0.21592e+01, 0.38949e+00,-0.23556e-02,-0.33022e-05,-0.52838e-02, &
     0.24513e-05, 0.55300e+04, 0.55400e+04/
     data ((acr(k,j),k=1,8),j= 49, 56) /                              &
     0.23029e+01, 0.15736e+00,-0.35795e-02, 0.21673e-05,-0.59680e-02, &
     0.60863e-05, 0.55400e+04, 0.55500e+04,                           &
     0.19540e+01, 0.11711e+00,-0.44087e-02, 0.51030e-05,-0.68665e-02, &
     0.91701e-05, 0.55500e+04, 0.55600e+04,                           &
     0.11950e+01, 0.29396e+00,-0.10618e-02, 0.34067e-05,-0.24231e-02, &
    -0.25820e-05, 0.55600e+04, 0.55700e+04,                           &
     0.48095e+01, 0.12465e+01, 0.19344e-02,-0.15456e-04,-0.68788e-03, &
    -0.12870e-04, 0.55700e+04, 0.55800e+04,                           &
     0.11674e+02, 0.15114e+01,-0.25504e-02, 0.82500e-06,-0.45912e-02, &
     0.22777e-05, 0.55800e+04, 0.55900e+04,                           &
     0.23702e+01, 0.84024e+00,-0.82688e-03, 0.29123e-05,-0.20134e-02, &
    -0.50547e-05, 0.55900e+04, 0.56000e+04,                           &
     0.34064e+01, 0.11326e+01,-0.25002e-02, 0.89997e-06,-0.48837e-02, &
     0.26082e-05, 0.56000e+04, 0.56100e+04,                           &
     0.17392e+01, 0.31991e+00,-0.45126e-02, 0.86603e-05,-0.61278e-02, &
     0.93747e-05, 0.56100e+04, 0.56200e+04/
     data ((acr(k,j),k=1,8),j= 57, 64) /                              &
     0.43480e+01, 0.58786e+00,-0.43137e-02, 0.46334e-05,-0.70133e-02, &
     0.95290e-05, 0.56200e+04, 0.56300e+04,                           &
     0.66586e+01, 0.18023e+01,-0.19866e-02,-0.15163e-05,-0.39533e-02, &
    -0.85993e-06, 0.56300e+04, 0.56400e+04,                           &
     0.43959e+01, 0.91267e+00,-0.18359e-02,-0.26379e-05,-0.45115e-02, &
     0.13137e-05, 0.56400e+04, 0.56500e+04,                           &
     0.29732e+01, 0.90097e+00,-0.45783e-03,-0.86540e-05,-0.35958e-02, &
    -0.28318e-05, 0.56500e+04, 0.56600e+04,                           &
     0.27758e+01, 0.82876e+00, 0.19221e-02,-0.14165e-04,-0.10447e-02, &
    -0.11375e-04, 0.56600e+04, 0.56700e+04,                           &
     0.14345e+01, 0.92016e+00, 0.46313e-02,-0.19911e-04, 0.18272e-02, &
    -0.20246e-04, 0.56700e+04, 0.56800e+04,                           &
     0.10486e+01, 0.44650e+00,-0.38086e-03, 0.56985e-05,-0.18929e-02, &
    -0.11857e-05, 0.56800e+04, 0.56900e+04,                           &
     0.94797e+00, 0.35596e+00,-0.22904e-03, 0.84588e-05,-0.16889e-02, &
    -0.37956e-06, 0.56900e+04, 0.57000e+04/
     data ((acr(k,j),k=1,8),j= 65, 72) /                              &
     0.93528e+00, 0.41431e+00, 0.12800e-02,-0.87918e-06,-0.14444e-02, &
    -0.47052e-05, 0.57000e+04, 0.57100e+04,                           &
     0.16454e+01, 0.57474e+00, 0.26463e-02,-0.13124e-04,-0.37027e-03, &
    -0.12452e-04, 0.57100e+04, 0.57200e+04,                           &
     0.20351e+01, 0.34637e+00, 0.26986e-02,-0.19613e-04,-0.46781e-03, &
    -0.14709e-04, 0.57200e+04, 0.57300e+04,                           &
     0.22445e+01, 0.69714e+00, 0.19373e-02,-0.16752e-04,-0.85242e-03, &
    -0.12956e-04, 0.57300e+04, 0.57400e+04,                           &
     0.32153e+01, 0.75956e+00, 0.13362e-02,-0.14189e-04,-0.10894e-02, &
    -0.11571e-04, 0.57400e+04, 0.57500e+04,                           &
     0.34596e+01, 0.73536e+00, 0.12876e-02,-0.12547e-04,-0.45370e-03, &
    -0.12912e-04, 0.57500e+04, 0.57600e+04,                           &
     0.49414e+01, 0.10526e+01, 0.64547e-03,-0.10024e-04,-0.96857e-03, &
    -0.11592e-04, 0.57600e+04, 0.57700e+04,                           &
     0.75119e+01, 0.14100e+01, 0.28819e-03,-0.74440e-05,-0.12786e-02, &
    -0.10307e-04, 0.57700e+04, 0.57800e+04/
     data ((acr(k,j),k=1,8),j= 73, 80) /                              &
     0.75235e+01, 0.13230e+01,-0.49407e-03,-0.68773e-05,-0.26290e-02, &
    -0.66380e-05, 0.57800e+04, 0.57900e+04,                           &
     0.10241e+02, 0.16323e+01,-0.13776e-02,-0.53214e-05,-0.37904e-02, &
    -0.23185e-05, 0.57900e+04, 0.58000e+04,                           &
     0.91281e+01, 0.93597e+00,-0.35876e-02, 0.23809e-05,-0.59879e-02, &
     0.59495e-05, 0.58000e+04, 0.58100e+04,                           &
     0.39872e+01, 0.67670e+00,-0.44503e-02, 0.68501e-05,-0.67178e-02, &
     0.97535e-05, 0.58100e+04, 0.58200e+04,                           &
     0.29318e+01, 0.24918e+00,-0.51233e-02, 0.77265e-05,-0.76555e-02, &
     0.11663e-04, 0.58200e+04, 0.58300e+04,                           &
     0.36596e+01, 0.91227e+00,-0.26988e-02, 0.23063e-05,-0.47124e-02, &
     0.25366e-05, 0.58300e+04, 0.58400e+04,                           &
     0.69904e+01, 0.93296e+00,-0.37966e-02, 0.61051e-05,-0.55025e-02, &
     0.73047e-05, 0.58400e+04, 0.58500e+04,                           &
     0.12815e+02, 0.19354e+01,-0.21432e-02, 0.81298e-06,-0.42246e-02, &
     0.22902e-05, 0.58500e+04, 0.58600e+04/
     data ((acr(k,j),k=1,8),j= 81, 88) /                              &
     0.10246e+02, 0.21238e+01, 0.46316e-04,-0.60726e-05,-0.22300e-02, &
    -0.53836e-05, 0.58600e+04, 0.58700e+04,                           &
     0.95647e+01, 0.19821e+01, 0.12796e-02,-0.74626e-05,-0.92616e-03, &
    -0.10886e-04, 0.58700e+04, 0.58800e+04,                           &
     0.48776e+01, 0.17308e+01, 0.12461e-02,-0.81397e-05,-0.56127e-04, &
    -0.13882e-04, 0.58800e+04, 0.58900e+04,                           &
     0.57337e+01, 0.87446e+00, 0.38286e-02,-0.17040e-04, 0.69410e-03, &
    -0.15730e-04, 0.58900e+04, 0.59000e+04,                           &
     0.74720e+01, 0.73641e+00, 0.28865e-02,-0.20558e-04,-0.66892e-03, &
    -0.13580e-04, 0.59000e+04, 0.59100e+04,                           &
     0.70944e+01, 0.66112e+00, 0.24095e-02,-0.21458e-04,-0.99935e-03, &
    -0.13636e-04, 0.59100e+04, 0.59200e+04,                           &
     0.81168e+01, 0.58477e+00, 0.13314e-02,-0.17142e-04,-0.10021e-02, &
    -0.13588e-04, 0.59200e+04, 0.59300e+04,                           &
     0.10972e+02, 0.47138e+00,-0.30541e-03,-0.96826e-05,-0.14478e-02, &
    -0.96570e-05, 0.59300e+04, 0.59400e+04/
     data ((acr(k,j),k=1,8),j= 89, 96) /                              &
     0.73055e+01, 0.38828e+00,-0.15833e-02,-0.39472e-05,-0.24534e-02, &
    -0.52105e-05, 0.59400e+04, 0.59500e+04,                           &
     0.18257e+00, 0.15141e+00, 0.32545e-02,-0.25286e-04, 0.32305e-03, &
    -0.19965e-04, 0.59500e+04, 0.59600e+04,                           &
     0.88030e+01, 0.28605e+00,-0.31163e-02, 0.13762e-05,-0.40606e-02, &
     0.10755e-05, 0.59600e+04, 0.59700e+04,                           &
     0.68729e+01, 0.20493e+00,-0.40682e-02, 0.54680e-05,-0.47670e-02, &
     0.43804e-05, 0.59700e+04, 0.59800e+04,                           &
     0.29955e+01, 0.47987e+00,-0.24740e-02, 0.99034e-05,-0.84087e-03, &
    -0.75046e-05, 0.59800e+04, 0.59900e+04,                           &
     0.14364e+02, 0.14135e+01, 0.66597e-02,-0.31277e-04, 0.45742e-02, &
    -0.31988e-04, 0.59900e+04, 0.60000e+04,                           &
     0.65167e+02, 0.97529e+00,-0.19616e-02,-0.22633e-05,-0.46527e-02, &
     0.10271e-05, 0.60000e+04, 0.60100e+04,                           &
     0.35647e+01, 0.20344e+00,-0.52252e-02, 0.10688e-04,-0.57674e-02, &
     0.80346e-05, 0.60100e+04, 0.60200e+04/
     data ((acr(k,j),k=1,8),j= 97,104) /                              &
     0.34141e+01, 0.14357e+00,-0.54304e-02, 0.95186e-05,-0.72736e-02, &
     0.11430e-04, 0.60200e+04, 0.60300e+04,                           &
     0.76350e+01, 0.43320e+00,-0.31157e-02, 0.13276e-04,-0.12770e-02, &
     0.17678e-06, 0.60300e+04, 0.60400e+04,                           &
     0.16118e+02, 0.28374e+00,-0.41664e-02, 0.49599e-05,-0.56931e-02, &
     0.57765e-05, 0.60400e+04, 0.60500e+04,                           &
     0.19177e+02, 0.33833e+00,-0.35167e-02, 0.17094e-05,-0.61198e-02, &
     0.66957e-05, 0.60500e+04, 0.60600e+04,                           &
     0.16850e+02, 0.88879e+00,-0.21460e-02,-0.32028e-05,-0.42847e-02, &
     0.27277e-06, 0.60600e+04, 0.60700e+04,                           &
     0.21077e+02, 0.25352e+00,-0.93983e-03,-0.90824e-05,-0.38178e-02, &
    -0.28130e-05, 0.60700e+04, 0.60800e+04,                           &
     0.14763e+02, 0.24370e+00, 0.58706e-03,-0.14814e-04,-0.30455e-02, &
    -0.58930e-05, 0.60800e+04, 0.60900e+04,                           &
     0.12160e+02, 0.31237e+00, 0.25210e-02,-0.21431e-04,-0.11996e-02, &
    -0.11964e-04, 0.60900e+04, 0.61000e+04/
     data ((acr(k,j),k=1,8),j=105,112) /                              &
     0.10601e+02, 0.42935e+00, 0.43443e-02,-0.27697e-04, 0.31712e-03, &
    -0.17641e-04, 0.61000e+04, 0.61100e+04,                           &
     0.77832e+01, 0.54373e+00, 0.57906e-02,-0.26532e-04, 0.14380e-02, &
    -0.19125e-04, 0.61100e+04, 0.61200e+04,                           &
     0.57902e+01, 0.11869e+01, 0.54008e-02,-0.18503e-04, 0.23381e-02, &
    -0.20741e-04, 0.61200e+04, 0.61300e+04,                           &
     0.40292e+01, 0.58128e+00, 0.79694e-02,-0.24981e-04, 0.32731e-02, &
    -0.22237e-04, 0.61300e+04, 0.61400e+04,                           &
     0.30446e+01, 0.82071e+00, 0.77456e-02,-0.13929e-04, 0.53212e-02, &
    -0.25965e-04, 0.61400e+04, 0.61500e+04,                           &
     0.19431e+01, 0.98831e+00, 0.60640e-02,-0.22019e-04, 0.32998e-02, &
    -0.24905e-04, 0.61500e+04, 0.61600e+04,                           &
     0.97862e+00, 0.77724e+00, 0.73198e-02,-0.33144e-04, 0.43341e-02, &
    -0.32300e-04, 0.61600e+04, 0.61700e+04,                           &
     0.41035e+00, 0.31003e+00, 0.44331e-02,-0.27381e-04, 0.19194e-02, &
    -0.24710e-04, 0.61700e+04, 0.61800e+04/
     data ((acr(k,j),k=1,8),j=113,120) /                              &
     0.19122e+00, 0.12211e+00, 0.46368e-02,-0.30314e-04, 0.17315e-02, &
    -0.25376e-04, 0.61800e+04, 0.61900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.61900e+04, 0.62000e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.62000e+04, 0.62100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.62100e+04, 0.62200e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.62200e+04, 0.62300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.62300e+04, 0.62400e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.62400e+04, 0.62500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.62500e+04, 0.62600e+04/
     data ((acr(k,j),k=1,8),j=121,128) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.62600e+04, 0.62700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.62700e+04, 0.62800e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.62800e+04, 0.62900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.62900e+04, 0.63000e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.63000e+04, 0.63100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.63100e+04, 0.63200e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.63200e+04, 0.63300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.63300e+04, 0.63400e+04/
     data ((acr(k,j),k=1,8),j=129,136) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.63400e+04, 0.63500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.63500e+04, 0.63600e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.63600e+04, 0.63700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.63700e+04, 0.63800e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.63800e+04, 0.63900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.63900e+04, 0.64000e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.64000e+04, 0.64100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.64100e+04, 0.64200e+04/
     data ((acr(k,j),k=1,8),j=137,144) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.64200e+04, 0.64300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.64300e+04, 0.64400e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.64400e+04, 0.64500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.64500e+04, 0.64600e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.64600e+04, 0.64700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.64700e+04, 0.64800e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.64800e+04, 0.64900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.64900e+04, 0.65000e+04/
     data ((acr(k,j),k=1,8),j=145,152) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.65000e+04, 0.65100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.65100e+04, 0.65200e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.65200e+04, 0.65300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.65300e+04, 0.65400e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.65400e+04, 0.65500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.65500e+04, 0.65600e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.65600e+04, 0.65700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.65700e+04, 0.65800e+04/
     data ((acr(k,j),k=1,8),j=153,160) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.65800e+04, 0.65900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.65900e+04, 0.66000e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.66000e+04, 0.66100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.66100e+04, 0.66200e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.66200e+04, 0.66300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.66300e+04, 0.66400e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.66400e+04, 0.66500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.66500e+04, 0.66600e+04/
     data ((acr(k,j),k=1,8),j=161,168) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.66600e+04, 0.66700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.66700e+04, 0.66800e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.66800e+04, 0.66900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.66900e+04, 0.67000e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.67000e+04, 0.67100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.67100e+04, 0.67200e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.67200e+04, 0.67300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.67300e+04, 0.67400e+04/
     data ((acr(k,j),k=1,8),j=169,176) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.67400e+04, 0.67500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.67500e+04, 0.67600e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.67600e+04, 0.67700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.67700e+04, 0.67800e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.67800e+04, 0.67900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.67900e+04, 0.68000e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.68000e+04, 0.68100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.68100e+04, 0.68200e+04/
     data ((acr(k,j),k=1,8),j=177,184) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.68200e+04, 0.68300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.68300e+04, 0.68400e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.68400e+04, 0.68500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.68500e+04, 0.68600e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.68600e+04, 0.68700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.68700e+04, 0.68800e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.68800e+04, 0.68900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.68900e+04, 0.69000e+04/
     data ((acr(k,j),k=1,8),j=185,192) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.69000e+04, 0.69100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.69100e+04, 0.69200e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.69200e+04, 0.69300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.69300e+04, 0.69400e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.69400e+04, 0.69500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.69500e+04, 0.69600e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.69600e+04, 0.69700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.69700e+04, 0.69800e+04/
     data ((acr(k,j),k=1,8),j=193,200) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.69800e+04, 0.69900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.69900e+04, 0.70000e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.70000e+04, 0.70100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.70100e+04, 0.70200e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.70200e+04, 0.70300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.70300e+04, 0.70400e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.70400e+04, 0.70500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.70500e+04, 0.70600e+04/
     data ((acr(k,j),k=1,8),j=201,208) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.70600e+04, 0.70700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.70700e+04, 0.70800e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.70800e+04, 0.70900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.70900e+04, 0.71000e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.71000e+04, 0.71100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.71100e+04, 0.71200e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.71200e+04, 0.71300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.71300e+04, 0.71400e+04/
     data ((acr(k,j),k=1,8),j=209,216) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.71400e+04, 0.71500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.71500e+04, 0.71600e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.71600e+04, 0.71700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.71700e+04, 0.71800e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.71800e+04, 0.71900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.71900e+04, 0.72000e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.72000e+04, 0.72100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.72100e+04, 0.72200e+04/
     data ((acr(k,j),k=1,8),j=217,224) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.72200e+04, 0.72300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.72300e+04, 0.72400e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.72400e+04, 0.72500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.72500e+04, 0.72600e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.72600e+04, 0.72700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.72700e+04, 0.72800e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.72800e+04, 0.72900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.72900e+04, 0.73000e+04/
     data ((acr(k,j),k=1,8),j=225,232) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.73000e+04, 0.73100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.73100e+04, 0.73200e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.73200e+04, 0.73300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.73300e+04, 0.73400e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.73400e+04, 0.73500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.73500e+04, 0.73600e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.73600e+04, 0.73700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.73700e+04, 0.73800e+04/
     data ((acr(k,j),k=1,8),j=233,240) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.73800e+04, 0.73900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.73900e+04, 0.74000e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.74000e+04, 0.74100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.74100e+04, 0.74200e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.74200e+04, 0.74300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.74300e+04, 0.74400e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.74400e+04, 0.74500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.74500e+04, 0.74600e+04/
     data ((acr(k,j),k=1,8),j=241,248) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.74600e+04, 0.74700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.74700e+04, 0.74800e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.74800e+04, 0.74900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.74900e+04, 0.75000e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.75000e+04, 0.75100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.75100e+04, 0.75200e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.75200e+04, 0.75300e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.75300e+04, 0.75400e+04/
     data ((acr(k,j),k=1,8),j=249,256) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.75400e+04, 0.75500e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.75500e+04, 0.75600e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.75600e+04, 0.75700e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.75700e+04, 0.75800e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.75800e+04, 0.75900e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.75900e+04, 0.76000e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.76000e+04, 0.76100e+04,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.76100e+04, 0.76200e+04/
!
    do i=1,8
        a(i)=acr(i,inu)
    enddo
!
    return
end
