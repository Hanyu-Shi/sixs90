subroutine oxyg5(a,inu)
    implicit none
    real(8) :: a(8)
    real(8) :: acr(8,256)
    integer :: inu,j,k,i
!   oxygen (12740 - 15290 cm-1)
!
     data ((acr(k,j),k=1,8),j=  1,  8) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.12740e+05, 0.12750e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.12750e+05, 0.12760e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.12760e+05, 0.12770e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.12770e+05, 0.12780e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.12780e+05, 0.12790e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.12790e+05, 0.12800e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.12800e+05, 0.12810e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.12810e+05, 0.12820e+05/
     data ((acr(k,j),k=1,8),j=  9, 16) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.12820e+05, 0.12830e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.12830e+05, 0.12840e+05,                           &
     0.14615e-07, 0.36385e-01, 0.66900e-01,-0.26760e-03, 0.64873e-01, &
    -0.26352e-03, 0.12840e+05, 0.12850e+05,                           &
     0.17551e-07, 0.18194e-01, 0.63056e-01,-0.25222e-03, 0.61029e-01, &
    -0.24814e-03, 0.12850e+05, 0.12860e+05,                           &
     0.55571e-07, 0.34936e-01, 0.60516e-01,-0.24079e-03, 0.58832e-01, &
    -0.23863e-03, 0.12860e+05, 0.12870e+05,                           &
     0.11652e-06, 0.36823e-01, 0.57157e-01,-0.22751e-03, 0.55426e-01, &
    -0.22511e-03, 0.12870e+05, 0.12880e+05,                           &
     0.22758e-06, 0.37524e-01, 0.54060e-01,-0.21528e-03, 0.52277e-01, &
    -0.21261e-03, 0.12880e+05, 0.12890e+05,                           &
     0.41844e-06, 0.44301e-01, 0.51439e-01,-0.20300e-03, 0.50656e-01, &
    -0.20073e-03, 0.12890e+05, 0.12900e+05/
     data ((acr(k,j),k=1,8),j= 17, 24) /                              &
     0.11196e-05, 0.64417e-01, 0.48398e-01,-0.19212e-03, 0.47119e-01, &
    -0.18843e-03, 0.12900e+05, 0.12910e+05,                           &
     0.13443e-05, 0.56296e-01, 0.46024e-01,-0.17943e-03, 0.45942e-01, &
    -0.17773e-03, 0.12910e+05, 0.12920e+05,                           &
     0.19719e-05, 0.64603e-01, 0.44605e-01,-0.17094e-03, 0.44515e-01, &
    -0.17298e-03, 0.12920e+05, 0.12930e+05,                           &
     0.26027e-05, 0.60089e-01, 0.42943e-01,-0.16604e-03, 0.41770e-01, &
    -0.16646e-03, 0.12930e+05, 0.12940e+05,                           &
     0.57954e-05, 0.10762e+00, 0.40374e-01,-0.15770e-03, 0.38909e-01, &
    -0.15614e-03, 0.12940e+05, 0.12950e+05,                           &
     0.43014e-05, 0.84704e-01, 0.42188e-01,-0.15762e-03, 0.40357e-01, &
    -0.15977e-03, 0.12950e+05, 0.12960e+05,                           &
     0.79614e-05, 0.97284e-01, 0.41471e-01,-0.16264e-03, 0.38436e-01, &
    -0.15588e-03, 0.12960e+05, 0.12970e+05,                           &
     0.17998e-04, 0.71534e-01, 0.38440e-01,-0.15372e-03, 0.36186e-01, &
    -0.14871e-03, 0.12970e+05, 0.12980e+05/
     data ((acr(k,j),k=1,8),j= 25, 32) /                              &
     0.27733e-04, 0.65530e-01, 0.34900e-01,-0.13823e-03, 0.34326e-01, &
    -0.14018e-03, 0.12980e+05, 0.12990e+05,                           &
     0.87111e-04, 0.82142e-01, 0.31863e-01,-0.12227e-03, 0.33009e-01, &
    -0.13017e-03, 0.12990e+05, 0.13000e+05,                           &
     0.59912e-04, 0.64608e-01, 0.31064e-01,-0.11475e-03, 0.36751e-01, &
    -0.12672e-03, 0.13000e+05, 0.13010e+05,                           &
     0.26987e-03, 0.37736e-01, 0.26012e-01,-0.10390e-03, 0.24588e-01, &
    -0.97152e-04, 0.13010e+05, 0.13020e+05,                           &
     0.59914e-03, 0.37261e-01, 0.22365e-01,-0.89458e-04, 0.20337e-01, &
    -0.85375e-04, 0.13020e+05, 0.13030e+05,                           &
     0.12384e-02, 0.38141e-01, 0.19000e-01,-0.75998e-04, 0.16972e-01, &
    -0.71916e-04, 0.13030e+05, 0.13040e+05,                           &
     0.23802e-02, 0.39440e-01, 0.15907e-01,-0.63629e-04, 0.13880e-01, &
    -0.59546e-04, 0.13040e+05, 0.13050e+05,                           &
     0.78692e-02, 0.59793e-01, 0.11910e-01,-0.46882e-04, 0.10043e-01, &
    -0.43832e-04, 0.13050e+05, 0.13060e+05/
     data ((acr(k,j),k=1,8),j= 33, 40) /                              &
     0.14112e-01, 0.63200e-01, 0.88039e-02,-0.34775e-04, 0.68655e-02, &
    -0.31246e-04, 0.13060e+05, 0.13070e+05,                           &
     0.15060e-01, 0.44198e-01, 0.62748e-02,-0.25099e-04, 0.42475e-02, &
    -0.21017e-04, 0.13070e+05, 0.13080e+05,                           &
     0.19365e-01, 0.45501e-01, 0.45529e-02,-0.18211e-04, 0.25253e-02, &
    -0.14128e-04, 0.13080e+05, 0.13090e+05,                           &
     0.35332e-01, 0.70599e-01, 0.26881e-02,-0.10606e-04, 0.67201e-03, &
    -0.66426e-05, 0.13090e+05, 0.13100e+05,                           &
     0.32021e-01, 0.71785e-01, 0.13302e-02,-0.52395e-05,-0.70181e-03, &
    -0.11794e-05, 0.13100e+05, 0.13110e+05,                           &
     0.20424e-01, 0.76507e-01, 0.34001e-03,-0.13505e-05,-0.17047e-02, &
     0.27972e-05, 0.13110e+05, 0.13120e+05,                           &
     0.83211e-02, 0.50488e-01, 0.63467e-04,-0.25366e-06,-0.19594e-02, &
     0.38105e-05, 0.13120e+05, 0.13130e+05,                           &
     0.39472e-01, 0.95766e-01, 0.76617e-03,-0.30190e-05,-0.12832e-02, &
     0.11287e-05, 0.13130e+05, 0.13140e+05/
     data ((acr(k,j),k=1,8),j= 41, 48) /                              &
     0.57473e-01, 0.11569e+00, 0.28314e-02,-0.10907e-04, 0.85516e-03, &
    -0.71029e-05, 0.13140e+05, 0.13150e+05,                           &
     0.46463e-01, 0.16239e+00, 0.72656e-02,-0.26424e-04, 0.59616e-02, &
    -0.24138e-04, 0.13150e+05, 0.13160e+05,                           &
     0.70486e-02, 0.18365e+00, 0.16807e-01,-0.58784e-04, 0.17894e-01, &
    -0.64531e-04, 0.13160e+05, 0.13170e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13170e+05, 0.13180e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13180e+05, 0.13190e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13190e+05, 0.13200e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13200e+05, 0.13210e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13210e+05, 0.13220e+05/
     data ((acr(k,j),k=1,8),j= 49, 56) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13220e+05, 0.13230e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13230e+05, 0.13240e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13240e+05, 0.13250e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13250e+05, 0.13260e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13260e+05, 0.13270e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13270e+05, 0.13280e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13280e+05, 0.13290e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13290e+05, 0.13300e+05/
     data ((acr(k,j),k=1,8),j= 57, 64) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13300e+05, 0.13310e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13310e+05, 0.13320e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13320e+05, 0.13330e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13330e+05, 0.13340e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13340e+05, 0.13350e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13350e+05, 0.13360e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13360e+05, 0.13370e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13370e+05, 0.13380e+05/
     data ((acr(k,j),k=1,8),j= 65, 72) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13380e+05, 0.13390e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13390e+05, 0.13400e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13400e+05, 0.13410e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13410e+05, 0.13420e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13420e+05, 0.13430e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13430e+05, 0.13440e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13440e+05, 0.13450e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13450e+05, 0.13460e+05/
     data ((acr(k,j),k=1,8),j= 73, 80) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13460e+05, 0.13470e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13470e+05, 0.13480e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13480e+05, 0.13490e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13490e+05, 0.13500e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13500e+05, 0.13510e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13510e+05, 0.13520e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13520e+05, 0.13530e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13530e+05, 0.13540e+05/
     data ((acr(k,j),k=1,8),j= 81, 88) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13540e+05, 0.13550e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13550e+05, 0.13560e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13560e+05, 0.13570e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13570e+05, 0.13580e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13580e+05, 0.13590e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13590e+05, 0.13600e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13600e+05, 0.13610e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13610e+05, 0.13620e+05/
     data ((acr(k,j),k=1,8),j= 89, 96) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13620e+05, 0.13630e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13630e+05, 0.13640e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13640e+05, 0.13650e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13650e+05, 0.13660e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13660e+05, 0.13670e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13670e+05, 0.13680e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13680e+05, 0.13690e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13690e+05, 0.13700e+05/
     data ((acr(k,j),k=1,8),j= 97,104) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13700e+05, 0.13710e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13710e+05, 0.13720e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13720e+05, 0.13730e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13730e+05, 0.13740e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13740e+05, 0.13750e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13750e+05, 0.13760e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13760e+05, 0.13770e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13770e+05, 0.13780e+05/
     data ((acr(k,j),k=1,8),j=105,112) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13780e+05, 0.13790e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13790e+05, 0.13800e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13800e+05, 0.13810e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13810e+05, 0.13820e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13820e+05, 0.13830e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13830e+05, 0.13840e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13840e+05, 0.13850e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13850e+05, 0.13860e+05/
     data ((acr(k,j),k=1,8),j=113,120) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13860e+05, 0.13870e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13870e+05, 0.13880e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13880e+05, 0.13890e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13890e+05, 0.13900e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13900e+05, 0.13910e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13910e+05, 0.13920e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13920e+05, 0.13930e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13930e+05, 0.13940e+05/
     data ((acr(k,j),k=1,8),j=121,128) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13940e+05, 0.13950e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13950e+05, 0.13960e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13960e+05, 0.13970e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13970e+05, 0.13980e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13980e+05, 0.13990e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.13990e+05, 0.14000e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14000e+05, 0.14010e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14010e+05, 0.14020e+05/
     data ((acr(k,j),k=1,8),j=129,136) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14020e+05, 0.14030e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14030e+05, 0.14040e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14040e+05, 0.14050e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14050e+05, 0.14060e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14060e+05, 0.14070e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14070e+05, 0.14080e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14080e+05, 0.14090e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14090e+05, 0.14100e+05/
     data ((acr(k,j),k=1,8),j=137,144) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14100e+05, 0.14110e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14110e+05, 0.14120e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14120e+05, 0.14130e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14130e+05, 0.14140e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14140e+05, 0.14150e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14150e+05, 0.14160e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14160e+05, 0.14170e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14170e+05, 0.14180e+05/
     data ((acr(k,j),k=1,8),j=145,152) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14180e+05, 0.14190e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14190e+05, 0.14200e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14200e+05, 0.14210e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14210e+05, 0.14220e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14220e+05, 0.14230e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14230e+05, 0.14240e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14240e+05, 0.14250e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14250e+05, 0.14260e+05/
     data ((acr(k,j),k=1,8),j=153,160) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14260e+05, 0.14270e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14270e+05, 0.14280e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14280e+05, 0.14290e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14290e+05, 0.14300e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14300e+05, 0.14310e+05,                           &
     0.32848e-07, 0.36386e-01, 0.53505e-01,-0.21402e-03, 0.51477e-01, &
    -0.20994e-03, 0.14310e+05, 0.14320e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14320e+05, 0.14330e+05,                           &
     0.11021e-06, 0.36386e-01, 0.48248e-01,-0.19299e-03, 0.46221e-01, &
    -0.18891e-03, 0.14330e+05, 0.14340e+05/
     data ((acr(k,j),k=1,8),j=161,168) /                              &
     0.34571e-06, 0.36386e-01, 0.43260e-01,-0.17304e-03, 0.41233e-01, &
    -0.16896e-03, 0.14340e+05, 0.14350e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14350e+05, 0.14360e+05,                           &
     0.10131e-05, 0.36386e-01, 0.38540e-01,-0.15416e-03, 0.36513e-01, &
    -0.15008e-03, 0.14360e+05, 0.14370e+05,                           &
     0.27746e-05, 0.36442e-01, 0.34089e-01,-0.13635e-03, 0.32054e-01, &
    -0.13222e-03, 0.14370e+05, 0.14380e+05,                           &
     0.36059e-05, 0.18240e-01, 0.29930e-01,-0.11971e-03, 0.27892e-01, &
    -0.11557e-03, 0.14380e+05, 0.14390e+05,                           &
     0.34852e-05, 0.18279e-01, 0.29884e-01,-0.11951e-03, 0.27832e-01, &
    -0.11527e-03, 0.14390e+05, 0.14400e+05,                           &
     0.16898e-04, 0.36461e-01, 0.26000e-01,-0.10399e-03, 0.23962e-01, &
    -0.99845e-04, 0.14400e+05, 0.14410e+05,                           &
     0.37525e-04, 0.37306e-01, 0.22363e-01,-0.89448e-04, 0.20327e-01, &
    -0.85313e-04, 0.14410e+05, 0.14420e+05/
     data ((acr(k,j),k=1,8),j=169,176) /                              &
     0.77568e-04, 0.38179e-01, 0.18999e-01,-0.75991e-04, 0.16964e-01, &
    -0.71868e-04, 0.14420e+05, 0.14430e+05,                           &
     0.76440e-04, 0.19568e-01, 0.15927e-01,-0.63702e-04, 0.13885e-01, &
    -0.59534e-04, 0.14430e+05, 0.14440e+05,                           &
     0.20947e-03, 0.38957e-01, 0.14055e-01,-0.55384e-04, 0.12236e-01, &
    -0.52506e-04, 0.14440e+05, 0.14450e+05,                           &
     0.35608e-03, 0.40125e-01, 0.11458e-01,-0.45137e-04, 0.95851e-02, &
    -0.41989e-04, 0.14450e+05, 0.14460e+05,                           &
     0.56079e-03, 0.41622e-01, 0.91259e-02,-0.35944e-04, 0.72131e-02, &
    -0.32578e-04, 0.14460e+05, 0.14470e+05,                           &
     0.81523e-03, 0.43384e-01, 0.70608e-02,-0.27807e-04, 0.51286e-02, &
    -0.24309e-04, 0.14470e+05, 0.14480e+05,                           &
     0.16639e-02, 0.67344e-01, 0.50076e-02,-0.19758e-04, 0.30273e-02, &
    -0.15990e-04, 0.14480e+05, 0.14490e+05,                           &
     0.14146e-02, 0.46358e-01, 0.31060e-02,-0.12424e-04, 0.10767e-02, &
    -0.83329e-05, 0.14490e+05, 0.14500e+05/
     data ((acr(k,j),k=1,8),j=177,184) /                              &
     0.22165e-02, 0.71957e-01, 0.16398e-02,-0.64770e-05,-0.38835e-03, &
    -0.24326e-05, 0.14500e+05, 0.14510e+05,                           &
     0.15542e-02, 0.74780e-01, 0.64291e-03,-0.25306e-05,-0.13922e-02, &
     0.15886e-05, 0.14510e+05, 0.14520e+05,                           &
     0.31245e-03, 0.26878e-01, 0.95018e-04,-0.38003e-06,-0.19316e-02, &
     0.36996e-05, 0.14520e+05, 0.14530e+05,                           &
     0.15805e-02, 0.97644e-01, 0.30344e-03,-0.10773e-05,-0.13616e-02, &
     0.59540e-05, 0.14530e+05, 0.14540e+05,                           &
     0.44256e-02, 0.14544e+00, 0.20504e-02,-0.74913e-05, 0.62397e-03, &
    -0.95059e-06, 0.14540e+05, 0.14550e+05,                           &
     0.39335e-02, 0.29035e+00, 0.78502e-02,-0.23749e-04, 0.80984e-02, &
    -0.27655e-04, 0.14550e+05, 0.14560e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14560e+05, 0.14570e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14570e+05, 0.14580e+05/
     data ((acr(k,j),k=1,8),j=185,192) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14580e+05, 0.14590e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14590e+05, 0.14600e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14600e+05, 0.14610e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14610e+05, 0.14620e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14620e+05, 0.14630e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14630e+05, 0.14640e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14640e+05, 0.14650e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14650e+05, 0.14660e+05/
     data ((acr(k,j),k=1,8),j=193,200) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14660e+05, 0.14670e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14670e+05, 0.14680e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14680e+05, 0.14690e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14690e+05, 0.14700e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14700e+05, 0.14710e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14710e+05, 0.14720e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14720e+05, 0.14730e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14730e+05, 0.14740e+05/
     data ((acr(k,j),k=1,8),j=201,208) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14740e+05, 0.14750e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14750e+05, 0.14760e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14760e+05, 0.14770e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14770e+05, 0.14780e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14780e+05, 0.14790e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14790e+05, 0.14800e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14800e+05, 0.14810e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14810e+05, 0.14820e+05/
     data ((acr(k,j),k=1,8),j=209,216) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14820e+05, 0.14830e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14830e+05, 0.14840e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14840e+05, 0.14850e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14850e+05, 0.14860e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14860e+05, 0.14870e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14870e+05, 0.14880e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14880e+05, 0.14890e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14890e+05, 0.14900e+05/
     data ((acr(k,j),k=1,8),j=217,224) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14900e+05, 0.14910e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14910e+05, 0.14920e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14920e+05, 0.14930e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14930e+05, 0.14940e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14940e+05, 0.14950e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14950e+05, 0.14960e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14960e+05, 0.14970e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14970e+05, 0.14980e+05/
     data ((acr(k,j),k=1,8),j=225,232) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14980e+05, 0.14990e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.14990e+05, 0.15000e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15000e+05, 0.15010e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15010e+05, 0.15020e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15020e+05, 0.15030e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15030e+05, 0.15040e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15040e+05, 0.15050e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15050e+05, 0.15060e+05/
     data ((acr(k,j),k=1,8),j=233,240) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15060e+05, 0.15070e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15070e+05, 0.15080e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15080e+05, 0.15090e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15090e+05, 0.15100e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15100e+05, 0.15110e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15110e+05, 0.15120e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15120e+05, 0.15130e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15130e+05, 0.15140e+05/
     data ((acr(k,j),k=1,8),j=241,248) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15140e+05, 0.15150e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15150e+05, 0.15160e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15160e+05, 0.15170e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15170e+05, 0.15180e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15180e+05, 0.15190e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15190e+05, 0.15200e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15200e+05, 0.15210e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15210e+05, 0.15220e+05/
     data ((acr(k,j),k=1,8),j=249,256) /                              &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15220e+05, 0.15230e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15230e+05, 0.15240e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15240e+05, 0.15250e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15250e+05, 0.15260e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15260e+05, 0.15270e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15270e+05, 0.15280e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15280e+05, 0.15290e+05,                           &
     0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
     0.00000e+00, 0.15290e+05, 0.15300e+05/
!
    do i=1,8
        a(i)=acr(i,inu)
    enddo
!
    return
end