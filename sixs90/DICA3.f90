subroutine dica3(a,inu)
    implicit none
    real(8) :: a(8)
    real(8) :: acr(8,256)
    integer :: inu,j,k,i
!   carbon dioxide (7620 - 10170 cm-1)
!
      data ((acr(k,j),k=1,8),j=  1,  8) /                              &
      0.41135e-04, 0.13491e+00, 0.19511e-01,-0.88592e-04, 0.17169e-01, &
     -0.86383e-04, 0.76200e+04, 0.76300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.76300e+04, 0.76400e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.76400e+04, 0.76500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.76500e+04, 0.76600e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.76600e+04, 0.76700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.76700e+04, 0.76800e+04,                           &
      0.69843e-05, 0.58690e-01, 0.17996e-01,-0.84950e-04, 0.14986e-01, &
     -0.79255e-04, 0.76800e+04, 0.76900e+04,                           &
      0.44856e-04, 0.11610e+00, 0.12313e-01,-0.61208e-04, 0.94917e-02, &
     -0.56817e-04, 0.76900e+04, 0.77000e+04/
      data ((acr(k,j),k=1,8),j=  9, 16) /                              &
      0.21119e-03, 0.14823e+00, 0.58288e-02,-0.35255e-04, 0.29752e-02, &
     -0.30772e-04, 0.77000e+04, 0.77100e+04,                           &
      0.68368e-03, 0.18822e+00, 0.21812e-03,-0.13229e-04,-0.27425e-02, &
     -0.80771e-05, 0.77100e+04, 0.77200e+04,                           &
      0.80401e-03, 0.20648e+00,-0.32887e-02, 0.50708e-07,-0.62117e-02, &
      0.59400e-05, 0.77200e+04, 0.77300e+04,                           &
      0.36897e-03, 0.20612e+00,-0.45166e-02, 0.47173e-05,-0.74494e-02, &
      0.10697e-04, 0.77300e+04, 0.77400e+04,                           &
      0.11094e-02, 0.31021e+00,-0.22536e-02,-0.10224e-05,-0.33444e-02, &
      0.21129e-05, 0.77400e+04, 0.77500e+04,                           &
      0.65848e-03, 0.26193e+00, 0.27594e-02,-0.21278e-04, 0.25217e-03, &
     -0.18776e-04, 0.77500e+04, 0.77600e+04,                           &
      0.73155e-04, 0.30739e+00, 0.13041e-01,-0.63190e-04, 0.10499e-01, &
     -0.60136e-04, 0.77600e+04, 0.77700e+04,                           &
      0.19363e-04, 0.19417e+00, 0.14647e-01,-0.71772e-04, 0.11659e-01, &
     -0.65892e-04, 0.77700e+04, 0.77800e+04/
      data ((acr(k,j),k=1,8),j= 17, 24) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.77800e+04, 0.77900e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.77900e+04, 0.78000e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.78000e+04, 0.78100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.78100e+04, 0.78200e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.78200e+04, 0.78300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.78300e+04, 0.78400e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.78400e+04, 0.78500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.78500e+04, 0.78600e+04/
      data ((acr(k,j),k=1,8),j= 25, 32) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.78600e+04, 0.78700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.78700e+04, 0.78800e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.78800e+04, 0.78900e+04,                           &
      0.37190e-04, 0.18910e+00, 0.32484e-03,-0.13576e-04,-0.26862e-02, &
     -0.82867e-05, 0.78900e+04, 0.79000e+04,                           &
      0.94656e-04, 0.36300e+00,-0.20395e-02,-0.42752e-05,-0.49939e-02, &
      0.12049e-05, 0.79000e+04, 0.79100e+04,                           &
      0.73888e-04, 0.33612e+00,-0.27854e-02,-0.14358e-05,-0.56280e-02, &
      0.37427e-05, 0.79100e+04, 0.79200e+04,                           &
      0.35986e-04, 0.14439e+00,-0.40916e-02, 0.30556e-05,-0.70423e-02, &
      0.90581e-05, 0.79200e+04, 0.79300e+04,                           &
      0.77290e-04, 0.22754e+00,-0.15562e-02,-0.62734e-05,-0.44592e-02, &
     -0.10627e-05, 0.79300e+04, 0.79400e+04/
      data ((acr(k,j),k=1,8),j= 33, 40) /                              &
      0.18388e-04, 0.92491e-01, 0.27526e-02,-0.24014e-04,-0.34204e-03, &
     -0.18092e-04, 0.79400e+04, 0.79500e+04,                           &
      0.19936e-06, 0.10968e-02, 0.17030e-02,-0.20156e-04,-0.13872e-02, &
     -0.14222e-04, 0.79500e+04, 0.79600e+04,                           &
      0.49455e-06, 0.19615e-02,-0.15846e-02,-0.68674e-05,-0.42920e-02, &
     -0.17132e-05, 0.79600e+04, 0.79700e+04,                           &
      0.27828e-06, 0.13177e-02,-0.36989e-02, 0.11926e-05,-0.65515e-02, &
      0.69563e-05, 0.79700e+04, 0.79800e+04,                           &
      0.38372e-06, 0.17475e-02,-0.38411e-02, 0.17875e-05,-0.67043e-02, &
      0.75942e-05, 0.79800e+04, 0.79900e+04,                           &
      0.73276e-06, 0.30110e-02,-0.67794e-03,-0.94857e-05,-0.34111e-02, &
     -0.51641e-05, 0.79900e+04, 0.80000e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.80000e+04, 0.80100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.80100e+04, 0.80200e+04/
      data ((acr(k,j),k=1,8),j= 41, 48) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.80200e+04, 0.80300e+04,                           &
      0.37029e-07, 0.36096e-03, 0.12109e-01,-0.62089e-04, 0.89059e-02, &
     -0.55639e-04, 0.80300e+04, 0.80400e+04,                           &
      0.21335e-06, 0.10751e-02, 0.88149e-02,-0.48284e-04, 0.57133e-02, &
     -0.42623e-04, 0.80400e+04, 0.80500e+04,                           &
      0.46462e-06, 0.10860e-02, 0.49284e-02,-0.32910e-04, 0.17555e-02, &
     -0.26942e-04, 0.80500e+04, 0.80600e+04,                           &
      0.15991e-05, 0.18385e-02, 0.64633e-03,-0.15477e-04,-0.22392e-02, &
     -0.10372e-04, 0.80600e+04, 0.80700e+04,                           &
      0.17752e-05, 0.16397e-02,-0.24966e-02,-0.34666e-05,-0.51909e-02, &
      0.18847e-05, 0.80700e+04, 0.80800e+04,                           &
      0.24423e-04, 0.26254e+00, 0.25790e-01,-0.77107e-04, 0.25403e-01, &
     -0.11405e-03, 0.80800e+04, 0.80900e+04,                           &
      0.61857e-04, 0.23309e+00, 0.21977e-01,-0.77807e-04, 0.20482e-01, &
     -0.95402e-04, 0.80900e+04, 0.81000e+04/
      data ((acr(k,j),k=1,8),j= 49, 56) /                              &
      0.13473e-03, 0.23293e+00, 0.19374e-01,-0.84779e-04, 0.16814e-01, &
     -0.82971e-04, 0.81000e+04, 0.81100e+04,                           &
      0.33293e-03, 0.42794e+00, 0.16857e-01,-0.74779e-04, 0.15630e-01, &
     -0.74735e-04, 0.81100e+04, 0.81200e+04,                           &
      0.45545e-03, 0.42084e+00, 0.14131e-01,-0.64833e-04, 0.11587e-01, &
     -0.62246e-04, 0.81200e+04, 0.81300e+04,                           &
      0.39267e-03, 0.44590e+00, 0.14220e-01,-0.66328e-04, 0.10179e-01, &
     -0.58346e-04, 0.81300e+04, 0.81400e+04,                           &
      0.14095e-02, 0.66221e+00, 0.12879e-01,-0.63481e-04, 0.10259e-01, &
     -0.59679e-04, 0.81400e+04, 0.81500e+04,                           &
      0.25744e-02, 0.34109e+00, 0.81434e-02,-0.35365e-04, 0.10242e-01, &
     -0.44763e-04, 0.81500e+04, 0.81600e+04,                           &
      0.53482e-02, 0.12345e+00, 0.21618e-02,-0.21398e-04,-0.89174e-03, &
     -0.15761e-04, 0.81600e+04, 0.81700e+04,                           &
      0.85974e-02, 0.12902e+00,-0.11939e-02,-0.82690e-05,-0.41809e-02, &
     -0.25012e-05, 0.81700e+04, 0.81800e+04/
      data ((acr(k,j),k=1,8),j= 57, 64) /                              &
      0.11093e-01, 0.20887e+00,-0.36354e-02, 0.13570e-05,-0.65940e-02, &
      0.74093e-05, 0.81800e+04, 0.81900e+04,                           &
      0.90124e-02, 0.22951e+00,-0.41975e-02, 0.34952e-05,-0.71581e-02, &
      0.96182e-05, 0.81900e+04, 0.82000e+04,                           &
      0.22977e-01, 0.54634e+00, 0.54679e-03,-0.74452e-05, 0.32715e-04, &
     -0.77515e-05, 0.82000e+04, 0.82100e+04,                           &
      0.40967e-04, 0.33171e+00, 0.29828e-01,-0.12884e-03, 0.27502e-01, &
     -0.12751e-03, 0.82100e+04, 0.82200e+04,                           &
      0.13117e-03, 0.36661e+00, 0.24705e-01,-0.10516e-03, 0.23688e-01, &
     -0.11029e-03, 0.82200e+04, 0.82300e+04,                           &
      0.41427e-03, 0.33236e+00, 0.18964e-01,-0.81986e-04, 0.18612e-01, &
     -0.88547e-04, 0.82300e+04, 0.82400e+04,                           &
      0.11268e-02, 0.30228e+00, 0.13786e-01,-0.60455e-04, 0.14251e-01, &
     -0.66968e-04, 0.82400e+04, 0.82500e+04,                           &
      0.38631e-02, 0.31150e+00, 0.80823e-02,-0.38679e-04, 0.84273e-02, &
     -0.42578e-04, 0.82500e+04, 0.82600e+04/
      data ((acr(k,j),k=1,8),j= 65, 72) /                              &
      0.60039e-02, 0.26991e+00, 0.40499e-02,-0.23093e-04, 0.51036e-02, &
     -0.27510e-04, 0.82600e+04, 0.82700e+04,                           &
      0.14968e-01, 0.24443e+00,-0.60624e-03,-0.89134e-05,-0.14888e-02, &
     -0.75779e-05, 0.82700e+04, 0.82800e+04,                           &
      0.15831e-01, 0.42226e+00,-0.20274e-02, 0.62203e-05,-0.16360e-03, &
     -0.21524e-05, 0.82800e+04, 0.82900e+04,                           &
      0.86272e-02, 0.42161e+00,-0.28717e-02, 0.16391e-04, 0.58996e-03, &
      0.15552e-04, 0.82900e+04, 0.83000e+04,                           &
      0.38978e-01, 0.60506e+00,-0.18963e-03,-0.46658e-05,-0.11243e-02, &
     -0.48784e-05, 0.83000e+04, 0.83100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.83100e+04, 0.83200e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.83200e+04, 0.83300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.83300e+04, 0.83400e+04/
      data ((acr(k,j),k=1,8),j= 73, 80) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.83400e+04, 0.83500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.83500e+04, 0.83600e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.83600e+04, 0.83700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.83700e+04, 0.83800e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.83800e+04, 0.83900e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.83900e+04, 0.84000e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.84000e+04, 0.84100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.84100e+04, 0.84200e+04/
      data ((acr(k,j),k=1,8),j= 81, 88) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.84200e+04, 0.84300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.84300e+04, 0.84400e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.84400e+04, 0.84500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.84500e+04, 0.84600e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.84600e+04, 0.84700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.84700e+04, 0.84800e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.84800e+04, 0.84900e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.84900e+04, 0.85000e+04/
      data ((acr(k,j),k=1,8),j= 89, 96) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.85000e+04, 0.85100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.85100e+04, 0.85200e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.85200e+04, 0.85300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.85300e+04, 0.85400e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.85400e+04, 0.85500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.85500e+04, 0.85600e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.85600e+04, 0.85700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.85700e+04, 0.85800e+04/
      data ((acr(k,j),k=1,8),j= 97,104) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.85800e+04, 0.85900e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.85900e+04, 0.86000e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.86000e+04, 0.86100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.86100e+04, 0.86200e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.86200e+04, 0.86300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.86300e+04, 0.86400e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.86400e+04, 0.86500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.86500e+04, 0.86600e+04/
      data ((acr(k,j),k=1,8),j=105,112) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.86600e+04, 0.86700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.86700e+04, 0.86800e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.86800e+04, 0.86900e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.86900e+04, 0.87000e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.87000e+04, 0.87100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.87100e+04, 0.87200e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.87200e+04, 0.87300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.87300e+04, 0.87400e+04/
      data ((acr(k,j),k=1,8),j=113,120) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.87400e+04, 0.87500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.87500e+04, 0.87600e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.87600e+04, 0.87700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.87700e+04, 0.87800e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.87800e+04, 0.87900e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.87900e+04, 0.88000e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.88000e+04, 0.88100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.88100e+04, 0.88200e+04/
      data ((acr(k,j),k=1,8),j=121,128) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.88200e+04, 0.88300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.88300e+04, 0.88400e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.88400e+04, 0.88500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.88500e+04, 0.88600e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.88600e+04, 0.88700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.88700e+04, 0.88800e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.88800e+04, 0.88900e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.88900e+04, 0.89000e+04/
      data ((acr(k,j),k=1,8),j=129,136) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.89000e+04, 0.89100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.89100e+04, 0.89200e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.89200e+04, 0.89300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.89300e+04, 0.89400e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.89400e+04, 0.89500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.89500e+04, 0.89600e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.89600e+04, 0.89700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.89700e+04, 0.89800e+04/
      data ((acr(k,j),k=1,8),j=137,144) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.89800e+04, 0.89900e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.89900e+04, 0.90000e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.90000e+04, 0.90100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.90100e+04, 0.90200e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.90200e+04, 0.90300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.90300e+04, 0.90400e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.90400e+04, 0.90500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.90500e+04, 0.90600e+04/
      data ((acr(k,j),k=1,8),j=145,152) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.90600e+04, 0.90700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.90700e+04, 0.90800e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.90800e+04, 0.90900e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.90900e+04, 0.91000e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.91000e+04, 0.91100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.91100e+04, 0.91200e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.91200e+04, 0.91300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.91300e+04, 0.91400e+04/
      data ((acr(k,j),k=1,8),j=153,160) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.91400e+04, 0.91500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.91500e+04, 0.91600e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.91600e+04, 0.91700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.91700e+04, 0.91800e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.91800e+04, 0.91900e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.91900e+04, 0.92000e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.92000e+04, 0.92100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.92100e+04, 0.92200e+04/
      data ((acr(k,j),k=1,8),j=161,168) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.92200e+04, 0.92300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.92300e+04, 0.92400e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.92400e+04, 0.92500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.92500e+04, 0.92600e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.92600e+04, 0.92700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.92700e+04, 0.92800e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.92800e+04, 0.92900e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.92900e+04, 0.93000e+04/
      data ((acr(k,j),k=1,8),j=169,176) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.93000e+04, 0.93100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.93100e+04, 0.93200e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.93200e+04, 0.93300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.93300e+04, 0.93400e+04,                           &
      0.99593e-05, 0.60237e-01, 0.97616e-02,-0.52148e-04, 0.66534e-02, &
     -0.46124e-04, 0.93400e+04, 0.93500e+04,                           &
      0.43567e-04, 0.12051e+00, 0.54505e-02,-0.34277e-04, 0.24633e-02, &
     -0.29032e-04, 0.93500e+04, 0.93600e+04,                           &
      0.88924e-04, 0.12477e+00, 0.12273e-02,-0.17740e-04,-0.18429e-02, &
     -0.11984e-04, 0.93600e+04, 0.93700e+04,                           &
      0.15573e-03, 0.16489e+00,-0.21342e-02,-0.44631e-05,-0.50460e-02, &
      0.11731e-05, 0.93700e+04, 0.93800e+04/
      data ((acr(k,j),k=1,8),j=177,184) /                              &
      0.94382e-04, 0.17579e+00,-0.41243e-02, 0.31944e-05,-0.71235e-02, &
      0.93601e-05, 0.93800e+04, 0.93900e+04,                           &
      0.21829e-03, 0.27491e+00,-0.34678e-02, 0.82722e-06,-0.64858e-02, &
      0.70242e-05, 0.93900e+04, 0.94000e+04,                           &
      0.22700e-03, 0.36616e+00, 0.19245e-02,-0.15790e-04,-0.15673e-03, &
     -0.15436e-04, 0.94000e+04, 0.94100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.94100e+04, 0.94200e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.94200e+04, 0.94300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.94300e+04, 0.94400e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.94400e+04, 0.94500e+04,                           &
      0.40618e-04, 0.33672e+00, 0.16922e-01,-0.79759e-04, 0.13854e-01, &
     -0.74271e-04, 0.94500e+04, 0.94600e+04/
      data ((acr(k,j),k=1,8),j=185,192) /                              &
      0.80260e-04, 0.35567e+00, 0.13366e-01,-0.66285e-04, 0.10448e-01, &
     -0.60815e-04, 0.94600e+04, 0.94700e+04,                           &
      0.17445e-03, 0.15326e+00, 0.81302e-02,-0.44538e-04, 0.56112e-02, &
     -0.41022e-04, 0.94700e+04, 0.94800e+04,                           &
      0.33041e-03, 0.26053e+00, 0.49592e-02,-0.28310e-04, 0.49032e-02, &
     -0.33415e-04, 0.94800e+04, 0.94900e+04,                           &
      0.89723e-03, 0.15820e+00,-0.11481e-03,-0.12250e-04,-0.31268e-02, &
     -0.67086e-05, 0.94900e+04, 0.95000e+04,                           &
      0.87248e-03, 0.13522e+00,-0.29599e-02,-0.13666e-05,-0.58130e-02, &
      0.43278e-05, 0.95000e+04, 0.95100e+04,                           &
      0.57391e-03, 0.20153e+00,-0.43771e-02, 0.41806e-05,-0.73125e-02, &
      0.10212e-04, 0.95100e+04, 0.95200e+04,                           &
      0.21060e-02, 0.33852e+00,-0.26269e-02,-0.19082e-05,-0.56393e-02, &
      0.38240e-05, 0.95200e+04, 0.95300e+04,                           &
      0.87766e-03, 0.35363e+00, 0.44041e-02,-0.24109e-04, 0.31707e-02, &
     -0.26057e-04, 0.95300e+04, 0.95400e+04/
      data ((acr(k,j),k=1,8),j=193,200) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.95400e+04, 0.95500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.95500e+04, 0.95600e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.95600e+04, 0.95700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.95700e+04, 0.95800e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.95800e+04, 0.95900e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.95900e+04, 0.96000e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.96000e+04, 0.96100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.96100e+04, 0.96200e+04/
      data ((acr(k,j),k=1,8),j=201,208) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.96200e+04, 0.96300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.96300e+04, 0.96400e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.96400e+04, 0.96500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.96500e+04, 0.96600e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.96600e+04, 0.96700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.96700e+04, 0.96800e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.96800e+04, 0.96900e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.96900e+04, 0.97000e+04/
      data ((acr(k,j),k=1,8),j=209,216) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.97000e+04, 0.97100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.97100e+04, 0.97200e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.97200e+04, 0.97300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.97300e+04, 0.97400e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.97400e+04, 0.97500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.97500e+04, 0.97600e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.97600e+04, 0.97700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.97700e+04, 0.97800e+04/
      data ((acr(k,j),k=1,8),j=217,224) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.97800e+04, 0.97900e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.97900e+04, 0.98000e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.98000e+04, 0.98100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.98100e+04, 0.98200e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.98200e+04, 0.98300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.98300e+04, 0.98400e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.98400e+04, 0.98500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.98500e+04, 0.98600e+04/
      data ((acr(k,j),k=1,8),j=225,232) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.98600e+04, 0.98700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.98700e+04, 0.98800e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.98800e+04, 0.98900e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.98900e+04, 0.99000e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.99000e+04, 0.99100e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.99100e+04, 0.99200e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.99200e+04, 0.99300e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.99300e+04, 0.99400e+04/
      data ((acr(k,j),k=1,8),j=233,240) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.99400e+04, 0.99500e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.99500e+04, 0.99600e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.99600e+04, 0.99700e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.99700e+04, 0.99800e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.99800e+04, 0.99900e+04,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.99900e+04, 0.10000e+05,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10000e+05, 0.10010e+05,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10010e+05, 0.10020e+05/
      data ((acr(k,j),k=1,8),j=241,248) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10020e+05, 0.10030e+05,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10030e+05, 0.10040e+05,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10040e+05, 0.10050e+05,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10050e+05, 0.10060e+05,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10060e+05, 0.10070e+05,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10070e+05, 0.10080e+05,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10080e+05, 0.10090e+05,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10090e+05, 0.10100e+05/
      data ((acr(k,j),k=1,8),j=249,256) /                              &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10100e+05, 0.10110e+05,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10110e+05, 0.10120e+05,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10120e+05, 0.10130e+05,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10130e+05, 0.10140e+05,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10140e+05, 0.10150e+05,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10150e+05, 0.10160e+05,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10160e+05, 0.10170e+05,                           &
      0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, 0.00000e+00, &
      0.00000e+00, 0.10170e+05, 0.10180e+05/
!
    do i=1,8
        a(i)=acr(i,inu)
    enddo
!
    return
end
