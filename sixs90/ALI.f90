subroutine ali(iwa)
    implicit none
    common /sixs_ffu/ s(1501),wlinf,wlsup
    real(8) :: sr(9,1501),wli(9),wls(9)
    real(8) :: wlinf,wlsup,s
    integer :: iwa,l,i

    data (sr(1,l),l=1,1501)/ 60*0.,                             &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0001,0.0011,0.0084,0.1313,0.7159,0.8481,0.9010,   &
     0.9452,0.9771,0.9889,0.9906,0.8846,0.1411,0.0208,0.0038,   &
     0.0006,0.0003,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,                                                    &
     1400*0./

    data (sr(2,l),l=1,1501) /  65*0.,                           &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0001,0.0004,0.0014,0.0044,0.0098,0.0226,0.0763,0.3166,   &
     0.5795,0.6084,0.6413,0.6608,0.6771,0.7032,0.7233,0.7249,   &
     0.7253,0.7512,0.7960,0.8345,0.8558,0.8720,0.8795,0.8775,   &
     0.8799,0.8978,0.9296,0.9687,0.9861,0.9854,0.9962,0.9680,   &
     0.6394,0.2610,0.1008,0.0442,0.0219,0.0113,0.0058,0.0029,   &
     0.0014,0.0008,0.0004,0.0002,0.0001,0.0001,0.0001,0.0001,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,                                             &
     1370*0./

    data (sr(3,l),l=1,1501)/ 95*0.,                             &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0001,0.0001,0.0002,   &
     0.0005,0.0011,0.0021,0.0038,0.0081,0.0217,0.0658,0.1825,   &
     0.3777,0.6131,0.7679,0.8196,0.8427,0.8554,0.8567,0.8502,   &
     0.8490,0.8573,0.8732,0.8969,0.9042,0.9062,0.9089,0.9126,   &
     0.9170,0.9170,0.9187,0.9329,0.9514,0.9537,0.9664,0.9717,   &
     0.9651,0.9676,0.9866,0.9995,0.9876,0.8677,0.5832,0.3143,   &
     0.1601,0.0824,0.0413,0.0182,0.0075,0.0029,0.0012,0.0006,   &
     0.0003,0.0002,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,                                             &
     1340*0./

    data (sr(4,l),l=1,1501)/ 125*0.,                            &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0001,0.0001,0.0003,0.0005,   &
     0.0009,0.0012,0.0015,0.0022,0.0036,0.0070,0.0150,0.0334,   &
     0.0692,0.1318,0.2436,0.4438,0.6893,0.8373,0.8585,0.8435,   &
     0.8591,0.8960,0.9181,0.9271,0.9341,0.9209,0.8920,0.8939,   &
     0.9305,0.9463,0.9286,0.9264,0.9547,0.9759,0.9617,0.9657,   &
     0.9958,0.9733,0.7396,0.3597,0.1325,0.0483,0.0195,0.0089,   &
     0.0044,0.0025,0.0014,0.0008,0.0006,0.0004,0.0002,0.0002,   &
     0.0001,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,                               &
     1300*0./

    data (sr(5,l),l=1,1501)/ 195*0.,                            &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0002,0.0003,   &
     0.0006,0.0015,0.0040,0.0114,0.0300,0.0717,0.1817,0.5186,   &
     0.9563,0.9936,0.9853,0.9929,0.9897,0.9949,0.9829,0.9638,   &
     0.9497,0.9633,0.9018,0.6182,0.2122,0.0549,0.0186,0.0086,   &
     0.0053,0.0035,0.0021,0.0012,0.0006,0.0005,0.0003,0.0001,   &
     0.0001,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,                                      &
     1255*0./

    data (sr(6,l),l=1,1501)/ 215*0.,                            &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0001,0.0001,0.0001,0.0002,0.0003,0.0005,   &
     0.0009,0.0019,0.0041,0.0103,0.0298,0.0997,0.3584,0.8497,   &
     1.0001,0.9855,0.9949,0.9972,0.9853,0.9723,0.9730,0.9684,   &
     0.9588,0.9316,0.9186,0.9158,0.9197,0.8857,0.8538,0.7765,   &
     0.6116,0.3680,0.1625,0.0610,0.0237,0.0101,0.0049,0.0026,   &
     0.0015,0.0008,0.0005,0.0005,0.0004,0.0003,0.0002,0.0002,   &
     0.0001,0.0003,0.0002,0.0005,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,                                             &
     1220*0./

    data (sr(7,l),l=1,1501)/ 345*0.,                            &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0001,   &
     0.0001,0.0001,0.0002,0.0003,0.0002,0.0002,0.0002,0.0003,   &
     0.0003,0.0004,0.0004,0.0008,0.0009,0.0015,0.0022,0.0034,   &
     0.0058,0.0092,0.0140,0.0192,0.0267,0.0355,0.0498,0.0723,   &
     0.1116,0.1858,0.3199,0.5328,0.7562,0.9017,0.9402,0.9297,   &
     0.9222,0.9260,0.9362,0.9446,0.9478,0.9472,0.9457,0.9467,   &
     0.9516,0.9584,0.9682,0.9755,0.9810,0.9824,0.9807,0.9787,   &
     0.9726,0.9715,0.9696,0.9765,0.9842,0.9951,0.9996,0.9910,   &
     0.9682,0.9357,0.9049,0.8833,0.8706,0.8215,0.6849,0.4661,   &
     0.2654,0.1402,0.0764,0.0440,0.0265,0.0169,0.0113,0.0082,   &
     0.0059,0.0049,0.0039,0.0031,0.0026,0.0022,0.0018,0.0014,   &
     0.0011,0.0007,0.0005,0.0003,0.0002,0.0002,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,                        &
     1055*0./

    data (sr(8,l),l=1,1501)/ 480*0.,                            &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0001,0.0001,0.0001,0.0001,0.0001,0.0002,0.0002,0.0003,   &
     0.0004,0.0005,0.0006,0.0007,0.0009,0.0011,0.0013,0.0020,   &
     0.0026,0.0037,0.0057,0.0086,0.0135,0.0214,0.0334,0.0518,   &
     0.0782,0.1146,0.1613,0.2193,0.2934,0.3864,0.4983,0.6170,   &
     0.7217,0.7857,0.8018,0.7861,0.7596,0.7366,0.7249,0.7270,   &
     0.7409,0.7616,0.7886,0.8175,0.8443,0.8656,0.8828,0.8966,   &
     0.9056,0.9125,0.9196,0.9252,0.9278,0.9311,0.9341,0.9360,   &
     0.9368,0.9382,0.9396,0.9385,0.9388,0.9372,0.9338,0.9264,   &
     0.9186,0.9089,0.8990,0.8885,0.8803,0.8738,0.8699,0.8706,   &
     0.8741,0.8822,0.8928,0.9058,0.9197,0.9353,0.9493,0.9631,   &
     0.9740,0.9831,0.9906,0.9959,0.9996,0.9996,0.9997,0.9969,   &
     0.9945,0.9895,0.9838,0.9744,0.9619,0.9451,0.9220,0.8901,   &
     0.8512,0.8073,0.7585,0.7115,0.6652,0.6221,0.5803,0.5416,   &
     0.5010,0.4546,0.4014,0.3432,0.2823,0.2229,0.1714,0.1290,   &
     0.0972,0.0733,0.0557,0.0427,0.0332,0.0263,0.0209,0.0168,   &
     0.0135,0.0109,0.0087,0.0070,0.0056,0.0043,0.0034,0.0027,   &
     0.0021,0.0017,0.0014,0.0011,0.0008,0.0007,0.0007,0.0004,   &
     0.0003,0.0003,0.0001,0.0003,0.0002,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,                                                    &
     860*0./

    data (sr(9,l),l=1,1501)/ 685*0.,                            &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0001,   &
     0.0002,0.0003,0.0003,0.0002,0.0004,0.0004,0.0004,0.0006,   &
     0.0006,0.0007,0.0008,0.0009,0.0010,0.0012,0.0016,0.0018,   &
     0.0023,0.0027,0.0033,0.0039,0.0047,0.0058,0.0069,0.0084,   &
     0.0104,0.0126,0.0158,0.0197,0.0244,0.0307,0.0386,0.0493,   &
     0.0625,0.0793,0.1010,0.1283,0.1642,0.2080,0.2619,0.3250,   &
     0.3987,0.4798,0.5613,0.6402,0.7113,0.7732,0.8233,0.8595,   &
     0.8875,0.9056,0.9190,0.9281,0.9355,0.9401,0.9440,0.9459,   &
     0.9476,0.9500,0.9503,0.9509,0.9513,0.9487,0.9486,0.9485,   &
     0.9469,0.9481,0.9464,0.9445,0.9443,0.9444,0.9397,0.9380,   &
     0.9410,0.9445,0.9410,0.9407,0.9432,0.9454,0.9459,0.9466,   &
     0.9506,0.9525,0.9571,0.9572,0.9580,0.9594,0.9627,0.9642,   &
     0.9669,0.9677,0.9685,0.9710,0.9724,0.9726,0.9760,0.9763,   &
     0.9768,0.9779,0.9773,0.9775,0.9796,0.9790,0.9802,0.9816,   &
     0.9806,0.9813,0.9820,0.9820,0.9840,0.9838,0.9831,0.9845,   &
     0.9846,0.9847,0.9874,0.9885,0.9903,0.9939,0.9941,0.9965,   &
     0.9974,0.9977,0.9989,0.9996,0.9982,0.9992,0.9979,0.9963,   &
     0.9947,0.9917,0.9887,0.9871,0.9839,0.9805,0.9787,0.9764,   &
     0.9744,0.9733,0.9682,0.9627,0.9567,0.9485,0.9397,0.9281,   &
     0.9115,0.8898,0.8607,0.8228,0.7726,0.7155,0.6492,0.5808,   &
     0.5112,0.4415,0.3763,0.3172,0.2647,0.2208,0.1824,0.1500,   &
     0.1233,0.1013,0.0832,0.0687,0.0567,0.0469,0.0390,0.0324,   &
     0.0268,0.0223,0.0185,0.0155,0.0130,0.0109,0.0093,0.0078,   &
     0.0067,0.0056,0.0048,0.0041,0.0035,0.0030,0.0026,0.0022,   &
     0.0019,0.0016,0.0014,0.0012,0.0011,0.0009,0.0008,0.0007,   &
     0.0006,0.0005,0.0005,0.0004,0.0004,0.0003,0.0003,0.0003,   &
     0.0002,0.0002,0.0002,0.0002,0.0001,0.0001,0.0001,0.0001,   &
     0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,   &
     0.0001,0.0001,0.0001,0.0001,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,                 &
     570*0./

    wli(1)=0.4225
    wls(1)=0.4625
    wli(2)=0.4325
    wls(2)=0.550
    wli(3)=0.500
    wls(3)=0.630
    wli(4)=0.5775
    wls(4)=0.730
    wli(5)=0.7525
    wls(5)=0.8375
    wli(6)=0.8025
    wls(6)=0.935
    wli(7)=1.130
    wls(7)=1.345
    wli(8)=1.470
    wls(8)=1.820
    wli(9)=1.980
    wls(9)=2.530
    do i=1,1501
        s(i)=sr(iwa,i)
    enddo
    wlinf=wli(iwa)
    wlsup=wls(iwa)
    return
end