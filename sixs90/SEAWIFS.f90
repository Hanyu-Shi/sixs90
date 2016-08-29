subroutine seawifs(iwa)
    implicit none
    common /sixs_ffu/ s(1501),wlinf,wlsup
    real(8) :: sr(8,1501),wli(8),wls(8)
    real(8) :: wlinf,wlsup,s
    integer :: iwa,l,i

    data (sr(1,l),l=1,1501) /  53*0.,                           &
     0.0000,0.0000,0.0009,0.0019,0.0038,0.0105,0.0301,0.1429,   &
     0.4118,0.6546,0.8044,0.8504,0.9117,0.9738,0.9964,0.9076,   &
     0.6166,0.3124,0.1076,0.0321,0.0128,0.0076,0.0048,0.0027,   &
     0.0014,0.0013,0.0010,0.0006,0.0005,0.0004,0.0002,0.0001,   &
     0.0001,0.0001,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0001,0.0004,0.0016,0.0035,0.0020,   &
     0.0002,0.0003,0.0013,0.0018,0.0016,0.0005,0.0004,0.0010,   &
     0.0015,0.0007,0.0002,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0001,   &
     0.0001,0.0002,0.0001,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     1320*0./

    data (sr(2,l),l=1,1501)/ 52*0.,                             &
     0.0020,0.0030,0.0011,0.0001,0.0000,0.0001,0.0004,0.0002,   &
     0.0005,0.0008,0.0001,0.0000,0.0001,0.0002,0.0007,0.0021,   &
     0.0032,0.0059,0.0157,0.0438,0.1324,0.3255,0.5775,0.7513,   &
     0.8901,0.9444,0.9849,0.9801,0.8595,0.6550,0.3439,0.1200,   &
     0.0435,0.0160,0.0062,0.0025,0.0006,0.0002,0.0002,0.0001,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0001,0.0006,0.0005,0.0001,0.0000,0.0000,0.0000,   &
     0.0001,0.0003,0.0008,0.0007,0.0002,0.0001,0.0001,0.0001,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,                                                    &
     1368*0./

    data (sr(3,l),l=1,1501)/ 52*0.,                             &
     0.0006,0.0005,0.0003,0.0003,0.0004,0.0001,0.0000,0.0000,   &
     0.0000,0.0000,0.0003,0.0003,0.0002,0.0008,0.0010,0.0004,   &
     0.0004,0.0005,0.0002,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0002,0.0007,0.0021,0.0075,0.0214,0.0620,0.1682,   &
     0.4015,0.7129,0.9347,0.9984,0.9760,0.9615,0.9843,0.9579,   &
     0.7082,0.3494,0.1304,0.0438,0.0180,0.0080,0.0044,0.0026,   &
     0.0019,0.0013,0.0009,0.0006,0.0003,0.0002,0.0001,0.0000,   &
     0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,   &
     0.0002,0.0002,0.0004,0.0012,0.0027,0.0050,0.0057,0.0041,   &
     0.0028,0.0019,0.0017,0.0016,0.0015,0.0015,0.0014,0.0011,   &
     0.0008,0.0005,0.0003,0.0002,0.0001,0.0001,0.0001,0.0001,   &
     0.0002,0.0001,0.0001,0.0001,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0001,0.0001,0.0001,0.0001,   &
     0.0001,0.0001,0.0001,0.0001,0.0001,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0001,0.0002,0.0002,0.0002,0.0001,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0001,0.0001,0.0002,0.0002,0.0001,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,                                                    &
     1192*0./

    data (sr(4,l),l=1,1501)/ 52*0.,                             &
     0.0004,0.0003,0.0004,0.0004,0.0003,0.0003,0.0003,0.0001,   &
     0.0001,0.0002,0.0001,0.0000,0.0000,0.0000,0.0001,0.0002,   &
     0.0007,0.0011,0.0005,0.0006,0.0010,0.0006,0.0003,0.0007,   &
     0.0031,0.0004,0.0002,0.0001,0.0001,0.0001,0.0000,0.0001,   &
     0.0001,0.0002,0.0004,0.0006,0.0009,0.0011,0.0012,0.0013,   &
     0.0017,0.0024,0.0048,0.0084,0.0187,0.0446,0.1264,0.3214,   &
     0.6542,0.9315,0.9909,0.9685,0.9747,0.9953,0.9943,0.9131,   &
     0.6639,0.3417,0.1349,0.0442,0.0147,0.0052,0.0021,0.0005,   &
     0.0003,0.0002,0.0002,0.0001,0.0001,0.0001,0.0001,0.0001,   &
     0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,   &
     0.0001,0.0001,0.0001,0.0002,0.0003,0.0006,0.0007,0.0007,   &
     0.0005,0.0003,0.0001,0.0001,0.0001,0.0001,0.0001,0.0002,   &
     0.0002,0.0001,0.0001,0.0001,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0001,0.0001,0.0001,0.0001,0.0001,   &
     0.0001,0.0001,0.0001,0.0001,0.0001,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0001,0.0001,0.0001,0.0001,0.0001,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0001,0.0001,0.0002,0.0002,0.0002,0.0001,0.0001,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,                                                    &
     1192*0./

    data (sr(5,l),l=1,1501)/ 53*0.,                             &
     0.0001,0.0002,0.0003,0.0002,0.0006,0.0008,0.0004,0.0008,   &
     0.0012,0.0010,0.0007,0.0006,0.0012,0.0016,0.0012,0.0015,   &
     0.0021,0.0028,0.0032,0.0033,0.0034,0.0035,0.0035,0.0030,   &
     0.0024,0.0023,0.0024,0.0025,0.0017,0.0012,0.0016,0.0027,   &
     0.0034,0.0024,0.0017,0.0013,0.0022,0.0022,0.0006,0.0002,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0001,0.0001,0.0003,0.0011,0.0032,0.0065,0.0148,0.0438,   &
     0.1398,0.4291,0.7892,0.9311,0.9446,0.9682,0.9986,0.9200,   &
     0.6663,0.3553,0.1425,0.0592,0.0270,0.0151,0.0095,0.0054,   &
     0.0029,0.0018,0.0013,0.0013,0.0010,0.0007,0.0005,0.0003,   &
     0.0002,0.0003,0.0004,0.0004,0.0002,0.0002,0.0003,0.0002,   &
     0.0001,0.0002,0.0002,0.0002,0.0002,0.0002,0.0004,0.0006,   &
     0.0008,0.0009,0.0008,0.0006,0.0008,0.0013,0.0020,0.0033,   &
     0.0044,0.0042,0.0035,0.0027,0.0022,0.0019,0.0016,0.0020,   &
     0.0027,0.0040,0.0057,0.0073,0.0069,0.0059,0.0046,0.0037,   &
     0.0030,0.0024,0.0018,0.0012,0.0009,0.0007,0.0003,0.0002,   &
     0.0002,0.0002,0.0001,0.0001,0.0001,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0001,0.0002,0.0002,0.0002,   &
     0.0001,0.0001,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0001,0.0001,0.0001,0.0001,   &
     0.0001,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0001,   &
     0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,   &
     0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,   &
     0.0001,0.0001,0.0001,0.0000,0.0000,0.0000,0.0000,0.0001,   &
     0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,   &
     0.0000,0.0000,0.0000,0.0001,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,                               &
     1140*0./

    data (sr(6,l),l=1,1501)/ 53*0.,                             &
     0.0000,0.0001,0.0001,0.0001,0.0002,0.0002,0.0002,0.0002,   &
     0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0002,   &
     0.0003,0.0003,0.0003,0.0002,0.0002,0.0003,0.0005,0.0006,   &
     0.0006,0.0004,0.0002,0.0002,0.0002,0.0003,0.0005,0.0005,   &
     0.0006,0.0008,0.0009,0.0009,0.0011,0.0012,0.0011,0.0011,   &
     0.0011,0.0014,0.0015,0.0013,0.0006,0.0001,0.0000,0.0001,   &
     0.0001,0.0003,0.0003,0.0001,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0001,0.0004,0.0004,0.0001,0.0001,   &
     0.0002,0.0001,0.0003,0.0007,0.0007,0.0006,0.0004,0.0005,   &
     0.0005,0.0002,0.0003,0.0005,0.0004,0.0004,0.0005,0.0005,   &
     0.0003,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,   &
     0.0002,0.0003,0.0004,0.0002,0.0001,0.0001,0.0001,0.0001,   &
     0.0001,0.0003,0.0006,0.0007,0.0005,0.0003,0.0003,0.0010,   &
     0.0020,0.0040,0.0122,0.0174,0.0349,0.1116,0.3640,0.7775,   &
     0.9903,0.9545,0.9641,0.9947,0.9807,0.8514,0.5827,0.3029,   &
     0.1355,0.0653,0.0330,0.0183,0.0105,0.0072,0.0047,0.0035,   &
     0.0024,0.0015,0.0010,0.0007,0.0007,0.0006,0.0005,0.0004,   &
     0.0004,0.0002,0.0001,0.0001,0.0001,0.0000,0.0000,0.0000,   &
     0.0001,0.0000,0.0001,0.0001,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0001,0.0001,0.0001,0.0001,   &
     0.0001,0.0001,0.0001,0.0001,0.0000,0.0000,0.0000,0.0001,   &
     0.0001,0.0001,0.0002,0.0002,0.0002,0.0001,0.0001,0.0001,   &
     0.0001,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0001,0.0001,0.0001,0.0001,0.0001,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0001,0.0001,0.0001,0.0001,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0001,0.0001,   &
     0.0001,0.0001,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0001,0.0001,0.0001,0.0001,0.0001,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,   &
     0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,                               &
     1180*0./

    data (sr(7,l),l=1,1501)/ 52*0.,                             &
     0.0003,0.0001,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0001,0.0001,0.0002,0.0004,0.0005,   &
     0.0006,0.0005,0.0003,0.0003,0.0002,0.0002,0.0001,0.0001,   &
     0.0001,0.0001,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0002,0.0008,0.0013,0.0019,   &
     0.0020,0.0021,0.0024,0.0020,0.0012,0.0008,0.0008,0.0011,   &
     0.0018,0.0024,0.0018,0.0012,0.0005,0.0003,0.0001,0.0001,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0001,0.0001,0.0001,0.0002,0.0003,0.0005,0.0011,   &
     0.0019,0.0033,0.0056,0.0093,0.0169,0.0306,0.0572,0.1043,   &
     0.1895,0.3276,0.5220,0.7086,0.8383,0.9120,0.9363,0.9353,   &
     0.9380,0.9356,0.9262,0.9210,0.9387,0.9655,0.9922,0.9609,   &
     0.8382,0.6797,0.4984,0.3299,0.2107,0.1310,0.0838,0.0556,   &
     0.0393,0.0278,0.0200,0.0164,0.0143,0.0118,0.0097,0.0081,   &
     0.0071,0.0063,0.0054,0.0048,0.0043,0.0040,0.0036,0.0032,   &
     0.0028,0.0026,0.0025,0.0024,0.0020,0.0019,0.0019,0.0020,   &
     0.0019,0.0019,0.0018,0.0019,0.0020,0.0020,0.0019,0.0018,   &
     0.0019,0.0019,0.0020,0.0021,0.0021,0.0021,0.0020,0.0018,   &
     0.0018,0.0019,0.0019,0.0019,0.0020,0.0021,0.0022,0.0022,   &
     0.0024,0.0026,0.0029,0.0035,0.0042,0.0049,0.0055,0.0059,   &
     0.0057,0.0054,0.0050,0.0045,0.0039,0.0033,0.0029,0.0024,   &
     0.0021,0.0018,0.0017,0.0016,0.0015,0.0015,0.0015,0.0016,   &
     0.0016,0.0016,0.0018,0.0019,0.0020,0.0022,0.0022,0.0020,   &
     0.0016,0.0012,0.0010,0.0007,0.0005,0.0004,0.0003,0.0003,   &
     0.0002,0.0002,0.0001,0.0001,0.0001,0.0001,0.0001,0.0000,   &
     0.0000,0.0000,0.0001,0.0001,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0001,0.0001,0.0001,0.0001,0.0002,0.0002,   &
     0.0001,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0001,0.0002,0.0002,0.0001,0.0001,0.0001,0.0002,   &
     0.0001,0.0002,0.0005,0.0006,0.0006,0.0006,0.0006,0.0003,   &
     0.0001,0.0000,0.0000,0.0000,0.0001,                        &
     1140*0./

    data (sr(8,l),l=1,1501)/ 52*0.,                             &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0002,0.0008,0.0010,0.0010,0.0012,0.0010,0.0006,   &
     0.0006,0.0006,0.0005,0.0004,0.0003,0.0002,0.0002,0.0001,   &
     0.0001,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0001,0.0001,0.0002,0.0002,0.0002,0.0003,0.0003,0.0004,   &
     0.0005,0.0008,0.0015,0.0025,0.0035,0.0044,0.0059,0.0097,   &
     0.0085,0.0067,0.0072,0.0066,0.0052,0.0041,0.0035,0.0037,   &
     0.0033,0.0029,0.0037,0.0035,0.0027,0.0017,0.0011,0.0011,   &
     0.0016,0.0022,0.0028,0.0030,0.0026,0.0014,0.0007,0.0005,   &
     0.0004,0.0002,0.0002,0.0001,0.0001,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,   &
     0.0000,0.0000,0.0000,0.0001,0.0003,0.0008,0.0013,0.0015,   &
     0.0017,0.0018,0.0028,0.0075,0.0158,0.0233,0.0231,0.0107,   &
     0.0056,0.0030,0.0026,0.0027,0.0035,0.0064,0.0103,0.0152,   &
     0.0132,0.0092,0.0043,0.0018,0.0012,0.0008,0.0007,0.0009,   &
     0.0008,0.0006,0.0008,0.0007,0.0006,0.0007,0.0009,0.0008,   &
     0.0009,0.0011,0.0012,0.0027,0.0021,0.0025,0.0039,0.0037,   &
     0.0043,0.0057,0.0079,0.0102,0.0157,0.0245,0.0422,0.0739,   &
     0.1388,0.2544,0.4430,0.6626,0.8436,0.9476,0.9851,0.9964,   &
     0.9962,0.9827,0.9602,0.9206,0.8634,0.8359,0.8598,0.9266,   &
     0.9736,0.9164,0.7136,0.4528,0.2495,0.1301,0.0733,0.0435,   &
     0.0289,0.0199,0.0150,0.0112,0.0073,0.0057,0.0046,0.0041,   &
     0.0035,0.0029,0.0028,0.0025,0.0023,0.0020,0.0018,0.0017,   &
     0.0017,0.0018,0.0017,0.0017,0.0017,0.0016,0.0015,0.0016,   &
     0.0016,0.0016,0.0015,0.0016,0.0018,0.0017,0.0017,0.0019,   &
     0.0020,0.0021,0.0021,0.0022,0.0023,0.0025,0.0027,0.0027,   &
     0.0028,0.0028,0.0027,0.0026,0.0025,0.0024,0.0023,0.0021,   &
     0.0021,0.0021,0.0022,0.0022,0.0023,0.0025,0.0025,0.0025,   &
     0.0023,0.0020,0.0019,0.0015,0.0011,0.0009,0.0008,0.0009,   &
     0.0008,0.0006,0.0005,0.0005,0.0004,0.0005,0.0004,0.0003,   &
     0.0003,0.0003,0.0002,0.0001,0.0001,0.0001,0.0001,0.0002,   &
     0.0002,0.0001,0.0001,0.0001,0.0000,0.0000,0.0000,0.0002,   &
     0.0002,0.0002,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,   &
     0.0000,0.0000,0.0001,0.0001,0.0000,                        &
     1140*0./

    wli(1)=0.3825
    wls(1)=0.7
    wli(2)=0.38
    wls(2)=0.58
    wli(3)=0.38
    wls(3)=1.02
    wli(4)=0.38
    wls(4)=1.02
    wli(5)=0.3825
    wls(5)=1.15
    wli(6)=0.3825
    wls(6)=1.05
    wli(7)=0.38
    wls(7)=1.15
    wli(8)=0.38
    wls(8)=1.15
    do i=1,1501
        s(i)=sr(iwa,i)
    enddo
    wlinf=wli(iwa)
    wlsup=wls(iwa)
    return
end
