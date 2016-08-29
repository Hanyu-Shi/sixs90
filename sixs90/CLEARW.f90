subroutine clearw (r)
    implicit none
    real(8) :: sr(1501),r(1501)
    integer :: l,i

    data (sr(l),l=1,135)/  58*0.,                              &
     .00000, .02050, .04100, .04100, .04100, .04100, .04100,   &
     .04100, .04100, .04100, .04100, .04100, .04100, .04100,   &
     .04100, .04100, .04100, .04100, .04100, .04100, .04100,   &
     .04100, .04100, .04100, .04100, .04100, .04100, .04100,   &
     .04100, .04100, .04100, .04100, .04100, .04100, .04100,   &
     .04100, .04100, .04100, .04100, .04100, .04100, .04100,   &
     .04100, .04150, .04200, .04250, .04300, .04350, .04400,   &
     .04400, .04400, .04500, .04600, .04650, .04700, .04800,   &
     .04900, .04950, .05000, .05100, .05200, .05300, .05400,   &
     .05450, .05500, .05550, .05600, .05750, .05900, .05950,   &
     .06000, .06050, .06100, .06100, .06100, .06000, .05900/
    data (sr(l),l=136,1501)/                                   &
     .05800, .05700, .05550, .05400, .05350, .05300, .05200,   &
     .05100, .05050, .05000, .04950, .04900, .04800, .04700,   &
     .04650, .04600, .04600, .04600, .04550, .04500, .04450,   &
     .04400, .04350, .04300, .04300, .04300, .04200, .04100,   &
     .04050, .04000, .03900, .03800, .03750, .03700, .03700,   &
     .03700, .03650, .03600, .03450, .03300, .03250, .03200,   &
     .03150, .03100, .03000, .02900, .02800, .02700, .02550,   &
     .02400, .02350, .02300, .02200, .02100, .01950, .01800,   &
     .01650, .01500, .01350, .01200, .01050, .00900, .00850,   &
     .00800, .00700, .00600, .00500, .00400, .00300, .00200,   &
     .00150, .00100, .00050, .00000, .00000, .00000, .00000,   &
     .00000, .00000, .00000, .00000, .00000, .00000, .00000,   &
     .00000, .00000,                                           &
     1280*0./
    do i=1,1501
        r(i)=sr(i)
    enddo
    return
 end
