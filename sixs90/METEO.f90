subroutine meteo
    implicit none
    common /sixs_ffu/ s(1501),wlinf,wlsup
    real(8) :: sr(1501)
    real(8) :: s,wlinf,wlsup
    integer :: l,i

    data (sr(l),l=1,110)/  40*0.,                              &
        .00,    .00,    .00,    .01,    .01,    .01,    .02,   &
        .02,    .02,    .02,    .02,    .02,    .03,    .03,   &
        .04,    .04,    .04,    .05,    .05,    .05,    .06,   &
        .06,    .07,    .07,    .07,    .08,    .08,    .09,   &
        .09,    .10,    .10,    .10,    .11,    .11,    .12,   &
        .12,    .12,    .13,    .14,    .14,    .15,    .15,   &
        .16,    .16,    .17,    .17,    .18,    .18,    .19,   &
        .20,    .20,    .21,    .21,    .22,    .23,    .24,   &
        .24,    .25,    .26,    .27,    .28,    .28,    .29,   &
        .30,    .30,    .31,    .32,    .33,    .34,    .35/
    data (sr(l),l=111,180)/                                    &
        .35,    .36,    .37,    .38,    .39,    .40,    .40,   &
        .41,    .42,    .43,    .44,    .45,    .46,    .48,   &
        .49,    .50,    .51,    .52,    .53,    .55,    .56,   &
        .57,    .58,    .60,    .61,    .62,    .63,    .64,   &
        .65,    .65,    .66,    .67,    .67,    .68,    .69,   &
        .69,    .70,    .71,    .71,    .72,    .73,    .73,   &
        .74,    .76,    .77,    .78,    .78,    .79,    .80,   &
        .81,    .82,    .83,    .84,    .85,    .86,    .87,   &
        .88,    .89,    .89,    .91,    .92,    .93,    .94,   &
        .95,    .96,    .96,    .97,    .98,    .98,    .99/
    data (sr(l),l=181,250)/                                    &
        .99,    .99,    .99,   1.00,   1.00,   1.00,   1.00,   &
       1.00,   1.00,   1.00,   1.00,   1.00,   1.00,    .99,   &
        .99,    .99,    .99,    .98,    .98,    .98,    .98,   &
        .98,    .97,    .97,    .97,    .97,    .97,    .97,   &
        .97,    .96,    .96,    .96,    .96,    .96,    .96,   &
        .96,    .96,    .96,    .96,    .95,    .95,    .95,   &
        .94,    .93,    .93,    .92,    .92,    .91,    .90,   &
        .89,    .89,    .88,    .88,    .87,    .86,    .86,   &
        .85,    .85,    .84,    .84,    .83,    .82,    .82,   &
        .81,    .80,    .80,    .79,    .79,    .78,    .77/
    data (sr(l),l=251,1501)/                                   &
        .77,    .76,    .76,    .75,    .75,    .74,    .74,   &
        .74,    .73,    .73,    .72,    .71,    .70,    .68,   &
        .67,    .65,    .64,    .63,    .62,    .61,    .60,   &
        .59,    .58,    .57,    .56,    .55,    .54,    .53,   &
        .52,    .51,    .50,    .49,    .49,    .48,    .47,   &
        .46,    .45,    .43,    .42,    .41,    .40,    .39,   &
        .38,    .37,    .36,    .35,    .34,    .33,    .31,   &
        .30,    .29,    .28,    .28,    .27,    .25,    .24,   &
        .23,    .22,    .21,    .20,    .19,    .18,    .17,   &
        .16,    .15,    .14,    .13,    .12,    .11,    .11,   &
        .10,    .09,    .08,    .08,    .08,    .07,    .06,   &
        .06,    .05,    .05,    .05,    .04,    .04,    .03,   &
        .03,    .02,    .02,    .01,    .01,    .01,    .01,   &
        .01,    .00,    .00,    .00,                           &
      1156*0./
    wlinf=0.3499999
    wlsup=1.11
    do i=1,1501
        s(i)=sr(i)
    enddo
    return
 end
