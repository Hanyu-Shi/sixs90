subroutine posnoa (month,jday,tu,nc,xlonan,hna,campm,asol,phi0,avis,phiv,xlon,xlat)
    implicit none
    logical:: ier
    real(8) :: tu,xlonan,hna,asol,phi0,avis,phiv,xlon,xlat,campm
    real(8) :: pi,r,ai,an,ylonan,t,hnam,u,delt,d,y,z
    real(8) :: ylat,cosy,siny,ylon,ylo1,zlat,zlon,xnum,xden
    integer :: month,jday,nc,iwr
    common/sixs_ier/iwr,ier

    pi=3.1415926
    r=860./6378.155
    ai=98.96*pi/180.
    an=360.*pi/(6119.*180.)
    ylonan=xlonan*pi/180.
    t=tu*3600
    hnam=hna
    hnam=hnam*3600
    u=t-hnam
    u=campm*u*an
    delt=((nc-(2048+1)/2.)*55.385/((2048.-1)/2.))
    delt=campm*delt*pi/180.
    avis=asin((1+r)*sin(delt))
    d=avis-delt
    y=cos(d)*cos(ai)*sin(u)-sin(ai)*sin(d)
    z=cos(d)*sin(ai)*sin(u)+cos(ai)*sin(d)
    ylat=asin(z)
    cosy=cos(d)*cos(u)/cos(ylat)
    siny=y/cos(ylat)
    ylon=asin(siny)

    if(cosy.le.0.) then
        if(siny.gt.0) ylon=pi-ylon
        if(siny.le.0) ylon=-(pi+ylon)
    endif

    ylo1=ylon+ylonan-(t-hnam)*2.*pi/86400.
    xlat=ylat*180./pi
    xlon=ylo1*180./pi

    call possol(month,jday,tu,xlon,xlat,asol,phi0)
    if(ier)return

    zlat=asin(sin(ai)*sin(u))
    zlon=atan2(cos(ai)*sin(u),cos(u))
    if(nc.ne.1024)then
        xnum=sin(zlon-ylon)*cos(zlat)/sin(abs(d))
        xden=(sin(zlat)-sin(ylat)*cos(d))/cos(ylat)/sin(abs(d))
        phiv=atan2(xnum,xden)
    else
        phiv=0.
    endif
    phiv=phiv*180./pi
    avis=abs(avis)
    avis=avis*180./pi
    return
end
