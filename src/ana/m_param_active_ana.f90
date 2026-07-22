module param_active_ana
implicit none

type param_forward_ana
 character(50) :: outputfolder
 ! frequency info
 integer(4)    :: nfreq ! # of frequency for this
 real(8),allocatable,dimension(:)   :: freq ! frequencies
 ! observatory info
 integer(4)    :: nobs  ! # of observatories
 integer(4)    :: lonlatflag !# 1: lonlat [deg], 2: xy [km]
 real(8)       :: lonorigin,latorigin ! [deg]
 integer(4)    :: nobsr
 real(8),      allocatable,dimension(:,:) :: xyz_r
 real(8),      allocatable,dimension(:)   :: A_r,sigma_r
 character(3)                             :: UTM ! UTM zone
 real(8),      allocatable,dimension(:,:) :: lonlataltobs!lon[deg],lat[deg],alt[km]
 real(8),      allocatable,dimension(:,:) :: xyzobs ! x eastward[km],y northward[km],z[km]
 character(50),allocatable,dimension(:)   :: obsname
 !# output spatial distribution 2017.10.11
 integer(4)    :: ixyflag      ! 0: nothing, 1: (nx,ny) surface values, 2: triangle surface
 character(50) :: xyfilehead   ! 2017.10.11 (ixyflag = 1 or 2)
 integer(4)    :: nx,ny        ! 2017.10.11 (ixyflag = 1     )
 !#
 real(8)       :: zorigin=0.d0 ! added on 2017.02.21

 integer(4)    :: nlayer ! 0: (partitioned) homegeneous earth, 1: file is given
 real(8)       :: sigma_air=1.d-8 ! [S/m]
 integer(4),allocatable,dimension(:) :: index ! element id for nphys 2, 2017.05.10
 real(8),   allocatable,dimension(:) :: sigma ! [S/m]
 real(8),   allocatable,dimension(:) :: sigma_depth ! [km] 2021.12.19
 real(8),   allocatable,dimension(:) :: rho   ! [Ohm.m]
end type param_forward_ana

type param_source
 integer(4) :: lonlatflag !# 1:lonlat [deg], 2: xy [km]
 integer(4) :: nsource    ! # of sources  added on 2017.07.11
 character(50),allocatable,dimension(:) :: sourcename ! 2017.07.11
 real(8),allocatable,dimension(:,:) :: xs1  ! [km] xyz coordinate of start electrode
 real(8),allocatable,dimension(:,:) :: xs2  ! [km] xyz coordinate of end electrode
 real(8),allocatable,dimension(:,:) :: lonlats1 ! [deg] lonlat of start electrode
 real(8),allocatable,dimension(:,:) :: lonlats2 ! [deg] xyz coordinate of end electrode
 real(8)                            :: I ! source current through wire [A]
 real(8) :: ds ! [m] dipole divide len for analytical cal 2021.12.20
end type param_source

type respdata
integer(4) :: nobs
complex(8),allocatable,dimension(:) :: ftobs  ! 2021.09.15
real(8),   allocatable,dimension(:) :: fpobsamp
real(8),   allocatable,dimension(:) :: fpobsphase
real(8),   allocatable,dimension(:) :: fsobsamp
real(8),   allocatable,dimension(:) :: fsobsphase
real(8),   allocatable,dimension(:) :: ftobsamp
real(8),   allocatable,dimension(:) :: ftobsphase
end type

contains
!######################################### ALLOCATERESPDATA
subroutine ALLOCATERESPDATA(nobs,resp)
implicit none
integer(4),    intent(in) :: nobs
type(respdata),intent(out) :: resp

resp%nobs=nobs
allocate(resp%ftobs     (resp%nobs)) ! 2021.09.15
allocate(resp%fpobsamp  (resp%nobs))
allocate(resp%fpobsphase(resp%nobs))
allocate(resp%fsobsamp  (resp%nobs))
allocate(resp%fsobsphase(resp%nobs))
allocate(resp%ftobsamp  (resp%nobs))
allocate(resp%ftobsphase(resp%nobs))

return
end subroutine
!################################################## subroutine readparam
!# 2021.12.19
subroutine readparam_ana(c_param,sparam)
implicit none
type(param_forward_ana),          intent(out) :: c_param
type(param_source),               intent(out) :: sparam
integer(4)                  :: i,j,nobs,input=2,nsource
integer(4)                  :: l ! 2020.09.28
character(50)               :: site
real(8)                     :: lonorigin,latorigin,a
character(100)              :: paramfile
integer(4),    parameter    :: n = 1000 ! 2020.09.28
character(200),dimension(n) :: lines    ! 2020.09.28
!#
write(*,*) ""
write(*,*) "<Please input the forward parameter file>" ! 2020.09.28
read(*,'(a)') paramfile           ! 2020.09.28
call readcontrolfile(paramfile,n,lines) ! 2020.09.28

open(input,file="tmp.ctl")
!#[2]# read mesh_param file
write(*,*) ""

read(input,10) c_param%outputfolder
write(*,   41) " output folder    : ",trim(c_param%outputfolder)   ! 2020.09.29

!# frequency info
read(input,*) c_param%nfreq          ! 12->* 2021.09.02
write(*,'(a,i3)') " # of frequency : ",c_param%nfreq
allocate(c_param%freq(c_param%nfreq))
do i=1, c_param%nfreq
 read(input,*) c_param%freq(i)        ! 12->* 2021.09.02
 write(*,'(f9.4,a)')     c_param%freq(i)," [Hz]" ! 2020.09.29
end do

!# observatory info
 read(input,*) c_param%nobs ! 11->* 2021.09.02
 write(*,'(a,i3)') " # of observatories : ",c_param%nobs ! 2020.09.17
 allocate(c_param%lonlataltobs(3,c_param%nobs))
 allocate(c_param%xyzobs      (3,c_param%nobs))
 allocate(c_param%obsname(c_param%nobs))
! allocate(c_param%A_obs(c_param%nobs))       ! added on 2016.10.20
! allocate(c_param%sigma_obs(c_param%nobs))   ! added on 2016.10.20

!# lonlatflag
 read(input,*)  c_param%lonlatflag    ! 1: lonlatalt, 2:xyz ! 11->* 2021.09.02
 if ( c_param%lonlatflag .eq. 1) then ! lonlat
  read(input,*) c_param%UTM           ! 10->* 2021.09.02
  write(*,*) "UTM zone : ",c_param%UTM ! 2020.09.29
  read(input,*) c_param%lonorigin,c_param%latorigin ! 12->* 2021.09.02
  write(*,*) "lonorigin =",c_param%lonorigin ! 2020.09.17
  write(*,*) "latorigin =",c_param%latorigin  ! 2020.09.17
 end if
 lonorigin = c_param%lonorigin
 latorigin = c_param%latorigin
 nobs      = c_param%nobs

!# read site information
if (c_param%lonlatflag .eq. 1) write(*,'(a)') "< conversion of Lon Lat to UTM x (east), y (north) [km]>"!2021.09.29
write(*,*) "" !2021.09.29

 do i=1, nobs
   read(input,10) c_param%obsname(i)
   site=c_param%obsname(i) ! 2021.09.02

   ! lonlatalt
   if (c_param%lonlatflag .eq. 1 ) then
    read(input,*) (c_param%lonlataltobs(j,i),j=1,3)
    write(*,'(1x,a,a,3f15.7)') trim(site)," :",c_param%lonlataltobs(1:3,i) !2021.09.02
    call UTMXY(c_param%lonlataltobs(1:2,i),&
        & lonorigin,latorigin,c_param%xyzobs(1:2,i),c_param%UTM)
    c_param%xyzobs(3,i) = c_param%lonlataltobs(3,i)
    write(*,'(1x,a,2f15.7,a)') " UTM>",c_param%xyzobs(1:2,i)," [km]" ! 2021.09.29
    write(*,*) ""  ! 2021.09.29

   ! xyz
   else if (c_param%lonlatflag .eq. 2 ) then ! xyz
    read(input,*) (c_param%xyzobs(j,i),j=1,3)

   end if
  end do

!#[2]# xy reading  2017.10.11
 write(*,*) ""
 write(*,*) "< input ixyflag: 0 for nothing, 1 for xy plan view >"
 read(input,*) c_param%ixyflag    ! 12->* 2021.09.02
 write(*,'(a,i3)') " ixyflag =",c_param%ixyflag ! 2020.09.17
 if ( c_param%ixyflag .eq. 1) then ! 2017.10.11
  read(input,*) c_param%xyfilehead     ! 2021.09.02
  read(input,*) c_param%nx,c_param%ny  ! 2021.09.02
  write(*,*) "nx,ny =",c_param%nx,c_param%ny
  write(*,*) "xyfilehead = ",c_param%xyfilehead
 end if                                ! 2017.10.11
 if ( c_param%ixyflag .eq. 2) then     ! 2017.10.11
  read(input,*) c_param%xyfilehead     ! 2021.09.02
 end if                                ! 2017.10.11

!#[3]# read source param
 sparam%lonlatflag = c_param%lonlatflag

 ! # of sources 2017.07.11
 write(*,*) ""
 write(*,*) "< input # of source wires > "!2017.07.11
 read(input,*) sparam%nsource       !2017.07.11, ! 11->* 2021.09.02
 nsource = sparam%nsource
 write(*,'(a,i3)') " # of source wires (nsource) =",nsource ! 2020.09.17
 allocate(sparam%xs1(3,nsource),    sparam%xs2(3,nsource))!2017.07.11
 allocate(sparam%lonlats1(2,nsource),sparam%lonlats2(2,nsource))!2017.07.11
 allocate(sparam%sourcename(nsource))

 ! lonlat input
 do i=1,nsource ! 2017.07.11
  read(input,*) sparam%sourcename(i) ! 10->* 2021.09.02
  write(*,41) " source name: ",sparam%sourcename(i) ! 2020.09.29
  if (sparam%lonlatflag .eq. 1) then ! lonlat
   read(input,*) sparam%lonlats1(1:2,i),sparam%xs1(3,i)
   read(input,*) sparam%lonlats2(1:2,i),sparam%xs2(3,i)
   write(*,44) " Start point (lon,lat,z)=",sparam%lonlats1(1:2,i),sparam%xs1(3,i)
   write(*,44) " End   point (lon,lat,z)=",sparam%lonlats1(1:2,i),sparam%xs2(3,i)
   call UTMXY(sparam%lonlats1(1:2,i),lonorigin,latorigin,sparam%xs1(1:2,i),c_param%UTM)
   call UTMXY(sparam%lonlats2(1:2,i),lonorigin,latorigin,sparam%xs2(1:2,i),c_param%UTM)
   write(*,*) " converted with UTM to :"
   write(*,'(a,3f15.7)') " xs1 (x,y,z)=",sparam%xs1(1:3,i) !2020.09.29
   write(*,'(a,3f15.7)') " xs2 (x,y,z)=",sparam%xs2(1:3,i) !2020.09.29
 !# xy input
 else if (sparam%lonlatflag .eq. 2) then ! xyz
  read(input,*) sparam%xs1(1:3,i)
  read(input,*) sparam%xs2(1:3,i)
 else
  goto 99
 end if
 end do ! 2017.07.11

 read(input,*) sparam%I ! [A] ! 11->* 2021.09.02
 write(*,'(a,g15.7,a)') " Electric source current :",sparam%I," [A]" ! 2020.09.29

 read(input,*) sparam%ds ! [m] 2021.12.20
 write(*,'(a,f12.5,a)') " Dipole division length ",sparam%ds, "[m]"
!#[4]## lonlatflag : 1 -> 2 because xyzobs is already set
 c_param%lonlatflag = 2


!#[5]## read conductivity information
  read(input,*) c_param%sigma_air ! 11->* 2021.09.02
  write(*,'(a,g15.7,a)') " sigma_air =",c_param%sigma_air," [S/m]"
   write(*,*) "" ! 2020.09.29
   write(*,*) "<Input # of layers in land region>"
   read(input,*) c_param%nlayer ! # of physical volume in land, ! 11->* 2021.09.02
   write(*,*) "# of layer is ",c_param%nlayer
   allocate( c_param%sigma(      c_param%nlayer) )   ! 2021.12.20
   allocate( c_param%rho(        c_param%nlayer) )   ! 2021.12.20
   allocate( c_param%sigma_depth(0:c_param%nlayer-1) )   ! 2021.12.20
   c_param%sigma_depth(0) =0.d0 ! ground surface

   write(*,*) "" ! 2020.09.29

   if ( c_param%nlayer .lt. 1 ) goto 101
   if (c_param%nlayer .ge. 2 ) then
    write(*,*) "<Inuput base depth of each layer up to nlayer-1 [km]"  !
    do i=1,c_param%nlayer-1
     read(input,*) c_param%sigma_depth(i)
     write(*,'(a,f15.7,a)') "Sigma depth", c_param%sigma_depth(i),"[km]"
    end do
   end if

   write(*,*) "<Inuput subsurface conductivity [S/m] for each layer>"  ! 2017.09.28
   do i=1,c_param%nlayer
    read(input,*) c_param%sigma(i) !2021.12.20
    write(*,'(i5,a,f15.7,a)') i,"sigma_land=",c_param%sigma(i)," [S/m]"
    c_param%rho(i)=1./c_param%sigma(i)
   end do

 close(input)
 write(*,*) "### READ FORWARD PARAMETERS END!! ###"

return

99 continue
write(*,*) "GEGEGE!"
write(*,*) "c_param%lonlatflag=",c_param%lonlatflag
write(*,*) "sparam%lonlatflag=",sparam%lonlatflag
stop
101 continue
write(*,*) "GEGEGE!"
write(*,*) "nalyer should be >=1 c_param%nlayer",c_param%nlayer ! 2021.12.20
stop
! in 10 to 32, "20x" is omitted, initial 20 characters are omitted in control file in readcontrol 2021.09.02
10 format(a)      ! 2021.09.02
11 format(i10)    ! 2021.09.02
12 format(g15.7)  ! 2021.09.02
21 format(2i10)   ! 2021.09.02
22 format(2g12.5) ! 2021.09.02
32 format(3g12.5) ! 2021.09.02
40 format(a)
41 format(a,a)
42 format(a,i3)
43 format(a,f9.4)
44 format(a,3f15.7)
end subroutine readparam_ana

!######################################################### UTMXY
!# Coded on 2016.10.12 by T.MINAMI
subroutine UTMXY(lonlat,lonorigin,latorigin,xyout,zone)
implicit none
real(8),intent(in) :: lonlat(2),lonorigin,latorigin
character(3),intent(in) :: zone
real(8) ,intent(out) :: xyout(2)
real(8) :: xorigin,yorigin,x,y

call UTMGMT(lonlat(1),lonlat(2), x,       y,      zone,0)
call UTMGMT(lonorigin,latorigin, xorigin, yorigin,zone,0)
!write(*,*) "lon=",lon,"lonorigin=",lonorigin,"lat=",lat,"latorigin=",latorigin
!write(*,*) "x=",x,"xorigin=",xorigin,"y=",y,"yorigin=",yorigin

xyout(1) = (x - xorigin)/1.d3 ! [km]
xyout(2) = (y - yorigin)/1.d3 ! [km]

!write(*,*) "xout=",xout,"yout=",yout

return
end subroutine
!######################################################### UTMGMT
subroutine UTMGMT(xin,yin,xout,yout,zone,iflag)
implicit none
real(8),intent(in) :: xin,yin
real(8),intent(out) :: xout,yout
character(3),intent(in) :: zone
character(37) :: values
integer(4),intent(in) :: iflag ! 0: LONLAT2UTM, 1:UTM2LONLAT

integer(4) ::  izone          ! 20200728
read(zone(1:2),*) izone       ! 20200728

!#[1]## prepare values
!write(values,'(g18.10,1x,g18.10)') xin,yin ! commented out 20200728

!#[2]## use mapproject
if (iflag .eq. 0 ) then ! LONLAT2UTM
! CALL system("echo "//values//" | gmt mapproject -Ju"//zone(1:3)//"/1.0 -F -C > tmp.dat")! commented out 20200728
 CALL utm_geo_tm(xin,yin,xout,yout,izone,iflag) ! 20200728
else if (iflag .eq. 1 ) then ! UTM2LONLAT
! CALL system("echo "//values//" | gmt mapproject -Ju"//zone(1:3)//"/1.0 -F -C -I > tmp.dat")! commented out 20200728
 CALL utm_geo_tm(xout,yout,xin,yin,izone,iflag) ! 20200728
else
 write(*,*) "GEGEGE! iflag should be 0 (LONLAT2UTM) or 1(UTM2LONLAT), iflag=",iflag
 stop
end if

!#[3]## read xout and yout ! commented out 20200728
! open(12,file="tmp.dat")
!  read(12,*) xout,yout
! close(12)
! call system("rm tmp.dat")

return
end subroutine

!######################################################### UTMGMT_N
subroutine UTMGMT_N(n,xin,yin,xout,yout,zone,iflag)
implicit none
integer(4),intent(in) :: n
real(8),   intent(in) :: xin(n),yin(n)
real(8),   intent(out) :: xout(n),yout(n)
character(3),intent(in) :: zone
character(37) :: values
integer(4),intent(in) :: iflag ! 0: LONLAT2UTM, 1:UTM2LONLAT
integer(4) :: i

integer(4) ::  izone          ! 20200728
read(zone(1:2),*) izone       ! 20200728
write(*,*) "izone",izone

!#[1]## prepare input file ! commented out 20200728
!open(11,file="in.dat")
! write(11,'(2g18.10)') (xin(i),yin(i),i=1,n)
!close(11)

!#[2]## use mapproject
if (iflag .eq. 0 ) then ! LONLAT2UTM
 ! CALL system("cat in.dat | gmt mapproject -Ju"//zone(1:3)//"/1.0 -F -C > tmp.dat")! commented out 20200728
 do i=1,n ! 20200728
  CALL utm_geo_tm(xin(i),yin(i),xout(i),yout(i),izone,iflag) ! 20200728
 end do   ! 20200728
else if (iflag .eq. 1 ) then ! UTM2LONLAT
 !CALL system("cat in.dat | gmt mapproject -Ju"//zone(1:3)//"/1.0 -F -C -I > tmp.dat") ! commented out 20200728
 do i=1,n ! 20200728
  CALL utm_geo_tm(xout(i),yout(i),xin(i),yin(i),izone,iflag) ! 20200728
 end do   ! 20200728
else
 write(*,*) "GEGEGE! iflag should be 0 (LONLAT2UTM) or 1(UTM2LONLAT), iflag=",iflag
 stop
end if

!#[3]## read xout and yout   commented out 20200728
! open(12,file="tmp.dat")
!  do i=1,n
!   read(12,*) xout(i),yout(i)
!  end do
! close(12)
! call system("rm tmp.dat")

return
end subroutine


!##################################################################### utm_geo
!  Added by Takuto Minami on July 28, 2020
!================================================================
!
!  UTM (Universal Transverse Mercator) projection from the USGS
!
!=====================================================================

  subroutine utm_geo_tm(rlon,rlat,rx,ry,UTM_PROJECTION_ZONE,iway) ! _tm is added 20200728 Takuto Minami

! convert geodetic longitude and latitude to UTM, and back
! use iway = ILONGLAT2UTM for long/lat to UTM, IUTM2LONGLAT for UTM to lat/long
! a list of UTM zones of the world is available at www.dmap.co.uk/utmworld.htm

  implicit none

! include "constants.h"   ! commented out 20200728 Takuto Minami
 logical :: suppress_utm_projection=.false. ! inserted on Dec. 17,2015
                        !refer page 22 of specfem3d-manual.pdf

! flag for projection from latitude/longitude to UTM, and back ! 20200728 Takuto Minami from
  integer, parameter :: ILONGLAT2UTM = 0, IUTM2LONGLAT = 1     ! 20200728 Takuto Minami
  double precision, parameter :: PI = 3.141592653589793d0      ! 20200728 Takuto Minami

!
!-----CAMx v2.03
!
!     UTM_GEO performs UTM to geodetic (long/lat) translation, and back.
!
!     This is a Fortran version of the BASIC program "Transverse Mercator
!     Conversion", Copyright 1986, Norman J. Berls (Stefan Musarra, 2/94)
!     Based on algorithm taken from "Map Projections Used by the USGS"
!     by John P. Snyder, Geological Survey Bulletin 1532, USDI.
!
!     Input/Output arguments:
!
!        rlon                  Longitude (deg, negative for West)
!        rlat                  Latitude (deg)
!        rx                    UTM easting (m)
!        ry                    UTM northing (m)
!        UTM_PROJECTION_ZONE  UTM zone
!        iway                  Conversion type
!                              ILONGLAT2UTM = geodetic to UTM
!                              IUTM2LONGLAT = UTM to geodetic
!

  integer utm_projection_zone,iway
  double precision rx,ry,rlon,rlat

  double precision, parameter :: degrad=pi/180., raddeg=180./pi
  double precision, parameter :: semimaj=6378206.4d0, semimin=6356583.8d0
  double precision, parameter :: scfa=.9996d0
  double precision, parameter :: north=0.d0, east=500000.d0

  double precision e2,e4,e6,ep2,xx,yy,dlat,dlon,zone,cm,cmr,delam
  double precision f1,f2,f3,f4,rm,rn,t,c,a,e1,u,rlat1,dlat1,c1,t1,rn1,r1,d
  double precision rx_save,ry_save,rlon_save,rlat_save

  if(suppress_utm_projection) then
    if (iway == ilonglat2utm) then
      rx = rlon
      ry = rlat
    else
      rlon = rx
      rlat = ry
    endif
    return
  endif

! save original parameters
  rlon_save = rlon
  rlat_save = rlat
  rx_save = rx
  ry_save = ry

! define parameters of reference ellipsoid
  e2=1.0-(semimin/semimaj)**2.0
  e4=e2*e2
  e6=e2*e4
  ep2=e2/(1.-e2)

  if (iway == iutm2longlat) then
    xx = rx
    yy = ry
  else
    dlon = rlon
    dlat = rlat
  endif
!
!----- Set Zone parameters
!
  zone = dble(utm_projection_zone)
  cm = zone*6.0 - 183.0
  cmr = cm*degrad
!
!---- Lat/Lon to UTM conversion
!
  if (iway == ilonglat2utm) then

  rlon = degrad*dlon
  rlat = degrad*dlat

  delam = dlon - cm
  if (delam < -180.) delam = delam + 360.
  if (delam > 180.) delam = delam - 360.
  delam = delam*degrad

  f1 = (1. - e2/4. - 3.*e4/64. - 5.*e6/256)*rlat
  f2 = 3.*e2/8. + 3.*e4/32. + 45.*e6/1024.
  f2 = f2*sin(2.*rlat)
  f3 = 15.*e4/256.*45.*e6/1024.
  f3 = f3*sin(4.*rlat)
  f4 = 35.*e6/3072.
  f4 = f4*sin(6.*rlat)
  rm = semimaj*(f1 - f2 + f3 - f4)
  if (dlat == 90. .or. dlat == -90.) then
    xx = 0.
    yy = scfa*rm
  else
    rn = semimaj/sqrt(1. - e2*sin(rlat)**2)
    t = tan(rlat)**2
    c = ep2*cos(rlat)**2
    a = cos(rlat)*delam

    f1 = (1. - t + c)*a**3/6.
    f2 = 5. - 18.*t + t**2 + 72.*c - 58.*ep2
    f2 = f2*a**5/120.
    xx = scfa*rn*(a + f1 + f2)
    f1 = a**2/2.
    f2 = 5. - t + 9.*c + 4.*c**2
    f2 = f2*a**4/24.
    f3 = 61. - 58.*t + t**2 + 600.*c - 330.*ep2
    f3 = f3*a**6/720.
    yy = scfa*(rm + rn*tan(rlat)*(f1 + f2 + f3))
  endif
  xx = xx + east
  yy = yy + north

!
!---- UTM to Lat/Lon conversion
!
  else

  xx = xx - east
  yy = yy - north
  e1 = sqrt(1. - e2)
  e1 = (1. - e1)/(1. + e1)
  rm = yy/scfa
  u = 1. - e2/4. - 3.*e4/64. - 5.*e6/256.
  u = rm/(semimaj*u)

  f1 = 3.*e1/2. - 27.*e1**3./32.
  f1 = f1*sin(2.*u)
  f2 = 21.*e1**2/16. - 55.*e1**4/32.
  f2 = f2*sin(4.*u)
  f3 = 151.*e1**3./96.
  f3 = f3*sin(6.*u)
  rlat1 = u + f1 + f2 + f3
  dlat1 = rlat1*raddeg
  if (dlat1 >= 90. .or. dlat1 <= -90.) then
    dlat1 = dmin1(dlat1,dble(90.) )
    dlat1 = dmax1(dlat1,dble(-90.) )
    dlon = cm
  else
    c1 = ep2*cos(rlat1)**2.
    t1 = tan(rlat1)**2.
    f1 = 1. - e2*sin(rlat1)**2.
    rn1 = semimaj/sqrt(f1)
    r1 = semimaj*(1. - e2)/sqrt(f1**3)
    d = xx/(rn1*scfa)

    f1 = rn1*tan(rlat1)/r1
    f2 = d**2/2.
    f3 = 5.*3.*t1 + 10.*c1 - 4.*c1**2 - 9.*ep2
    f3 = f3*d**2*d**2/24.
    f4 = 61. + 90.*t1 + 298.*c1 + 45.*t1**2. - 252.*ep2 - 3.*c1**2
    f4 = f4*(d**2)**3./720.
    rlat = rlat1 - f1*(f2 - f3 + f4)
    dlat = rlat*raddeg

    f1 = 1. + 2.*t1 + c1
    f1 = f1*d**2*d/6.
    f2 = 5. - 2.*c1 + 28.*t1 - 3.*c1**2 + 8.*ep2 + 24.*t1**2.
    f2 = f2*(d**2)**2*d/120.
    rlon = cmr + (d - f1 + f2)/cos(rlat1)
    dlon = rlon*raddeg
    if (dlon < -180.) dlon = dlon + 360.
    if (dlon > 180.) dlon = dlon - 360.
  endif
  endif

  if (iway == iutm2longlat) then
    rlon = dlon
    rlat = dlat
    rx = rx_save
    ry = ry_save
  else
    rx = xx
    ry = yy
    rlon = rlon_save
    rlat = rlat_save
  endif

  end subroutine utm_geo_tm

!############################################ whoelread ! 2020.09.28
subroutine readcontrolfile(filename,n,lines,ikeep) ! 2021.10.04 ikeep is added
implicit none
character(100),intent(in)  :: filename
integer(4),    intent(in)  :: n
character(200),intent(out) :: lines(n)
integer(4) :: i
integer(4),optional,intent(in) ::  ikeep ! 1 means keep 20 columns 2021.10.04
i=0
write(*,'(a,a)') " file is ",trim(filename)
lines(:)=" " ! initialize 2020.12.10

open(1,file=filename)
do
i=i+1
read(1,'(a)',end=99) lines(i)
end do
99 continue

close(1)

!do i=1,20
!write(*,'(a)') trim(lines(i))
!end do

call subtractcommentout(n,lines)

open(1,file="tmp.ctl")
do i=1,n
if ( present(ikeep) .and. ikeep .eq. 1 ) then !2021.10.04
 write(1,'(a)')  lines(i)(1:len_trim(lines(i))) ! cut initial 20 characters 2021.10.04
else   ! 2021.10.04
 write(1,'(a)')  lines(i)(21:len_trim(lines(i))) ! cut initial 20 characters 2021.09.02
end if ! 2021.10.04
! write(1,'(a)') trim(lines(i))
end do
close(1)

write(*,'(a)') " ### readcontrolfile END!! ###"

return
end subroutine
!############################################# subtractcommentout ! 2020.09.28
subroutine subtractcommentout(n,lines)
implicit none
integer(4),    intent(in)    :: n
character(200),intent(inout) :: lines(n)
character(200)               :: lines_out(n)
integer(4) :: i,j
j=0
do i=1,n
if ( lines(i)(1:2) .ne. "##" ) then
  j=j+1
  lines_out(j) = lines(i)
!  write(*,'(a)')trim(lines_out(j))
 end if
end do

lines = lines_out

return
end subroutine

!######################################################### UTMGMT_N2 2021.09.29
! moved from meshgen2_bell.f90 on 2021.09.29 Name is changed from UTMGMT_N to UTMGMT_N2
subroutine UTMGMT_N2(n,xin,yin,xout,yout,zone,iflag)
implicit none
integer(4),  intent(in)  :: n
real(8),     intent(in)  :: xin(n),yin(n)
real(8),     intent(out) :: xout(n),yout(n)
character(3),intent(in)  :: zone
character(37)            :: values
integer(4),intent(in)    :: iflag ! 0: LONLAT2UTM, 1:UTM2LONLAT
integer(4)               :: i
!#[1]## prepare input file
open(11,file="in.dat")
 write(11,'(2g18.10)') (xin(i),yin(i),i=1,n)
close(11)

!#[2]## use mapproject
if (iflag .eq. 0 ) then ! LONLAT2UTM
 CALL system("cat in.dat | mapproject -Ju"//zone(1:3)//"/1.0 -F -C > tmp.dat")
else if (iflag .eq. 1 ) then ! UTM2LONLAT
 CALL system("cat in.dat | mapproject -Ju"//zone(1:3)//"/1.0 -F -C -I > tmp.dat")
else
 write(*,*) "GEGEGE! iflag should be 0 (LONLAT2UTM) or 1(UTM2LONLAT), iflag=",iflag
 stop
end if

!#[3]## read xout and yout
 open(12,file="tmp.dat")
  do i=1,n
   read(12,*) xout(i),yout(i)
  end do
 close(12)
 ! call system("rm tmp.dat")
return
end



end module param_active_ana
