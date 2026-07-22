module param_mt
use param
implicit none

type param_forward_mt ! 2021.12.14
 ! mesh info
 integer(4)    :: itopoflag    ! 0: no topo, 1: topo are from topofiles 2017.09.29
 integer(4)    :: nfile        ! 2017.09.27
 character(50),allocatable,dimension(:)   :: topofile ! 2017.09.27
 real(8),      allocatable,dimension(:,:) :: lonlatshift  ! 2017.09.27
 character(50) :: g_meshfile   ! global mesh file
 integer(4)    :: surface_id_ground=6 ! surface consisting of triangles for ground
 character(50) :: z_meshfile    ! 2d triangle mesh file
 character(50) :: g_lineinfofile ! added on April 21, 2016
 character(50) :: outputfolder
 character(50) :: header2d
 character(50) :: header3d
 ! frequency info
 integer(4)    :: nfreq ! # of frequency for this
 real(8),allocatable,dimension(:)   :: freq ! frequencies
 ! observatory info
 integer(4)    :: nobs  ! # of observatories
 integer(4)    :: lonlatflag !# 1: lonlat [deg], 2: xy [km]
 real(8)       :: wlon,elon,slat,nlat ! [deg]
 real(8)       :: lonorigin,latorigin ! [deg]
 real(8)       :: lenout    ! [km]
 real(8)       :: upzin     ! [km]
 real(8)       :: downzin   ! [km]
 real(8)       :: zmax      ! [km]
 real(8)       :: zmin      ! [km]
 real(8)       :: sizein    ! [km]
 real(8)       :: sizebo    ! [km]
 real(8)       :: sigma_obs ! [km] ! converted to dimension
 real(8)       :: A_obs     ! [km] ! on2016.10.20
 real(8)       :: dlen_source ! [km]  2017.10.12
 real(8)       :: sigma_src ! [km]
 real(8)       :: A_src     ! [km]
 integer(4)    :: nobsr
 real(8),      allocatable,dimension(:,:) :: xyz_r
 real(8),      allocatable,dimension(:)   :: A_r,sigma_r
 real(8)                                  :: xbound(4),ybound(4),zbound(4)
 character(3)                             :: UTM ! UTM zone
 real(8),      allocatable,dimension(:,:) :: lonlataltobs!lon[deg],lat[deg],alt[km]
 real(8),      allocatable,dimension(:,:) :: xyzobs ! x eastward[km],y northward[km],z[km]
 character(50),allocatable,dimension(:)   :: obsname
 !# output spatial distribution 2017.10.11
 integer(4)    :: ixyflag      ! 0: nothing, 1: (nx,ny) surface values, 2: triangle surface
 character(50) :: xyfilehead   ! 2017.10.11 (ixyflag = 1 or 2)
 integer(4)    :: nx,ny        ! 2017.10.11 (ixyflag = 1     )
 !#
 real(8)       :: xyzminmax(6) ! added on 2017.02.21
 real(8)       :: zorigin=0.d0 ! added on 2017.02.21
 !# conductivity structure
 character(50) :: condfile
 integer(4)    :: condflag
end type

!type param_bell ! inserted on 2016.10.12 for Kusatsushirane
! character(50) :: bellgeofile
! integer(4) :: ilonlatflag ! 1 for lonlat, 2 for xy [km]
! real(8)  :: lon_bell, lat_bell
! real(8)  :: radius![km]
! real(8)  :: width ![km]
! real(8)  :: ztop  ! < 0 [km] how deep is the top of bell
! real(8)  :: zbot  ! < 0 [km] how deep is the bottom of bell
! real(8)  :: zthic ! [km] thickness of top of the bell
! real(8)  :: reso_bell ! [km] mesh resolution for inside bell
!--
! real(8)  :: xyz_bell(3)
!end type

contains

!################################################## subroutine readparam
!# change the
subroutine readparam_mt(c_param,g_cond) ! 2021.12.15
implicit none
type(param_forward_mt),       intent(out) :: c_param
type(param_cond),   optional, intent(out) :: g_cond
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
write(*,*) "<input 0 for no topo, 1 for topofiles (itopoflag)>" ! 2018.10.02
read(input,*)    c_param%itopoflag      ! 11->* 2021.09.02
write(*,'(a,i3)') " itopoflag =", c_param%itopoflag    ! 2020.09.17

if ( c_param%itopoflag .eq. 0 ) then ! 2017.09.29
 goto 101 ! 2017.09.29
else if (c_param%itopoflag .eq. 1 ) then ! 2017.09.29
  read(input,*) c_param%nfile            ! 12->* 2021.09.02
  allocate(c_param%topofile(c_param%nfile) )            ! 2017.09.27
  allocate(c_param%lonlatshift(2,c_param%nfile) )       ! 2017.09.27
 do i=1,c_param%nfile                                  ! 2017.09.27
  read(input,10) c_param%topofile(i)
  read(input,*) c_param%lonlatshift(1:2,i)             ! 2021.09.28
 end do                                                ! 2017.09.27
else         ! 2017.09.29
 write(*,*) "GEGEGE! itopoflag should be 0 or 1: itopoflag",c_param%itopoflag!2017.09.29
 stop        ! 2017.09.29
end if       ! 2017.09.29
101 continue ! 2017.09.29
write(*,*) ""
write(*,   40) " <Please enter global meshfile>"
read(input,10) c_param%g_meshfile
write(*,   41) " global mesh file : ",trim(c_param%g_meshfile)     ! 2020.09.29
read(input,10) c_param%z_meshfile
write(*,   41) " 2dz    mesh file : ",trim(c_param%z_meshfile)     ! 2020.09.29
read(input,10) c_param%g_lineinfofile
write(*,   41) " line info   file : ",trim(c_param%g_lineinfofile) ! 2020.09.29
read(input,10) c_param%outputfolder
write(*,   41) " output folder    : ",trim(c_param%outputfolder)   ! 2020.09.29
read(input,10) c_param%header2d
write(*,   41) "header2d          : ", trim(c_param%header2d)      ! 2020.10.31
read(input,10) c_param%header3d
write(*,   41) "header3d          : ", trim(c_param%header3d)      ! 2020.10.31

!# frequency info
read(input,*) c_param%nfreq          ! 12->* 2021.09.02
write(*,'(a,i3)') " # of frequency : ",c_param%nfreq
allocate(c_param%freq(c_param%nfreq))
do i=1, c_param%nfreq
 read(input,*) c_param%freq(i)        ! 12->* 2021.09.02
 write(*,'(f9.4,a)')     c_param%freq(i)," [Hz]" ! 2020.09.29
end do

!# xy plane boundary by xy coordinate [km]
 write(*,*) ""
 write(*,*) "< input xbound and ybound >"    ! 2020.09.29
 read(input,*) c_param%xbound(2)             ! 12->* 2021.09.02
 read(input,*) c_param%xbound(3)             ! 12->* 2021.09.02
 read(input,*) c_param%ybound(2)             ! 12->* 2021.09.02
 read(input,*) c_param%ybound(3)             ! 12->* 2021.09.02
 write(*,43) " xbound(2)=",c_param%xbound(2) ! 2020.09.29
 write(*,43) " xbound(3)=",c_param%xbound(3) ! 2020.09.29
 write(*,43) " ybound(2)=",c_param%ybound(2) ! 2020.09.29
 write(*,43) " ybound(3)=",c_param%ybound(3) ! 2020.09.29

!# mesh info
 write(*,*) ""
 write(*,*) "< input meshinfo >" ! 2020.09.29
 read(input,*) c_param%lenout    ! 12->* 2021.09.02
 read(input,*) c_param%upzin     ! 12->* 2021.09.02
 read(input,*) c_param%downzin   ! 12->* 2021.09.02
 read(input,*) c_param%zmax      ! 12->* 2021.09.02
 read(input,*) c_param%zmin      ! 12->* 2021.09.02
 read(input,*) c_param%sizein    ! 12->* 2021.09.02
 read(input,*) c_param%sizebo    ! 12->* 2021.09.02
 read(input,*) c_param%sigma_obs ! commented out on 2016.10.20 ! 12->* 2021.09.02
 read(input,*) c_param%A_obs     ! commented out on 2016.10.20 ! 12->* 2021.09.02
 read(input,*) c_param%dlen_source ! [km] 2017.10.12 ! 12->* 2021.09.02
 read(input,*) c_param%sigma_src ! 12->* 2021.09.02
 read(input,*) c_param%A_src     ! 12->* 2021.09.02
 write(*,43) " lenout      =",c_param%lenout ! 2020.09.29
 write(*,43) " upzin       =",c_param%upzin
 write(*,43) " downzin     =",c_param%downzin
 write(*,43) " zmax        =",c_param%zmax
 write(*,43) " zmin        =",c_param%zmin
 write(*,43) " sizein      =",c_param%sizein
 write(*,43) " sizebo      =",c_param%sizebo
 write(*,43) " sigma_obs   =",c_param%sigma_obs ! commented out on 2016.10.20
 write(*,43) " A_obs       =",c_param%A_obs     ! commented out on 2016.10.20
 write(*,43) " dlen_source =",c_param%dlen_source ! [km] 2017.10.12
 write(*,43) " sigma_src   =",c_param%sigma_src
 write(*,43) " A_src       =",c_param%A_src

 !# calculate xbound(1:4), ybound(1:4), zbound(1:4)
 c_param%zbound(1:4) = (/c_param%zmin,c_param%downzin,c_param%upzin,c_param%zmax/)
 c_param%xbound(1) =   c_param%xbound(2) - c_param%lenout
 c_param%xbound(4) = - c_param%xbound(1)
 c_param%ybound(1) =   c_param%ybound(2) - c_param%lenout
 c_param%ybound(4) = - c_param%ybound(1)
 write(*,*) "" ! 2020.09.29
 write(*,'(a,4f9.4)') " xbound =",c_param%xbound(1:4) ! 2021.09.02
 write(*,'(a,4f9.4)') " ybound =",c_param%ybound(1:4) ! 2021.09.02
 write(*,'(a,4f9.4)') " zbound =",c_param%zbound(1:4) ! 2021.09.02
 write(*,*) "" ! 2021.09.02
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
 nobs = c_param%nobs

!# read site information
if (c_param%lonlatflag .eq. 1) write(*,'(a)') "< conversion of Lon Lat to UTM x (east), y (north) [km]>"!2021.09.29
write(*,*) "" !2021.09.29

 do i=1, nobs
   read(input,10) c_param%obsname(i)
   site=c_param%obsname(i) ! 2021.09.02
   write(*,*) i,c_param%obsname(i) ! 2021.12.17
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

!#[3]## for source parameters

!#[4]## lonlatflag : 1 -> 2 because xyzobs is laready set
 c_param%lonlatflag = 2

!#[5]## read conductivity information
 if ( present(g_cond) ) then     ! if the optional argument, g_cond, is present
  read(input,*) g_cond%sigma_air ! 11->* 2021.09.02
  write(*,'(a,g15.7,a)') " sigma_air =",g_cond%sigma_air," [S/m]"
  read(input,*) g_cond%condflag  ! 0:homogeneous, 1:file
  if (g_cond%condflag .eq. 0 ) then                         ! 2017.09.29
   write(*,*) "" ! 2020.09.29
   write(*,*) "<Input # of physical volumes in land region>"
   read(input,*) g_cond%nvolume ! # of physical volume in land, ! 11->* 2021.09.02
   write(*,*) "nvolume=",g_cond%nvolume
   allocate( g_cond%sigma_land(g_cond%nvolume) )                  ! 2017.09.28
   write(*,*) "" ! 2020.09.29
   write(*,*) "<Inuput land sigma [S/m] for each physical volume>"  ! 2017.09.28
   do i=1,g_cond%nvolume
    read(input,*) g_cond%sigma_land(i) ! conductivity in land region 2017.09.28, ! 12->* 2021.09.02
    write(*,*) i,"sigma_land=",g_cond%sigma_land(i),"[S/m]"
   end do
  else if (g_cond%condflag .eq. 1) then ! file
   read(input,'(a50)') g_cond%condfile  ! 2021.10.04
   write(*,*) "cond file is",g_cond%condfile
   CALL READCOND(g_cond)          ! read conductivity structure
  else
   write(*,*) "GEGEGE condflag should be 0 or 1 : condflag=",g_cond%condflag
   stop
  end if
 end if

 close(input)
 write(*,*) "### READ FORWAR PARAM END!! ###"

return

99 continue
write(*,*) "GEGEGE!"
write(*,*) "c_param%lonlatflag=",c_param%lonlatflag
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
44 format(a,3f9.4)
end subroutine readparam_mt

end module param_mt
