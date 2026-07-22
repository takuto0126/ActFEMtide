! Coded on 2015.08.14
! to generate simple 3d mesh
!# Copied from volcano/3Dfwd_minami/mesh_test/3dgeo.f90 on Nov. 12, 2015
!# Here is study/volcano/3Dtest/
!# meshed for 2 layer model is generated
!### Please change the parameter k, the number of layers, before running this code ###
!#   when k =3 (# of layers beneath the surface is 3),
!#   ----------    z= height
!#   |        |
!#   |--------| z=0
!#   |--------| z=  z(1)
!#   |--------| z=  z(2)
!#   ---------- z=  z(3) = -depth
!###############################################################
program gen3dgeo
use param
use gmsh_geo
implicit none
!include "meshpara.f90"
type(param_forward)   :: g_param
type(param_source)    :: sparam
type(geo_info_3d)     :: geo3d ! see m_gmsh_geo.f90
character(50)         :: outgeo="3dmesh.geo"
!character(50) :: outgeo1="3dmesh_obs.geo"
!## variables for geofile #####################
integer(4)            :: k ! number of layers beneath the ground surface
real(8),   dimension(:),   allocatable :: z ! z(i) is the lower bound of i-th layer beneath the ground
real(8),   dimension(:,:), allocatable :: xyz ! x : eastward, y: northward, z: upward
real(8),   dimension(:,:), allocatable :: xyz_obs ! x : eastward, y: northward, z: upward
integer(4),dimension(:,:), allocatable :: lines, lineloop, surfaceloop
integer(4),dimension(:,:), allocatable :: splines_obs
integer(4),dimension(:),   allocatable :: volume
integer(4)    :: node, nline, nlineloop, nsurfaceloop, nvolume
integer(4)    :: node_obs,i
integer(4)    :: nobs   ! 2017.09.08
real(8)       :: xout,yout,zmax,zmin,depth,height,distance ! added on 2016.10.17
!## variables for pos file#####################
character(50) :: outpos="bgmesh.pos"
real(8)       :: lc1, lc2, lc3, lc_obs, radius
real(8)       :: height1, height2, height3, height4, dlen, dlen2,shift
!
!#[0]## read parameters
  CALL READPARAM(g_param,sparam)

!#[1]## set
  distance = g_param%xbound(4)
  height   = g_param%zbound(4)
  zmax     = g_param%zbound(4)
  zmin     = g_param%zbound(1)
  xout     = g_param%xbound(4)
  yout     = g_param%ybound(4)
  depth    = - g_param%zbound(1)

!#[0]## give k and allocate
k=1       ! # of lyaers beneath the surface
!node=8+4*k
node     = 8+4*k
nobs     = g_param%nobs  ! 2017.09.08
node_obs = nobs*4        ! 2017.09.08
nline    = 12 + 8*k ! 4*3 means the lines above the ground surcface including the surface, 2*k indicates lower
nlineloop=6+5*k
nsurfaceloop=k+1
nvolume=k+1
allocate( xyz(3,node), xyz_obs(3,node_obs) )
allocate( lines(2,nline),splines_obs(5,g_param%nobs) )
allocate( lineloop(4,nlineloop))
allocate( surfaceloop(6,nsurfaceloop), volume(nvolume), z(k) )
!z(1)=-0.2d0   ! [km]
!z(2)=-1.d0   ! [km]
!z(3)=-depth ! [km]
z(1)=-depth
!#[1]## points
call GENPOINTS(xyz, height,g_param%xbound,g_param%ybound, k, z, node)

!#[2]## lines
call GENLINES(lines, nline, k)

!#[3]## lineloops
call GENLINELOOP(lines, nline, lineloop, nlineloop, k)

!#[4]## surfaceloops
call GENSURFACELOOP(surfaceloop, nsurfaceloop, k)

!#[5]## volume
call GENVOLUME(volume, nvolume, k)

!#[OBS]## OBS points, line, lineloop
radius=0.001d0 !   1 m
lc_obs=0.0005d0 ! 50 cm
call GENPOINTS_OBS  (g_param,xyz_obs,node_obs,radius)
call GENLINES_OBS   (nobs,g_param,splines_obs,node)   ! 2017.09.08
call GENGEOINFO(xyz,node,xyz_obs,node_obs,lines,nline,nobs,splines_obs,g_param,&
                    & lineloop,nlineloop,surfaceloop,nsurfaceloop,geo3d)

!#[6]## output geofile
write(*,*) "before outgeofile"
call OUTGEOFILE(outgeo, xyz,node,lines,nline,lineloop,nlineloop,surfaceloop,nsurfaceloop,volume,nvolume,g_param)
!call OUTGEOFILE_obs(outgeo, xyz,node,lines,nline,lineloop,nlineloop,surfaceloop,nsurfaceloop,volume,nvolume,g_param)　! using the sphere

!call OUTGEOINFO3D(outgeo,geo3d) ! see m_gmsh_geo.f90

!#[7]## output posfile ##################

call calobsr(sparam,g_param) ! cal g_param%nobsr, xyz_r,sigma_r,A_r (m_param.f90)
call outbgmesh3d(outpos,g_param)

stop

!#[7]## output posfile ##################
!#  |<- dlen2->|<-dlen2->|
!#  ----------------------
!#  |   |<dlen>|<dlen>|  |   y2   ^ y
!#  |   ---------------  |   yc2  |
!#  |   |             |  |
!#  |   ---------------  |   yc1
!#  ----------------------   y1    -> x

g_param%nobs = 0

call OUTPOSFILE2_obs(outpos,g_param)
!
end program gen3dgeo
!######################################### GEN3DGEO
subroutine GENGEOINFO(xyz,node,xyz_obs,node_obs,lines,nline,nobs,splines_obs,g_param,&
                    & lineloop,nlineloop,surfaceloop,nsurfaceloop,geo3d)!2017.09.08
use gmsh_geo
use param
implicit none
type(param_forward),intent(in)  :: g_param
integer(4),         intent(in)  :: nobs    ! 2017.09.08
integer(4),         intent(in)  :: node,node_obs,nline,nlineloop,nsurfaceloop
real(8),            intent(in)  :: xyz(3,node)
real(8),            intent(in)  :: xyz_obs(3,node_obs)
integer(4),         intent(in)  :: lines(2,nline)
integer(4),         intent(in)  :: splines_obs(5,nobs)
integer(4),         intent(in)  :: lineloop(4,nlineloop)
integer(4),         intent(in)  :: surfaceloop(6,nsurfaceloop)
type(geo_info_3d),  intent(out) :: geo3d
integer(4)                      :: i,j

!#[1]# point
CALL INITPOINT(geo3d%point, node+node_obs )
geo3d%point%xyz(1:3,1:node)               =     xyz(1:3,1:node)
geo3d%point%xyz(1:3,node+1:node+node_obs) = xyz_obs(1:3,1:node_obs)
geo3d%point%lc(1:node)                    = 20.d0    ! 20km
geo3d%point%lc(node+1:node+node_obs)      = 0.0001d0 ! 10cm

!#[2]# line
CALL INITLINE(geo3d%line,     nline+nobs    ,5) !  5 is closed 4-point spline
geo3d%line%linetype(1:nline)            = 1 ! normal line composed of 2 nodes
geo3d%line%linetype(nline+1:nline+nobs) = 2 ! spline composed of 5 nodes
geo3d%line%npoint(1:nline)              = 2
geo3d%line%npoint(nline+1:nline+nobs)   = 5
geo3d%line%node_id(1:2,1:nline)            = lines(1:2,1:nline)! spline from 5 nodes
geo3d%line%node_id(1:5,nline+1:nline+nobs) = splines_obs(1:5,1:nobs) ! spline composed of 5 nodes

!#[3]# lineloop
CALL INITLINELOOP(geo3d%lineloop, nlineloop+nobs,4) ! 4 is max # of lines composisng lineloop
geo3d%lineloop%nline(1:nlineloop) = 4 ! normal lineloop composed of 4 lines
geo3d%lineloop%nline(nlineloop+1:nlineloop+nobs) = 1 ! lineloop composed of 1 spline
geo3d%lineloop%line_id(1:4,1:nlineloop)          = lineloop(1:4,1:nlineloop)
geo3d%lineloop%line_id(1,nlineloop+1:nlineloop+nobs) = (/(i,i=nline+1,nline+nobs)/)

!#[4]# Plane
CALL INITPLANE(geo3d%plane, nlineloop+nobs,nobs+1) ! nobs + 1 is planes consisting surface
geo3d%plane%nlineloop(1:nlineloop+nobs)    = 1 ! normal plane composed of 1 lineloop
geo3d%plane%nlineloop(6)              = nobs+1 ! only the ground surface composed of nobs +1
geo3d%plane%lineloop_id(1,1:nlineloop+nobs)   = (/(i,i=1,nlineloop+nobs)/)
geo3d%plane%lineloop_id(1:nobs+1,6) = (/6,(nlineloop+i,i=1,nobs)/) ! 6th plane is ground

!#[5]# Planeloop
CALL INITPLANELOOP(geo3d%planeloop,nsurfaceloop,  6+nobs)
geo3d%planeloop%nplane(1:nsurfaceloop)   = 6 ! normal planeloop composed of 1 lineloop
geo3d%planeloop%nplane(1:2)              = 6+nobs ! air and the first layer
geo3d%planeloop%plane_id(1:6,1:nsurfaceloop)          = surfaceloop(1:6,1:nsurfaceloop)
geo3d%planeloop%plane_id(7:6+nobs,1) = (/(i,i=nlineloop+1,nlineloop+nobs)/) ! air
geo3d%planeloop%plane_id(7:6+nobs,2) = (/(i,i=nlineloop+1,nlineloop+nobs)/) ! first layer

!#[6]# volume
CALL INITVOLUME(geo3d%volume,   nsurfaceloop,       1)
geo3d%volume%nplaneloop(1:nsurfaceloop)   = 1 ! normal volume composed of 1 surfaceloop
geo3d%volume%planeloop_id(1,1:nsurfaceloop) = (/(i,i=1,nsurfaceloop)/)

write(*,*) "nlineloop=",nlineloop
write(*,*) "nobs=",nobs
write(*,*) "size=",size(geo3d%plane%lineloop_id(1,:))


return
end

!######################################### GENLINES_OBS
subroutine GENLINES_OBS(nobs,g_param,splines_obs,node)
use param
implicit none
integer(4),         intent(in)    :: node, nobs ! 2017.09.08
type(param_forward),intent(in)    :: g_param
integer(4),         intent(inout) :: splines_obs(5,nobs) ! 207.09.08
integer(4)                        :: i,ii

do i=1,nobs
 ii=(i-1)*4+node
 splines_obs(1:5,i)=(/ii+1,ii+2,ii+3,ii+4,ii+1/)
end do

write(*,*) "### GENLINES_OBS END!! ###"
return
end
!######################################### GENPOINTS_OBS
subroutine GENPOINTS_OBS(g_param,xyz_obs,node_obs,radius)
use param
implicit none
integer(4),         intent(in)    :: node_obs
type(param_forward),intent(in)    :: g_param
real(8),            intent(inout) :: xyz_obs(3,node_obs)
real(8),            intent(in)    :: radius
real(8) :: x0,y0,z0
integer(4) :: i,j

do i=1,g_param%nobs
x0=g_param%xyzobs(1,i)
y0=g_param%xyzobs(2,i)
z0=g_param%xyzobs(3,i)
 xyz_obs(1:3,(i-1)*4+1)=(/x0+radius,y0,       z0/) ! x component
 xyz_obs(1:3,(i-1)*4+2)=(/x0,       y0+radius,z0/) ! x component
 xyz_obs(1:3,(i-1)*4+3)=(/x0-radius,y0,       z0/) ! x component
 xyz_obs(1:3,(i-1)*4+4)=(/x0,       y0-radius,z0/) ! x component
!do j=1,4
!write(*,*) "xyz_obs(1:3,(i-1)*4+j)=",i,j,xyz_obs(1:3,(i-1)*4+j)
!end do
end do

write(*,*) "### GENPOINTS_OBS END!! ###"
return
end
!
!#########################################  GENPOINTS
subroutine GENPOINTS(xyz, height,xbound,ybound, k, z, node)
implicit none
integer(4),intent(in)  :: node, k
real(8),   intent(in)  :: xbound(4),ybound(4)
real(8),   intent(in)  ::  z(k), height
real(8),   intent(out) :: xyz(3,node)
integer(4) :: i
do i=-1,k ! k+2 is the # of node layers, when the k layers exist beneath the ground surface
if (i .eq. -1 ) xyz(3, 1:4 )=height
if (i .eq. 0 )  xyz(3, 5:8 )=0.d0
if (i .ge. 1 )  xyz(3, 8+(i-1)*4+1: 8+(i-1)*4+4)=z(i)
xyz(1:2,4*(i+1)+1)=(/ xbound(1), ybound(4)/) ! left rear
xyz(1:2,4*(i+1)+2)=(/ xbound(1), ybound(1) /) ! left front
xyz(1:2,4*(i+1)+3)=(/ xbound(4), ybound(1)/)   ! right fromt
xyz(1:2,4*(i+1)+4)=(/ xbound(4), ybound(4)/)      ! right rear
end do
write(*,*) "### GENPOINTS END!! ###"
return
end
!####################################### GENLINES
subroutine GENLINES(lines, nline, k)
implicit none
integer(4), intent(in)  :: k, nline ! k is # of layers beneath the ground surface
integer(4), intent(out) :: lines(2,nline)
integer(4)              :: nshift, lshift, i, j

!#[1-1]# initial 4 horizontal lines
  lines(1:2, 1 )=(/1, 2/)
  lines(1:2, 2 )=(/2, 3/)
  lines(1:2, 3 )=(/3, 4/)
  lines(1:2, 4 )=(/4, 1/)
!#[2]## i-th lines, where i=0 means vertical and horizontal lines are above and on the ground surface, respectively
  do i=0, k
!#[2-1]# i-th 4 vertical lines
    lshift=8*i  + 4
    nshift=4*i + 0
    do j=1,4
      lines(1:2, lshift+j )=(/nshift+j, nshift+j+4/)
    end do
!#[2-2]# i-th 4 horizontal lines
    lines(1:2, lshift+5 )=(/nshift+5, nshift+6/)
    lines(1:2, lshift+6 )=(/nshift+6, nshift+7/)
    lines(1:2, lshift+7 )=(/nshift+7, nshift+8/)
    lines(1:2, lshift+8 )=(/nshift+8, nshift+5/)
  end do
write(*,*) "### GENLINES END!! ###"
return
end
!######################################## GENLINELOOP
subroutine GENLINELOOP(lines, nline, lineloop, nlineloop, k)
implicit none
integer(4),intent(in)  :: lines(2, nline), nline, k, nlineloop
integer(4),intent(out) :: lineloop(4, nlineloop)
integer(4)             :: i, lshift, lineloopshift
!#   ^ y       1   4            5   8
!#   |-> x     2   3    (top    6   7  (surface
!#---------------------------
!#[1]# top horizontal loop
lineloop(1:4,1)=(/1,2,3,4/)    !       top loop
!#[2]# i-th layer loop
do i=0,k
!#[2-1] 4 vertical lineloop
  lineloopshift=5*i + 1
  lshift=8*i
  lineloop(1:4,lineloopshift+1)=(/lshift+5, lshift+9, -(lshift+6), -(lshift+1) /)
  lineloop(1:4,lineloopshift+2)=(/lshift+6, lshift+10, -(lshift+7), -(lshift+2) /)
  lineloop(1:4,lineloopshift+3)=(/lshift+7, lshift+11, -(lshift+8), -(lshift+3) /)
  lineloop(1:4,lineloopshift+4)=(/lshift+8, lshift+12, -(lshift+5), -(lshift+4) /)
!#[2-2] 1 horizontal lineloop
  lineloop(1:4,lineloopshift+5)=(/lshift+9, lshift+10, lshift+11, lshift+12 /)
end do
write(*,*) "### GENLINELOOP END!! ###"
return
end
!######################################## GENSURFACELOOP
subroutine GENSURFACELOOP(surfaceloop, nsurfaceloop, k)
implicit none
integer(4),intent(in)  :: nsurfaceloop, k
integer(4),intent(out) :: surfaceloop(6,nsurfaceloop)
integer(4)             :: i, llshift
do i=0,k
  llshift=5*i ! lineloopshift
  surfaceloop(1:6, i+1)=(/llshift+1,llshift+2,llshift+3,llshift+4,llshift+5,llshift+6/)
end do
write(*,*) "### GENSURFACELOOP END!! ###"
return
end
!######################################## GENSURFACELOOP
subroutine GENVOLUME(volume, nvolume, k)
implicit none
integer(4),intent(in) :: nvolume,k
integer(4),intent(out) :: volume(nvolume)
integer(4) :: i
do i=0,k
 volume(i+1)=i+1
end do
write(*,*) "### GENVOLUME END!! ###"
return
end
!######################################## OUTGEOFILE_OBS
! add sphere for obs positions
subroutine OUTGEOFILE_obs(outgeo, xyz,node,lines,nline,lineloop,nlineloop,surfaceloop,&
                        & nsurfaceloop,volume,nvolume,g_param)
use param
implicit none
type(param_forward),intent(in) :: g_param
integer(4),         intent(in) :: node, nline, nlineloop, nsurfaceloop, nvolume
integer(4),         intent(in) :: lines(2,nline), lineloop(4,nlineloop)
integer(4),         intent(in) :: surfaceloop(6,nsurfaceloop),volume(nvolume)
real(8),            intent(in) :: xyz(3,node)
character(50),      intent(in) :: outgeo
real(8)                        :: lc2, x,y,z,lc,lc0,r
integer(4)                     :: i,j,k,ishif,lishif,lloop,lsurface,lsurfaceloop,lvolume

!#
lc2=30.d0
open(1,file=outgeo)
!write(1,*) "lc1=10.0;"
!#[1]## Points
write(1,*) "lc2=",lc2,";"
do i=1,node
  write(1,*) "Point(",i,")={",xyz(1,i),",",xyz(2,i),",",xyz(3,i),",", lc2,"};"
end do
!#[2]## Lines
do i=1, nline
write(1,*) "Line(",i,")={",lines(1,i),",",lines(2,i),"};"
end do
!#[3]## Lineloops and Plane Surface
do i=1,nlineloop
  write(1,*) "Line Loop(",i,")={",lineloop(1,i),",",lineloop(2,i),",",lineloop(3,i),",",lineloop(4,i),"};"
end do
 ishif   =node ! node id shift
 lishif  =nline ! line id shift
 lloop   =nlineloop ! line loop id shift
 lsurface=nlineloop ! surface id shift
 lsurfaceloop=nsurfaceloop ! surfaceloop id
 lvolume=nvolume    ! volume id
do i=1,g_param%nobs
 x=g_param%xyzobs(1,i)
 y=g_param%xyzobs(2,i)
 z=g_param%xyzobs(3,i)
! lc=0.001
! r=0.005
 lc0=0.00020
 lc=0.002
 r=0.005
! lc=0.004*dsqrt(x**2.d0+y**2.d0+z**2.d0)
! r =0.010*dsqrt(x**2.d0+y**2.d0+z**2.d0) ! 0 - 1 [km]
 call OUTSPHERE(1,x,y,z,lc0,lc,r,ishif,lishif,lloop,lsurface,lsurfaceloop,lvolume)
end do
do i=1,nlineloop
 if ( i .ne. 6) write(1,*) "Plane Surface(",i,")={",i,"};"
 if ( i .eq. 6) then
  write(1,*) "Plane Surface(6)={ 6,",(nlineloop+9*j,",",j=1,g_param%nobs-1),nlineloop+9*g_param%nobs,"};"
  ! 9th line loop is horizontal circle lineloop for each sphere
 end if
end do
!#[4]## Surfaceloops
do i=1,nsurfaceloop
  if (i .ge. 3) write(1,*) "Surface loop(",i,")={",(surfaceloop(j,i),",",j=1,5),surfaceloop(6,i),"};"
  if (i .le. 2) then
   write(1,*) "Surface loop(",i,")={",(surfaceloop(j,i),",",j=1,6),&
   &((nlineloop+(j-1)*9+(i-1)*4+k,",",k=1,4),j=1,g_param%nobs-1),&
   & nlineloop+9*(g_param%nobs-1)+(i-1)*4 +1,",",&
   & nlineloop+9*(g_param%nobs-1)+(i-1)*4 +2,",",&
   & nlineloop+9*(g_param%nobs-1)+(i-1)*4 +3,",",&
   & nlineloop+9*(g_param%nobs-1)+(i-1)*4 +4,"};"
  end if
end do
!#[5]## Volumes
do i=1,nvolume
  write(1,*) "Volume(",i,")={",volume(i),"};"
end do
do i=1,nvolume
if ( i .le. 2) then
 write(1,*)"Physical Volume(",i,")={",i,(",",nvolume+(j-1)*2+i,j=1,g_param%nobs),"};"
else
write(1,*) "Physical Volume(",i,")={",i,"};"
end if
end do
close(1)
write(*,*) "### OUTGEOFILE_obs END!! ###"
return
end
!

!######################################## OUTGEOFILE
subroutine OUTGEOFILE(outgeo, xyz,node,lines,nline,lineloop,nlineloop,surfaceloop,nsurfaceloop,volume,nvolume,g_param)
use param
implicit none
type(param_forward),intent(in) :: g_param
integer(4),         intent(in) :: node, nline, nlineloop, nsurfaceloop, nvolume
integer(4),         intent(in) :: lines(2,nline), lineloop(4,nlineloop)
integer(4),         intent(in) :: surfaceloop(6,nsurfaceloop), volume(nvolume)
real(8),            intent(in) :: xyz(3,node)
character(50),      intent(in) :: outgeo
real(8)                        :: lc2
integer(4)                     :: i,j
!#
lc2=g_param%sizebo
open(1,file=outgeo)
!write(1,*) "lc1=10.0;"
!#[1]## Points
write(1,*) "lc2=",lc2,";"
do i=1,node
  write(1,*) "Point(",i,")={",xyz(1,i),",",xyz(2,i),",",xyz(3,i),",", lc2,"};"
end do
!#[2]## Lines
do i=1, nline
write(1,*) "Line(",i,")={",lines(1,i),",",lines(2,i),"};"
end do
!#[3]## Lineloops and Plane Surface
do i=1,nlineloop
  write(1,*) "Line Loop(",i,")={",lineloop(1,i),",",lineloop(2,i),",",lineloop(3,i),",",lineloop(4,i),"};"
  write(1,*) "Plane Surface(",i,")={",i,"};"
end do
!#[4]## Surfaceloops
do i=1,nsurfaceloop
  write(1,*) "Surface loop(",i,")={",(surfaceloop(j,i),",",j=1,5),surfaceloop(6,i),"};"
end do
!#[5]## Volumes
do i=1,nvolume
  write(1,*) "Volume(",i,")={",volume(i),"};"
  write(1,*) "Physical Volume(",i,")={",volume(i),"};"
end do

write(1,*) " " ! Inserted on 2016.10.11
100 continue
!# close
close(1)
write(*,*) "### OUTGEOFILE END!! ###"
return
end
!
!########################################### OUTPOSFILE2_OBS
subroutine OUTPOSFILE2_OBS(outpos,g_param)
use param
implicit none
type(param_forward),intent(in) :: g_param
character(50),      intent(in) :: outpos
real(8)                        :: lc1,lc2,lc3, dlen,dlen2,shift
real(8)                        :: height1, height2, height3, height4,xout,yout
integer(4)                     :: ifile=2
real(8)                        :: x0,y0,z0,x,y,z,dd,dx(-3:3),s(4)
integer(4)                     :: nx,ny,nz
integer(4),parameter           :: nxmax=50,nymax=400,nzmax=50
real(4)                        :: xgrd0(nxmax+1),ygrd0(nymax+1),zgrd0(nzmax+1)
real(8)                        :: xgrd(nxmax+1),ygrd(nymax+1),zgrd(nzmax+1)
integer(4)                     :: i,j,k,kk,jj
real(8)                        :: s1,s2,s3,s4,s5,s6,s7,s8,value_OBS,v(8)

!#[1]## set
xout = g_param%xbound(4)
yout = g_param%ybound(4)
dlen=0.1d0  ! 0.001
height1=g_param%zbound(4)
dlen2=min(xout,yout,height1)
height2=0.01d0  ! 0.01
height3=-0.01d0 ! 0.01
height4=g_param%zbound(1)
lc1=0.05d0 ! 0.0005
lc2=0.001d0
!lc3=20.d0  ! 20.d0
lc3=g_param%sizebo  ! 20.d0
shift=20.d0


!####### Initial grid ######
xgrd0(1:19)=(/-real(xout),-50.,-10.,-5.,-1.,-0.1,-0.05,-0.01,-0.003,&
          & 0.,0.003,0.01,0.05,0.1, 1.,5.,10.,50.,real(xout)/)
ygrd0(1:19)=(/-real(yout),-50.,-10.,-5.,-1.,-0.1,-0.05,-0.01,-0.003,&
          & 0.,0.003,0.01,0.05, 0.1,1.,5.,10.,50.,real(yout)/)
zgrd0(1:19)=(/real(height4),-50.,-10.,-5.,-1.,-0.1,-0.05,-0.01,-0.003,&
          & 0.,0.003,0.01,0.05, 0.1,1.,5.,10.,50.,real(height1)/)
do i=1,19
 xgrd(i) = dble(xgrd0(i))
 ygrd(i) = dble(ygrd0(i))
 zgrd(i) = dble(zgrd0(i))
end do

nx=19
ny=19
nz=19
!####### Initial grid ######
dd= 0.001d0 ! 1m for threshold
!dx(-3:3)= (/-0.02,-0.012,-0.005,0.,0.005,0.012,0.02/) ! 2 m
!dx(-3:3)= (/-0.025,-0.015,-0.005,0.,0.005,0.015,0.025/) ! 2 m
dx(-3:3)= (/-0.025,-0.025,-0.010,0.,0.010,0.025,0.025/) ! 2 m
jj=2

!#[1]## modify the grid
do i=1,g_param%nobs
x0=g_param%xyzobs(1,i)
y0=g_param%xyzobs(2,i)
z0=g_param%xyzobs(3,i)
write(*,*) "x0=",x0,"y0=",y0,"z0=",z0
!#[x]
do j=-jj,jj
 x = x0 + dx(j)
 call updategrd(xgrd,nxmax,nx,x,dd)
end do
!#[y]
do j=-jj,jj
 y = y0 + dx(j)
 call updategrd(ygrd,nymax,ny,y,dd)
end do
!#[z]
do j=-jj,jj
 z = z0 + dx(j)
 call updategrd(zgrd,nzmax,nz,z,dd)
end do
!
end do
do i=1,nx
write(*,*) "xgrd (",i,")=",xgrd(i)
end do
do i=1,ny
write(*,*) "ygrd (",i,")=",ygrd(i)
end do
do i=1,nz
write(*,*) "zgrd (",i,")=",zgrd(i)
end do

!#[2]## calculate values and output .pos file
open(ifile,file=outpos)
write(ifile,*) 'View "backgraound mesh"{'
do k=1,nz-1
 do j=1,ny-1
  do i=1,nx-1                                                 !     z   y       x
  s1=value_OBS(xgrd(i)  ,ygrd(j+1),zgrd(k+1),lc1,lc3,dlen,dlen2,shift,g_param) ! 1: up   top     left
  s2=value_OBS(xgrd(i)  ,ygrd(j)  ,zgrd(k+1),lc1,lc3,dlen,dlen2,shift,g_param) ! 2: up   bottom  left
  s3=value_OBS(xgrd(i+1),ygrd(j)  ,zgrd(k+1),lc1,lc3,dlen,dlen2,shift,g_param) ! 3: up   bottom  right
  s4=value_OBS(xgrd(i+1),ygrd(j+1),zgrd(k+1),lc1,lc3,dlen,dlen2,shift,g_param) ! 4: up   top     right
  s5=value_OBS(xgrd(i)  ,ygrd(j+1),zgrd(k)  ,lc1,lc3,dlen,dlen2,shift,g_param) ! 5: down top     left
  s6=value_OBS(xgrd(i)  ,ygrd(j)  ,zgrd(k)  ,lc1,lc3,dlen,dlen2,shift,g_param) ! 6: down bottom  left
  s7=value_OBS(xgrd(i+1),ygrd(j)  ,zgrd(k)  ,lc1,lc3,dlen,dlen2,shift,g_param) ! 7: down bottom  right
  s8=value_OBS(xgrd(i+1),ygrd(j+1),zgrd(k)  ,lc1,lc3,dlen,dlen2,shift,g_param) ! 8: down top     right
  v(1:8)=(/s1,s2,s3,s4,s5,s6,s7,s8/)
  call SHWRITE(ifile,xgrd(i),xgrd(i+1),ygrd(j),ygrd(j+1),zgrd(k+1), zgrd(k), v)
  end do
 end do
end do
write(2,*) "};"
close(2)
write(*,*) "### OUTPOSFILE END!! ###"

return
end
!#######################################################    value
function value_OBS(x,y,z,lc1,lc3,dlen,dlen2,shift,g_param)
use param
implicit none
type(param_forward),intent(in) :: g_param
real(8),            intent(in) :: x,y,z,shift ! [km]
real(8),            intent(in) :: lc1,lc3,dlen,dlen2
real(8)                        :: value_OBS, value_OBS_def
real(8)                        :: r,k,A,B,robs,robsmin,power
integer(4) :: i
r=dsqrt(x**2. + y**2. + z**2.)
!## for origin area
if ( r .le. dlen) then ! dlen = 0.1
  value_OBS=lc1 ! 0.03
 if ( r .lt. 0.040 ) value_OBS =0.010d0  ! 8 m
 if ( r .lt. 0.030 ) value_OBS =0.008d0  ! 7 m
 if ( r .lt. 0.020 ) value_OBS =0.006d0  ! 5 m
 if ( r .lt. 0.010 ) value_OBS =0.0040d0 ! 3 m
 if ( r .lt. 0.004 ) value_OBS =0.0030d0 ! 3 m
 if ( r .lt. 0.003 ) value_OBS =0.0010d0 ! 1 m
 if ( r .lt. 0.001 ) value_OBS =0.0005d0 ! 10 cm
end if
if ( r .ge. dlen .and. r .lt. dlen2 ) then
!# [1] A*exp(k*r)
! k=dlog(lc3/lc1)/(dlen2-dlen)
! A=lc1*dexp(-k*dlen)
! value=A*exp(k*r)./complot.sh
!# [2] A*r^3+B
!  shift=10.d0
!  power=3.d0  ! power should be around 2.5 ~ 2.7 when solved in local PC on May 15, 2016
  power=2.7d0
  A=(lc3-lc1)/((dlen2+shift)**power - (dlen+shift)**power)
  B=lc1-A*((dlen+shift)**power)
  value_OBS=A*((r+shift)**power)+B
!  if ( dabs(x) .lt. 1. .and. dabs(y) .lt. 1. ) then
!   if (dabs(z-(-0.2)) .lt. 0.051 ) value_OBS=0.5*value_OBS
!   if (dabs(z-(-0.2)) .lt. 0.02 ) value_OBS=0.02d0
!  end if
end if
!## for each observatory
value_OBS_def=value_OBS
do i=1,g_param%nobs
robs=dsqrt((g_param%xyzobs(1,i)-x)**2.d0 + (g_param%xyzobs(2,i)-y)**2.d0 + (g_param%xyzobs(3,i)-z)**2.d0)
 robsmin=0.001d0
 if ( robs .lt. 0.040 ) value_OBS =value_OBS_def*0.2 + robsmin*0.8 ! 7m
 if ( robs .lt. 0.030 ) value_OBS =value_OBS_def*0.2 + robsmin*0.8  ! 7m
 if ( robs .lt. 0.020 ) value_OBS =value_OBS_def*0.1 + robsmin*0.9  ! 5m
 if ( robs .lt. 0.010 ) value_OBS =value_OBS_def*0.05 + robsmin*0.95 ! 2 cm
 if ( robs .lt. 0.004 ) value_OBS =value_OBS_def*0.02 + robsmin*0.98 ! 2 cm
 if ( robs .lt. 0.003 ) value_OBS =value_OBS_def*0.01 + robsmin*0.99 ! 1 m
 if ( robs .lt. 0.001 ) value_OBS =value_OBS_def*0.0 + robsmin*1.0  ! 1 m
!
! if ( robs .lt. 0.040 ) value_OBS =0.005 ! 7m
! if ( robs .lt. 0.030 ) value_OBS =0.005 ! 7m
! if ( robs .lt. 0.020 ) value_OBS =0.004 ! 5m
! if ( robs .lt. 0.010 ) value_OBS =0.003 ! 2 cm
! if ( robs .lt. 0.004 ) value_OBS =0.003 ! 2 cm
! if ( robs .lt. 0.003 ) value_OBS =0.002 ! 1 m
! if ( robs .lt. 0.001 ) value_OBS =0.001 ! 10 cm
end do
if ( r .ge. dlen2 ) value_OBS=lc3
return
end

!########################################### upgradegrd
subroutine updategrd(xgrd,nxmax,nx,x,dd)
implicit none
integer(4), intent(in) :: nxmax
real(8),    intent(inout) :: xgrd(nxmax)
integer(4),intent(inout) :: nx
real(8),intent(in) :: dd,x
integer(4) :: k,kk
kk=0
!write(*,*) "input y=",x
 do k=1,nx
  if ( xgrd(k) .lt. x ) kk=k ! remember k
  if ( dabs(xgrd(k) - x) .lt. dd ) goto 10 ! do nothing
 end do
 xgrd(kk+2:nx+1)=xgrd(kk+1:nx) ! shift
 xgrd(kk+1)=x ! new x
 nx=nx+1
! write(*,*) "new y =",xgrd(kk+1)
 if (nxmax .lt. nx) then
  write(*,*) "GEGEGE nx is greater than nx! nxmax=",nxmax,"nx=",nx
  stop
 end if
10 continue
return
end
!
!######################################################### horizontal9
!# tetrahedron size, s
!# s1   s2
!#  ---------------  z2(1) ^
!#  |   |    |    |        |  z
!#  ---------------  z2(2)
!# s3 s4
!#------  v ------------------
!#    1 < 4
!#    v   ^   z=z2(1)
!#    2 > 3
!#  &
!#    5 < 8
!#    v   ^   z=z2(2)
!#    6 > 7
subroutine horizontal9(ifile,x2,y2,xc2,yc2,z2,s)
implicit none
integer(4), intent(in) :: ifile
real(8),    intent(in) :: x2(2),y2(2),xc2(2),yc2(2),z2(2)
real(8),    intent(in) :: s(4)
real(8)                :: v(8)
real(8)                :: s1,s2,s3,s4

s1=s(1)
s2=s(2)
s3=s(3)
s4=s(4)
!#[1]# west
v(1:8)=(/s1,s1,s2,s1,s3,s3,s4,s3/)         ! north west
call SHWRITE(ifile, x2(1), xc2(1), yc2(2), y2(2), z2(1), z2(2), v)
v(1:8)=(/s1,s1,s2,s2,s3,s3,s4,s4/)         ! west
call SHWRITE(ifile, x2(1), xc2(1), yc2(1), yc2(2), z2(1), z2(2), v)
v(1:8)=(/s1,s1,s1,s2,s3,s3,s3,s4/)         ! south west
call SHWRITE(ifile, x2(1), xc2(1), y2(1),  yc2(1), z2(1), z2(2), v)

!#[2]# center
v(1:8)=(/s1,s2,s2,s1,s3,s4,s4,s3/)
call SHWRITE(ifile, xc2(1), xc2(2), yc2(2), y2(2), z2(1), z2(2), v)  ! north center
v(1:8)=(/s2,s2,s2,s2,s4,s4,s4,s4/)
call SHWRITE(ifile, xc2(1), xc2(2), yc2(1), yc2(2), z2(1), z2(2), v) !  center
v(1:8)=(/s2,s1,s1,s2,s4,s3,s3,s4/)
call SHWRITE(ifile, xc2(1), xc2(2), y2(1), yc2(1), z2(1), z2(2), v) !  south center

!#[3]# east
v(1:8)=(/s1,s2,s1,s1,s3,s4,s3,s3/)
call SHWRITE(ifile, xc2(2), x2(2), yc2(2), y2(2), z2(1), z2(2), v)  ! north center
v(1:8)=(/s2,s2,s1,s1,s4,s4,s3,s3/)
call SHWRITE(ifile, xc2(2), x2(2), yc2(1), yc2(2),z2(1), z2(2), v) !  center
v(1:8)=(/s2,s1,s1,s1,s4,s3,s3,s3/)
call SHWRITE(ifile, xc2(2), x2(2), y2(1), yc2(1), z2(1), z2(2), v) !  south center

return
end

!######################################################## OUTSPHERE
subroutine OUTSPHERE(ifile,x,y,z,lc0,lc,r,ishif,lishif,lloop,lsurface,&
                   & lsurfaceloop,lvolume)
implicit none
integer(4),intent(in)    :: ifile
real(8),   intent(in)    :: x,y,z,lc0,lc,r
integer(4),intent(inout) :: ishif,lishif,lloop,lsurface,lsurfaceloop,lvolume

write(ifile,*)"x=",x,";"
write(ifile,*)"y=",y,";"
write(ifile,*)"z=",z,";"
write(ifile,*)"lc0=",lc0,";"
write(ifile,*)"lc=",lc,";"
write(ifile,*)"r=",r,";"
write(ifile,*)"ishif=",ishif,";"
write(ifile,*)"lishif=",lishif,";"
write(ifile,*)"lloop=",lloop,";"
write(ifile,*)"lsurface=",lsurface,";"
write(ifile,*)"lsurfaceloop=",lsurface,";"
write(ifile,*)"lvolume=",lvolume,";"
write(ifile,*)"Point(ishif+1) = {x, y, z,   lc0}; // origin"
write(ifile,*)"Point(ishif+2) = {x, y, z+r, lc}; // top"
write(ifile,*)"Point(ishif+3) = {x+r, y, z, lc}; // right"
write(ifile,*)"Point(ishif+4) = {x, y+r, z, lc}; // up"
write(ifile,*)"Point(ishif+5) = {x-r, y, z, lc}; // left"
write(ifile,*)"Point(ishif+6) = {x, y-r, z, lc}; // down"
write(ifile,*)"Point(ishif+7) = {x, y, z-r, lc}; // bottom"
write(ifile,*)"Circle(lishif + 1)={ishif +3,ishif +1,ishif +2};"
write(ifile,*)"Circle(lishif + 2)={ishif +4,ishif +1,ishif +2};"
write(ifile,*)"Circle(lishif + 3)={ishif +5,ishif +1,ishif +2};"
write(ifile,*)"Circle(lishif + 4)={ishif +6,ishif +1,ishif +2};"
write(ifile,*)"Circle(lishif + 5)={ishif +3,ishif +1,ishif +4};"
write(ifile,*)"Circle(lishif + 6)={ishif +4,ishif +1,ishif +5};"
write(ifile,*)"Circle(lishif + 7)={ishif +5,ishif +1,ishif +6};"
write(ifile,*)"Circle(lishif + 8)={ishif +6,ishif +1,ishif +3};"
write(ifile,*)"Circle(lishif + 9) ={ishif +3,ishif +1,ishif +7};"
write(ifile,*)"Circle(lishif + 10)={ishif +4,ishif +1,ishif +7};"
write(ifile,*)"Circle(lishif + 11)={ishif +5,ishif +1,ishif +7};"
write(ifile,*)"Circle(lishif + 12)={ishif +6,ishif +1,ishif +7};"
write(ifile,*)"Line loop(lloop +1)={-(lishif+1),lishif+ 5,lishif+ 2};"
write(ifile,*)"Line loop(lloop +2)={-(lishif+ 2),lishif+ 6,lishif+ 3};"
write(ifile,*)"Line loop(lloop +3)={-(lishif+ 3),lishif+ 7,lishif+ 4};"
write(ifile,*)"Line loop(lloop +4)={-(lishif+ 4),lishif+ 8,lishif+ 1};"
write(ifile,*)"Line loop(lloop +5)={-(lishif+ 9),lishif+ 5,lishif+ 10};"
write(ifile,*)"Line loop(lloop +6)={-(lishif+ 10),lishif+ 6,lishif+ 11};"
write(ifile,*)"Line loop(lloop +7)={-(lishif+ 11),lishif+ 7,lishif+ 12};"
write(ifile,*)"Line loop(lloop +8)={-(lishif+ 12),lishif+ 8,lishif+ 9};"
write(ifile,*)"Line loop(lloop +9)={lishif+5,lishif+6,lishif+7,lishif+8};"
write(ifile,*)"Ruled Surface(lsurface + 1)={lloop+1} In Sphere {ishif + 1};"
write(ifile,*)"Ruled Surface(lsurface + 2)={lloop+2} In Sphere {ishif + 1};"
write(ifile,*)"Ruled Surface(lsurface + 3)={lloop+3} In Sphere {ishif + 1};"
write(ifile,*)"Ruled Surface(lsurface + 4)={lloop+4} In Sphere {ishif + 1};"
write(ifile,*)"Ruled Surface(lsurface + 5)={lloop+5} In Sphere {ishif + 1};"
write(ifile,*)"Ruled Surface(lsurface + 6)={lloop+6} In Sphere {ishif + 1};"
write(ifile,*)"Ruled Surface(lsurface + 7)={lloop+7} In Sphere {ishif + 1};"
write(ifile,*)"Ruled Surface(lsurface + 8)={lloop+8} In Sphere {ishif + 1};"
write(ifile,*)"Plane Surface(lsurface + 9)={lloop+9};"
write(ifile,*)"Surface loop(lsurfaceloop + 1)={lsurface+1,lsurface+2,lsurface+3,lsurface+4,lsurface+9};"
write(ifile,*)"Surface loop(lsurfaceloop + 2)={lsurface+5,lsurface+6,lsurface+7,lsurface+8,lsurface+9};"
write(ifile,*)"Volume(lvolume + 1)={lsurfaceloop + 1};"
write(ifile,*)"Volume(lvolume + 2)={lsurfaceloop + 2};"
ishif=ishif+7
lishif=lishif+12
lloop=lloop+9
lsurface=lsurface+9
lsurfaceloop=lsurfaceloop+2
lvolume=lvolume+2
return
end

