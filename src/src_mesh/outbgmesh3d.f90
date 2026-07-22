!# Coded on 2016.09.22
!# Copied from tsunami/3D/ana_comp/src/mesh/outbgmesh3d.f90
subroutine outbgmesh3d(posfile,g_meshpara)
use horizontalresolution
use param
implicit none

type(param_forward),intent(in)     :: g_meshpara
character(50),      intent(in)     :: posfile
integer(4)                         :: ifile
real(8)                            :: si,sb
integer(4)                         :: nx,ny,nz,i,j,k,jj
integer(4),         parameter      :: nxmax=500,nymax=500,nzmax=500
real(8)                            :: xgrd(nxmax+1),ygrd(nymax+1),zgrd(nzmax+1)
real(8)                            :: dd,dx(-11:11)
real(8)                            :: v1,v2,v3,v4,v5,v6,v7,v8,x0,y0,z0
real(8)                            :: xx,yy,zz,r2,sigma,v(8)
integer(4)                         :: nobsr
real(8),allocatable,dimension(:,:) :: xyz_r        ! 2017.09.08
real(8),allocatable,dimension(:)   :: sigma_r, A_r ! 2017.09.08

!#[0]## set input
xgrd(1:4) = g_meshpara%xbound(1:4)
ygrd(1:4) = g_meshpara%ybound(1:4)
xgrd(1)   = xgrd(1) -20.0
xgrd(4)   = xgrd(4) +20.0
ygrd(1)   = ygrd(1) -20.0
ygrd(4)   = ygrd(4) +20.0
zgrd(4)   = g_meshpara%zmax + 20.0
zgrd(3)   = g_meshpara%upzin
zgrd(2)   = g_meshpara%downzin
zgrd(1)   = g_meshpara%zmin - 20.0
nobsr     = g_meshpara%nobsr
allocate(xyz_r(3,nobsr),sigma_r(nobsr),A_r(nobsr))
sigma_r   = g_meshpara%sigma_r ! [km]
xyz_r     = g_meshpara%xyz_r
A_r       = g_meshpara%A_r ! [km]
!sigma     = g_meshpara%sigma_obs ! [km]
!nobs      = g_meshpara%nobs
!xyzobs    = g_meshpara%xyzobs
! added on 2016.10.14

!####### Initial grid ######
r2=dsqrt(2.d0)
ifile=1
nx=4
ny=4
nz=4
!#
!dd= 0.001d0 ! 1m for threshold for distance between two grid lines
dd= 0.040d0 ! 1m for threshold for distance between two grid lines
jj=8

!#[1]## modify the grid
!goto 100
do i=1,nobsr
!dx(-11:11)= (/-3.0,-2.5,-2.0,-1.5,-1.0,-0.75,-0.5,-0.35,-0.20,-0.10,-0.05,0.0,&
!     &  0.05,  0.10, 0.20,0.35,0.50,0.75,1.0,1.5,2.0,2.5,3.0/) * sigma_r(i)
dx(-8:8)= (/-2.5,-2.0,-1.5,-1.0,-0.75,-0.5,-0.25,-0.10,0.0,&
     &    0.10, 0.25,0.50,0.75,1.0,1.5,2.0,2.5/) * sigma_r(i)
 x0=xyz_r(1,i) ! [km]
 y0=xyz_r(2,i) ! [km]
 z0=xyz_r(3,i) ! [km]
! write(*,*) "x0=",x0,"y0=",y0,"z0=",z0
 do j=-jj,jj
 !#[x]
  xx = x0 + dx(j)
  call updategrd(xgrd,nxmax,nx,xx,dd)
 !#[y]
  yy = y0 + dx(j)
  call updategrd(ygrd,nymax,ny,yy,dd)
 !#[z]
  zz = z0 + dx(j)
  call updategrd(zgrd,nzmax,nz,zz,dd)
 end do
!
end do
!do i=1,nx
! write(*,*) "xgrd (",i,")=",xgrd(i)
!end do
do i=1,ny
 call value_3d_r(0.d0,ygrd(i),g_meshpara%zorigin,g_meshpara,v1) ! 2017.09.08
! write(*,*) "ygrd (",i,")=",ygrd(i),"v=",v1
end do
!do i=1,nz
! write(*,*) "zgrd (",i,")=",zgrd(i)
!end do


open(ifile,file=posfile)
write(ifile,*)'View "backgraound mesh"{'
do k=1,nz-1
 do j=1,ny-1
  do i=1,nx-1
   call value_3d_r(xgrd(i  ),ygrd(j+1),zgrd(k  ),g_meshpara,v1)
   call value_3d_r(xgrd(i  ),ygrd(j  ),zgrd(k  ),g_meshpara,v2)
   call value_3d_r(xgrd(i+1),ygrd(j  ),zgrd(k  ),g_meshpara,v3)
   call value_3d_r(xgrd(i+1),ygrd(j+1),zgrd(k  ),g_meshpara,v4)
  !
   call value_3d_r(xgrd(i  ),ygrd(j+1),zgrd(k+1),g_meshpara,v5)
   call value_3d_r(xgrd(i  ),ygrd(j  ),zgrd(k+1),g_meshpara,v6)
   call value_3d_r(xgrd(i+1),ygrd(j  ),zgrd(k+1),g_meshpara,v7)
   call value_3d_r(xgrd(i+1),ygrd(j+1),zgrd(k+1),g_meshpara,v8)
  v(1:8)=(/v1,v2,v3,v4,v5,v6,v7,v8/)
  call SHWRITE(ifile,xgrd(i),xgrd(i+1),ygrd(j),ygrd(j+1),zgrd(k),zgrd(k+1), v)
  end do
 end do
end do
write(ifile,*)"};"
close(ifile)
write(*,*) "### outbgmesh3d end! ###"
return
end subroutine outbgmesh3d

!######################################################### SHWRITE
subroutine SHWRITE(ifile,x1,x2,y1,y2,z1,z2,v)
! x is the horizontal axis, y is the vertical axis
implicit none
integer(4),intent(in)   :: ifile
real(8),   intent(in)   :: x1,x2,y1,y2,z1,z2
real(8),   intent(in)   :: v(8) ! 2017.09.08
real(8),   dimension(8) :: x,y,z
integer(4)              :: i,j

x(1:8)=(/x1,x1,x2,x2,x1,x1,x2,x2/)
y(1:8)=(/y2,y1,y1,y2, y2,y1,y1,y2/)
z(1:8)=(/z1,z1,z1,z1,z2,z2,z2,z2/)
write(ifile,*) "SH(",x(1),",",y(1),",",z(1)
write(ifile,*)  (",",x(j),",",y(j),",",z(j),j=2,8)
write(ifile,*) "){",v(1)
write(ifile,*)  (",",v(i),i=2,8)
write(ifile,*)"};"
return
end subroutine SHWRITE


