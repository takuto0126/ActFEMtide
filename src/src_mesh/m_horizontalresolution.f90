module horizontalresolution
implicit none

contains
!########################################## value_r
! Coded on 2016.10.14
subroutine value_r(x,y,g_meshpara,v_r)
use param
implicit none
type(param_forward),intent(in)     :: g_meshpara
real(8),            intent(in)     :: x,y
real(8),            intent(out)    :: v_r     ! 2017.09.08
real(8),allocatable,dimension(:,:) :: xyz     ! 2017.09.08
real(8),allocatable,dimension(:)   :: v       ! 2017.09.08
real(8),allocatable,dimension(:,:) :: xyz_r      ! 2017.09.08
real(8),allocatable,dimension(:)   :: sigma_r,A_r! 2017.09.08
real(8)                            :: si,sb, x1,x2,y1,y2,sigma,A
real(8)                            :: robs
integer(4)                         :: i,nobs
integer(4)                         :: nobsr

!#[0]## set input
nobs = g_meshpara%nobs
allocate(xyz(3,nobs)) ! 2017.09.08
xyz  = g_meshpara%xyzobs
x1   = g_meshpara%xbound(2)
x2   = g_meshpara%xbound(3)
y1   = g_meshpara%ybound(2)
y2   = g_meshpara%ybound(3)
si   = g_meshpara%sizein
sb   = g_meshpara%sizebo
!A = g_meshpara%A_obs ! [km]
!sigma = g_meshpara%sigma_obs ! [km]
! added on 2016.10.14
nobsr   = g_meshpara%nobsr
allocate( xyz_r(3,nobsr)                     ) ! 2017.09.08
allocate( v(nobsr),sigma_r(nobsr),A_r(nobsr) ) ! 2017.09.08
A_r     = g_meshpara%A_r ! [km]
sigma_r = g_meshpara%sigma_r ! [km]
xyz_r   = g_meshpara%xyz_r

!#[1]## value default
 if ( x .lt. x1 .or. x2 .lt. x .or. y .lt. y1 .or. y2 .lt. y) then
  v_r=sb ! 2017.09.08
 else
!#[2]## value for focus area or near observatories
  do i=1,nobsr
   robs=dsqrt((xyz_r(1,i)-x)**2.d0 + (xyz_r(2,i)-y)**2.d0 )
   v(i) = si
   if ( robs .le. 3.*sigma_r(i) ) then
    v(i) = si - (si - A_r(i))*exp(-robs**2./2./sigma_r(i)**2.)
   end if
  end do
  v_r = si            ! 2017.09.08
  do i=1,nobsr
   v_r = min(v_r,v(i))! 2017.09.08
  end do
 end if

return
end subroutine
!########################################## value_3d_r
subroutine value_3d_r(x,y,z,g_meshpara,v_3d_r)
use param
implicit none
type(param_forward),intent(in)     :: g_meshpara
real(8),            intent(in)     :: x,y,z
real(8),            intent(out)    :: v_3d_r ! 2017.09.08
real(8)                            :: si,sb, x1,x2,y1,y2,z1,z2,sigma,A
real(8)                            :: robs
real(8),allocatable,dimension(:,:) :: xyz ! 2017.09.08
integer(4)                         :: i,nobs
real(8)                            :: x02,y02
integer(4)                         :: nobsr
real(8),allocatable,dimension(:,:) :: xyz_r      ! 2017.09.08
real(8),allocatable,dimension(:)   :: sigma_r,A_r! 2017.09.08

!#[0]## set input
nobs = g_meshpara%nobs
allocate(xyz(3,nobs))
xyz  = g_meshpara%xyzobs
x1   = g_meshpara%xbound(2)
x2   = g_meshpara%xbound(3)
x02  = g_meshpara%xbound(4)
y1   = g_meshpara%ybound(2)
y2   = g_meshpara%ybound(3)
y02  = g_meshpara%ybound(4)
z2   = g_meshpara%upzin
z1   = g_meshpara%downzin
si   = g_meshpara%sizein
sb   = g_meshpara%sizebo
!A     = g_meshpara%A_obs ! [km]
!sigma = g_meshpara%sigma_obs
! added on 2016.10.14
nobsr   = g_meshpara%nobsr
allocate( xyz_r(3,nobsr)            ) ! 2017.09.08
allocate( sigma_r(nobsr),A_r(nobsr) ) ! 2017.09.08
A_r     = g_meshpara%A_r ! [km]
sigma_r = g_meshpara%sigma_r ! [km]
xyz_r   = g_meshpara%xyz_r

!#[1]## value default
 if ( x .lt. x1 .or. x2 .lt. x .or. &
   &  y .lt. y1 .or. y2 .lt. y .or. &
   &  z .lt. z1 .or. z2 .lt. z ) then
!  robs=dsqrt(x**2.d0 + y**2.d0 + z**2.d0)
  v_3d_r = sb ! 2017.09.08
!  value_3d_r = g_meshpara%sizein/(2.*g_meshpara%sigma_src + 0.25)*(robs + 0.25)
!  if (  x2 .lt. x .and. y1 .lt. y .and. y .lt. y2) value_3d_r= si + (sb-si)/(x02-x2)*(x-x2)
!  if (  x .lt. x1 .and. y1 .lt. y .and. y .lt. y2) value_3d_r= si + (sb-si)/(x02-x2)*(x1-x)
!  if (  y2 .lt. y .and. x1 .lt. x .and. x .lt. x2) value_3d_r= si + (sb-si)/(y02-y2)*(y-y2)
!  if (  y .lt. y1 .and. x1 .lt. x .and. x .lt. x2) value_3d_r= si + (sb-si)/(y02-y2)*(y1-y)
 else
!#[2]## value for focus area or near observatories
  v_3d_r=si   ! 2017.09.08
  do i=1,nobsr
   robs=dsqrt((xyz_r(1,i)-x)**2.d0 + (xyz_r(2,i)-y)**2.d0 + (xyz_r(3,i)-z)**2.d0)
   if ( robs .le. 3.*sigma_r(i) ) then
    v_3d_r = min(v_3d_r, si - (si - A_r(i))*exp(-robs**2./2./sigma_r(i)**2.) )
!     v_3d_r = min(value_3d_r, si/7.389*exp(-robs**2./2./0.33/0.33))
   else
    v_3d_r = min(v_3d_r,si) ! 2017.09.08
   end if
  end do
 end if

return
end subroutine


!########################################## value
function value(x,y,g_meshpara)
use param
implicit none
type(param_forward),intent(in) :: g_meshpara
real(8),            intent(in) :: x,y
real(8)                        :: si,sb, x1,x2,y1,y2,sigma,A
real(8)      :: value, robs
real(8),allocatable,dimension(:,:) :: xyz ! 2017.09.08
real(8),allocatable,dimension(:)   :: v   ! 2017.09.08
integer(4) :: i,nobs

!#[0]## set input
nobs = g_meshpara%nobs
allocate(xyz(3,nobs),v(nobs)) ! 2017.09.08
xyz  = g_meshpara%xyzobs
x1 = g_meshpara%xbound(2)
x2 = g_meshpara%xbound(3)
y1 = g_meshpara%ybound(2)
y2 = g_meshpara%ybound(3)
si = g_meshpara%sizein
sb = g_meshpara%sizebo
A  = g_meshpara%A_obs ! [km]
sigma = g_meshpara%sigma_obs ! [km]

!#[1]## value default
 if ( x .lt. x1 .or. x2 .lt. x .or. y .lt. y1 .or. y2 .lt. y) then
  value=sb
 else
!#[2]## value for focus area or near observatories
  do i=1,nobs
   robs=dsqrt((xyz(1,i)-x)**2.d0 + (xyz(2,i)-y)**2.d0 )
   v(i) = si
   if ( robs .le. 3.*sigma ) then
    v(i) = si - (si - A)*exp(-robs**2./2./sigma**2.)
   end if
  end do
  value = si
  do i=1,nobs
   value = min(value,v(i))
  end do
 end if

return
end function


!########################################## value_3d
function value_3d(x,y,z,g_meshpara)
use param
implicit none
type(param_forward),intent(in)     :: g_meshpara
real(8),            intent(in)     :: x,y,z
real(8)                            :: si,sb, x1,x2,y1,y2,z1,z2,sigma,A
real(8)                            :: value_3d, robs
real(8),allocatable,dimension(:,:) :: xyz
integer(4)                         :: i,nobs

!#[0]## set input
nobs = g_meshpara%nobs
allocate(xyz(3,nobs)) ! 2017.09.08
xyz  = g_meshpara%xyzobs
x1   = g_meshpara%xbound(2)
x2   = g_meshpara%xbound(3)
y1   = g_meshpara%ybound(2)
y2   = g_meshpara%ybound(3)
z2   = g_meshpara%upzin
z1   = g_meshpara%downzin
si   = g_meshpara%sizein
sb   = g_meshpara%sizebo
A     = g_meshpara%A_obs ! [km]
sigma = g_meshpara%sigma_obs

!#[1]## value default
 if ( x .lt. x1 .or. x2 .lt. x .or. &
   &  y .lt. y1 .or. y2 .lt. y .or. &
   &  z .lt. z1 .or. z2 .lt. z ) then
  value_3d=sb
 else
!#[2]## value for focus area or near observatories
  value_3d=si
  do i=1,nobs
   robs=dsqrt((xyz(1,i)-x)**2.d0 + (xyz(2,i)-y)**2.d0 + (xyz(3,i)-z)**2.d0)
   if ( robs .le. 3.*sigma ) then
    value_3d = min(value_3d, si - (si - A)*exp(-robs**2./2./sigma**2.) )
   else
    value_3d = min(value_3d,si)
   end if
  end do
 end if

return
end function

!########################################### upgradegrd
subroutine updategrd(xgrd,nxmax,nx,x,dd)
implicit none
integer(4),intent(in)    :: nxmax
real(8),   intent(in)    :: dd,x
real(8),   intent(inout) :: xgrd(nxmax)
integer(4),intent(inout) :: nx
integer(4)               :: k,kk
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
end subroutine

end module horizontalresolution


