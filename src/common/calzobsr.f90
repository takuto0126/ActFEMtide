! Coded on 2017.02.21
subroutine calzobsr(h_mesh,g_param)
use param
use mesh_type
use triangle
implicit none
type(mesh),            intent(in)     :: h_mesh
type(param_forward),   intent(inout)  :: g_param
type(grid_list_type)                  :: glist
integer(4)                            :: nobs,nobsr,nx,ny ! 2021.09.29 nobs added
real(8),   allocatable,dimension(:,:) :: xyz_r,xyz
integer(4),allocatable,dimension(:,:) :: n3k
real(8),   allocatable,dimension(:)   :: znew
real(8)                               :: a3(3)
integer(4)                            :: iele,n1,n2,n3,j,ntri
real(8)                               :: xyzminmax(6),zorigin

write(*,'(a)') " ### CALZOBSR START!! ###"!2021.09.29

!#[1]## set
allocate(xyz(3,h_mesh%node),n3k(h_mesh%ntri,3))
allocate(xyz_r(3,g_param%nobsr))
allocate(znew(g_param%nobsr))
nobs      = g_param%nobs  ! 2021.09.29
nobsr     = g_param%nobsr
xyz       = h_mesh%xyz    ! triangle mesh
n3k       = h_mesh%n3
ntri      = h_mesh%ntri
xyz_r     = g_param%xyz_r
xyzminmax = g_param%xyzminmax
!write(*,*) "xyzminmax"

!#[2]## cal z for nobsr
nx=300;ny=300
CALL allocate_2Dgrid_list(nx,ny,ntri,glist)   ! see m_triangle.f90
CALL gen2Dgridforlist(xyzminmax,glist) ! see m_mesh_type.f90
CALL classifytri2grd(h_mesh,glist)   ! classify ele to glist,see


!#[3] search for the triangle including (x1,y1)
write(*,*) ""
do j=1,nobsr
if ( j == 1) write(*,'(a)') "< z is calculated for receiver points>" ! 2021.09.29
if ( j == nobs + 1) write(*,'(a)') "< z is calculated for virtual observatories along source wire>" ! 2021.09.29

call findtriwithgrid(h_mesh,glist,xyz_r(1:2,j),iele,a3)
! do i=1,ntri
    n1 = n3k(iele,1)
    n2 = n3k(iele,2)
    n3 = n3k(iele,3)
    znew(j) = a3(1)*xyz(3,n1)+a3(2)*xyz(3,n2)+a3(3)*xyz(3,n3) + xyz_r(3,j)
!
if ( j .le. nobs ) then ! 2021.09.29
   write(*,'(i3,1x,a3,a,3f14.7,a)') j,g_param%obsname(j)," (x,y,z)=",xyz_r(1:2,j),znew(j)," [km]" ! 2021.09.29
else
   write(*,'(i3,1x,a,3f14.7,a)') j,"    (x,y,z)=",xyz_r(1:2,j),znew(j)," [km]" ! 2021.09.29
end if ! 2021.09.29

end do
write(*,*) "" !2021.09.29

!#[4]## set znew to xyz_r
 g_param%xyz_r(3,1:nobsr) = znew(1:nobsr)

!#[5]## reflect topo
 call findtriwithgrid(h_mesh,glist,(/0.d0,0.d0/),iele,a3)
 zorigin = a3(1)*xyz(3,n1)+a3(2)*xyz(3,n2)+a3(3)*xyz(3,n3)
 g_param%zorigin = zorigin
 g_param%upzin   = g_param%upzin   + zorigin
 g_param%downzin = g_param%downzin + zorigin

write(*,'(a)') " ### CALZOBSR END!! ###"
return
end
