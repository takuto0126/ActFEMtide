!# coded on 2021.06.01
subroutine cond3dto2d(g_mesh,g_surface,g_cond)
use surface_type
use param
use mesh_type
implicit none
type(surface),    intent(inout) :: g_surface(6)
type(mesh),       intent(in)    :: g_mesh
type(param_cond), intent(in)    :: g_cond
integer(4),allocatable,dimension(:,:) :: n4flag
integer(4) :: ntri,ntet,nphys2,iele_tet,i,j

ntet   = g_mesh%ntet
n4flag = g_mesh%n4flag

!# allocate g_surface(j)%cond(:)
do j=2,5
 ntri = g_surface(j)%ntri
 allocate( g_surface(j)%cond(ntri) )
end do

do j=2,5 ! 5 faces except top surface
 
 ntri = g_surface(j)%ntri
 if ( .not. allocated(g_surface(j)%cond)) allocate(g_surface(j)%cond(ntri))

 do i=1,ntri

 iele_tet = g_surface(j)%ifacetri_to_tet(i)

 if (      n4flag(iele_tet,1) .eq. 1 ) then ! in the air
  g_surface(j)%cond(i) = g_cond%sigma_air
 else if ( n4flag(iele_tet,1) .eq. 2 ) then ! under ground
  if ( g_cond%condflag .eq. 0)g_surface(j)%cond(i)=g_cond%sigma_land(1) ! 2017.09.29
  if ( g_cond%condflag .eq. 1)g_surface(j)%cond(i)=g_cond%sigma(iele_tet - g_cond%nphys1)
 else
  write(*,*) "GEGEGE n4flag=",n4flag(iele_tet,1)
  stop
 end if
 end do ! surface triangle loop

end do ! surface loop

return
end
