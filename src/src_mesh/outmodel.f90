!# Coded on 2017.02.22
program outmodel
use mesh_type
use param
implicit none
type(param_forward) :: g_param
type(param_source)  :: s_param
type(param_cond)    :: g_cond
type(mesh)          :: h_mesh
character(50)  :: name3d,mshfile,outfile="cond.msh"

!#[1]## read param
 CALL READPARAM(g_param,s_param,g_cond)
 if (g_cond%condflag .ne. 0 ) then
  write(*,*) "g_cond%condflag should be 0 in this program"
  write(*,*) "condflag=",g_cond%condflag
  stop
 end if

!#[2]## read 3dmesh
 name3d=g_param%header3d
 mshfile=name3d(1:len_trim(name3d))//".msh"
 CALL READMESH_TOTAL(h_mesh,mshfile)

!#[3]## OUT condmodel
 CALL OUTHOMCOND(h_mesh,g_cond,outfile)

end program

!########################################### OUTHOMCOND
subroutine OUTHOMCOND(h_mesh,g_cond,outfile)
use mesh_type
use param
implicit none
type(mesh),intent(in)       :: h_mesh
type(param_cond),intent(in) :: g_cond
character(50),   intent(in) :: outfile
integer(4) :: ntet,nlin,npoi,ntri
integer(4) :: j, ishift,icount
real(8),allocatable,dimension(:) :: rho_land ! 2017.09.28
integer(4) :: nvolume ! 2017.09.28
character(150) :: a

!#[0]## set
ntet = h_mesh%ntet
nlin = h_mesh%nlin
ntri = h_mesh%ntri
npoi = h_mesh%npoi
nvolume = g_cond%nvolume      ! 2017.09.28
allocate(rho_land(nvolume))   ! 2017.09.28
rho_land(:) = 1.d0/g_cond%sigma_land(:) ! 2017.09.28

!#[1]## create vector fields
open(1,file=outfile)
icount=0
do j=1,ntet
 if (h_mesh%n4flag(j,1) .eq. 2) icount=icount+1
end do
write(1,'(a)') "$ElementData"
write(1,'(a)') "1"
write(1,'(a)') '"A rho model view"'
write(1,'(a)') "1"
write(1,'(a)') "0.0"
write(1,'(a)') "3"
write(1,'(a)') "0"
write(1,'(a)') "1" ! means only one (scalar) value is assigned to element
write(1,'(i10)') icount
ishift=nlin + ntri + npoi
do j=1,ntet
 if ( h_mesh%n4flag(j,1) .ge. 2 ) then  ! only on land 2017.09.28
  write(1,*) ishift+j,rho_land(h_mesh%n4flag(j,1)-1) ! 2017.09.28
 end if
end do
write(1,'(a)') "$EndElementData"
close(1)

return
end
