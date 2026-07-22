! Coded on 2017.03.06
program makeface
use mesh_type  ! see m_mesh_type.f90
use face_type  ! see m_face_type.f90
use param      ! see m_param.f90
use horizontalresolution ! see m_horizontalresolution.f90
implicit none
type(param_forward) :: g_param
type(param_source)  :: s_param
character(50) :: name,filename
type(mesh)            :: g_mesh     ! see m_mesh_type.f90
type(face_info)       :: g_face     ! see m_face_type.f90
integer(4) :: nface, ntet,i,j,node
integer(4),allocatable,dimension(:,:) :: n4

!#[0]##
  CALL readparam(g_param,s_param)
  name = g_param%header3d
  filename=name(1:len_trim(name))//".msh"

!#[1]## Mesh READ
  CALL READMESH_TOTAL(g_mesh,filename)

!#[2]## Line information
  ntet = g_mesh%ntet  ! 2107.09.08
  node = g_mesh%node  ! 2017.09.08
  allocate(n4(ntet,4))! 2017.09.08
  n4   = g_mesh%n4    ! 2107.09.08

  CALL MKFACE(  g_face,node,ntet,4, n4) ! make g_line
  CALL MKN4FACE(g_face,node,ntet,   n4) ! make g_line%n6line
  nface= g_face%nface
  ntet = g_face%ntet
  CALL FACE2ELEMENT(g_face)

!#[3]## Output lines and n6line
open(10,file="faceinfo.dat")
 write(10,*) nface
 write(10,'(3i10)') ((g_face%face(j,i),j=1,3),i=1,nface)
 write(10,*) ntet
 write(10,'(4i10)') ((g_face%n4face(i,j),j=1,4),i=1,ntet)
close(10)

end program makeface
