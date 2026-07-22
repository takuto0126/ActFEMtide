! Takuto Minami
! coded on Jan 5, 2021
program changemodel
use param
use mesh_type
use cond_type ! see m_cond_type.f90
implicit none
type(param_cmod) :: m_param
type(mesh)       :: h_mesh
integer(4)       :: nmodel,nphys2
integer(4)       :: i
integer(4),allocatable,dimension(:) :: ele2model,id,nmodelfile
real(8),   allocatable,dimension(:) :: rho_model,rho_model_new

!#[1]##
call readctl(m_param) ! see below

!#[2]## read 3d msh tetrahedral mesh
CALL READMESH_TOTAL(h_mesh,m_param%mshfile)

!#[3]## read model connectfile and modelfile
!#[3-1]## read connection
 open(1,file=m_param%connectfile,status='old',err=90)
 read(1,*) nphys2,nmodel
 allocate(id(nphys2),ele2model(nphys2))
 do i=1,nphys2
  read(1,'(2i10)') id(i),ele2model(i)
 end do
 close(1)
 write(*,*) "### READ CONNECTIONFILE END!! ###"

!#[3-2]## read model
open(1,file=m_param%inmodelfile,status='old',err=91) ! 20127.12.21
 read(1,*) nmodel
 write(*,*) " # of model parameters (nmodel)=", nmodel
 allocate(rho_model(nmodel))
 do i=1,nmodel
  read(1,*) rho_model(i)
 end do
close(1)
write(*,*) "### READ MODELFILE END!! ###"


!#[4]## change the model
allocate(rho_model_new(nmodel))
CALL ASSIGNRHO(h_mesh,nphys2,ele2model,nmodel,rho_model,rho_model_new,m_param)

!#[8]# output new model
open(1,file=m_param%outmodelfile)
 write(1,*) nmodel
 do i=1,nmodel
  write(1,*) rho_model_new(i)
 end do
close(1)
write(*,*) "OUTPUT NEW MODEL END",m_param%outmodelfile

stop
90 continue
write(*,*) "No connect file",m_param%connectfile
stop
91 continue
write(*,*) "No model file",m_param%inmodelfile
stop
end program

!####################################
subroutine readctl(m_param)
use cond_type
implicit none
type(param_cmod),intent(out) :: m_param
integer(4)       :: input=5,i

!#[1]# skip header
 read(input,*)

!#[2]# read
 read(input,10) m_param%mshfile
 read(input,10) m_param%connectfile
 read(input,10) m_param%inmodelfile
 read(input,10) m_param%outmodelfile

!#[2]# ratio
 read(input,12) m_param%ncuboid
 allocate(m_param%g_cuboid(m_param%ncuboid))

!#[5]#
 do i=1,m_param%ncuboid
   read(input,13) m_param%g_cuboid(i)%xminmax(1:2)
   read(input,13) m_param%g_cuboid(i)%yminmax(1:2)
   read(input,13) m_param%g_cuboid(i)%zminmax(1:2)
   read(input,11) m_param%g_cuboid(i)%rho
 end do

10 format(20x,a)
11 format(20x,g15.7)
12 format(20x,i10)
13 format(20x,2g15.7)
return
end

!####################################  ASSIGNCOND
!# copied from gridmdl2msh.f90
subroutine ASSIGNRHO(h_mesh,nphys2,ele2model,nmodel,rho_model,rho_model_new,m_param)
use mesh_type
use cond_type
implicit none
type(mesh),      intent(in)  :: h_mesh
type(param_cmod),intent(in)  :: m_param
integer(4),      intent(in)  :: nphys2,nmodel
integer(4),      intent(in)  :: ele2model(nphys2)
real(8),         intent(in)  ::rho_model(nmodel)
real(8),         intent(out) :: rho_model_new(nmodel)
integer(4)                            :: node,ntet
integer(4),allocatable,dimension(:,:) :: n4
integer(4),allocatable,dimension(:)   :: n4flag
real(8),   allocatable,dimension(:,:) :: xyz,cxyz
integer(4) :: i,j,k,l,nphys1
real(8) :: x2(2),y2(2),z2(2)

!
ntet            = h_mesh%ntet
node            = h_mesh%node
allocate(n4(ntet,4),xyz(3,node),n4flag(ntet))
allocate(cxyz(3,ntet))
n4         = h_mesh%n4
n4flag(:)  = h_mesh%n4flag(:,1)
xyz        = h_mesh%xyz
nphys1     = ntet - nphys2
!write(*,*) "check0"

!#[2]## assign center of gravity
do i=1,nphys2
 ![2-1]# center of gravity of this tetrahedron
 l = nphys1 + i
 do j=1,3
 cxyz(j,l)=(xyz(j,n4(l,1))+xyz(j,n4(l,2))+xyz(j,n4(l,3))+xyz(j,n4(l,4)))/4.d0
 end do
end do
!write(*,*) "check1"

!#[3]## change model
  rho_model_new = rho_model ! initial
  do i=1,m_param%ncuboid
   x2=m_param%g_cuboid(i)%xminmax
   y2=m_param%g_cuboid(i)%yminmax
   z2=m_param%g_cuboid(i)%zminmax
   
   do j=1,nphys2 ! tetrahedron loop
    l=nphys1 + j
    if (z2(1) .lt. cxyz(3,l) .and. cxyz(3,l) .lt. z2(2) ) then
    if (y2(1) .lt. cxyz(2,l) .and. cxyz(2,l) .lt. y2(2) ) then
    if (x2(1) .lt. cxyz(1,l) .and. cxyz(1,l) .lt. x2(2) ) then
     rho_model_new(ele2model(j)) = m_param%g_cuboid(i)%rho
    end if
    end if
    end if
   end do
  end do


write(*,*) "### ASSIGNRHO END!! ###"
return
end
