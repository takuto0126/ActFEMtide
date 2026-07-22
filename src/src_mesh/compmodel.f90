!# coded on 2017.06.14
program compmodel
use mesh_type
implicit none
type(mesh) :: g_mesh,ref_mesh
character(50) :: modelfile != "model06.msh"
character(50) :: refmodelfile != "model01.msh"
character(50) :: outfile="model_d.msh"
integer(4)    :: ifile,n1,n2,j
real(8)       :: rho1,rho2

write(*,*) "input reference model file"
read(5,'(a50)') refmodelfile
write(*,*) "input model file to be compared"
read(5,'(a50)') modelfile

call READMESH_TOTAL(g_mesh,modelfile)
call READMESH_TOTAL(ref_mesh,refmodelfile)

n1=ref_mesh%nmodel
n2=  g_mesh%nmodel
if ( n1 .ne. n2) then
 write(*,*) "GEGEGE n1=",n1,"n2=",n2
 stop
end if

ifile=1
open(ifile,file=outfile)
 call MESHOUT(ifile,g_mesh)
 write(1,'(a)') "$ElementData"
 write(1,'(a)') "1"
 write(1,'(a)') '"A rho model view"'
 write(1,'(a)') "1"
 write(1,'(a)') "0.0"
 write(1,'(a)') "3"
 write(1,'(a)') "0"
 write(1,'(a)') "1" ! means only one (scalar) value is assigned to element
 write(1,'(i10)') g_mesh%nmodel
 do j=1,g_mesh%nmodel
  rho1 = ref_mesh%cmodel(j)
  rho2 = g_mesh%cmodel(j)
  write(1,'(i10,e15.7)') g_mesh%index(j),log10(rho2) -log10(rho1)
 end do
 write(1,'(a)') "$EndElementData"

close(ifile)

end program compmodel
