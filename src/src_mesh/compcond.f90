program compcond
implicit none
character(50) :: fcond(2)
integer(4)    :: i,j,n1,n2
real(8),allocatable,dimension(:,:) :: rho
integer(4),allocatable,dimension(:,:) :: index
real(8)                          :: rho1,rho2

write(*,*) "Please input condfile1"
read(*,'(a)') fcond(1)
write(*,*) "Please input condfile2"
read(*,'(a)') fcond(2)

 open(1,file=fcond(1))
 open(2,file=fcond(2))
 do j=1,11
  read(1,*)
  read(2,*)
 end do
 read(1,*) n1
 read(2,*) n2

 if ( n1 .ne. n2 ) then
  write(*,*) "GEGEGE!", n1,"not equal to",n2
 stop
 end if

 allocate(rho(n1,2),index(n1,2))
 do j=1,n1
  read(1,*) index(j,1),rho(j,1)
  read(2,*) index(j,2),rho(j,2)
 end do
 close(1)
 close(2)

open(1,file="./cond_d.msh")
 write(1,'(a)') "$MeshFormat"
 write(1,'(a)') "2.2 0 8"
 write(1,'(a)') "$EndMeshFormat"
 write(1,'(a)') "$ElementData"
 write(1,'(a)') "1"
 write(1,'(a)') '"A rho model view"'
 write(1,'(a)') "1"
 write(1,'(a)') "0.0"
 write(1,'(a)') "3"
 write(1,'(a)') "0"
 write(1,'(a)') "1" ! means only one (scalar) value is assigned to element
 write(1,'(i10)') n1
 do j=1,n1
  rho1 = rho(j,1)
  rho2 = rho(j,2)
  write(1,'(i10,e15.7)') index(j,1),log10(rho2) -log10(rho1)
 end do
 write(1,'(a)') "$EndElementData"

close(1)

end program
