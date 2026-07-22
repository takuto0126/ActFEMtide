!# coded on 2018.11.08
program n_model2cond
use modelpart
use param
use mesh_type
implicit none

type(model)            :: g_model      ! see m_modelpart.f90
type(param_cond)       :: g_cond       ! see m_param.f90 : goal conductivity
type(mesh)             :: g_mesh
character(50)          :: mshfile,connectfile,modelfile,outfile
integer(4)             :: nphys1,nphys2, nmodel,i,j,nlin,npoi,ntri,ishift
integer(4),allocatable,dimension(:) :: id,ele2model,index
real(8),   allocatable,dimension(:) :: value
!# for mask
real(8)               :: threshold
character(50)         :: sensfile
integer(4)            :: iunit ! 0 : Ohm.m, 1: log10(Ohm.m) 2021.01.08
integer(4)            :: imask ! 0:nomask, 1:mask
real(8),   allocatable,dimension(:) :: value_mask

!#[0]## read parameters
  call readfilename(mshfile,connectfile,modelfile,&
&                       iunit,imask,threshold,sensfile,outfile)!2021.01.08

!#[1]## read mesh
 write(*,*) "mshfile",mshfile
  CALL READMESH_TOTAL(g_mesh,mshfile)

!#[2]## read model
!#[2-1]## read connectfile
  open(1,file=connectfile)
  read(1,*) nphys2,nmodel
  allocate(id(nphys2),ele2model(nphys2))
  do i=1,nphys2
   read(1,'(2i10)') id(i),ele2model(i)
  end do
  close(1)
  write(*,*) "### READ CONNECTIONFILE END!! ###"

 !#[2-2]## read model
 allocate(value(nmodel),value_mask(nmodel))
 open(1,file=modelfile,status='old',err=90) ! 2017.12.21
  read(1,*,err=80) nmodel
  do i=1,nmodel
   read(1,*) value(i)
  end do
 close(1)
 write(*,*) "### READ MODELFILE END!! ###"
 !#
 if (imask .eq. 1 ) then
  open(1,file=sensfile)
   read(1,*) nmodel
   do i=1,nmodel
    read(1,*) value_mask(i)
   end do
  close(1)
 end if
 !#
 goto 100
  80 continue
  close(1)
  write(*,*) "Cannot read modelfile"
  stop
  goto 100
  90 continue
  write(*,*) "File is not exist",modelfile
  stop
 100 continue ! 2018.01.17
 nphys1 = g_mesh%ntet - nphys2
 write(*,*) "nphys1",nphys1
 write(*,*) "nphys2",nphys2 ! 2020.12.02
 allocate(index(nphys2))
 do i=1,nphys2
  index(i) = nphys1+i
 end do

 !#[2-3]## set model
  g_model%nmodel       = nmodel
  g_model%nphys1       = g_mesh%ntet -nphys2
  g_model%nphys2       = nphys2
  g_model%ele2model    = ele2model
  if ( iunit .eq. 0 ) then             ! 2021.01.08
   g_model%logrho_model = log10(value) ! 2021.01.08
  elseif (iunit .eq. 1) then           ! 2021.01.08
   g_model%logrho_model = value        ! 2021.01.08
  else                                 ! 2021.01.08
   write(*,*) "GEGEGE iunit",iunit     ! 2021.01.08
   write(*,*) "input should be 0: [Ohm.m] or 1:(log10(Ohm.m))" ! 2021.01.08
   stop  ! 2021.01.08
  end if ! 2021.01.08
  g_model%index        = index

!#[3]## generate cond
  !#[3-1]## initial g_cond
  g_cond%condflag  = 0 ! 0:homogeneous
  g_cond%sigma_air = 1.d-8 ! cond for air
  g_cond%nvolume   = 1  ! # of physical volume in land
  allocate(g_cond%sigma_land(1))
  g_cond%sigma_land(1)=0.01

!#[4]## generate cond
  if (iunit .eq. 0 ) then ! unit [Ohm.m] 2021.01.08
   call model2cond(g_model,g_cond,0) ! 0 means g_cond%rho = 10**g_model%logrho_model
  elseif (iunit .eq. 1) then ! unit log10([Ohm.m]) 2021.01.08
   call model2cond(g_model,g_cond,1) ! 1 means g_cond%rho = g_model%logrho_model
  end if

!#[5]## output
 npoi = g_mesh%npoi
 nlin = g_mesh%nlin
 ntri = g_mesh%ntri
 open(1,file=outfile)
 write(1,'(a)') "$MeshFormat"     ! 2017.09.13
 write(1,'(a)') "2.2 0 8"        ! 2017.09.13
 write(1,'(a)') "$EndMeshFormat"  ! 2017.09.13
 write(1,'(a)') "$ElementData"
 write(1,'(a)') "1"
 write(1,'(a)') '"A rho model view"'
 write(1,'(a)') "1"
 write(1,'(a)') "0.0"
 write(1,'(a)') "3"
 write(1,'(a)') "0"
 write(1,'(a)') "1" ! means only one (scalar) value is assigned to element
 write(1,'(i10)') nphys2
 ishift = npoi + nlin + ntri
 do j=1,nphys2
!  write(*,*) "j=",j,"nphys2=",nphys2,"ele2model(j)=",ele2model(j),"nmodel=",nmodel
!   write(1,*) id(j),g_cond%rho(j)
  if ( imask .ne. 1 ) then
   write(1,*) ishift+index(j),g_cond%rho(j)
  else
   if ( value_mask(ele2model(j)) .gt. threshold ) then
    write(1,*) ishift+index(j),g_cond%rho(j)
   else
    write(1,*) ishift+index(j),10.0 ! 2019.05.26
   end if
  end if
 end do
 write(1,'(a)') "$EndElementData"
close(1)

end program

subroutine readfilename(mshfile,connectfile,modelfile,&
&                       iunit,imask,threshold,sensfile,outfile)!2021.01.08 iunit is added
implicit none
character(50)          :: mshfile,connectfile,modelfile,outfile,sensfile
integer(4)             :: iunit ! 2021.01.08
integer(4)             :: imask
real(8)                :: threshold

write(*,*) "input mshfile"
read(5,'(a)') mshfile
write(*,*) "connect connectfile"
read(5,'(a)') connectfile
write(*,*) "input modelfile"
read(5,'(a)') modelfile
write(*,*) "parameter input : 0 for Ohm.m , 1 for log10(Ohm.m)" ! 2021.01.08
read(5,'(i5)') iunit
write(*,*) "0 for no mask, 1 for mask model"
read(5,'(i5)') imask ! 2020.12.01
if ( imask .eq. 1 ) then
 write(*,*) "input sensitivity model"
 read(5,'(a)') sensfile
 write(*,*) " input threshold"
 read(5,*)  threshold
end if


write(*,*) "input output condfile"
read(5,'(a)') outfile

return
end

