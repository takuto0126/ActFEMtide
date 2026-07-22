!# coded on 2022.10.24
module cond_model_type
use param
use mesh_type
implicit none

type(info_cond_model)
integer(4) :: iflag ! 0: condfile, 1:model file
character(100) :: condfile
character(100) :: modelfile
character(100) :: connectfile
end type

contains

!##################################################################################
subroutine READGENCOND(g_info_cond_model,h_mesh,g_cond)
 implicit none
 type(info_cond_model),intent(in)  :: g_info_cond_model
 type(mesh),           intent(in)  :: h_mesh
 type(param_cond),     intent(out) :: g_cond
 integer(4) :: iflag   ! 0: condfile, 1:model file
 integer(4) :: iexist

  !#[1]## set
   iflag = g_info_cond_model%iflag ! 0: condfile, 1:model file
    !# flag = 0 requires:
     sigmahomo = g_info_cond_model%sigmahomo            ! 2018.10.04
    !# flag = 1 requires:
     g_cond%condfile = g_info_cond_model%condfile      ! 2017.07.19
    !# flag = 2 requires:
     connectfile = g_info_cond_model%connectfile
     modelfile   = g_info_cond_model%modelfile

  !#[2]## generate cond from either of cond or model
   if     ( iflag .eq. 0 ) then                ! 2018.10.04
     sigmahomo = g_info_cond_model%sigmahomo            ! 2018.10.04
     CALL SETCOND(g_cond,h_mesh,sigmahomo)              ! 2018.10.04

   elseif ( iflag .eq. 1 ) then ! condfile       2018.03.18
     CALL READCOND(g_cond)                              ! 2017.07.19

   elseif ( iflag .eq. 2 ) then ! modelfile
     CALL READMODEL2COND(g_cond,connectfile,modelfile)

   else                                                ! 2018.10.04
     write(*,*) "GEGEGE! icondflag_ref",icondflag_ref   ! 2018.10.04
     stop                                               ! 2018.10.04

 end if



return
end

!##################################################### READMODEL2COND
!# coded on 2018.06.21
subroutine readmodel2cond(r_cond,connectfile,modelfile)
 use modelpart
 use param
 implicit none
 type(param_cond), intent(inout)     :: r_cond
 character(50),    intent(in)        :: connectfile
 character(50),    intent(in)        :: modelfile
 integer(4)                          :: nphys2,nmodel,nmodel2
 integer(4)                          :: i
 integer(4),allocatable,dimension(:) :: id,ele2model
 real(8),   allocatable,dimension(:) :: rho

 !#[1]## read connectfile
  open(1,file=connectfile,status='old',err=90)
  read(1,*) nphys2,nmodel
  allocate(id(nphys2),ele2model(nphys2))
  do i=1,nphys2
   read(1,'(2i10)') id(i),ele2model(i)
  end do
  close(1)

 !#[2]## read model
 allocate(rho(nmodel))
 open(1,file=modelfile,status='old',err=80) ! 20127.12.21
  read(1,*,err=81) nmodel2
  if ( nmodel .ne. nmodel2 ) then
   write(*,*) "GEGEGE nmodel",nmodel,"nmodel2",nmodel2
   stop
  end if
  do i=1,nmodel
   read(1,*) rho(i)
  end do
 close(1)

 !#[3]## gen r_cond
 allocate(r_cond%rho(  nphys2) )
 allocate(r_cond%sigma(nphys2) )
 allocate(r_cond%index(nphys2)) ! 2018.03.20
 r_cond%nphys2 = nphys2
 r_cond%sigma  = -9999
 do i=1,nphys2
  r_cond%rho(i)=rho(ele2model(i))
  if (abs(r_cond%rho(i)) .gt. 1.d-10 ) r_cond%sigma(i) = 1.d0/r_cond%rho(i)
 end do

 return
  90 continue
  write(*,*) "File is not exist",connectfile
  stop
  80 continue
  write(*,*) "File is not exist",modelfile
  stop
  81 continue
  write(*,*) "File is not exist",modelfile,"line",i
  stop
 end




end module
