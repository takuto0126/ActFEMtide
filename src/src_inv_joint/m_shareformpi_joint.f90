! Modified on 2017.08.31 for multisource amp phase inversion
! Coded on 2017.05.50
! to share structures between processes for MPI parallelization
module shareformpi_joint ! ap -> joint 2021.12.25
use param
use mesh_type
use line_type
use outresp     ! 2017.06.02
use face_type   ! 2017.06.02
use matrix      ! 2017.06.05
use param_jointinv ! 2021.12.25
use modelpart   ! 2017.06.05
use caltime     ! 2017.12.22
implicit none
include 'mpif.h'

contains
!############################################################
!# modified on 2017.08.31 for multisource inversion
!# Coded on 2017.06.07
subroutine shareapinv(g_param,sparam,g_cond,g_mesh,g_line,g_param_joint,g_model,ip)
implicit none
integer(4),         intent(in)      :: ip
type(param_forward),intent(inout)   :: g_param
type(param_source), intent(inout)   :: sparam
type(param_cond),   intent(inout)   :: g_cond
type(mesh),         intent(inout)   :: g_mesh
type(line_info),    intent(inout)   :: g_line
!type(face_info),    intent(inout) :: g_face
type(param_joint), intent(inout) :: g_param_joint ! 2017.08.31
type(model),          intent(inout) :: g_model

!  write(*,*) "### SHAREAPINV START ###"
  CALL SHAREFORWARD(g_param,sparam,g_cond,ip) ! see m_shareformpi.f90
  CALL SHAREMESHLINE(g_mesh,g_line,ip)        ! see m_shareformpi.f90
!  CALL SHAREFACE(g_face,ip)                  ! see m_shareformpi.f90 2017.08.31
  CALL SHAREINVPARAJOINT(g_param_joint,ip)       ! 2017.08.31 see below
  CALL SHAREMODEL(g_model,ip)                 ! see m_shareformpi.f90
return
end

!############################################################
! coded on 2021.12.30
subroutine sharemt(g_param_mt,g_surface,ip)
  use surface_type
  use param_mt
  implicit none
  integer(4),             intent(in)   :: ip
  type(param_forward_mt),intent(inout) :: g_param_mt
  type(surface),         intent(inout) :: g_surface(6)
  integer(4)                           :: i

  do i=1,6
   call sharesurface(g_surface(i),ip)
  end do
  write(*,*) "### SHARE SURFACE END!! ###"

  call shareparamforwardmt(g_param_mt,ip)
  
  if ( ip .ne. 0 ) then ! 2021.12.31
   g_surface(1)%facetype="xy" ! top
   g_surface(2)%facetype="xz" ! north
   g_surface(3)%facetype="yz" ! west
   g_surface(4)%facetype="xz" ! south
   g_surface(5)%facetype="yz" ! east
   g_surface(6)%facetype="xy" ! bottom
  end if

  return
  end
!############################################################
  !# coded on 2021.12.30
  subroutine shareparamforwardmt(g_param,ip)
    use param_mt
    implicit none
    type(param_forward_mt),intent(inout) :: g_param
    integer(4),            intent(in)    :: ip
    integer(4) :: errno,nfreq,nobsr,nobs
    integer(4) :: nfile ! 2017.12.13

!#[3]## broadcast
 call MPI_BCAST(g_param%itopoflag,      1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno) ! 2017.12.13
 if ( g_param%itopoflag .eq. 1) then      ! 2017.12.13 when topography is considered
  call MPI_BCAST(g_param%nfile,         1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno) ! 2017.12.13
  nfile = g_param%nfile                   ! 2017.12.13
  if ( ip .ne. 0 ) then                   ! 2017.12.13
   allocate(g_param%topofile(nfile))      ! 2017.12.13
   allocate(g_param%lonlatshift(2,nfile)) ! 2017.12.13
  end if                                  ! 2017.12.13
  call MPI_BCAST(g_param%topofile,  nfile*50, MPI_CHAR,   0,MPI_COMM_WORLD,errno)! 2017.12.13
  call MPI_BCAST(g_param%lonlatshift,nfile*2, MPI_REAL8,  0,MPI_COMM_WORLD,errno)! 2017.12.13
 end if                                  ! 2017.09.29
 call MPI_BCAST(g_param%g_meshfile,    50, MPI_CHAR,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%surface_id_ground, 1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%z_meshfile,    50, MPI_CHAR,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%g_lineinfofile,50, MPI_CHAR,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%outputfolder,  50, MPI_CHAR,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%header2d,      50, MPI_CHAR,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%header3d,      50, MPI_CHAR,0,MPI_COMM_WORLD,errno)

 call MPI_BCAST(g_param%nfreq,          1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 nfreq = g_param%nfreq
 if ( ip .ne. 0 ) allocate(g_param%freq(nfreq))
 call MPI_BCAST(g_param%freq(1),     nfreq, MPI_REAL8,0,MPI_COMM_WORLD,errno)

 call MPI_BCAST(g_param%nobs,          1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%lonlatflag,    1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%wlon,     1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%elon,     1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%slat,     1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%nlat,     1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%lonorigin,1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%latorigin,1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%lenout,   1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%upzin,    1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%downzin,  1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%zmax,     1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%zmin,     1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%sizein,   1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%sizebo,   1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%sigma_obs,1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%A_obs,    1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%sigma_src,1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%A_src,    1, MPI_REAL8,0,MPI_COMM_WORLD,errno)

 call MPI_BCAST(g_param%nobsr,    1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 nobsr = g_param%nobsr
 if ( ip .ne. 0) then
  allocate(g_param%xyz_r(3,nobsr))
  allocate(g_param%A_r(nobsr)    )
  allocate(g_param%sigma_r(nobsr))
 end if

 call MPI_BCAST(g_param%xyz_r,nobsr*3, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%A_r,    nobsr, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%sigma_r,nobsr, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%xbound,     4, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%ybound,     4, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%zbound,     4, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%UTM,        3, MPI_CHAR, 0,MPI_COMM_WORLD,errno)

 nobs = g_param%nobs
 if ( ip .ne. 0) then
  allocate(g_param%lonlataltobs(3,nobs))
  allocate(g_param%xyzobs(3,nobs)      )
  allocate(g_param%obsname(nobs)       )
 end if

 call MPI_BCAST(g_param%lonlataltobs,3*nobs, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%xyzobs,      3*nobs, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%obsname,    50*nobs, MPI_CHAR, 0,MPI_COMM_WORLD,errno)

 call MPI_BCAST(g_param%xyzminmax,   6, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%zorigin,     1, MPI_REAL8,0,MPI_COMM_WORLD,errno)

 call MPI_BCAST(g_param%condfile,   50, MPI_CHAR,    0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%condflag,    1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)


 if (ip .eq. 0 ) write(*,'(a)') " ### SHARE FPARAMFORWARDMT   END!! ###" ! 2020.09.17
    return
    end

  !############################################################
  !# coded on 2021.12.30
subroutine sharesurface(g_surface,ip)
use surface_type
  implicit none
  integer(4),   intent(in)    :: ip
  type(surface),intent(inout) :: g_surface
  integer(4) :: errno
  integer(4) :: node3d,nline,ntri,node

  call MPI_BCAST(g_surface%nline,     1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_surface%ntri,      1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_surface%node,      1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_surface%node3d,    1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  node3d = g_surface%node3d
  ntri   = g_surface%ntri
  nline  = g_surface%nline
  node   = g_surface%node

  if ( ip .ne. 0 ) then
   allocate( g_surface%x1x2_face(    2,node)   )
   allocate( g_surface%inodeface_to_3D(node)   ) 
   allocate( g_surface%inode3D_to_face(node3d) )
   allocate( g_surface%n3line(ntri,3)          )
   allocate( g_surface%n3(ntri,3)              )
   allocate( g_surface%ilineface_to_3D(nline)  )
   allocate( g_surface%line(2,nline)           )
   allocate( g_surface%ifacetri_to_face(ntri)  )
   allocate( g_surface%ifacetri_to_tet( ntri)  )
   !
   allocate( g_surface%iflag_bound(nline))
  end if

  call MPI_BCAST(g_surface%x1x2_face(1,1),     2*node, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_surface%n3(1,1),            3*ntri, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_surface%n3line(1,1),        3*ntri, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_surface%inodeface_to_3D(1),   node, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_surface%inode3D_to_face(1), node3d, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_surface%line(1,1),         2*nline, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_surface%ilineface_to_3D(1),  nline, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_surface%ifacetri_to_face(1),  ntri, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_surface%ifacetri_to_tet(1),   ntri, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)

  call MPI_BCAST(g_surface%iflag_bound(1),      nline, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  return
  end 
!############################################################
!# Coded on 2017.06.05
subroutine SHAREMODEL(h_model,ip)
implicit none
integer(4), intent(in)    :: ip
type(model),intent(inout) :: h_model
integer(4) :: nmodel,nphys1,nphys2,errno, ip_from
type(watch) :: t_watch ! 2017.12.22

 call watchstart(t_watch) ! 2017.12.22
 call MPI_BCAST(h_model%nmodel,      1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(h_model%nmodelactive,1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno) ! 2018.06.21
 call MPI_BCAST(h_model%nphys1,      1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(h_model%nphys2,      1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 nmodel = h_model%nmodel
 nphys1 = h_model%nphys1
 nphys2 = h_model%nphys2

 if ( ip .ne. 0) then
  if ( .not. allocated(h_model%index)       ) allocate( h_model%index(nphys2)       )
  if ( .not. allocated(h_model%ele2model)   ) allocate( h_model%ele2model(nphys2)   )
  if ( .not. allocated(h_model%rho_model)   ) allocate( h_model%rho_model(nmodel)   )
  if ( .not. allocated(h_model%logrho_model)) allocate( h_model%logrho_model(nmodel))
  if ( .not. allocated(h_model%iactive)     ) allocate( h_model%iactive(     nmodel)) ! 2018.06.21
 end if

 call MPI_BCAST(h_model%index(1),     nphys2, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(h_model%ele2model(1), nphys2, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(h_model%iactive(1)  , nmodel, MPI_INTEGER4,0,MPI_COMM_WORLD,errno) ! 2018.06.21
 call MPI_BCAST(h_model%rho_model,    nmodel, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(h_model%logrho_model, nmodel, MPI_REAL8,   0,MPI_COMM_WORLD,errno)

 ip_from = 0
 call sharerealcrsmatrix(h_model%model2ele,ip_from,ip) ! see below

 call watchstop(t_watch) ! 2017.12.22
 if (ip .eq. 0) write(*,'(a,f7.3,a)') " ### SHAREMODEL  END!! ###   Time =",t_watch%time," [min]"!2020.09.17

return
end

!############################################################
!# Modified on 2017.08.31 for multisource inversion
!# Coded on 2017.06.08
subroutine SHAREINVPARAJOINT(g_param_joint,ip) ! 2017.08.31
implicit none
integer(4),             intent(in)    :: ip
type(param_joint),intent(inout) :: g_param_joint   ! 2017.08.31
integer(4)                            :: i,j,nobs,nfreq,errno
integer(4)                            :: nsr_inv, nobs_s ! 2017.07.13
integer(4)                            :: icomp           ! 2018.10.05

 call MPI_BCAST(g_param_joint%ijoint,  1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)    ! 2022.10.14
 call MPI_BCAST(g_param_joint%nobs,    1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)    ! 2017.08.31
 call MPI_BCAST(g_param_joint%nobs_mt, 1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)    ! 2022.01.04
 call MPI_BCAST(g_param_joint%nfreq,   1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)    ! 2017.08.31
 call MPI_BCAST(g_param_joint%nfreq_mt,1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)    ! 2022.01.04
 call MPI_BCAST(g_param_joint%g_faceinfofile, 50, MPI_CHAR,0,MPI_COMM_WORLD,errno)! 2017.08.31
 call MPI_BCAST(g_param_joint%g_initcondfile, 50, MPI_CHAR,0,MPI_COMM_WORLD,errno)! 2017.08.31
 call MPI_BCAST(g_param_joint%nsr_inv, 1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)    ! 2017.08.31
 nobs     = g_param_joint%nobs    ! 2017.08.31
 nfreq    = g_param_joint%nfreq   ! 2017.08.31
 nsr_inv  = g_param_joint%nsr_inv ! 2017.07.13

 if ( ip .ne. 0) then
  allocate(g_param_joint%srcindex(nsr_inv) ) ! 2017.08.31
  allocate(g_param_joint%obsinfo( nsr_inv) ) ! 2017.08.31
 end if
  call MPI_BCAST(g_param_joint%srcindex,nsr_inv,MPI_INTEGER4, 0,MPI_COMM_WORLD,errno) !2017.08.31

 do j=1,nsr_inv       ! 2017.07.13
  call MPI_BCAST(g_param_joint%obsinfo(j)%nobs_s,1,MPI_INTEGER4, 0,MPI_COMM_WORLD,errno) !2017.08.31
  nobs_s = g_param_joint%obsinfo(j)%nobs_s                ! 2017.08.31
  if ( ip .ne. 0) then                                    ! 2017.08.31
   allocate(g_param_joint%obsinfo(j)%ampfile(5,nobs_s) )  ! 2018.10.04
   allocate(g_param_joint%obsinfo(j)%phafile(5,nobs_s) )  ! 2018.10.04
   allocate(g_param_joint%obsinfo(j)%obsindex(2,5,nobs_s))! 2018.10.04
  end if                                                  ! 2017.08.31
  do icomp = 1,5                                          ! 2018.10.04
  do i=1,nobs_s                                           ! 2017.08.31
  call MPI_BCAST(g_param_joint%obsinfo(j)%ampfile(icomp,i),50,MPI_CHAR, 0,MPI_COMM_WORLD,errno)     !2018.10.04
  call MPI_BCAST(g_param_joint%obsinfo(j)%phafile(icomp,i),50,MPI_CHAR, 0,MPI_COMM_WORLD,errno)     !2018.10.04
  call MPI_BCAST(g_param_joint%obsinfo(j)%obsindex(1,icomp,i),2,MPI_INTEGER4,0,MPI_COMM_WORLD,errno)!2018.10.04
  end do         ! 2017.07.13
  end do         ! 2018.10.04
 end do          ! 2017.07.13

 call MPI_BCAST(g_param_joint%iflag_comp(1),5,MPI_INTEGER4,0,MPI_COMM_WORLD,errno) ! 2018.10.04

 call MPI_BCAST(g_param_joint%errorfloor_act, 1,MPI_REAL8,0,MPI_COMM_WORLD,errno)      ! 2021.12.27
 call MPI_BCAST(g_param_joint%nalpha,     1,MPI_INTEGER4,0,MPI_COMM_WORLD,errno)   ! 2017.09.11
 call MPI_BCAST(g_param_joint%ialphaflag, 1,MPI_INTEGER4,0,MPI_COMM_WORLD,errno)   ! 2017.09.11

 if ( ip .eq. 0 ) write(*,'(a)') " ### SHAREINVPARAJOINT  END!! ###" ! 2022.01.04
return
end


!############################################################
! Coded on 2017.06.05
subroutine sharerealcrsmatrix(rmat,ip_from,ip)
implicit none
integer(4),           intent(in)    :: ip_from,ip
type(real_crs_matrix),intent(inout) :: rmat
integer(4)                          :: nrow,ncolm,ntot,errno

 call MPI_BCAST(rmat%nrow,  1, MPI_INTEGER4,ip_from,MPI_COMM_WORLD,errno)
 call MPI_BCAST(rmat%ncolm, 1, MPI_INTEGER4,ip_from,MPI_COMM_WORLD,errno)
 call MPI_BCAST(rmat%ntot,  1, MPI_INTEGER4,ip_from,MPI_COMM_WORLD,errno)
 nrow  = rmat%nrow
 ncolm = rmat%ncolm
 ntot  = rmat%ntot

if ( ip .ne. ip_from ) then
 if (.not. allocated(rmat%stack)) allocate(rmat%stack(0:nrow))
 if (.not. allocated(rmat%item) ) allocate(rmat%item(ntot))
 if (.not. allocated(rmat%val)  ) allocate(rmat%val(ntot))
 if ( allocated(rmat%stack) .and. (nrow + 1 .ne. size(rmat%stack))) then
  deallocate(rmat%stack) ; allocate(rmat%stack(0:nrow))
 end if
 if ( allocated(rmat%item) .and. (ntot .ne. size(rmat%item))) then
  deallocate(rmat%item) ; allocate(rmat%item(ntot))
 end if
 if ( allocated(rmat%val) .and. (ntot .ne. size(rmat%val))) then
  deallocate(rmat%val) ; allocate(rmat%val(ntot))
 end if
end if

 call MPI_BCAST(rmat%stack(0),nrow+1, MPI_INTEGER4,ip_from,MPI_COMM_WORLD,errno)
 call MPI_BCAST(rmat%item(1),   ntot, MPI_INTEGER4,ip_from,MPI_COMM_WORLD,errno)
 call MPI_BCAST(rmat%val(1),    ntot, MPI_REAL8,   ip_from,MPI_COMM_WORLD,errno)

return
end

!############################################################
! Coded on 2022.01.05
subroutine sharecomplexcrsmatrix(rmat,ip_from,ip)
  implicit none
  integer(4),           intent(in)    :: ip_from,ip
  type(complex_crs_matrix),intent(inout) :: rmat
  integer(4)                          :: nrow,ncolm,ntot,errno
  
   call MPI_BCAST(rmat%nrow,  1, MPI_INTEGER4,ip_from,MPI_COMM_WORLD,errno)
   call MPI_BCAST(rmat%ncolm, 1, MPI_INTEGER4,ip_from,MPI_COMM_WORLD,errno)
   call MPI_BCAST(rmat%ntot,  1, MPI_INTEGER4,ip_from,MPI_COMM_WORLD,errno)
   nrow  = rmat%nrow
   ncolm = rmat%ncolm
   ntot  = rmat%ntot
  
  if ( ip .ne. ip_from ) then
   if (.not. allocated(rmat%stack)) allocate(rmat%stack(0:nrow))
   if (.not. allocated(rmat%item) ) allocate(rmat%item(ntot))
   if (.not. allocated(rmat%val)  ) allocate(rmat%val(ntot))
   if ( allocated(rmat%stack) .and. (nrow + 1 .ne. size(rmat%stack))) then
    deallocate(rmat%stack) ; allocate(rmat%stack(0:nrow))
   end if
   if ( allocated(rmat%item) .and. (ntot .ne. size(rmat%item))) then
    deallocate(rmat%item) ; allocate(rmat%item(ntot))
   end if
   if ( allocated(rmat%val) .and. (ntot .ne. size(rmat%val))) then
    deallocate(rmat%val) ; allocate(rmat%val(ntot))
   end if
  end if
  
   call MPI_BCAST(rmat%stack(0),nrow+1, MPI_INTEGER4,ip_from,MPI_COMM_WORLD,errno)
   call MPI_BCAST(rmat%item(1),   ntot, MPI_INTEGER4,ip_from,MPI_COMM_WORLD,errno)
   call MPI_BCAST(rmat%val(1),    ntot, MPI_COMPLEX8,ip_from,MPI_COMM_WORLD,errno)!2022.01.05
  
  return
  end
!############################################################
! Coded on 2017.06.02
subroutine shareface(g_face,ip)
implicit none
type(face_info),intent(inout) :: g_face
integer(4),intent(in) :: ip
integer(4) :: errno,nface,node,ntri,ntet

 call MPI_BCAST(g_face%nface, 1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_face%node,  1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_face%ntet,  1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_face%ntri,  1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)

nface = g_face%nface
node  = g_face%node
ntri  = g_face%ntri
ntet  = g_face%ntet

if ( ip .ne. 0) then
 allocate(g_face%face(3,nface))
 allocate(g_face%n4face(ntet,4))
 allocate(g_face%face2ele(2,nface))
end if

 call MPI_BCAST(g_face%face(1,1),    3*nface, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_face%n4face(1,1),   4*ntet, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_face%face2ele(1,1),2*nface, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)

return
end

!############################################################
! Coded on 2017.06.02
subroutine sharerespdata(resp,ip_from)
implicit none
integer(4),    intent(in)    :: ip_from
type(respdata),intent(inout) :: resp
integer(4)                   :: nobs,errno

nobs = resp%nobs
call MPI_BCAST(resp%fpobsamp(1),  nobs, MPI_REAL8, ip_from, MPI_COMM_WORLD,errno)
call MPI_BCAST(resp%fpobsphase(1),nobs, MPI_REAL8, ip_from, MPI_COMM_WORLD,errno)
call MPI_BCAST(resp%fsobsamp(1),  nobs, MPI_REAL8, ip_from, MPI_COMM_WORLD,errno)
call MPI_BCAST(resp%fsobsphase(1),nobs, MPI_REAL8, ip_from, MPI_COMM_WORLD,errno)
call MPI_BCAST(resp%ftobsamp(1),  nobs, MPI_REAL8, ip_from, MPI_COMM_WORLD,errno)
call MPI_BCAST(resp%ftobsphase(1),nobs, MPI_REAL8, ip_from, MPI_COMM_WORLD,errno)

return
end
!############################################################
! Coded on 2022.01.02
subroutine shareimpdata(timp,ip_from)
  implicit none
  integer(4),    intent(in)    :: ip_from
  type(respmt),  intent(inout) :: timp
  integer(4)                   :: nobs,errno
  
  nobs = timp%nobs
  call MPI_BCAST(timp%zxx(1),  nobs, MPI_COMPLEX16, ip_from, MPI_COMM_WORLD,errno)
  call MPI_BCAST(timp%zxy(1),  nobs, MPI_COMPLEX16, ip_from, MPI_COMM_WORLD,errno)
  call MPI_BCAST(timp%zyx(1),  nobs, MPI_COMPLEX16, ip_from, MPI_COMM_WORLD,errno)
  call MPI_BCAST(timp%zyy(1),  nobs, MPI_COMPLEX16, ip_from, MPI_COMM_WORLD,errno)
  call MPI_BCAST(timp%rhoxx(1),nobs, MPI_REAL8,     ip_from, MPI_COMM_WORLD,errno)
  call MPI_BCAST(timp%rhoxy(1),nobs, MPI_REAL8,     ip_from, MPI_COMM_WORLD,errno)
  call MPI_BCAST(timp%rhoyx(1),nobs, MPI_REAL8,     ip_from, MPI_COMM_WORLD,errno)
  call MPI_BCAST(timp%rhoyy(1),nobs, MPI_REAL8,     ip_from, MPI_COMM_WORLD,errno)
  call MPI_BCAST(timp%phaxx(1),nobs, MPI_REAL8,     ip_from, MPI_COMM_WORLD,errno)
  call MPI_BCAST(timp%phaxy(1),nobs, MPI_REAL8,     ip_from, MPI_COMM_WORLD,errno)
  call MPI_BCAST(timp%phayx(1),nobs, MPI_REAL8,     ip_from, MPI_COMM_WORLD,errno)
  call MPI_BCAST(timp%phayy(1),nobs, MPI_REAL8,     ip_from, MPI_COMM_WORLD,errno)
  
  return
  end
  !############################################################
! Coded on 2017.05.31
subroutine sharemeshline(g_mesh,g_line,ip)
implicit none
integer(4),     intent(in)    :: ip
type(mesh),     intent(inout) :: g_mesh
type(line_info),intent(inout) :: g_line

call sharemesh(g_mesh,ip)
call shareline(g_line,ip)

if (ip .eq. 0 ) write(*,'(a)') " ### SHAREMESHLINE END!! ###" ! 2020.09.17

return
end
!############################################################
! Coded on 2017.05.31
subroutine shareline(g_line,ip)
implicit none
integer(4),intent(in)    :: ip
type(line_info),intent(inout) :: g_line
integer(4) :: nline,node,ntet,errno

 call MPI_BCAST(g_line%nline, 1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_line%node,  1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_line%ntet,  1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_line%ntri,  1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)

 node  = g_line%node
 nline = g_line%nline
 ntet  = g_line%ntet
 if (ip .ne. 0) then
  allocate(g_line%line(2,nline))
  allocate(g_line%n6line(ntet,6))
 end if

!# Default
 call MPI_BCAST(g_line%line(1,1),    2*nline, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_line%n6line(1,1),   ntet*6, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)

!# after MKLINE is conducted
! if ( node .gt. 0 ) then
! if ( ip .ne. 0 ) then
!    allocate(g_line%line_stack(0:node))
!    allocate(g_line%line_item(nline))
!  end if
!  call MPI_BCAST(g_line%line_stack(0), node+1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
!  call MPI_BCAST(g_line%line_item(1),   nline, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
! end if

if ( ip .eq. 0) write(*,'(a)') " ### SHARELINE     END!! ###" ! 2020.09.18

return
end
!############################################################
! Coded on 2017.05.31
subroutine sharemesh(g_mesh,ip)
implicit none
integer(4),intent(in)    :: ip
type(mesh),intent(inout) :: g_mesh
integer(4) :: errno,node,ntet,ntri,nlin,npoi

 call MPI_BCAST(g_mesh%meshname, 50, MPI_CHAR,0,MPI_COMM_WORLD,errno)

 call MPI_BCAST(g_mesh%node, 1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_mesh%ntet, 1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_mesh%ntri, 1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_mesh%nlin, 1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_mesh%npoi, 1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_mesh%icoordinateflag, 1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 node = g_mesh%node
 ntet = g_mesh%ntet
 ntri = g_mesh%ntri
 nlin = g_mesh%nlin
 npoi = g_mesh%npoi

 if ( ip .ne. 0) then
  allocate(g_mesh%xyz(3,node)) ! default
  if (g_mesh%icoordinateflag .eq. 2 )  allocate(g_mesh%lonlatalt(3,node))
  if (g_mesh%icoordinateflag .eq. 3 )  allocate(g_mesh%xyzspherical(3,node))
  allocate(g_mesh%n4(ntet,4),g_mesh%n4flag(ntet,2))
  allocate(g_mesh%n3(ntri,3),g_mesh%n3flag(ntri,2))
  allocate(g_mesh%n2(nlin,2),g_mesh%n2flag(nlin,2))
  allocate(g_mesh%n1(npoi,1),g_mesh%n1flag(npoi,2))
 end if

 !#
 call MPI_BCAST(g_mesh%lonorigin,     1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_mesh%latorigin,     1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_mesh%xyz(1,1),          node*3, MPI_REAL8,0,MPI_COMM_WORLD,errno)

 if  (g_mesh%icoordinateflag .eq. 2 ) then
  call MPI_BCAST(g_mesh%lonlatalt(1,1),   node*3, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 end if
 if  (g_mesh%icoordinateflag .eq. 3 ) then
  call MPI_BCAST(g_mesh%xyzspherical(1,1),node*3, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 end if

 !#
 if ( ntet .gt. 0) then
  call MPI_BCAST(g_mesh%n4(1,1),    ntet*4, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_mesh%n4flag(1,1),ntet*2, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 end if
 if ( ntri .gt. 0 ) then
  call MPI_BCAST(g_mesh%n3(1,1),    ntri*3, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_mesh%n3flag(1,1),ntri*2, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 end if
 if ( nlin .gt. 0 ) then
  call MPI_BCAST(g_mesh%n2(1,1),    nlin*2, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_mesh%n2flag(1,1),nlin*2, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 end if
 !#
 if ( npoi .gt. 0 ) then
  call MPI_BCAST(g_mesh%n1(1,1),    npoi*1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_mesh%n1flag(1,1),npoi*2, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 end if
 !# the followings are not shred for the present on 2017.05.31
 !integer(4) :: nmodel
 !integer(4),allocatable,dimension(:) :: iele2model ! i-th ele corresponds to iele2model(i)-th model field
 !real(8),allocatable,dimension(:) :: cmodel ! conductivity for i-th field

if ( ip .eq. 0) write(*,'(a)') " ### SHAREMESH     END!! ###" ! 2020.09.17

return
end
!############################################################
subroutine shareforward(g_param,sparam,g_cond,ip)
implicit none
integer(4),intent(in) :: ip
type(param_source), intent(inout) :: sparam
type(param_forward),intent(inout) :: g_param
type(param_cond),   intent(inout) :: g_cond

call sharecond(g_cond,ip)
call sharesource(sparam,ip) ! 2017.09.04
call sharefparam(g_param,ip)

if (ip .eq. 0 ) write(*,'(a)') " ### SHAREFORWARD  END!! ###" ! 2020.09.17

return
end

!##################################################### shareparam
subroutine sharefparam(g_param,ip)
implicit none
type(param_forward),intent(inout) :: g_param
integer(4),intent(in) :: ip
integer(4) :: errno,nfreq,nobsr,nobs
integer(4) :: nfile ! 2017.12.13

!#[3]## broadcast
 call MPI_BCAST(g_param%itopoflag,      1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno) ! 2017.12.13
 if ( g_param%itopoflag .eq. 1) then      ! 2017.12.13 when topography is considered
  call MPI_BCAST(g_param%nfile,         1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno) ! 2017.12.13
  nfile = g_param%nfile                   ! 2017.12.13
  if ( ip .ne. 0 ) then                   ! 2017.12.13
   allocate(g_param%topofile(nfile))      ! 2017.12.13
   allocate(g_param%lonlatshift(2,nfile)) ! 2017.12.13
  end if                                  ! 2017.12.13
  call MPI_BCAST(g_param%topofile,  nfile*50, MPI_CHAR,   0,MPI_COMM_WORLD,errno)! 2017.12.13
  call MPI_BCAST(g_param%lonlatshift,nfile*2, MPI_REAL8,  0,MPI_COMM_WORLD,errno)! 2017.12.13
 end if                                  ! 2017.09.29
 call MPI_BCAST(g_param%g_meshfile,    50, MPI_CHAR,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%surface_id_ground, 1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%z_meshfile,    50, MPI_CHAR,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%g_lineinfofile,50, MPI_CHAR,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%outputfolder,  50, MPI_CHAR,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%header2d,      50, MPI_CHAR,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%header3d,      50, MPI_CHAR,0,MPI_COMM_WORLD,errno)

 call MPI_BCAST(g_param%nfreq,          1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 nfreq = g_param%nfreq
 if ( ip .ne. 0 ) allocate(g_param%freq(nfreq))
 call MPI_BCAST(g_param%freq(1),     nfreq, MPI_REAL8,0,MPI_COMM_WORLD,errno)

 call MPI_BCAST(g_param%nobs,          1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%lonlatflag,    1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%wlon,     1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%elon,     1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%slat,     1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%nlat,     1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%lonorigin,1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%latorigin,1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%lenout,   1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%upzin,    1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%downzin,  1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%zmax,     1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%zmin,     1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%sizein,   1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%sizebo,   1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%sigma_obs,1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%A_obs,    1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%sigma_src,1, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%A_src,    1, MPI_REAL8,0,MPI_COMM_WORLD,errno)

 call MPI_BCAST(g_param%nobsr,    1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 nobsr = g_param%nobsr
 if ( ip .ne. 0) then
  allocate(g_param%xyz_r(3,nobsr))
  allocate(g_param%A_r(nobsr)    )
  allocate(g_param%sigma_r(nobsr))
 end if

 call MPI_BCAST(g_param%xyz_r,nobsr*3, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%A_r,    nobsr, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%sigma_r,nobsr, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%xbound,     4, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%ybound,     4, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%zbound,     4, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%UTM,        3, MPI_CHAR, 0,MPI_COMM_WORLD,errno)

 nobs = g_param%nobs
 if ( ip .ne. 0) then
  allocate(g_param%lonlataltobs(3,nobs))
  allocate(g_param%xyzobs(3,nobs)      )
  allocate(g_param%obsname(nobs)       )
 end if
 call MPI_BCAST(g_param%lonlataltobs,3*nobs, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%xyzobs,      3*nobs, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%obsname,    50*nobs, MPI_CHAR, 0,MPI_COMM_WORLD,errno)

 call MPI_BCAST(g_param%xyzminmax,   6, MPI_REAL8,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%zorigin,     1, MPI_REAL8,0,MPI_COMM_WORLD,errno)

 call MPI_BCAST(g_param%condfile,   50, MPI_CHAR,    0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%condflag,    1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)


if (ip .eq. 0 ) write(*,'(a)') " ### SHAREFPARAM   END!! ###" ! 2020.09.17
return
end
!##################################################### sharesource
!# modified on 2017.09.04 for multiple source
subroutine sharesource(sparam,ip)
implicit none
type(param_source),intent(inout) :: sparam
integer(4),intent(in) :: ip
integer(4) :: errno,nsource

 call MPI_BCAST(sparam%lonlatflag, 1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(sparam%nsource,    1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)!2017.07.11
 nsource = sparam%nsource ! 2017.07.12

 if ( ip .ne. 0) then      ! 2017.07.11
  allocate(sparam%sourcename(nsource) ) ! 2017.07.11
  allocate(sparam%xs1(3,nsource)      ) ! 2017.07.11
  allocate(sparam%xs2(3,nsource)      ) ! 2017.07.11
  allocate(sparam%lonlats1(2,nsource) ) ! 2017.07.11
  allocate(sparam%lonlats2(2,nsource) ) ! 2017.07.11
 end if

 call MPI_BCAST(sparam%xs1,        3*nsource, MPI_REAL8,   0,MPI_COMM_WORLD,errno)! 2017.07.11
 call MPI_BCAST(sparam%xs2,        3*nsource, MPI_REAL8,   0,MPI_COMM_WORLD,errno)! 2017.07.11
 call MPI_BCAST(sparam%lonlats1,   2*nsource, MPI_REAL8,   0,MPI_COMM_WORLD,errno)! 2017.07.11
 call MPI_BCAST(sparam%lonlats2,   2*nsource, MPI_REAL8,   0,MPI_COMM_WORLD,errno)! 2017.07.11
 call MPI_BCAST(sparam%I,          1,         MPI_REAL8,   0,MPI_COMM_WORLD,errno)

!if (ip .eq. 1 ) then ! commented out on 2020.09.17
! write(*,*) "sparam%xs1",sparam%xs1
! write(*,*) "sparam%xs2",sparam%xs2
! write(*,*) "sparam%lonlats1",sparam%lonlats1
! write(*,*) "sparam%lonlats2",sparam%lonlats2
! write(*,*) "sparam%I",sparam%I
! end if

return
end


!##################################################### sharecond
subroutine sharecond(g_cond,ip)
implicit none
type(param_cond),intent(inout) :: g_cond
integer(4),intent(in) :: ip
integer(4) :: nphys2,errno,ibyte
integer(4) :: nvolume ! 2017.12.13

!#[1]## share nphys2
!write(*,*) "sharecond start!"

if ( ip .eq. 0) then
 nphys2 = g_cond%nphys2
 write(*,*) "ip=0 nphys2=",nphys2
end if
call MPI_BCAST(nphys2,1,MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
g_cond%nphys2 = nphys2

!write(*,*) "nphys2=",nphys2,"ip=",ip

!#[2]## allocate except ip = 0
if ( ip .ne. 0 ) then
 allocate( g_cond%index(nphys2))
 allocate( g_cond%sigma(nphys2))
 allocate( g_cond%rho(nphys2))
end if

!#[3]## broadcast
call MPI_BCAST(g_cond%index(1), nphys2, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
call MPI_BCAST(g_cond%rho(1),   nphys2, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
call MPI_BCAST(g_cond%sigma(1), nphys2, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
call MPI_BCAST(g_cond%nphys1,        1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
call MPI_BCAST(g_cond%ntet,          1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
call MPI_BCAST(g_cond%condflag,      1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
call MPI_BCAST(g_cond%sigma_air,     1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
if ( g_cond%condflag .eq. 0 ) then                     ! 2017.12.13
 call MPI_BCAST(g_cond%nvolume,      1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno) ! 2017.12.13
 nvolume = g_cond%nvolume                              ! 2017.12.13
 if ( ip .ne. 0 ) allocate(g_cond%sigma_land(nvolume)) ! 2017.12.13
 call MPI_BCAST(g_cond%sigma_land, nvolume, MPI_REAL8,   0,MPI_COMM_WORLD,errno) ! 2017.12.13
end if                                                 ! 2017.12.13
if ( g_cond%condflag .eq. 1 ) then                     ! 2017.12.13
 call MPI_BCAST(g_cond%condfile,     50, MPI_CHAR ,   0,MPI_COMM_WORLD,errno)
end if                                                 ! 2017.12.13

return
end

end module
