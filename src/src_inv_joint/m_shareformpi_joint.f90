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
use shareformpi_mt ! 2025.09.17
use mpi
implicit none

contains
!############################################################
subroutine sharejointinv(g_param,sparam,g_cond,g_mesh,g_line,g_param_joint,g_model,ip)
 !# modified on 2017.08.31 for multisource inversion
 !# Coded on 2017.06.07
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
  CALL SHAREINVPARAJOINT(g_param_joint,ip)    ! 2017.08.31 see below
  CALL SHAREMODEL(g_model,ip)                 ! see m_shareformpi.f90
  if (ip .eq. 0) write(*,'(a)') " ### SHAREJOINTINV  END!! ###" ! 2026.03.03
 return
 end

!############################################################
subroutine SHAREMODEL(h_model,ip)
 !# Coded on 2017.06.05
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
subroutine SHAREINVPARAJOINT(g_param_joint,ip) ! 2017.08.31
 !# Modified on 2017.08.31 for multisource inversion
 !# Coded on 2017.06.08
 implicit none
 integer(4),             intent(in)    :: ip
 type(param_joint),      intent(inout) :: g_param_joint   ! 2017.08.31
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
 !# Tipper
 call MPI_BCAST(g_param_joint%iflag_tipper,  1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)! 2023.12.23
 
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
subroutine sharerealcrsmatrix(rmat,ip_from,ip)
 ! Coded on 2017.06.05
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
subroutine sharecomplexcrsmatrix(rmat,ip_from,ip)
 ! Coded on 2022.01.05
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
subroutine shareface(g_face,ip)
 ! Coded on 2017.06.02
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
subroutine sharesource(sparam,ip)
 !# modified on 2017.09.04 for multiple source
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



end module
