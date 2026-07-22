! Coded on 2017.05.50
! to share structures between processes for MPI parallelization
module shareformpi_mdm
use param
use mesh_type
use line_type
use outresp      ! 2017.06.02
use face_type    ! 2017.06.02
use matrix       ! 2017.06.05
use param_mdminv ! 2017.10.31
use modelpart    ! 2017.06.05
implicit none
include 'mpif.h'

contains
!############################################################
!# Coded on 2017.06.07
subroutine shareinv(g_param,sparam,g_cond,g_mesh,g_line,g_face,g_param_inv,g_model,ip)
implicit none
integer(4),           intent(in)    :: ip
type(param_forward),  intent(inout) :: g_param
type(param_source),   intent(inout) :: sparam
type(param_cond),     intent(inout) :: g_cond
type(mesh),           intent(inout) :: g_mesh
type(line_info),      intent(inout) :: g_line
type(face_info),      intent(inout) :: g_face
type(param_inversion),intent(inout) :: g_param_inv
type(model),          intent(inout) :: g_model

  CALL SHAREFORWARD(g_param,sparam,g_cond,ip) ! see m_shareformpi.f90
  CALL SHAREMESHLINE(g_mesh,g_line,ip)        ! see m_shareformpi.f90
  CALL SHAREFACE(g_face,ip)                   ! see m_shareformpi.f90
  CALL SHAREINVPARA(g_param_inv,ip)           ! see m_shareformpi.f90, 2017.07.13 for multisrc
  CALL SHAREMODEL(g_model,ip)                 ! see m_shareformpi.f90

return
end

!############################################################
!# Coded on 2017.06.05
subroutine SHAREMODEL(h_model,ip)
implicit none
integer(4), intent(in)    :: ip
type(model),intent(inout) :: h_model
integer(4) :: nmodel,nphys1,nphys2,errno, ip_from

 call MPI_BCAST(h_model%nmodel,  1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(h_model%nphys1,  1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(h_model%nphys2,  1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 nmodel = h_model%nmodel
 nphys1 = h_model%nphys1
 nphys2 = h_model%nphys2

 if ( ip .ne. 0) then
  if ( .not. allocated(h_model%index)       ) allocate( h_model%index(nphys2)       )
  if ( .not. allocated(h_model%ele2model)   ) allocate( h_model%ele2model(nphys2)   )
  if ( .not. allocated(h_model%rho_model)   ) allocate( h_model%rho_model(nmodel)   )
  if ( .not. allocated(h_model%logrho_model)) allocate( h_model%logrho_model(nmodel))
 end if

 call MPI_BCAST(h_model%index(1),     nphys2, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(h_model%ele2model(1), nphys2, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(h_model%rho_model,    nmodel, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(h_model%logrho_model, nmodel, MPI_REAL8,   0,MPI_COMM_WORLD,errno)

 ip_from = 0
 call sharerealcrsmatrix(h_model%model2ele,ip_from,ip) ! see below

 if (ip .eq. 0) write(*,*) "### SHAREMODEL END!! ###"

return
end

!############################################################
!# modified on 2017.07.13 for multiple sources
!# Coded on 2017.06.05
subroutine SHAREINVPARA(g_param_inv,ip)
implicit none
integer(4),           intent(in)    :: ip
type(param_inversion),intent(inout) :: g_param_inv
integer(4)                          :: i,j,nobs,nfreq,errno
integer(4)                          :: nsr_inv, nobs_s1,nobs_s2 ! 2017.07.13

 call MPI_BCAST(g_param_inv%nobs,    1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param_inv%nfreq,   1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param_inv%g_faceinfofile, 50, MPI_CHAR,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param_inv%g_initcondfile, 50, MPI_CHAR,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param_inv%nsr_inv, 1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno) !2017.07.13
 nobs     = g_param_inv%nobs
 nfreq    = g_param_inv%nfreq
 nsr_inv  = g_param_inv%nsr_inv ! 2017.07.13

 if ( ip .ne. 0) then
  allocate(g_param_inv%srcindex(nsr_inv) ) ! 2017.07.13
  allocate(g_param_inv%obsinfo1( nsr_inv) ) ! 2017.10.31
  allocate(g_param_inv%obsinfo2( nsr_inv) ) ! 2017.10.31
 end if
  call MPI_BCAST(g_param_inv%srcindex,nsr_inv,MPI_INTEGER4, 0,MPI_COMM_WORLD,errno) !2017.07.14

 do j=1,nsr_inv       ! 2017.07.13
  call MPI_BCAST(g_param_inv%obsinfo1(j)%nobs_s,1,MPI_INTEGER4, 0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_param_inv%obsinfo2(j)%nobs_s,1,MPI_INTEGER4, 0,MPI_COMM_WORLD,errno)
  nobs_s1 = g_param_inv%obsinfo1(j)%nobs_s ! 2017.07.13
  nobs_s2 = g_param_inv%obsinfo2(j)%nobs_s ! 2017.07.13
  if ( ip .ne. 0) then ! 2017.07.13
   allocate(g_param_inv%obsinfo1(j)%obsfile( nobs_s1) ) ! 2017.10.31
   allocate(g_param_inv%obsinfo1(j)%obsindex(nobs_s1))  ! 2017.10.31
   allocate(g_param_inv%obsinfo2(j)%obsfile( nobs_s2) ) ! 2017.10.31
   allocate(g_param_inv%obsinfo2(j)%obsindex(nobs_s2))  ! 2017.10.31
  end if                ! 2017.07.13
  do i=1,nobs_s1        ! 2017.10.31
   call MPI_BCAST(g_param_inv%obsinfo1(j)%obsfile(i),50,MPI_CHAR, 0,MPI_COMM_WORLD,errno)
   call MPI_BCAST(g_param_inv%obsinfo1(j)%obsindex(i),1,MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  end do                ! 2017.07.13
  do i=1,nobs_s2        ! 2017.10.31
   call MPI_BCAST(g_param_inv%obsinfo2(j)%obsfile(i),50,MPI_CHAR, 0,MPI_COMM_WORLD,errno)
   call MPI_BCAST(g_param_inv%obsinfo2(j)%obsindex(i),1,MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  end do         ! 2017.07.13
 end do          ! 2017.07.13

 call MPI_BCAST(g_param_inv%errorfloor,     1,MPI_REAL8,0,MPI_COMM_WORLD,errno)

 if ( ip .eq. 0 ) write(*,*) "### SHAREINVPARA END!! ###"
return
end

!############################################################
! Coded on 2017.06.05
subroutine sharerealcrsmatrix(rmat,ip_from,ip)
implicit none
integer(4),           intent(in)    :: ip_from,ip
type(real_crs_matrix),intent(inout) :: rmat
integer(4) :: nrow,ncolm,ntot,errno
type(real_crs_matrix) :: rmat2

 call MPI_BCAST(rmat%nrow,  1, MPI_INTEGER4,ip_from,MPI_COMM_WORLD,errno)
 call MPI_BCAST(rmat%ncolm, 1, MPI_INTEGER4,ip_from,MPI_COMM_WORLD,errno)
 call MPI_BCAST(rmat%ntot,  1, MPI_INTEGER4,ip_from,MPI_COMM_WORLD,errno)
 nrow  = rmat%nrow
 ncolm = rmat%ncolm
 ntot  = rmat%ntot

 allocate(rmat2%stack(0:nrow))
 allocate(rmat2%item(ntot)   )
 allocate(rmat2%val(ntot)    )

 if ( ip .eq. ip_from ) rmat2 = rmat
 if ( ip .ne. ip_from ) then         ! 2017.07.17
  rmat2%nrow  = nrow                 ! 2017.07.17
  rmat2%ncolm = ncolm                ! 2017.07.17
  rmat2%ntot  = ntot                 ! 2017.07.17
 end if                              ! 2017.07.17

 call MPI_BCAST(rmat2%stack(0),nrow+1, MPI_INTEGER4,ip_from,MPI_COMM_WORLD,errno)
 call MPI_BCAST(rmat2%item(1),   ntot, MPI_INTEGER4,ip_from,MPI_COMM_WORLD,errno)
 call MPI_BCAST(rmat2%val(1),    ntot, MPI_REAL8,   ip_from,MPI_COMM_WORLD,errno)

 if ( ip .ne. ip_from ) rmat = rmat2
! if ( ip .eq. 0 ) write(*,*) "rmat%ntot=",rmat%ntot

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

if ( ip .eq. 0 ) write(*,*) "### SHAREFACE END!! ###"
return
end

!############################################################
! Coded on 2017.06.02
subroutine sharerespdata(resp,ip_from)
implicit none
integer(4),intent(in) :: ip_from
type(respdata),intent(inout) :: resp
integer(4) :: nobs,errno

nobs = resp%nobs
call MPI_BCAST(resp%fpobsamp(1),  nobs,MPI_REAL8,ip_from,MPI_COMM_WORLD,errno)
call MPI_BCAST(resp%fpobsphase(1),nobs,MPI_REAL8,ip_from,MPI_COMM_WORLD,errno)
call MPI_BCAST(resp%fsobsamp(1),  nobs,MPI_REAL8,ip_from,MPI_COMM_WORLD,errno)
call MPI_BCAST(resp%fsobsphase(1), nobs,MPI_REAL8,ip_from,MPI_COMM_WORLD,errno)
call MPI_BCAST(resp%ftobsamp(1),  nobs,MPI_REAL8,ip_from,MPI_COMM_WORLD,errno)
call MPI_BCAST(resp%ftobsphase(1),nobs,MPI_REAL8,ip_from,MPI_COMM_WORLD,errno)

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

if (ip .eq. 0 ) write(*,*) "### SHAREMESHLINE END!! ###"

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
 if ( node .gt. 0 ) then
  if ( ip .ne. 0 ) then
    allocate(g_line%line_stack(0:node))
    allocate(g_line%line_item(nline))
  end if
  call MPI_BCAST(g_line%line_stack(0), node+1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
  call MPI_BCAST(g_line%line_item(1),   nline, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 end if

if ( ip .eq. 0) write(*,*) "### SHARELINE END!! ###"

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

if ( ip .eq. 0) write(*,*) "### SHAREMESH END!! ###"

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
call sharesource(sparam,ip) ! 2017.07.12
call sharefparam(g_param,ip)

if (ip .eq. 0 ) write(*,*) "### SHAREFORWARD END!! ###"

return
end

!##################################################### shareparam
subroutine sharefparam(g_param,ip)
implicit none
type(param_forward),intent(inout) :: g_param
integer(4),intent(in) :: ip
integer(4) :: errno,nfreq,nobsr,nobs
integer(4) :: nfile  ! 2017.09.29

!#[3]## broadcast
 call MPI_BCAST(g_param%itopoflag,      1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno) ! 2017.09.29
 if ( g_param%itopoflag .eq. 1) then      ! 2017.09.29 when topography is considered
  call MPI_BCAST(g_param%nfile,          1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno) ! 2017.09.29
  nfile = g_param%nfile                   ! 2017.09.29
  if ( ip .ne. 0 ) then                   ! 2017.09.29
   allocate(g_param%topofile(nfile))      ! 2017.09.29
   allocate(g_param%lonlatshift(2,nfile)) ! 2017.09.29
  end if                                  ! 2017.09.29
  call MPI_BCAST(g_param%topofile,  nfile*50, MPI_CHAR,   0,MPI_COMM_WORLD,errno)! 2017.09.29
  call MPI_BCAST(g_param%lonlatshift,nfile*2, MPI_REAL8,  0,MPI_COMM_WORLD,errno)! 2017.09.29
 end if                                  ! 2017.09.29
 call MPI_BCAST(g_param%g_meshfile,    50, MPI_CHAR,    0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%surface_id_ground, 1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%z_meshfile,    50, MPI_CHAR,    0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%g_lineinfofile,50, MPI_CHAR,    0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%outputfolder,  50, MPI_CHAR,    0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%header2d,      50, MPI_CHAR,    0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%header3d,      50, MPI_CHAR,    0,MPI_COMM_WORLD,errno)

 call MPI_BCAST(g_param%nfreq,          1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 nfreq = g_param%nfreq
 if ( ip .ne. 0 ) allocate(g_param%freq(nfreq))
 call MPI_BCAST(g_param%freq(1),     nfreq, MPI_REAL8,0,MPI_COMM_WORLD,errno)

 call MPI_BCAST(g_param%nobs,        1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%lonlatflag,  1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%wlon,        1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%elon,        1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%slat,        1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%nlat,        1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%lonorigin,   1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%latorigin,   1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%lenout,      1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%upzin,       1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%downzin,     1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%zmax,        1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%zmin,        1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%sizein,      1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%sizebo,      1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%sigma_obs,   1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%A_obs,       1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%sigma_src,   1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)
 call MPI_BCAST(g_param%A_src,       1, MPI_REAL8,   0,MPI_COMM_WORLD,errno)

 call MPI_BCAST(g_param%ixyflag,     1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno) ! 2017.10.12
 call MPI_BCAST(g_param%nx,          1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno) ! 2017.10.12
 call MPI_BCAST(g_param%ny,          1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno) ! 2017.10.12
 call MPI_BCAST(g_param%xyfilehead, 50, MPI_CHAR,    0,MPI_COMM_WORLD,errno) ! 2017.10.12

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


if (ip .eq. 0 ) write(*,*) "### SHAREFPARAM END!! ###"
return
end
!##################################################### sharesource
!# modified on 2017.07.11 for multiple source
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

return
end

!##################################################### sharecond
subroutine sharecond(g_cond,ip)
implicit none
type(param_cond),intent(inout) :: g_cond
integer(4),intent(in) :: ip
integer(4) :: nphys2,errno,ibyte
integer(4) :: nvolume ! 2017.09.29

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
if ( g_cond%condflag .eq. 0 ) then                     ! 2017.09.29
 call MPI_BCAST(g_cond%nvolume,      1, MPI_INTEGER4,0,MPI_COMM_WORLD,errno) ! 2017.09.29
 nvolume = g_cond%nvolume                              ! 2017.09.29
 if ( ip .ne. 0 ) allocate(g_cond%sigma_land(nvolume)) ! 2017.09.29
 call MPI_BCAST(g_cond%sigma_land, nvolume, MPI_REAL8,   0,MPI_COMM_WORLD,errno) ! 2017.09.29
end if                                                 ! 2017.09.29
if ( g_cond%condflag .eq. 1 ) then                     ! 2017.09.29
 call MPI_BCAST(g_cond%condfile,     50, MPI_CHAR ,   0,MPI_COMM_WORLD,errno)
end if                                                 ! 2017.09.29

return
end

end module
