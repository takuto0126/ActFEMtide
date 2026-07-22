! coded on March 3, 2016
! to solve the induction equation in terms of the magnetic field
! using edge-FEM
!
program n_ebfem_bxyz
use mesh_type
use param
use matrix
use line_type  ! see m_line_type.f90
use iccg_var_takuto
use outresp
use obs_type     ! see m_obs_type.f90 added on 2016.10.27
implicit none
type(param_forward)                         :: g_param
type(param_source)                          :: sparam
type(param_cond)                            :: g_cond     ! see m_param.f90
type(mesh)                                  :: g_mesh     ! see m_mesh_type.f90
type(mesh)                                  :: h_mesh     ! topography file
type(line_info)                             :: g_line     ! see m_line_type.f90
type(global_matrix)                         :: A          ! see m_iccg_var_takuto.f90
type(real_crs_matrix)                       :: coeffobs(2,3) !m_matrix.f90; 1:edge, 2:face
type(obs_info)                              :: obs_xy     ! added on 2016.10.27
type(respdata),allocatable,dimension(:,:,:) :: resp_xy    ! 2018.02.22
type(respdata),allocatable,dimension(:,:,:) :: resp5      ! 2018.02.22
integer(4),    allocatable,dimension(:,:)   :: table_dof
complex(8),    allocatable,dimension(:,:)   :: fp !(nline,nsr)[nT * km] 2018.02.22
complex(8),    allocatable,dimension(:,:)   :: fs !(nline,nsr)[nT * km] 2018.02.22
integer(4)                                  :: nline, ntet, nsr,nfreq ! 2018.02.22
integer(4),parameter                        :: icomp = 2 ! 1 for b field, 2 for e field
integer(4)                                  :: i,j,ip,iresfile,dofn
real(8)                                     :: omega, freq ! 2018.02.22
integer(4)                                  :: ixyflag     ! 2017.10.11
!##

!#[0]## read parameters
  CALL READPARAM(g_param,sparam,g_cond)

!#[1]## Mesh READ
   CALL READMESH_TOTAL(g_mesh,g_param%g_meshfile)
   CALL READMESH_TOTAL(h_mesh,g_param%z_meshfile)
   CALL GENXYZMINMAX(g_mesh,g_param)             ! generate xyzminmax 2017.10.12
   CALL PREPZSRCOBS(h_mesh,g_param,sparam)
   if (g_cond%condflag .eq. 1) then !  "1" means conductivity file is given
    g_cond%ntet   = g_mesh%ntet
    g_cond%nphys1 = g_mesh%ntet - g_cond%nphys2
   end if
   CALL SETNPHYS1INDEX2COND(g_mesh,g_cond) ! added on 2017.05.31; see below

!#[2]## Line information
  CALL READLINE(g_param%g_lineinfofile,g_line) ! see m_line_type.f90

!#[3]## Prepare primary magnetic field
!  CALL PREPXYS(g_param,sparam)         ! commented out on 2017.05.14
!  CALL PREPZOBS(g_param,g_mesh,sparam) ! calculate z of obs and source
  nline   = g_line%nline
  ntet    = g_line%ntet
  nsr     = sparam%nsource     ! 2018.02.22
  ixyflag = g_param%ixyflag    ! 2018.02.22
  nfreq   = g_param%nfreq      ! 2018.02.22

!#[4]## allocate respdata and open output files for each observatory

!#[5]## allocate global matrix
  dofn=1; ip=0
  allocate( table_dof(nline,dofn))
  CALL SET_TABLE_DOF(dofn,nline,table_dof,nline)
  CALL set_iccg_var7_dofn(dofn,6,nline,nline,nline,ntet,g_line%n6line,table_dof,A,ip)

!#[6]## prepare coefficients
  CALL PREPOBSCOEFF(g_param,g_mesh,g_line,coeffobs) ! for x,y,z component

!#[7]## set resp5
  allocate( resp5(5,nsr,nfreq)) ! 2018.02.22
  CALL ALLOCATERESP(g_param%nobs,nsr,resp5,ip,nfreq) ! 2018.02.22, see below

!#[3]## Prepare coefficients for values at observatories
  if ( ixyflag .eq. 1 .and. ip .lt. nfreq ) then
   CALL PREPOBSCOEFF_XY(g_param,g_mesh,h_mesh,g_line,obs_xy) ! 2018.02.22
   allocate( resp_xy(5,nsr,nfreq) )                          ! 2018.02.22
   CALL ALLOCATERESP(obs_xy%nobs,nsr,resp_xy,ip,nfreq)       ! 2018.02.22
  end if
  if ( ixyflag .eq. 2) then ! 2017.10.11
   write(*,*) "GEGEGE not suported in the current version"   ! 2017.10.12
   stop
  end if

 allocate( fs(nline,nsr) ) ! 2018.02.22 for multiple source, fp deleted on 2021.09.15

!==================================================================  freq loop
do i=1,nfreq

  freq  = g_param%freq(i)   ! 2018.02.22
  omega = 2.d0*(4.d0*datan(1.d0))*freq
  write(*,*) "frequency =",g_param%freq(i),"[Hz]"

 !#[7]## conduct forward calculation with model
 CALL forward_bxyz(A,g_mesh,g_line,nline,nsr,fs,g_param%freq(i),sparam,g_param,g_cond,ip)

 !#[8]## calculate response at every observation point
    !# calculate bx,by,bz
    CALL CALOBSEBCOMP(fs,nline,nsr,omega,coeffobs,resp5(:,:,i))!fp deleted 21.09.15
    ! E field at xy plane
    if ( ixyflag .eq. 1 ) then ! 2017.10.12
     CALL CALOBSEBCOMP(fs,nline,nsr,omega,obs_xy%coeff,resp_xy(:,:,i))!fp deleted 21.09.15
     CALL OUTFREQFILES2(freq,nsr,resp_xy(:,:,i),g_param,obs_xy)!m_outresp.f90
    end if

!#[9]## output resp to frequency file
    j=1
    if (.false.) CALL OUTFREQ(freq,g_param,resp5(:,j,i)) !2017.07.11, see below

end do ! freq loop end

  !#[8]## output resp to obs file
   CALL OUTOBSFILESFWD(g_param,sparam,nsr,resp5,nfreq) !2017.07.11 m_outresp.f90

end program n_ebfem_bxyz
!#############################################
!# copied from n_ebfem_bxyz_mpi.f90 on 2018.02.22
!# modified on 2017.07.11 to include multiple sources
!# coded on 2017.05.31
subroutine CALOBSEBCOMP(fs,nline,nsr,omega,coeffobs,resp5) ! fp,fs -> ft 2021.09.14
use matrix
use outresp
implicit none
real(8),              intent(in)    :: omega
integer(4),           intent(in)    :: nline, nsr
complex(8),           intent(inout) :: fs(nline,nsr) ! 2017.07.11
type(real_crs_matrix),intent(in)    :: coeffobs(2,3)
type(respdata),       intent(inout) :: resp5(5,nsr)                 !2017.07.11
integer(4)                          :: i  ! 2017.07.11

do i=1, nsr ! 2017.07.11
 CALL CALOBSRESP(fs(:,i),nline,coeffobs(2,1),resp5(1,i)  ) !bx,fp deleted 2021.09.15
 CALL CALOBSRESP(fs(:,i),nline,coeffobs(2,2),resp5(2,i)  ) !by,fp deleted 2021.09.15
 CALL CALOBSRESP(fs(:,i),nline,coeffobs(2,3),resp5(3,i)  ) !bz,fp deleted 2021.09.15
 fs(:,i)= - (0.d0,1.d0)*omega*fs(:,i) !E= -i*omega*A
 CALL CALOBSRESP(fs(:,i),nline,coeffobs(1,1),resp5(4,i)  ) ! ex,fp deleted 2021.09.15
 CALL CALOBSRESP(fs(:,i),nline,coeffobs(1,2),resp5(5,i)  ) ! ey,fp deleted 2021.09.15
end do      ! 2017.07.11

write(*,*) "### CALOBSEBCOMP END!! ###" ! 2017.07.12

return
end

!#############################################
!# copied from n_ebfem_bxyz_mpi.f90 on 2018.02.22
!# coded on 2017.05.31
subroutine OUTFREQ(freq,g_param,resp5)
use param
use outresp
implicit none
real(8),intent(in) :: freq
type(param_forward),intent(in) :: g_param
type(respdata),dimension(5) :: resp5 ! added on 2016.10.27 for plan view of ex

     CALL OUTFREQFILES(freq,resp5(1),g_param,"Hx") ! see m_outresp.f90
     CALL OUTFREQFILES(freq,resp5(2),g_param,"Hy")
     CALL OUTFREQFILES(freq,resp5(3),g_param,"Hz")
     CALL OUTFREQFILES(freq,resp5(4),g_param,"Ex") ! see m_outresp.f90
     CALL OUTFREQFILES(freq,resp5(5),g_param,"Ey")

return
end
!############################################## subroutine setNPHYS1INDEX2COND
! Modified on 2017.05.31
! Coded on 2017.05.12
subroutine SETNPHYS1INDEX2COND(g_mesh,g_cond)
use param
use mesh_type
implicit none
type(mesh),      intent(in)    :: g_mesh
type(param_cond),intent(inout) :: g_cond
integer(4) :: nphys1,nphys2,i,ntet,icount
integer(4),allocatable,dimension(:,:) :: n4flag

!#[0]## set
  ntet   = g_mesh%ntet
  allocate(n4flag(ntet,2))
  n4flag = g_mesh%n4flag

!#[1]## calculate nphys1 and nphys2
  nphys2 = 0
  do i=1,ntet
   if ( n4flag(i,1) .ge. 2 ) nphys2 = nphys2 + 1 ! count land elements 2017.09.29
  end do

  if (g_cond%condflag .eq. 1 ) then ! check when cond file is given
    if ( g_cond%nphys2 .ne. nphys2) then
     write(*,*) "GEGEGE nphys2=",nphys2,"g_cond%nphys2=",g_cond%nphys2
     stop
    end if
  end if

!#[2]## set nphys1 and nphys2
  nphys1 = ntet - nphys2
  g_cond%nphys1 = nphys1
  g_cond%nphys2 = nphys2 ! # of elements in 2nd physical volume (land)
  g_cond%ntet   = ntet
  write(*,*) "nphys1=",nphys1,"nphys2=",nphys2,"ntet=",g_mesh%ntet

!#[3]## prepare the rho, sigma in the case where file is not given
  ! If the file is given the following is allocated in READCOND in m_param.f90
  if (g_cond%condflag .eq. 0 ) then ! file is not given
   allocate( g_cond%sigma(nphys2) )
   allocate( g_cond%rho(  nphys2) )
   allocate( g_cond%index(nphys2) )
   do i=1,nphys2                          ! 2017.09.29
    if ( n4flag(nphys1+i,1) .le. 1 ) then ! 2017.09.29
     write(*,*) "GEGEGE! i=",i,"n4flag(nphys1+i,1)=",n4flag(nphys1+i,1),"nphys1=",nphys1
     stop                                 ! 2017.09.29
    end if                                ! 2017.09.29
    g_cond%sigma(i) = g_cond%sigma_land(n4flag(nphys1+i,1)-1) ! 2017.09.29
    g_cond%rho(i)   = 1.d0/g_cond%sigma(i)
   end do
  end if

!#[3]## set index
  do i=1,nphys2
   g_cond%index(i) = nphys1 +i
  end do

write(*,*) "### SETNPYS1INDEX2COND ###"
return
end
!###################################################################
! modified for spherical on 2016.11.20
! iflag = 0 for xyz
! iflag = 1 for xyzspherical
subroutine GENXYZMINMAX(em_mesh,g_param)
use param ! 2016.11.20
use mesh_type
implicit none
type(mesh),intent(in) :: em_mesh
type(param_forward),intent(inout) :: g_param
real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
real(8) :: xyz(3,em_mesh%node),xyzminmax(6)
integer(4) :: i
xyz = em_mesh%xyz ! normal
xmin=xyz(1,1) ; xmax=xyz(1,1)
ymin=xyz(2,1) ; ymax=xyz(2,1)
zmin=xyz(3,1) ; zmax=xyz(3,1)

do i=1,em_mesh%node
 xmin=min(xmin,xyz(1,i))
 xmax=max(xmax,xyz(1,i))
 ymin=min(ymin,xyz(2,i))
 ymax=max(ymax,xyz(2,i))
 zmin=min(zmin,xyz(3,i))
 zmax=max(zmax,xyz(3,i))
end do

xyzminmax(1:6)=(/xmin,xmax,ymin,ymax,zmin,zmax/)

!# set output
g_param%xyzminmax = xyzminmax

write(*,*) "### GENXYZMINMAX END!! ###"
return
end

!#############################################
!# copied from 2018.02.22
!# modified on 2017.10.12
!# modified on 2017.07.11 to include nsr
!# coded    on 2017.05.31
subroutine ALLOCATERESP(nobs,nsr,resp,ip,nfreq)
use outresp
use param
implicit none
integer(4),         intent(in)    :: nobs
integer(4),         intent(in)    :: nsr ! 2017.07.11
integer(4),         intent(in)    :: nfreq,ip
type(respdata),     intent(inout) :: resp(5,nsr,nfreq) !2017.07.11
integer(4)                        :: i,j,k

do j=1,nfreq
 do i=1,5
  do k=1,nsr ! 2017.07.11
   CALL ALLOCATERESPDATA(nobs,resp(i,k,j)) ! 2017.07.11
  end do     ! 2017.07.11
 end do
end do

if( ip .eq. 0) write(*,*) "### ALLOCATERESP END!! ###"
return
end
!######################################## OBSPRIMARYANA
subroutine OBSPRIMARYANA(sparam,g_param,freq,icomp)
use param
implicit none
integer(4),          intent(in) :: icomp
real(8),             intent(in) :: freq
type(param_forward), intent(in) :: g_param
type(param_source),  intent(in) :: sparam
real(8)    :: amp,phase,a(3),p(3),x3(3),x1(3),x2(3),I,ab(3),pb(3)
complex(8) :: one=(1.d0,0.d0),iunit=(0.d0,1.d0),e3(3),fp(3),b3(3),e1(3),e2(3),phi
integer(4) :: j,k,ntimes=20
real(8)    :: omega,pi, dd,ds
pi=4.d0*datan(1.d0)
omega=2.d0*pi*freq

x1(1:3)=sparam%xs1(1:3,1)*1.d3 ! [m]
x2(1:3)=sparam%xs2(1:3,1)*1.d3 ! [m]
ds=x2(1)-x1(1) ! [m]
I      =sparam%I        ! [A]

open(51,file="bpobs_ana.dat")
open(52,file="epobs_ana.dat")

do j=1,g_param%nobs
 x3(1:3)=(/g_param%xyzobs(1,j),g_param%xyzobs(2,j),g_param%xyzobs(3,j)/)*1.d3 ! [m]

!# magnetic field
!  CALL magfromwirecurrent(x3,x1,x2,I,ntimes,fp(1:3),freq) ! primary magnetic field [nT]
  CALL anahalfspace(x3,I,ds,freq,3,1,fp) ! 3 for total, 1 for magnetic field
  a(1:3)=(/  amp(fp(1)),  amp(fp(2)),  amp(fp(3))/)
  p(1:3)=(/phase(fp(1)),phase(fp(2)),phase(fp(3))/)
  write(51,*) g_param%xyzobs(1:2,j),(a(k),p(k),k=1,3)
!# electric field
!  CALL elefromwirecurrent(x3,x1,x2,I,freq,ntimes,fp(1:3),e3,phi) ! primary electric field [uV/m]
  CALL anahalfspace(x3,I,ds,freq,3,2,fp) ! 3 for total, 2 for electric field
  a(1:3)=(/  amp(fp(1)),  amp(fp(2)),  amp(fp(3))/)
  p(1:3)=(/phase(fp(1)),phase(fp(2)),phase(fp(3))/)
  !# Note the bx calculated from e1, e2 are consistent with those from magfromwirecurrent !!
  ! Then the magfromwire and elefromwirecurrent are correct !! on 2016.03.25
  dd=0.5d0 ! [m]
!  CALL elefromwirecurrent(x3+(/0.d0,dd,0.d0/),x1,x2,I,freq,ntimes,e1,b3) ! primary electric field [uV/m]
!  CALL elefromwirecurrent(x3-(/0.d0,dd,0.d0/),x1,x2,I,freq,ntimes,e2,b3) ! primary electric field [uV/m]
!  b3(1)=(0.d0,0.d0)
!  b3(2)=(0.d0,0.d0)
!  b3(3)=iunit/omega*(e1(1)-e2(1))/(2.*dd)*1.d3 ! [nT]
!  ab(1:3)=(/  amp(b3(1)),  amp(b3(2)),  amp(b3(3))/)
!  pb(1:3)=(/phase(b3(1)),phase(b3(2)),phase(b3(3))/)
  write(52,*) g_param%xyzobs(1:2,j),(a(k),p(k),k=1,3) !,(ab(k),pb(k),k=1,3)

end do

close(51) ; close(52)

write(*,*) "### OBSPRIMARYANA END!! ###"
return
end subroutine OBSPRIMARYANA
!######################################## OBSPRIMARYOUT
subroutine OBSPRIMARYOUT(fp,doftot,coeffobs,g_param,omega,icomp)
use matrix
use outresp
use param
implicit none
integer(4),           intent(in) :: doftot, icomp
complex(8),           intent(in) :: fp(doftot)
real(8),              intent(in) :: omega
type(real_crs_matrix),intent(in) :: coeffobs(2,3) ! see m_matrix.f90
type(param_forward),  intent(in) :: g_param
complex(8) :: fp_x(g_param%nobs), fp_new(doftot),iunit=(0.d0,1.d0)
complex(8) :: fp_y(g_param%nobs)
complex(8) :: fp_z(g_param%nobs)
real(8) :: amp,phase,a(3),p(3)
integer(4) :: i,j

!#[2]## calculate bp,bs,bt at observation points

if ( icomp .eq. 1) open(51,file="bpobs.dat")
if ( icomp .eq. 2) open(51,file="epobs.dat")
 CALL mul_matcrs_cv(coeffobs(1,1),fp,doftot,fp_x) ! see m_matrix.f90
 CALL mul_matcrs_cv(coeffobs(1,2),fp,doftot,fp_y) ! see m_matrix.f90
 CALL mul_matcrs_cv(coeffobs(1,3),fp,doftot,fp_z) ! see m_matrix.f90
 do i=1,g_param%nobs
  a(1:3)=(/   amp(fp_x(i)),    amp(fp_y(i)),    amp(fp_z(i))/)
  p(1:3)=(/ phase(fp_x(i)),  phase(fp_y(i)),  phase(fp_z(i))/)
  write(51,*) g_param%xyzobs(1:2,i),(a(j),p(j),j=1,3) ! fp_x, fp_y,fp_z
 end do
 close(51)

if ( icomp .eq. 2) then ! also bpobs
 open(51,file="bpobs.dat")
 fp_new=iunit/omega*fp ! [mV/km] -> [nT] ! B=i/omega*rotE, int rotE dl=Bds
 CALL mul_matcrs_cv(coeffobs(2,1),fp_new,doftot,fp_x) ! see m_matrix.f90
 CALL mul_matcrs_cv(coeffobs(2,2),fp_new,doftot,fp_y) ! see m_matrix.f90
 CALL mul_matcrs_cv(coeffobs(2,3),fp_new,doftot,fp_z) ! see m_matrix.f90
 do i=1,g_param%nobs
  a(1:3)=(/    amp(fp_x(i)),    amp(fp_y(i)),    amp(fp_z(i))/)
  p(1:3)=(/  phase(fp_x(i)),  phase(fp_y(i)),  phase(fp_z(i))/)
  write(51,*) g_param%xyzobs(1:2,i),(a(j),p(j),j=1,3) ! fp_x, fp_y,fp_z
 end do
 close(51)
end if

write(*,*) "### OBSPRIMARYOUT END!! ###"
return
end subroutine OBSPRIMARYOUT

!######################################## CALOBSRESP
!# Coded on Nov. 21, 2015
!# This calculates the output b fields and output results
subroutine CALOBSRESP(ft,doftot,coeffobs,resp) ! fp, fs -> ft 2021.09.15
use matrix
use outresp
implicit none
integer(4),           intent(in)    :: doftot
complex(8),           intent(in)    :: ft(doftot) ! 2021.09.15
type(real_crs_matrix),intent(in)    :: coeffobs ! see m_matrix.f90
type(respdata),       intent(inout) :: resp   ! see m_outresp.f90
complex(8),allocatable,dimension(:) :: ftobs
real(8) :: amp,phase
integer(4) :: i
allocate(ftobs(resp%nobs)) ! 2021.09.15

!#[1]## generate btotal

!#[2]## calculate bp,bs,bt at observation points
!CALL mul_matcrs_cv(coeffobs,fp(1:doftot),doftot,fpobs) ! see m_matrix.f90
!CALL mul_matcrs_cv(coeffobs,fs(1:doftot),doftot,fsobs) ! see m_matrix.f90
CALL mul_matcrs_cv(coeffobs,ft(1:doftot),doftot,ftobs) ! see m_matrix.f90
!write(*,*) "bpobs,bsobs,btobs is created!"

!#[3]## cal b3 comp and output
do i=1,resp%nobs
! resp%fpobsamp(i)  =amp  (fpobs(i)) ! amplitude of bz primary
! resp%fpobsphase(i)=phase(fpobs(i)) ! phase of primary fields
 resp%ftobsamp(i)  =amp  (ftobs(i)) ! amp of bz
 resp%ftobsphase(i)=phase(ftobs(i)) ! phase of bz
 resp%ftobs(i)     =ftobs(i)        ! 2021.09.15
! resp%fsobsamp(i)  =amp  (fsobs(i)) ! amp of bz   ! commented out on 2021.09.15
! resp%fsobsphase(i)=phase(fsobs(i)) ! phase of bz ! commented out on 2021.09.15
end do

write(*,*) "### CALOBSRESP END!! ###"
return
end
!######################################## function phase
function phase(c) ! [deg]
implicit none
complex(8),intent(in) :: c
real(8) :: phase
real(8),parameter :: pi=4.d0*datan(1.d0), r2d=180.d0/pi
 phase=datan2(dimag(c),dreal(c))*r2d
 return
end
!######################################## function amp
function amp(c)
implicit none
complex(8),intent(in) :: c
real(8) :: amp
 amp=dsqrt(dreal(c)**2.d0 + dimag(c)**2.d0)
 return
end

!###################################################################
! Coded on Feb 15, 2016
subroutine set_table_dof(dofn,nodtot,table_dof,doftot_ip)
implicit none
	integer(4),intent(in) :: dofn, nodtot, doftot_ip
	integer(4),intent(out) :: table_dof(nodtot,dofn)
	integer(4) :: i,j,dof_id
	dof_id=0
	do i=1,nodtot
       do j=1,dofn
	  dof_id = dof_id + 1
	  table_dof(i,j)= dof_id
       end do
	end do
	if (dof_id .ne. doftot_ip ) then
       write(*,*) "GEGEGE dof_id=",dof_id,"doftot_ip=",doftot_ip
	 stop
	end if
	write(*,*) "### SET TABLE_DOF END!! ###"
return
end


!####################################################### PREPOBSCOEFF
! copied from ../FEM_node/n_bzfem.f90
! adjusted to edge-FEM code
subroutine PREPOBSCOEFF(g_param,h_mesh,l_line,coeffobs)
use mesh_type
use line_type
use param
use matrix
use fem_edge_util
implicit none
type(mesh),           intent(in)  :: h_mesh
type(line_info),      intent(in)  :: l_line
type(param_forward),  intent(in)  :: g_param
type(real_crs_matrix),intent(out) :: coeffobs(2,3)     ! (1,(1,2,3)) for edge (x,y,z)
real(8) :: x3(3),a4(4),r6(6),len(6),w(3,6),elm_xyz(3,4),v
real(8) :: coeff(3,6),coeff_rot(3,6)
real(8) :: wfe(3,6) ! face basis function * rotation matrix
integer(4) :: iele,i,ii,j,k,l,jj,n6(6),ierr

!#[1]# allocate coeffobs
do i=1,2 ; do j=1,3
 coeffobs(i,j)%nrow=g_param%nobs
 coeffobs(i,j)%ntot=g_param%nobs * 6 ! for lines of tetrahedral mesh
 allocate(coeffobs(i,j)%stack(0:coeffobs(i,j)%nrow))
 allocate( coeffobs(i,j)%item(  coeffobs(i,j)%ntot))
 allocate(  coeffobs(i,j)%val(  coeffobs(i,j)%ntot))
 coeffobs(i,j)%stack(0)=0
end do ; end do

!#[2]# find element and set values to coeffobs
ii=0
do i=1,g_param%nobs
 if (g_param%lonlatflag .eq. 2 ) then ! xyobs is already set
  x3(1:3)=(/g_param%xyzobs(1,i),g_param%xyzobs(2,i),g_param%xyzobs(3,i)/) ! [km]
  call FINDELEMENT0(x3,h_mesh,iele,a4) ! see m_mesh_type.f90
!  write(*,*) "x3(1:3,i)=",x3(1:3)
!  write(*,*) "ieleobs(i)=",iele
!  write(*,*) "coeff(i,1:4)=",coeff(1:4)
  do j=1,2
   do k=1,3
   coeffobs(j,k)%stack(i)=coeffobs(j,k)%stack(i-1) + 6
   end do
  end do
  do j=1,4
    elm_xyz(1:3,j)=h_mesh%xyz(1:3,h_mesh%n4(iele,j))
  end do
  CALL EDGEBASISFUN(elm_xyz,x3,w  ,len,v) !v[km]^3,w[1/km],len[km],see m_fem_edge_util.f90
  CALL FACEBASISFUN(elm_xyz,x3,wfe,v) !v[km]^3,w[1/km],len[km],see m_fem_edge_util.f90
  do j=1,6
   ! coeff is values for line-integrated value
   coeff(1:3,j)     =   w(1:3,j) * isign(1,l_line%n6line(iele,j)) ! for icom-th component
   coeff_rot(1:3,j) = wfe(1:3,j) * isign(1,l_line%n6line(iele,j)) ! for icom-th component
   r6(j)=abs(l_line%n6line(iele,j))*1.d0
  end do
  n6(1:6)=(/1,2,3,4,5,6/)
  call SORT_INDEX(6,n6,r6) !sort n4 index by r4 : see sort_index.f90
  do j=1,6
   jj=abs(l_line%n6line(iele,n6(j))) ! jj is global node id
   ii=ii+1                  ! ii is entry id for coeffobs matrix
   do l=1,3 ! l for component
    coeffobs(1,l)%item(ii)=jj
    coeffobs(2,l)%item(ii)=jj
    coeffobs(1,l)%val(ii)=    coeff(l,n6(j)) ! for x component
    coeffobs(2,l)%val(ii)=coeff_rot(l,n6(j)) ! for x component
   end do
  end do
 else
  write(*,*) "GEGEGE! in PREPOBSCOEFF"
  write(*,*) "g_param%lonlatflag",g_param%lonlatflag,"should be 2 here."
  stop
 end if
end do

write(*,*) "### PREPOBSCOEFF END!! ###"
return

end

!####################################################### PREPOBSCOEFF_XY
!# copied from n_ebfem_bxyz_mpi.f90 on 2018.02.22
!# 2017.10.12
!# coefficients for obs based on nx*ny grid
subroutine PREPOBSCOEFF_XY(g_param,em_mesh,h_mesh,l_line,obs_xy)
use param
use mesh_type
use line_type
use matrix
use obs_type
use fem_edge_util
!use spherical
implicit none
type(param_forward),             intent(in)    :: g_param
type(mesh),                      intent(in)    :: em_mesh ! 3-D whole mesh
type(mesh),                      intent(in)    :: h_mesh  ! 2-D triangle topo mesh
type(line_info),                 intent(in)    :: l_line
type(obs_info),                  intent(out)   :: obs_xy
type(grid_list_type)                           :: glist
real(8)                                        :: xyzminmax(6)
integer(4)                                     :: nx,ny,nz !lower if element was not found
integer(4)                                     :: ntet,i,j,ii,nobs,ndat,nodek
real(8)                                        :: x3(3)
real(8)                                        :: xmin,xmax,ymin,ymax,dx,dy
integer(4)                                     :: nnx,nny,nnz=1 ! vertical plane

!#[0]## set input
 nnx            = g_param%nx         ! 2017.10.11
 nny            = g_param%ny         ! 2017.10.11
 nobs           = nnx*nny
 obs_xy%nobs    = nobs
 obs_xy%name    = g_param%xyfilehead ! 2017.10.12
 xyzminmax      = g_param%xyzminmax  ! 2016.11.20
 !# focus only on the center area
 xmin           = g_param%xbound(2)  ! 2017.10.12
 xmax           = g_param%xbound(3)  ! 2017.10.12
 ymin           = g_param%ybound(2)  ! 2017.10.12
 ymax           = g_param%ybound(3)  ! 2017.10.12
 ntet           = em_mesh%ntet

!#[1]# prepare coordinate of observatories for seafloor
 allocate(obs_xy%xyz_obs(3,nobs))
 write(*,*) "nobs=",nobs,"PREPOBSCOEFF_XY"
 dx=(xmax - xmin)/float(nnx)
 dy=(ymax - ymin)/float(nny)
 do j=1,nnx
  do i=1,nny
   ii=(j-1)*nny + i
   obs_xy%xyz_obs(1,ii)=(j-1)*dx + xmin + dx/2.d0
   obs_xy%xyz_obs(2,ii)=(i-1)*dy + ymin + dy/2.d0
  end do
 end do
 obs_xy%xyz_obs(3,:)=0.d0

!#[2]## read ki_mesh and set input added on 2016.11.20
 CALL PREPZOBSGRD(g_param,h_mesh,obs_xy) ! see below 2017.10.11

!#[3]## generate horizontal grid and classify all the elements to gridded boxes
 nx=1000;ny=1000;nz=10
 CALL allocate_grid_list(nx,ny,nz,ntet,glist)  ! see m_mesh_type.f90
 CALL gengridforlist(xyzminmax,glist)          ! see m_mesh_type.f90
 CALL classifyelement2grd(em_mesh,glist)       ! classify ele to glist,see m_mesh_type.f90

!#[4]# allocate coeffobs and calculate
 CALL CALCOEFF_LINE(em_mesh,l_line,glist,obs_xy) ! 2017.10.11

write(*,*) "### PREPOBSCOEFF_XY END!! ###"
return

99 continue
 write(*,*) "GEGEGE! ndat=",ndat,"is not equal to nodek",nodek,"PREPOBSCOEFF_XY"
 stop
end
!########################################################### PREPZOBS
! Copied from tsunami/3D/ana_comp/src/prepobscoeff.f90
! on 2017.10.11
! calcualte h_mesh%zobs from h_mesh%xyz
subroutine PREPZOBSGRD(g_param,h_mesh,obs)
use mesh_type
use obs_type  ! 2016.11.20
use matrix    ! 2016.11.20
use triangle  ! 2016.11.23 see m_triangle
use param     ! 2016.11.23
use triangle  ! 2016.11.23
implicit none
type(mesh),                     intent(in)     :: h_mesh  ! 2d mesh with topography
type(param_forward),            intent(in)     :: g_param
type(obs_info),                 intent(inout)  :: obs
real(8),            allocatable,dimension(:,:) :: xyzk ! triangle,tetra nodes
real(8),            allocatable,dimension(:,:) :: xyzobs
integer(4),         allocatable,dimension(:,:) :: n3k
type(grid_list_type) :: glist
integer(4)           :: ntri,nobs
integer(4)           :: i,j,k,n1,n2,n3,m1,m2,m3,nx,ny,iele
real(8),dimension(2) :: x12,x13,x23,v1,v2,v3
real(8)              :: a3(3),a,xyzminmax(6)

!#[0]## set
 nobs      = obs%nobs
 write(*,*) "nobs=",nobs,"PREPZOBSGRD"
 allocate( xyzk(3,h_mesh%node), xyzobs(3,nobs), n3k(h_mesh%ntri,3)  )
 xyzk      = h_mesh%xyz    ! triangle mesh
 n3k       = h_mesh%n3
 ntri      = h_mesh%ntri
 xyzobs    = obs%xyz_obs
 xyzminmax = g_param%xyzminmax
 write(*,*) "xyzminmax=",xyzminmax

!#[0]## generate element list
 nx=100;ny=100
 CALL allocate_2Dgrid_list(nx,ny,ntri,glist)   ! see m_mesh_type.f90
 CALL gen2Dgridforlist(xyzminmax,glist) ! see m_mesh_type.f90
 write(*,*) "h_mesh%ntri=",h_mesh%ntri
 CALL classifytri2grd(h_mesh,glist)   ! classify ele to glist,see

!#[1] search for the triangle including (x1,y1)
 write(*,*) "nobs=",nobs
 do j=1,nobs
 call findtriwithgrid(h_mesh,glist,xyzobs(1:2,j),iele,a3)
    n1 = n3k(iele,1)
    n2 = n3k(iele,2)
    n3 = n3k(iele,3)
    obs%xyz_obs(3,j) = a3(1)*xyzk(3,n1)+a3(2)*xyzk(3,n2)+a3(3)*xyzk(3,n3) - 0.001! 1m below
!  write(*,*) "j=",j,"xyzobs(1:2,j)=",xyzobs(1:2,j),"determined z=",obs%xyz_obs(3,j)
 end do

 write(*,*) "### PREPZOBSGRD END!! ###"
 return
end

!####################################################### cal coeff for observatories
! copied from tsunami/3D_ana_comp/src/solver/n_ebfem_tsunamiEM.f90 on 2016.10.27
!# coded on May 20, 2016
!# coeff_line(1,1:3) for A = [edge_basis_fun]{Al}
!# coeff_line(2,1:3) for B = [face_basis_fun]{Al}
subroutine CALCOEFF_LINE(h_mesh,l_line,glist,obs)
use mesh_type
use matrix
use line_type
use obs_type
use fem_edge_util
implicit none
type(mesh),intent(in) :: h_mesh
type(line_info),intent(in) :: l_line
type(grid_list_type),intent(in) :: glist
type(obs_info),intent(inout) :: obs
!# internal variables
integer(4) :: i,j,l,nrow,ncolm,ii,jj,iele,n6(6)
real(8) :: x3(3),elm_xyz(3,4),r6(6),w(3,6),wfe(3,6),len(6),v,a4(4),coeff(3,6),coeff_rot(3,6)

!#[1]## allocate coeff_line
nrow=obs%nobs
ncolm=6
do i=1,2 ; do j=1,3
 CALL allocate_real_crs_with_steady_ncolm(obs%coeff(i,j),nrow,ncolm) ! see m_matrix.f90
end do ; end do

!#[2]## find element with glist
do i=1,obs%nobs
  if (mod(i,10000) .eq. 0 ) write(*,*) "i=",i
!  write(*,*) "i-th triangle",i,"ntrik=",ntrik
  x3(1:3)=obs%xyz_obs(1:3,i)
  CALL FINDELEWITHGRID(h_mesh,glist,x3,iele,a4) ! see m_mesh_type.f90

  !#[2-3]##
  do j=1,4
    elm_xyz(1:3,j)=h_mesh%xyz(1:3,h_mesh%n4(iele,j))
  end do
  CALL EDGEBASISFUN(elm_xyz,x3,w,len,v) !v[km]^3,w[1/km],len[km],see m_fem_edge_util.f90
  CALL FACEBASISFUN(elm_xyz,x3,wfe,v) !v[km]^3,w[1/km],len[km],see m_fem_edge_util.f90
  do j=1,6
   ! coeff is values for line-integrated value
   coeff(1:3,j)     =   w(1:3,j) * isign(1,l_line%n6line(iele,j)) ! for icom-th component
   coeff_rot(1:3,j) = wfe(1:3,j) * isign(1,l_line%n6line(iele,j)) ! for icom-th component
   r6(j)=abs(l_line%n6line(iele,j))*1.d0
  end do
  n6(1:6)=(/1,2,3,4,5,6/)
  call SORT_INDEX(6,n6,r6) !sort n4 index by r4 : see sort_index.f90
  do j=1,6
   jj=abs(l_line%n6line(iele,n6(j))) ! jj is global node id
   ii=(i-1)*6+j                  ! ii is entry id for coeffobs matrix
   do l=1,3 ! l-th component
    obs%coeff(1,l)%item(ii)=jj
    obs%coeff(2,l)%item(ii)=jj
    obs%coeff(1,l)%val(ii)=    coeff(l,n6(j)) ! for x component
    obs%coeff(2,l)%val(ii)=coeff_rot(l,n6(j)) ! for x component
   end do
  end do
end do


return
end
!################################################################# PREPZSRCOBS
! modified on 2017.07.11 to include multi source
! Coded    on 2017.02.21
subroutine PREPZSRCOBS(h_mesh,g_param,s_param)
use param
use mesh_type
use triangle
implicit none
type(mesh),          intent(in)       :: h_mesh
type(param_forward), intent(inout)    :: g_param
type(param_source),  intent(inout)    :: s_param
type(grid_list_type)                  :: glist
integer(4)                            :: nobs,nx,ny
real(8),   allocatable,dimension(:,:) :: xyzobs,xyz
integer(4),allocatable,dimension(:,:) :: n3k
real(8),   allocatable,dimension(:)   :: znew
real(8)    :: a3(3)
integer(4) :: iele,n1,n2,n3,j,k,ntri
integer(4) :: nsr                                ! 2017.07.11
real(8),allocatable,dimension(:,:)    :: xs1,xs2 ! 2017.07.18
real(8)    :: xyzminmax(6)                       ! 2017.07.18

!#[0]## cal xyzminmax of h_mesh
!  CALL GENXYZMINMAX(h_mesh,g_param)  ! commented out  2017.10.12

!#[1]## set
nsr       = s_param%nsource     ! 2017.07.11
allocate(xs1(3,nsr),xs2(3,nsr)) ! 2017.07.11
allocate(xyz(3,h_mesh%node),n3k(h_mesh%ntri,3))
allocate(xyzobs(3,g_param%nobs))
allocate(znew(g_param%nobs))
nobs      = g_param%nobs
xyz       = h_mesh%xyz    ! triangle mesh
n3k       = h_mesh%n3
ntri      = h_mesh%ntri
xyzobs    = g_param%xyzobs
xs1       = s_param%xs1
xs2       = s_param%xs2
xyzminmax = g_param%xyzminmax


!#[2]## cal z for nobsr
nx=1000;ny=1000
CALL allocate_2Dgrid_list(nx,ny,ntri,glist)   ! see m_mesh_type.f90
CALL gen2Dgridforlist(xyzminmax,glist) ! see m_mesh_type.f90
CALL classifytri2grd(h_mesh,glist)   ! classify ele to glist,see


!#[3] search for the triangle including (x1,y1)
do j=1,nobs
    call findtriwithgrid(h_mesh,glist,xyzobs(1:2,j),iele,a3)
    n1 = n3k(iele,1); n2 = n3k(iele,2) ; n3 = n3k(iele,3)
    znew(j) = a3(1)*xyz(3,n1)+a3(2)*xyz(3,n2)+a3(3)*xyz(3,n3) + xyzobs(3,j)
!
    write(*,*) "xyobs(1:2,j)=",xyzobs(1:2,j)
    write(*,*) j,"/nobs",xyzobs(3,j),"->",znew(j),"[km]"
end do

!#[3-2]## source z
 do k=1,nsr                                                               ! 2017.07.11
!#  start point
    call findtriwithgrid(h_mesh,glist,xs1(1:2,k),iele,a3)                 ! 2017.07.11
    n1 = n3k(iele,1); n2 = n3k(iele,2) ; n3 = n3k(iele,3)
    xs1(3,k) = a3(1)*xyz(3,n1)+a3(2)*xyz(3,n2)+a3(3)*xyz(3,n3) + xs1(3,k) ! 2017.07.11
    !
    write(*,*) "k",k,"xs1(1:2,k)=",xs1(1:2,k)                             ! 2017.07.11
    write(*,*) "z",s_param%xs1(3,k),"->",xs1(3,k),"[km]"                  ! 2017.07.11

!#  end point
    call findtriwithgrid(h_mesh,glist,xs2(1:2,k),iele,a3)                 ! 2017.07.11
    n1 = n3k(iele,1); n2 = n3k(iele,2) ; n3 = n3k(iele,3)
    xs2(3,k) = a3(1)*xyz(3,n1)+a3(2)*xyz(3,n2)+a3(3)*xyz(3,n3) + xs2(3,k) ! 2017.07.11
    !
    write(*,*) "k",k,"xs2(1:2,k)=",xs2(1:2,k)                             ! 2017.07.11
    write(*,*) "z",s_param%xs2(3,k),"->",xs2(3,k),"[km]"                  ! 2017.07.11
 end do                                                                   ! 2017.07.11


!#[4]## set znew to xyz_r
    g_param%xyzobs(3,1:nobs) = znew(1:nobs)
    s_param%xs1(3,:) = xs1(3,:)   ! 2017.07.11
    s_param%xs2(3,:) = xs2(3,:)   ! 2017.07.11


write(*,*) "### PREPZSRCOBS END ###"
return
end
!########################################################### PREPZOBS
! calcualte h_mesh%zobs from h_mesh%xyz
subroutine PREPZOBS(g_param,h_mesh,sparam)
use mesh_type
use param
implicit none
type(mesh),        intent(in)    :: h_mesh
type(param_forward),intent(inout) :: g_param
type(param_source),intent(inout)  :: sparam
integer(4) :: isurface
integer(4) :: i,j,n1,n2,n3
real(8),dimension(2) :: x12,x13,x23,v1,v2,v3
real(8) :: a1,a2,a3,a
isurface = g_param%surface_id_ground
!##

!#[1] search for the triangle including (x1,y1)
do j=1,g_param%nobs
 do i=1,h_mesh%ntri
  if ( h_mesh%n3flag(i,1) .eq. g_param%surface_id_ground ) then !usually when n3flag(i,1)=6
   n1=h_mesh%n3(i,1)
   n2=h_mesh%n3(i,2)
   n3=h_mesh%n3(i,3)
   x12(1:2)= h_mesh%xyz(1:2,n2) - h_mesh%xyz(1:2,n1) ! [km]
   x13(1:2)= h_mesh%xyz(1:2,n3) - h_mesh%xyz(1:2,n1) ! [km]
   x23(1:2)= h_mesh%xyz(1:2,n3) - h_mesh%xyz(1:2,n2) ! [km]
   v1=g_param%xyzobs(1:2,j) - h_mesh%xyz(1:2,n1)
   v2=g_param%xyzobs(1:2,j) - h_mesh%xyz(1:2,n2)
   v3=g_param%xyzobs(1:2,j) - h_mesh%xyz(1:2,n3)
   a =( x13(1)*x12(2) - x13(2)*x12(1))
   a2=( x13(1)* v1(2) - x13(2)* v1(1))/a
   a3=(  v1(1)*x12(2) -  v1(2)*x12(1))/a
   a1=(  v2(1)*x23(2) -  v2(2)*x23(1))/a
   if (a1 .ge. 0. .and. a2 .ge. 0. .and. a3 .ge. 0. ) then
    g_param%xyzobs(3,j)=a1*h_mesh%xyz(3,n1)+a2*h_mesh%xyz(3,n2)+a3*h_mesh%xyz(3,n3) ! [km]
    write(*,*) "element # =",i
!   write(*,*) "a1,a2,a3=",a1,a2,a3
!   write(*,*) "xyzg(3,n1 - n3)=",xyzg(3,n1),xyzg(3,n2),xyzg(3,n3)
    goto 100
   end if
  end if
 end do
 100 continue
end do

write(*,*) "### PREPZOBS END!! ###"
return
end

!########################################################### CALXYUTM
! Coded on Dec. 24, 2015
subroutine CALXYUTM(n,lon,lat,lonorigin,latorigin,x,y,zone)
implicit none
integer(4),intent(in) :: n
character(3),intent(in) :: zone
real(8),intent(in)  :: lonorigin, latorigin
real(8),intent(in)  :: lon(n),lat(n)
real(8),intent(out) :: x(n),y(n)
real(8),dimension(n+1) :: lon1,lat1,x1,y1
integer(4) :: i

!#[1]# prepare
lon1(1)=lonorigin
lat1(1)=latorigin
lon1(2:n+1)=lon(1:n)
lat1(2:n+1)=lat(1:n)

!#[2]# cal UTM
CALL UTMGMT(n+1,lon1,lat1,x1,y1,zone(1:3),0) ! 0: LONLAT2XY

!#[3]# cal x, y
do i=1,n
 x(i)=x1(i+1)-x1(1)
 y(i)=y1(i+1)-y1(1)
end do

!#[4]# output
do i=1,n
 write(*,*)"lonlat=",lon1(i+1),lat1(i+1),"x,y=",x(i),y(i)
end do
return
end
!######################################################### UTMGMT
subroutine UTMGMT(n,xin,yin,xout,yout,zone,iflag)
implicit none
integer(4),intent(in) :: n
real(8),intent(in) :: xin(n),yin(n)
real(8),intent(out) :: xout(n),yout(n)
character(3),intent(in) :: zone
character(37) :: values
integer(4),intent(in) :: iflag ! 0: LONLAT2UTM, 1:UTM2LONLAT
integer(4) :: i
!#[1]## prepare input file
open(11,file="in.dat")
 write(11,'(2g18.10)') (xin(i),yin(i),i=1,n)
close(11)

!#[2]## use mapproject
if (iflag .eq. 0 ) then ! LONLAT2UTM
 CALL system("cat in.dat | mapproject -Ju"//zone(1:3)//"/1.0 -F -C > tmp.dat")
else if (iflag .eq. 1 ) then ! UTM2LONLAT
 CALL system("cat in.dat | mapproject -Ju"//zone(1:3)//"/1.0 -F -C -I > tmp.dat")
else
 write(*,*) "GEGEGE! iflag should be 0 (LONLAT2UTM) or 1(UTM2LONLAT), iflag=",iflag
 stop
end if

!#[3]## read xout and yout
 open(12,file="tmp.dat")
  do i=1,n
   read(12,*) xout(i),yout(i)
  end do
 close(12)
 call system("rm tmp.dat")
return
end
