! coded on 2021.06.08
program ebfem_3DMT
use mesh_type
use param    ! 2021.12.15
use param_mt ! 2021.12.15
use matrix
use line_type    ! see m_line_tyoe.f90
use iccg_var_takuto
use outresp      ! see ../common/m_outresp.f90 for respdata,respmt
use surface_type ! see src_2D/m_surface_type.f90
use obs_type     ! see m_obs_type.f90
use constants, only: pi,dmu
use face_type
implicit none
type(param_forward_mt)  :: g_param_mt ! 2021.12.15
type(param_cond)        :: g_cond  ! see src/common/m_param.f90
type(mesh)              :: g_mesh  ! see src/common/m_mesh_type.f90
type(mesh)              :: h_mesh  ! topography file
type(mesh)              :: ki_mesh ! 2026.07.29 read polygon_ki.msh
type(line_info)         :: g_line  ! m_line_type.f90
type(global_matrix)     :: A       ! see m_iccg_var_takuto.f90
type(real_crs_matrix)   :: coeffobs(2,3) !m_matrix.f90; 1:edge,
type(surface)           :: g_surface(6)
type(face_info)         :: g_face  ! see m_face_type.f90
integer(4)              :: nface,node,ntri,doftot,ntet
integer(4)              :: nline,nsr,nfreq,nfreq_ip
integer(4)              :: i,j,iresfile,i_surface,k
integer(4)              :: ixyflag
real(8)                 :: omega,freq
type(respdata),allocatable,dimension(:,:,:) :: resp5  ! 2018.02.22
type(respmt),  allocatable,dimension(:)     :: resp_mt !see m_outresp.f90
type(resptip), allocatable,dimension(:)     :: resp_tip! m_outresp.f90 2025.04.21
integer(4),    allocatable,dimension(:,:)   :: n4
complex(8),    allocatable,dimension(:,:)   :: al_MT ! 2021.09.14
integer(4)                                  :: ip=0,iele
character(1) ::num

!#[0]## read param
 CALL READPARAM_MT(g_param_mt,g_cond)!get parameter info and conductivity info

!#[1]## Mesh READ
 CALL READMESH_TOTAL(g_mesh,g_param_mt%g_meshfile) ! 3D mesh
 CALL READMESH_TOTAL(h_mesh,g_param_mt%z_meshfile) ! topography file
 CALL PREPZOBSMT(h_mesh,g_param_mt) ! 2026.07.29 modification for ActFEMtide
 CALL GENXYZMINMAX_MT(g_mesh,g_param_mt) ! generate xyzminmax 2021.12.15
  if (g_cond%condflag .eq. 1) then !"1" : conductivity file is given
   g_cond%ntet   = g_mesh%ntet
   g_cond%nphys1 = g_mesh%ntet - g_cond%nphys2
  end if
 CALL SETNPHYS1INDEX2COND(g_mesh,g_cond) ! see below

!#[2]## Line information
 CALL READLINE(g_param_mt%g_lineinfofile,g_line) ! see m_line_type.f90

 nline   = g_line%nline      ! 2021.09.14
 ntet    = g_mesh%ntet       ! 2107.09.08
 nsr     = 2                 ! 1 for ex, 2 for ey polarization 2021.09.14
 node    = g_mesh%node       ! 2017.09.08
 ixyflag = g_param_mt%ixyflag! 2018.02.22
 nfreq   = g_param_mt%nfreq  ! 2018.02.22
 n4      = g_mesh%n4         ! allocate n4(ntet,4) here 2020.10.04

!#[3]## Make face of 3D mesh
 CALL MKFACE(  g_face,node,ntet,4, n4) ! make g_line
 CALL MKN4FACE(g_face,node,ntet,   n4) ! make g_line%n6line
 CALL FACE2ELEMENT(g_face) ! see common/m_face_type.f90

!#[4]## Extract Faceinfo and fine boundary line of 3D mesh for TM mode calculation
!        Face2        ^ y
!        ----         |
! Face3 |    | Face5  |
!        ----         ----> x
!        Face4
 CALL EXTRACT6SURFACES(g_mesh,g_line,g_face,g_surface) ! ../src_2D/m_surface_type.f90
 CALL FINDBOUNDARYLINE(g_mesh,g_surface)  ! Find boudnary line m_surface_type.f90
 CALL PREPAOFSURFACE(g_surface(2:5),4,ip) ! allocate A and table_dof for 2DMT

!#[5]## allocate global matrix
 CALL set_size_of_A(6,nline,nline,nline,ntet,g_line%n6line,A,ip) ! m_iccg_var_takuto.f90 2021.09.16

!#[6]## prepare coefficients
do i=1,g_param_mt%nobs
 write(*,'(a,i3,a,3f15.7)') "i",i," xyzobs",g_param_mt%xyzobs(1:3,i)
end do
 CALL PREPOBSCOEFF(g_param_mt,g_mesh,g_line,coeffobs) ! for bx,by,bz,ex,ey component

!#[3]## Prepare coefficients for values at xyplane !2025.07.17 <under construction>
!  if ( ixyflag .eq. 1 ) then
!   CALL PREPOBSCOEFF_XY(g_param,g_mesh,h_mesh,g_line,obs_xy) ! 2018.02.22
!   allocate( resp_xy(5,nsr,nfreq),resp_xy_mt(nfreq) )        ! 225.07.17
!   CALL ALLOCATERESP(obs_xy%nobs,nsr,resp_xy,ip,nfreq)       ! 2018.02.22
!  end if

!#[7]## set resp5
allocate( resp5(5,nsr,nfreq), resp_mt(nfreq)) ! 2018.02.22
allocate( resp_tip(nfreq))                    ! 2025.04.21 resp_tip is added
CALL ALLOCATERESP_MT(g_param_mt%nobs,nsr,resp5,resp_mt,resp_tip,ip,nfreq)!2025.04.21,see below

allocate( al_MT(nline,2) ) ! 2021.09.14 for ex and ey source

!#[8]## set 2D conductivity at four side surfaces
CALL COND3DTO2D(g_mesh,g_surface,g_cond) ! see m_surface_type.f90 2021.06.01

!#[9]## frequency loop
do i = 1,nfreq
 freq  = g_param_mt%freq(i)
 write(*,*) "Freq",freq,"Hz" ! 2021.12.15
 omega = 2.*pi*freq

!#[9-1]## prepare boundary condition by 2DMT calculation at each freq
 do j=2,5 ! surface 2DTM loop, 2:north,3:west,4:south,5:east surface
  CALL forward_2DTM(g_surface(j),freq,g_param_mt,g_cond,ip,j)!## solve 2DTM for BC
 write(num,'(i1)') j
 !open(1,file="2D"//num//".dat")
 ! write(1,'(i4,2g15.7)') (k,g_surface(j)%bs(k),k=1,g_surface(j)%nline)
 !close(1)
 end do ! end surface loop


!#[9-2]## forward for both Ex and Ey polarization
 !write(*,*) "ip",ip
 !write(*,*) "g_param",g_param_mt%xbound(4)
 call forward_3DMT(A,g_param_mt,g_mesh,g_line,nline,al_MT,freq,g_cond,g_surface,ip)
 if (ip .eq. 0 .and. .false. ) then ! fs_mt01.dat
   open(1,file="fs_mt01.dat")
   write(1,'(i6,4f15.7)') (j,al_mt(j,1:2),j=1,nline)
   close(1)
   end if

!#[9-3]## calculate E,B at obs
 CALL CALOBSEBCOMP(al_MT,nline,nsr,omega,coeffobs,resp5(:,:,i)) !see below

!#[9-4]## cal MT impedance
 call CALRESPMT( resp5(:,1:2,i),resp_mt(i),omega )
 call CALRESPTIP(resp5(:,1:2,i),resp_tip(i),omega,ip) ! 2025.04.21

end do ! end frequency loop

 call OUTOBSRESPMT(g_param_mt,resp_mt,nfreq)
 call OUTOBSRESP_TIP(g_param_mt,resp_tip,nfreq) ! 2025.04.21

end program ebfem_3DMT

!#############################################
! coded on 2021.09.15
subroutine OUTOBSRESPMT(g_param_mt,resp_mt,nfreq)
use param_mt
use outresp
implicit none
type(param_forward_mt),intent(in) :: g_param_mt
integer(4),            intent(in) :: nfreq
type(respmt),          intent(in) :: resp_mt(nfreq)
character(50)  :: head, sour, site
character(70)  :: filename,filename2 ! 2021.12.15
integer(4)     :: i,j,k,l,nhead,nsite,nobs
real(8)        :: freq

!#[1]## set
 nobs  = g_param_mt%nobs
 head  = g_param_mt%outputfolder
 nhead = len_trim(head)

!#[2]##
 do l=1,nobs
  site  = g_param_mt%obsname(l)
  nsite = len_trim(site)
  filename  = head(1:nhead)//site(1:nsite)//"_MT.dat"
  filename2 = head(1:nhead)//site(1:nsite)//"_MT_imp.dat" ! 2021.12.15
   open(31,file=filename )
   open(32,file=filename2) ! 2021.12.15

   do i=1,nfreq
    freq = g_param_mt%freq(i)
    ! In homogeneous structure, arg(rho_yx)=pi/4, arg(rho_xy)=-3/4*pi, 
    ! because the coordinate system wehre z is upward positive
    write(31,'(9g15.7)') freq,resp_mt(i)%rhoxx(l),resp_mt(i)%phaxx(l),&
    &                         resp_mt(i)%rhoxy(l),resp_mt(i)%phaxy(l),&
    &                         resp_mt(i)%rhoyx(l),resp_mt(i)%phayx(l),&
    &                         resp_mt(i)%rhoyy(l),resp_mt(i)%phayy(l)
    write(32,'(9g15.7)') freq,real(resp_mt(i)%zxx(l)),imag(resp_mt(i)%zxx(l)),& ! 2021.12.15
    &                         real(resp_mt(i)%zxy(l)),imag(resp_mt(i)%zxy(l)),& ! 2021.12.15
    &                         real(resp_mt(i)%zyx(l)),imag(resp_mt(i)%zyx(l)),& ! 2021.12.15
    &                         real(resp_mt(i)%zyy(l)),imag(resp_mt(i)%zyy(l))    ! 2021.12.15
   end do

   close(31)
   close(32) ! 2021.12.15
  end do

 write(*,*) "### OUTOBSRESPMT END!! ###"

return
end

!########################################### OUTOBSRESP_TIP 2023.12.25
!# copied from OUTOBSFILESINV_TIP and modified
subroutine OUTOBSRESP_TIP(g_param_mt,resp_tip,nfreq)
 !# coded on 2023.12.25
 !# declaration
   use param_mt
   use outresp
   implicit none
   integer(4),                intent(in)    :: nfreq! 2017.09.11
   type(resptip),             intent(in)   :: resp_tip(nfreq) ! 2023.12.25
   type(param_forward_mt),    intent(in)    :: g_param_mt         ! 2022.01.02
   real(8)                                  :: freq            ! 2017.07.14
   integer(4)                               :: i,j,k,l,nobs    ! 2017.07.14
   character(100)                           :: filename1       ! 2022.01.22
   character(50)                            :: head,site       ! 2017.07.14
   integer(4)                               :: idat,ialphaflag ! 2017.09.11
   integer(4)                               :: icomp
  
 !#[1]## set
   nobs       = g_param_mt%nobs
   head       = g_param_mt%outputfolder     ! 2017.07.25
  
 !#[2]## output impedance and rho and phi to obs files
   do l=1,nobs
     site  = g_param_mt%obsname(l)
     filename1 = trim(head)//trim(site)//"_TIP.dat"     ! 2022.01.02
     open(31,file=filename1)

     do i=1,nfreq
       freq = g_param_mt%freq(i)
       write(31,'(5g15.7)') freq,resp_tip(i)%tx(l),resp_tip(i)%ty(l)
     end do

     close(31)
   end do

 write(*,*) "### OUTOBSRESP_TIP END!! ###"  
 return
 end

!#############################################
! coded on 2021.09.15
subroutine calrespmt(resp5,resp_mt,omega)
use outresp
use constants ! dmu,pi
implicit none
real(8),       intent(in)    :: omega
type(respdata),intent(in)    :: resp5(5,2) ! 1 for ex, 2 for ey polarization
type(respmt),  intent(inout) :: resp_mt
complex(8),allocatable,dimension(:,:) :: be5_ex,be5_ey ! be5 = bx,by,bz,ex,ey
complex(8) :: a,b,c,d,iunit=(0.d0,1.d0)
complex(8) :: det,z(2,2),bi(2,2),e(2,2)
integer(4) :: i,j,nobs
real(8)    :: coef,amp,phase

nobs = resp_mt%nobs
allocate(be5_ex(5,nobs),be5_ey(5,nobs))

! calculate impedance Z = E/B [mV/km]/[nT]
! (Ex_ex Ex_ey) = (Zxx Zxy)(Bx_ex Bx_ey)
! (Ey_ex Ey_ey) = (Zyx Zyy)(By_ex By_ey)
!  [E]=Z[B]
!  Z = [E][B]^-1

!# set bxyzexy_ex and bxyzexy_ey
do j=1,nobs
 do i=1,5
  be5_ex(i,j)=resp5(i,1)%ftobs(j) ! ex polarization
  be5_ey(i,j)=resp5(i,2)%ftobs(j) ! ey plarization
 end do
end do

!# calculate impedance
do j=1,nobs
 a = be5_ex(1,j) ! Bx_ex
 c = be5_ex(2,j) ! By_ex ! 2022.01.12
 b = be5_ey(1,j) ! Bx_ey ! 2022.01.12
 d = be5_ey(2,j) ! By_ey
 det = a*d - b*c
 bi(1,1:2)=(/ d, -b/)
 bi(2,1:2)=(/ -c, a/)
 bi = bi/det
 e(1,1:2)=(/be5_ex(4,j), be5_ey(4,j)/)
 e(2,1:2)=(/be5_ex(5,j), be5_ey(5,j)/)
 z = matmul(e,bi) ! 2022.01.12
 resp_mt%zxx(j) = z(1,1)
 resp_mt%zxy(j) = z(1,2)
 resp_mt%zyx(j) = z(2,1)
 resp_mt%zyy(j) = z(2,2)
! rho and pha, rhoa = mu/omega*|Z|**2.
coef = dmu/omega*1.d+6
 resp_mt%rhoxx(j) = coef*amp(z(1,1))**2. ! [Ohm.m]
 resp_mt%rhoxy(j) = coef*amp(z(1,2))**2. ! [Ohm.m]
 resp_mt%rhoyx(j) = coef*amp(z(2,1))**2. ! [Ohm.m]
 resp_mt%rhoyy(j) = coef*amp(z(2,2))**2. ! [Ohm.m]
!
 resp_mt%phaxx(j) = phase(z(1,1))
 resp_mt%phaxy(j) = phase(z(1,2))
 resp_mt%phayx(j) = phase(z(2,1))
 resp_mt%phayy(j) = phase(z(2,2))
end do

write(*,*) "### CALRESPMT END!! ###" ! 2024.10.06
return
end
!#############################################
! 2025.04.21 copied from n_inv_joint.f90
subroutine calresptip(resp5,resp_tip,omega,ip) ! 2023.12.23
  ! coded on 2023.12.23
  use outresp
  use constants ! dmu,pi
  implicit none
  real(8),       intent(in)     :: omega
  integer(4),    intent(in)     :: ip     ! 2022.12.05
  type(respdata),intent(in)     :: resp5(5,2) ! 1 for ex, 2 for ey polarization
  type(resptip),  intent(inout) :: resp_tip
  complex(8),    allocatable    :: be5_ex(:,:),be5_ey(:,:) ! be5 = bx,by,bz,ex,ey
  complex(8)                    :: a,b,c,d,iunit=(0.d0,1.d0)
  complex(8)                    :: det,txy(2,1),bi(2,2),e(2,1)
  integer(4)                    :: i,j,nobs
  real(8)                       :: coef,amp,phase
  
  nobs = resp_tip%nobs
  allocate(be5_ex(5,nobs),be5_ey(5,nobs))
  ! calculate tipper Txy = Bz/(Bx, By)
  ! (Bz_ex) = (Bx_ex By_ex) (Tx)
  ! (Bz_ey) = (Bx_ex By_ex) (Ty)
  ! [Bz] =[Bxy][T] 
  ! Then, Tx and Ty are obtained by
  ! [T] = [Bxy]^-1 [Bz]
  !
  !# set bxyzexy_ex and bxyzexy_ey
  !write(*,*) "nobs",nobs,"ip",ip
  do i=1,5
   do j=1,nobs
     be5_ex(i,j)=resp5(i,1)%ftobs(j) ! ex polarization
     be5_ey(i,j)=resp5(i,2)%ftobs(j) ! ey plarization
   end do
  end do
  
  write(*,*)
  !# calculate tipper
  do j=1,nobs
   a = be5_ex(1,j) ! Bx_ex
   b = be5_ex(2,j) ! By_ex ! 2022.12.25
   c = be5_ey(1,j) ! Bx_ey ! 2023.12.25
   d = be5_ey(2,j) ! By_ey
   det = a*d - b*c
   bi(1,1:2)=(/ d, -b/)
   bi(2,1:2)=(/ -c, a/)
   bi = bi/det
   e(1,1)= be5_ex(3,j) ! (bz_ex)
   e(2,1)= be5_ey(3,j) ! (bz_ey)
   txy = matmul(bi,e)
   resp_tip%tx(j) = txy(1,1) ! [nT]/[nT]
   resp_tip%ty(j) = txy(2,1)
  end do
  
  write(*,'(a,i2)') " ### CALRESPTIP     END !! ###  ip =",ip  ! 2023.12.23
  return
  end


!######################################## function phase
real(8) function phase(c) ! [deg]
implicit none
complex(8),intent(in) :: c
real(8),parameter :: pi=4.d0*datan(1.d0), r2d=180.d0/pi
 phase=datan2(dimag(c),dreal(c))*r2d
 return
end function phase
!######################################## function amp
real(8) function amp(c)
implicit none
complex(8),intent(in) :: c
!real(8) :: amp
 amp=dsqrt(dreal(c)**2.d0 + dimag(c)**2.d0)
 return
end function


!#############################################
!# modified on 2025.04.21 based on ALLOCATERESP_MT in n_inv_joint.f90
!# copied from ../solver/n_ebfem_bxyz.f90 on 2021.09.14
subroutine ALLOCATERESP_MT(nobs,nsr,resp,resp_mt,resp_tip,ip,nfreq)
use outresp
use param
implicit none
integer(4),         intent(in)    :: nobs
integer(4),         intent(in)    :: nsr ! 2017.07.11
integer(4),         intent(in)    :: nfreq,ip
type(respdata),     intent(inout) :: resp(5,nsr,nfreq) !2017.07.11
type(respmt),       intent(inout) :: resp_mt(nfreq)     !2021.09.14
type(resptip),      intent(inout) :: resp_tip(nfreq)   ! tippers for each freq
integer(4)                        :: i,j,k

do j=1,nfreq
 CALL ALLOCATERESPMT(  nobs,resp_mt(j)    ) ! 2021.09.14 m_outresp.f90
 CALL ALLOCATERESPTIP( nobs,resp_tip(j)   ) ! 2023.12.23 m_outresp.f90
 do i=1,5
   do k=1,nsr ! 2017.07.11
     CALL ALLOCATERESPDATA(nobs,resp(  i,k,j)) ! 2017.07.11
   end do     ! 2017.07.11
 end do
end do

if( ip .eq. 0) write(*,*) "### ALLOCATERESP_MT END!! ###"
return
end

!################################## subroutine setNPHYS1INDEX2COND
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
if ( n4flag(i,1) .ge. 3 ) nphys2 = nphys2 + 1 ! Modified on 2026.07.29 to adjust TMTGEM mesh
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
 CALL CALOBSRESP(fs(:,i),nline,coeffobs(1,1),resp5(4,i)  ) !ex,fp deleted 2021.09.15
 CALL CALOBSRESP(fs(:,i),nline,coeffobs(1,2),resp5(5,i)  ) !ey,fp deleted 2021.09.15
end do      ! 2017.07.11

!write(*,*) "### CALOBSEBCOMP END!! ###" ! 2017.07.12

return
end
!######################################## CALOBSRESP
!# Coded on Nov. 21, 2015
!# This calculates the output b fields and output results
subroutine CALOBSRESP(ft,nline,coeffobs,resp) ! fp, fs -> ft 2021.09.15
use matrix
use outresp
implicit none
integer(4),           intent(in)    :: nline
complex(8),           intent(in)    :: ft(nline) ! 2021.09.15
type(real_crs_matrix),intent(in)    :: coeffobs ! see m_matrix.f90
type(respdata),       intent(inout) :: resp   ! see m_outresp.f90
complex(8),allocatable,dimension(:) :: ftobs
real(8) :: amp,phase
integer(4) :: i
allocate(ftobs(resp%nobs)) ! 2021.09.15

!#[1]## generate btotal

!#[2]## calculate bp,bs,bt at observation points
CALL mul_matcrs_cv(coeffobs,ft(1:nline),nline,ftobs) ! see m_matrix.f90

!#[3]## cal b3 comp and output
do i=1,resp%nobs
 resp%ftobsamp(i)  =amp  (ftobs(i)) ! amp of bz
 resp%ftobsphase(i)=phase(ftobs(i)) ! phase of bz
 resp%ftobs(i)     =ftobs(i)        ! 2021.09.15
end do

!write(*,*) "### CALOBSRESP END!! ###"
return
end
!####################################################### PREPOBSCOEFF
! copied from ../solver/n_ebfem_bxyz.f90 on 2021.09.15
subroutine PREPOBSCOEFF(g_param,h_mesh,l_line,coeffobs) ! 2021.12.15
use mesh_type
use line_type
use param_mt ! 2021.12.15
use matrix
use fem_edge_util
implicit none
type(mesh),              intent(in)  :: h_mesh
type(line_info),         intent(in)  :: l_line
type(param_forward_mt),  intent(in)  :: g_param ! 2021.12.15
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

!###################################################################
! modified for spherical on 2016.11.20
! iflag = 0 for xyz
! iflag = 1 for xyzspherical
subroutine GENXYZMINMAX_MT(em_mesh,g_param_mt)
use param_mt ! 2016.11.20
use mesh_type
implicit none
type(mesh),            intent(inout) :: em_mesh ! 2021.10.13
type(param_forward_mt),intent(inout) :: g_param_mt ! 2021.12.15
real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
real(8) :: xyzminmax(6)
real(8), allocatable, dimension(:,:) :: xyz ! 2025.07.15
integer(4) :: i
allocate(xyz(3,em_mesh%node))   ! 2025.07.15
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

write(*,*) "xmin,xmax",xmin,xmax ! 2021.10.13
write(*,*) "ymin,ymax",ymin,ymax ! 2021.10.13
write(*,*) "zmin,zmax",zmin,zmax ! 2021.10.13

xyzminmax(1:6)=(/xmin,xmax,ymin,ymax,zmin,zmax/)

!# set output
g_param_mt%xyzminmax = xyzminmax ! 2021.12.15
em_mesh%xyzminmax = xyzminmax    ! 2021.10.13

write(*,*) "### GENXYZMINMAX END!! ###"
return
end

!################################################################# PREPZOBSMT
! modified based on PREPZOBSSRC in solver_mpi/n_ebfem_bxyz_mpi.f90 on 2026.07.29
subroutine PREPZOBSMT(h_mesh,g_param_mt)
use param_mt
use mesh_type
use triangle
use outresp
implicit none
type(mesh),             intent(in)    :: h_mesh
type(param_forward_mt), intent(inout) :: g_param_mt
type(grid_list_type)                  :: glist
integer(4)                            :: nobs,nx,ny
real(8),   allocatable,dimension(:,:) :: xyzobs,xyz
integer(4),allocatable,dimension(:,:) :: n3k
real(8),   allocatable,dimension(:)   :: znew
real(8)    :: a3(3)
integer(4) :: iele,n1,n2,n3,j,k,ntri,idev
integer(4) :: nsr                                ! 2017.07.11
real(8),allocatable,dimension(:,:)    :: xs1,xs2 ! 2017.07.18
real(8)    :: xyzminmax(6)                       ! 2017.07.18

!#[0]## cal xyzminmax of h_mesh
!  CALL GENXYZMINMAX(h_mesh,g_param)  ! commented out  2017.10.12

!#[1]## set
nsr       = s_param%nsource     ! 2017.07.11
allocate(xs1(3,nsr),xs2(3,nsr)) ! 2017.07.11
allocate(xyz(3,h_mesh%node),n3k(h_mesh%ntri,3))
allocate(xyzobs(3,g_param_mt%nobs))
allocate(znew(g_param_mt%nobs))
nobs      = g_param_mt%nobs
xyz       = h_mesh%xyz    ! triangle mesh
n3k       = h_mesh%n3
ntri      = h_mesh%ntri
xyzobs    = g_param_mt%xyzobs
xs1       = s_param%xs1
xs2       = s_param%xs2
xyzminmax = g_param_mt%xyzminmax


!#[2]## cal z for nobsr
nx=1000;ny=1000
CALL allocate_2Dgrid_list(nx,ny,ntri,glist)   ! see m_mesh_type.f90
write(*,*) "before gen2Dgridforlist in PREPZSRCOBS"
CALL gen2Dgridforlist(xyzminmax,glist) ! see m_mesh_type.f90
CALL classifytri2grd(h_mesh,glist)   ! classify ele to glist,see


!#[3] search for the triangle including (x1,y1)

idev = free_unit()
open(idev,file=trim(g_param%outputfolder)//"site.dat")! 2026.07.29
do j=1,nobs
    call findtriwithgrid(h_mesh,glist,xyzobs(1:2,j),iele,a3)
    n1 = n3k(iele,1); n2 = n3k(iele,2) ; n3 = n3k(iele,3)
    znew(j) = a3(1)*xyz(3,n1)+a3(2)*xyz(3,n2)+a3(3)*xyz(3,n3) + xyzobs(3,j)
!
    write(*,*) "------------------------------------------"
    write(*,'(i5,a,a)') j," site ID : ",trim(g_param%obsname(j))
    write(*,'(a,2f12.5,a,f12.5,a)') "lon  lat",g_param%lonlataltobs(1:2,j), " alt",g_param%lonlataltobs(3,j)," [km]"
    write(*,'(a,2f12.5,2(a,f12.5),a)') "  x    y",xyzobs(1:2,j), "   z",xyzobs(3,j)," -> ",znew(j)," [km]"
    write(idev,'(i5,5f12.5,a)') j,g_param%lonlataltobs(1:2,j),xyzobs(1:3,j),trim(g_param%obsname(j))
end do
close(idev)
    write(*,*)"-------------------------------------------"
 

!#[4]## set znew to xyz_r
    g_param_mt%xyzobs(3,1:nobs) = znew(1:nobs)
 

write(*,*) "### PREPZOBSMT END ###"
return
end