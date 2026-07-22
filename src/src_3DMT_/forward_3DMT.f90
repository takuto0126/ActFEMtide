!# modified for 3DMT for 2021.09.14
!# modified for multiple source on 2017.07.11
!# Coded on Sep 14, 2021
subroutine forward_3DMT(A,g_param_mt,h_mesh,l_line,nline,al_MT,freq,g_cond,g_surface,ip)! 2021.12.15
use mesh_type
use iccg_var_takuto
use line_type
use param        ! 2021.12.15
use param_mt     ! 2021.12.14
use constants    ! added on 2016.10.17
use surface_type ! 2021.09.14
implicit none
!include "../mesh_for_FEM/meshpara.f90" ! zmin, zmax, xout, yout! commented out on201610.17
!--------------- input and output variants ------------------
integer(4),            intent(in)        :: ip
integer(4),            intent(in)        :: nline
real(8),               intent(in)        :: freq
type(mesh),            intent(in)        :: h_mesh       ! see m_mesh_type.f90
type(line_info),       intent(in)        :: l_line       ! see m_line_type.f90
type(param_forward_mt),intent(in)        :: g_param_mt   ! 2021.12.15
type(param_cond)  ,    intent(in)        :: g_cond
type(surface),         intent(in)        :: g_surface(6) ! 2021.09.14
complex(8),            intent(out)       :: al_MT(nline,2) ! [mV/km * km] 1:ex, 2:ey polari 
type(global_matrix),   intent(inout)     :: A ! see m_iccg_var_takuto.f90
!--------------- internal variants
complex(8),allocatable,dimension(:,:) :: Avalue_bc ! 2021.09.14
complex(8),allocatable,dimension(:,:) :: b_vec     ! right hand side vector
logical,   allocatable,dimension(:,:) :: line_bc   ! 2021.09.14
integer(4) :: i,j,k,nsr=2
real(8)    :: omega,zmin,zmax,xout,yout
allocate( b_vec(nline,2),Avalue_bc(nline,2),line_bc(  nline,2) )

!#[1]## set
zmin = g_param_mt%zbound(1) ! 2021.12.15
zmax = g_param_mt%zbound(4) ! 2021.12.15
xout = g_param_mt%xbound(4) ! 2021.12.15
yout = g_param_mt%ybound(4) ! 2021.12.15
!write(*,*) "### forward_3DMT start nline"
write(*,*) zmin,zmax,xout,yout
!write(*,*) "ip",ip

!#[3]## SET Coefficient matrix and Generate Matrix

!-------- start frequency loop ----------- start frequency loop--------
omega=2.d0*pi*freq
!write(*,*) "freq=",freq,"[Hz]"

!#[4]## Initialize the matrices and vectors
write(*,*) "nline",nline
!write(*,*) "size of A%D",size(A%D)
!write(*,*) "size of A%AU",size(A%AU)
!write(*,*) "size of A%AL",size(A%AL)
!write(*,*) "size of al_MT",size(al_MT)
!write(*,*) "size of b_vec",size(b_vec)
CALL INITIALIZE(A,b_vec,al_MT,nline,2)! nsr=2 for MT, see m_iccg_var_takuto.f90 2021.09.16

!#[4-2]## Generate Matrix for CRS format : b_vec is set to zero
CALL GENMAT_MT(h_mesh,l_line,A,b_vec,omega,g_cond) ! ok 2021.12.15

!#[4-3]## Copy upper triangle to lower driangle
CALL COPY_UL_ICCG12(A,ip) ! ok 2021.09.14

!#[5]## Dirichlet boundary at calculation boundaries
!CALL GENBCCSEM(zmin,zmax,xout,yout,l_line,h_mesh,Avalue_bc,line_bc)
CALL GENBCMT(g_surface,omega,nline,Avalue_bc,line_bc) ! ok 2021.09.14

!#[6]## Set Boundary Condition for dirichlet boundary; Avalue_bc -> A and b_vec
CALL SET_BC_3Djoint(A, nline, nsr, b_vec, Avalue_bc, line_bc(:,1), ip) ! ok 2021.09.14

open(1,file="bc.dat")
do i=1,nline
 if (line_bc(i,1) ) write(1,'(i6,4g15.7)'),i,Avalue_bc(i,1:2)
end do
close(1)

!#[7]## Solve
!call solveMUMPS(doftot,A,b_vec,bs,ip)  ! for MacbookPro 15inch
call solvePARDISO(nline,nsr,A,b_vec,al_MT,ip) !　2017.07.11

write(*,*) "### forward_bxyz END !! ip=",ip,"freq=",freq,"###"
return
end subroutine forward_3DMT

!#######################################
!coded on 2021.09.14
! BC is homogeneous 1 mV/km in x (Ex polarization) and y (Ey polarization) at top surface
subroutine GENBCMT(g_surface,omega,nline,Avalue_bc,line_bc)
use surface_type
use line_type
use outerinnerproduct
implicit none
! 1:top,2:north,3:west,4:south,5:east,6:bottom
integer(4),   parameter     :: nsr = 2 ! ex and ey polarization
real(8),      intent(in)    :: omega
type(surface),intent(in)    :: g_surface(6)
integer(4),   intent(in)    :: nline
complex(8),   intent(inout) :: Avalue_bc(nline,nsr)
logical,      intent(inout) :: line_bc(  nline,nsr)
integer(4)                  :: iline,isr,j,l
real(8)                     :: xy1(3),xy2(3),dxy(3),xhat(3),yhat(3),dl
integer(4)                  :: np1,np2,nline_surface
real(8),     allocatable,dimension(:,:) :: x1x2_face
integer(4),  allocatable,dimension(:)   :: ilineface_to_3D
complex(8)                  :: iunit=(0.d0,1.d0)

!#[0]## initialization
Avalue_bc(1:nline,1:2) = 0.d0
line_bc(  1:nline,1:2) = .false.
xy1=0.d0; xy2=0.d0; dxy=0.d0
xhat=(/1.,0.,0./)
yhat=(/0.,1.,0./)

!#[1]# set Avalue for top
j=1
x1x2_face       = g_surface(1)%x1x2_face
nline_surface   = g_surface(1)%nline
ilineface_to_3D = g_surface(1)%ilineface_to_3D ! allocate and fill

do l=1,nline_surface
 iline = ilineface_to_3D(l)
 np1 = g_surface(1)%line(1,l)
 np2 = g_surface(1)%line(2,l)
 xy1(1:2) = g_surface(1)%x1x2_face(1:2,np1)
 xy2(1:2) = g_surface(1)%x1x2_face(1:2,np2)
 dxy = xy2 - xy1
 dl=dsqrt( dxy(1)**2. + dxy(2)**2. )
 ! E = - i*omega*A -> A = i /omega * E
 ! Al =  i/omega * (x_2 - x_1) *E
 Avalue_bc(iline,1)=iunit/omega * inner(dxy,xhat)*1. ! ex polarization
 Avalue_bc(iline,2)=iunit/omega * inner(dxy,yhat)*1. ! ey polarization
 line_bc(iline,1:2) = .true.
end do
deallocate(ilineface_to_3D)

! ex -> only north(surface 2) and south (surface 4)
! ey -> only west (surface 3) and east  (surface 5)

!#[2]# set Avalue_bc at 2~5 surfaces for MT source

 do j=2,5 ! side surfaces
 nline_surface   = g_surface(j)%nline
 ilineface_to_3D = g_surface(j)%ilineface_to_3D ! allocate and fill

 if ( j == 2 .or. j == 4 ) isr = 1 ! ex polarization
 if ( j == 3 .or. j == 5 ) isr = 2 ! ey polarization

 do l=1,nline_surface
  iline = ilineface_to_3D(l)
  ! Bs is El [mV/km*km]=[mV] -> Al = i/omega * El   2021.09.14
  Avalue_bc(iline,isr) = iunit/omega*g_surface(j)%bs(l) ! only for eather of ex or ey
  line_bc(iline,:)     = .true.                         ! both for ex and ey polarization
 end do !

 deallocate(ilineface_to_3D)
 end do ! j surface loop

!#[3]## bottom
j=6
nline_surface   = g_surface(j)%nline
ilineface_to_3D = g_surface(j)%ilineface_to_3D ! allocate and fill

do l=1,nline_surface
 iline = ilineface_to_3D(l)
 Avalue_bc(iline,isr) = 0.d0    ! only for eather of ex or ey
 line_bc(iline,:)     = .true.  ! both for ex and ey polarization
end do !


return
end

!#################################################################### GENMAT_MT
! modified on 2021.09.14 for 3DMT
! Coded on March 4, 2016
! replace integration by tetrai_table by analytical integration
subroutine GENMAT_MT(h_mesh,l_line,A,b_vec,omega,g_cond) ! 2021.12.14
use  outerinnerproduct
use  iccg_var_takuto ! b_vec is not included, see m_iccg_var_takuto.f90
use  mesh_type       ! see m_mesh_type.f90
use  line_type       ! see m_line_type.f90
use  fem_util        ! for volume, intv, (see m_fem_utiil.f90 )
use  fem_edge_util   ! see fem_edge_util.f90
use  param
!use  m_param_ana, only:cond,istructure ! see m_param_ana.f90 commented out 2021.07.17
use  constants,   only:pi,dmu          ! see m_constants.f90, 2017.07.11
implicit none
type(mesh),         intent(in)      :: h_mesh
type(line_info),    intent(in)      :: l_line
type(param_cond),   intent(in)      :: g_cond
type(global_matrix),intent(inout)   :: A
real(8),            intent(in)      :: omega
complex(8),         intent(out)     :: b_vec(l_line%nline,2)!Ex, Ey 2021.12.14
real(8)                             :: elm_xyz(3,4),xx(3,6), gn(3,4)
real(8)                             :: w(6,3), S(6,6), v, sigma,sigma0, yy
complex(8)                          :: iunit=(0.d0, 1.d0), rhs1
complex(8),         dimension(6,6)  :: elm_k, S1
integer(4),         dimension(6)    :: table_dof_elm, idirection
real(8),allocatable,dimension(:,:)  :: x3s,x3e
integer(4) :: iele, i, j, k, l, m, n, ii, jj, id_group
!---------------  scales ------------------------------------------------------
real(8), parameter                  :: L0=1.d+3  ! [m]  scale length
real(8)                             :: AA, a1(4),a2(4),sigma_bell,localPQ(3),x3p2(3),x3p1(3)
complex(8)                          :: BB
complex(8)                          :: b3(3,4),bl(6)
logical                             :: itrue
!real(8), allocatable,dimension(:,:) :: x3s,x3e  ! 2017.07.11
!real(8),             dimension(3)   :: x3p1,x3p2,localPQ,x1,x2,x3

!#[1] ## left-hand side matrix
do iele=1, h_mesh%ntet  ! start elemetn loop
!#
  !# [1] ## ! check the direction of edge, compared to the defined lines
  idirection(1:6)=1
  do j=1,6
    if ( l_line%n6line(iele, j) .lt. 0 ) idirection(j)=-1
  end do

  !# [2] ## Prepare the coordinates for 4 nodes of elements
  do j=1,4
   elm_xyz(1:3,j)=h_mesh%xyz(1:3,h_mesh%n4(iele,j)) ! [km]
  end do
  ! [x_mn]^T=L{x'34 x'14 x'42 x'23 x'31 x'12}=L[x'_lm]^T
  call calxmn(elm_xyz,xx)         ! see fem_util.f90
  call gradnodebasisfun(elm_xyz,gn,v) ! see fem_util.f90

  !# [3] ## First term from the rot rot, S
  ! [ int{ (rot w) (rot w)^T }dv ]{Bsl}
  ! Since rot w =1/3/v*x_mn, 
  !  int (rot w) cdot (rot w) dv = 1/9/v*(x_mn cdot x_m'n')
  AA=1.d0/9.d0/v
  S(:,:)=0.d0
  do j=1,6
    do k=1,6
	S(j,k)=inner(xx(1:3,j),xx(1:3,k))*idirection(j)*idirection(k)*AA  ! S is real
    end do
  end do   ! S [km*rad/s*S/m]

  !# [4] ## Second term from i * omega * mu * sigma * int{ sigma w cdot w }dv {Bsl}
  !# [4-1] ## assemble coefficient for i * omega*
  j = h_mesh%n4flag(iele,1)                 ! 2017.09.29
  if ( j .eq. 1 ) sigma=g_cond%sigma_air    ! 2017.09.29
  if ( j .ge. 2 ) then                      ! 2017.09.29
   if ( g_cond%condflag .eq. 0  )     sigma = g_cond%sigma_land(j-1) ! 2017.09.29
   if ( g_cond%condflag .eq. 1  ) then ! condflag = 1 -> file conductivity
    sigma = g_cond%sigma(iele - g_cond%nphys1) ! sigma store only nphys=2 element
   end if
  else if ( h_mesh%n4flag(iele,1) .ge. 3 ) then
    write(*,*) "GEGEGE h_mesh%n4flag(iele,1) = ",h_mesh%n4flag(iele,1)
    stop
  end if

!  if ( h_mesh%n4flag(iele,2) .ge. 4 ) sigma=0.01d0   ! ocean
  BB =iunit*omega*dmu*sigma *L0**2.d0 ! BB is complex in ww, dn_dx twice

  !# [4-2] ## assemble scheme No.2 ( analytical assembly)
  ! Since \nabla lambda_k =1/6/V*( x_ln \times x_lm )  is constant,
  ! int { w_i cdot w_j }dv can be calculated analytically
  S1(:,:)=(0.d0,0.d0)
  ! yy = int { w_i cdot w_j }dv, where w is vector shape function
  !      = int [n_k*gn(:,l) - n_l*gn(:,k) ] cdot [n_m*gn(:,n) - n_n*gn(:,m) ] dv
  !      = int ( n_k*n_m ) dv [ gn(:,l) cdot gn (:,n) ]      first
  !      - int ( n_k*n_n ) dv [ gn(:,l) cdot gn (:,m) ]       second
  !      - int ( n_l*n_m ) dv [ gn(:,k) cdot gn (:,n) ]       third
  !      + int ( n_l*n_n ) dv [ gn(:,k) cdot gn (:,m) ]      forth
  do i=1,6
    do j=1,6
      k=kl(i,1);l=kl(i,2) ; m=kl(j,1) ; n=kl(j,2) ! gn*gn [km^-2], intv [km^3], yy[km]
      yy =     intv(k,m,v)*inner(gn(:,l), gn(:,n))   & ! first term
     &	-  intv(k,n,v)*inner(gn(:,l), gn(:,m))   & ! second term
     &      -  intv(l,m,v)*inner(gn(:,k), gn(:,n))   & ! third term
     &      +  intv(l,n,v)*inner(gn(:,k), gn(:,m))     ! forth term
     S1(i,j)= yy*idirection(i)*idirection(j)*BB  ! S1 [km*rad/s*S/m]
    end do
  end do

  !# [5] ## Construct elemnt matrix, elm_k
  elm_k(:,:)=S(:,:)+S1(:,:) ! elm_k (complex), S(real)

  !# [6] ## Set global matrix from elm_k
  do i=1,6
    table_dof_elm(i)=l_line%n6line(iele,i)*idirection(i) ! make n6line positive
  end do
  CALL sup_iccg(elm_k,table_dof_elm,6,A%D,A%INU,A%IAU,A%AU,l_line%nline,A%iau_tot)! see m_iccg_var_takuto.f90 2021.09.16

end do    ! element loop end

b_vec(:,:)=0.d0 ! 2021.12.15

write(*,*) "### GENMAT END !! ###"

!  do i=1,l_line%nline
!   if (b_vec(i) .ne. 0.d0) write(*,*) i,"b=",b_vec(i)
!  end do
!  stop

return
end subroutine GENMAT_MT
!################################################### checkvalues
subroutine checkcoeff(elm_xyz,xx,gn,elm_k,S,dof1,S1,dof2,dof3,v)
implicit none
integer(4),intent(in) :: dof1,dof2,dof3
real(8),intent(in) :: v
real(8),intent(in) :: elm_xyz(3,4)
real(8),intent(in) :: xx(3,6)
real(8),intent(in) :: gn(3,4)
real(8),intent(in)    ::     S(dof1,dof1)
complex(8),intent(in) :: elm_k(dof1,dof1)
complex(8),intent(in) ::    S1(dof2,dof3)
integer(4) :: i,j,k
! elm_xyz(3,4)
write(*,*) "elm_xyz="
write(*,'(3g15.7)') ((elm_xyz(i,j),i=1,3),j=1,4)
! xx(3,6)
write(*,*)"xx"
write(*,'(3g15.7)') ((xx(i,j),i=1,3),j=1,6)
write(*,*) "volume=",v
! gn(3,4) (real)
write(*,*) "gn"
write(*,'(3g15.7)')((gn(i,j),i=1,3),j=1,4)
! S(4,4) (real)
write(*,*) "S"
do j=1,dof1
 write(*,*) (S(j,k),k=1,dof1)
end do
! elm_k(4,4) (complex)
write(*,*) "elm_k"
do j=1,dof1
 write(*,*) (elm_k(j,k),k=1,dof1)
end do
! S1 (complex)
write(*,*) "S1"
do i=1,dof2
 write(*,*) (S1(i,j),j=1,dof3)
end do
end

!############################################### checksourceelement
! Coded by T. MINAMI on May 10, 2016
! confirmed the sobroutine works well.
subroutine checksourceelement(x3s,x3e,elm_xyz,itrue,x3p1,x3p2)
use fem_edge_util
use fem_util
use outerinnerproduct
implicit none
real(8),intent(in) :: x3s(3),x3e(3) ! P: start point, Q: end point
real(8),intent(in) :: elm_xyz(3,4)
logical,intent(out) :: itrue ! the wire penetrating/touching source wire
real(8),intent(out) :: x3p1(3),x3p2(3)
real(8) :: PQ(3),OP(3),AP(3),OS(3),ABAC(3),uABAC(3),AS(3),BS(3),CS(3)
real(8) :: AB(3),AC(3),BC(3)
real(8) :: t,tdenom,lambda_A,lambda_B,lambda_C,x3p(3,2)
real(8) :: a(4),x3(3)
integer(4) :: iface,nA,nB,nC,i,j
logical,dimension(2) :: innerflag,onedgeflag

itrue=.false.
innerflag(1:2)=.false.
onedgeflag(1:2)=.false.
OP(1:3)=x3s(1:3)
PQ(1:3)=x3e(1:3) - x3s(1:3) ! start to end vector

!#[0]## inner or outer
call nodebasisfun(elm_xyz,x3s,a)
if ( 0.d0 .lt. a(1) .and. 0.d0 .lt. a(2) .and. 0.d0 .lt. a(3) .and. 0.d0 .lt. a(4)) then
 innerflag(1)=.true.
 itrue=.true.
 x3p1(1:3)=x3s(1:3)
end if
call nodebasisfun(elm_xyz,x3e,a)
if ( 0.d0 .lt. a(1) .and. 0.d0 .lt. a(2) .and. 0.d0 .lt. a(3) .and. 0.d0 .lt. a(4)) then
 innerflag(2)=.true.
 itrue=.true.
 x3p2(1:3)=x3e(1:3)
end if

!#[1]# check whether each face touch/penetrate PQ or not

do iface=1,4
 nA=lmn(iface,1) ! see fem_edge_util for lmn
 nB=lmn(iface,2)
 nC=lmn(iface,3)

 AP(1:3)=OP(1:3) - elm_xyz(1:3,nA)
 AB(1:3)=elm_xyz(1:3,nB) - elm_xyz(1:3,nA)
 AC(1:3)=elm_xyz(1:3,nC) - elm_xyz(1:3,nA)
 BC(1:3)=elm_xyz(1:3,nC) - elm_xyz(1:3,nB)
 ABAC(1:3)=outer(AB,AC) ! outer vector from the cell
 uABAC(1:3)=ABAC(1:3)/dsqrt(inner(ABAC,ABAC)) ! unit vector in ABAC direction

 tdenom=inner(PQ, ABAC) ! minus -> P is outer; plus ->  Q is outer
 if ( tdenom .eq. 0.d0 ) cycle ! PQ and iface is parallel
 t=inner(-AP, outer(AB,AC))/tdenom
 if ( t .lt. 0.d0 .or. 1.d0 .lt. t ) cycle ! not penetrate or touch the iface
 if ( tdenom .lt. 0.d0 ) onedgeflag(1) = .true. ! P is outer
 if ( tdenom .gt. 0.d0 ) onedgeflag(2) = .true. ! Q is outer
 OS(1:3)=OP(1:3)+t*PQ(1:3)     ! S is the point on plain, ABC, and 0=< t =<1
 AS(1:3)=OS(1:3)- elm_xyz(1:3,nA)
 BS(1:3)=OS(1:3)- elm_xyz(1:3,nB)
 CS(1:3)=OS(1:3)- elm_xyz(1:3,nC)
 lambda_A = inner(outer(BC,BS),uABAC)/inner(ABAC,uABAC)
 lambda_B = inner(outer(AS,AC),uABAC)/inner(ABAC,uABAC)
 lambda_C = inner(outer(AB,AS),uABAC)/inner(ABAC,uABAC)

 if (0.d0 .lt. lambda_A .and. 0.d0 .lt. lambda_B .and. 0.d0 .lt. lambda_C) then
  itrue=.true.
  x3(1:3)=    elm_xyz(1:3,nA)*lambda_A &
	   & +  elm_xyz(1:3,nB)*lambda_B &
	   & +  elm_xyz(1:3,nC)*lambda_C
  if ( tdenom .lt. 0.d0 ) x3p1(1:3)= x3(1:3)   ! P side point
  if ( tdenom .gt. 0.d0 ) x3p2(1:3)= x3(1:3)   ! Q side point
 end if
end do

!#[2]# check validity
 if (itrue) then
  if ( innerflag(1) .and. .not. innerflag(2))      then ! only P is in cell
   if (  .not. onedgeflag(1) .and. onedgeflag(2) ) goto 101
  else if ( innerflag(2) .and. .not. innerflag(1)) then ! only Q is in cell
   if ( .not. onedgeflag(2) .and. onedgeflag(1) ) goto 101
  else if ( innerflag(1) .and.  innerflag(2) )     then ! both P and Q are in cell
   if ( .not. onedgeflag(1) .and. .not. onedgeflag(2)) goto 101
  else if ( .not. innerflag(1) .and. .not. innerflag(2)) then ! both P and Q are out of cell
   if ( onedgeflag(1) .and. onedgeflag(2)  )            goto 101
  end if
  goto 100
 end if
 101 continue

!#[3] if itrue is ".true."
!if (itrue) then ! commented out on 2017.02.20
if (.false.) then
 write(*,*) "iture=",itrue
 write(*,*) "innerflag(1:2)=",innerflag(1:2)
 write(*,*) "onedgeflag(1:2)=",onedgeflag(1:2)
 write(*,'(a13,3g15.7)') ("elm_xyz(1:3)=",elm_xyz(1:3,j),j=1,4)
 write(*,*) "x3p1(1:3)",x3p1(1:3)
 write(*,*) "x3p2(1:3)",x3p2(1:3)
end if

!write(*,*) "### checksource element END!! ###"

return
100 continue
write(*,*) "GEGEGE innerflag(1:2)=",innerflag(1:2),"onedgeflag(1:2)=",onedgeflag(1:2)
write(*,'(a8,3g15.7)') ("elm_xyz=",elm_xyz(1:3,j),j=1,4)
write(*,*) "x3s=",x3s(1:3)
write(*,*) "x3e=",x3e(1:3)
stop
end
!########################################  GENBCCSEM
! Coded on October 19, 2015
! THis subroutine set the boundary condition for calculation boundaries
subroutine  GENBCCSEM(zmin,zmax,xout,yout,l_line,h_mesh,Avalue_bc,line_bc)
use mesh_type
use line_type
implicit none
real(8),   intent(in)  :: zmin,zmax,xout,yout ! These parameters are from meshpara.f90
type(mesh),     intent(in) :: h_mesh
type(line_info),intent(in) :: l_line
logical(4),intent(out) :: line_bc(l_line%nline)  ! if .true., the dirichlet bound set
complex(8),intent(out) :: Avalue_bc(l_line%nline) ! Dirichlet boundary value
integer(4),allocatable,dimension(:) :: line_group! 1:top,2:north,3:west,4:south, 5:east, 6:bottom
integer(4) :: i, j, ncount, n1, n2
real(8) :: x1, y1, z1, x2, y2, z2, xout1, yout1, zmin1, zmax1
!# For direct solvers
allocate(line_group(l_line%nline))

!#[0]## set boudary coordinates
xout1=xout-1.d0; yout1=yout-1.d0; zmin1=zmin+1.d0;  zmax1=zmax-1.d0

!#[1]## Initialize Avalue_bc and line_bc
line_bc(:)=.false.
Avalue_bc(:)=(0.d0, 0.d0)
line_group(:)=0  ! 1:top, 2:north, 3:west, 4:south, 5:east, 6:bottom

!#[2]## Initialize Avalue_bc and line_bc
do i=1,l_line%nline
     n1=l_line%line(1,i) ; n2=l_line%line(2,i)
     x1=h_mesh%xyz(1,n1) ; y1=h_mesh%xyz(2,n1) ; z1=h_mesh%xyz(3,n1)
     x2=h_mesh%xyz(1,n2) ; y2=h_mesh%xyz(2,n2) ; z2=h_mesh%xyz(3,n2)
     if         ( x1 .le. -xout1 .and.  x2 .le. -xout1) then
        line_group(i)=3
     else if  (x1 .ge. xout1 .and.  x2 .ge. xout1 ) then
        line_group(i)=5
     else if  ( y1 .le. -yout1 .and. y2 .le. -yout1 )  then
        line_group(i)=4
     else if  ( y1 .ge. yout1 .and. y2 .ge. yout1 ) then
        line_group(i)=2
     else if  ( z1 .le. zmin1 .and. z2 .le. zmin1 )  then
        line_group(i)=6
     else if  ( z1 .ge. zmax1 .and. z2 .ge. zmax1 ) then
	  line_group(i)=1
     end if
     if ( line_group(i) .ne. 0 ) then
	   line_bc(i) = .true.
	   Avalue_bc(i)=(0.d0, 0.d0)
     end if
end do    ! line loop

!#[3]## Output .msh file to confirm which lines are selected
CALL OUTBCLINES(line_bc,l_line%line,l_line%nline,h_mesh%xyz,h_mesh%node,line_group)

write(*,*) "### GENBCCSEM END!! ###"
return
end subroutine
!################################################  OUTBCLINES
! Coded on October 19,
subroutine OUTBCLINES(line_bc, line, nline, xyzg, nodeg, line_group)
implicit none
integer(4),intent(in) :: nline, nodeg, line(2,nline), line_group(nline)
real(8),intent(in) :: xyzg(3, nodeg)
logical, intent(in) :: line_bc(nline)
integer(4) :: i, nline_bc, icount
!#[1]## count the number of boundary lines
nline_bc=0
do i=1,nline
if ( line_bc(i) ) nline_bc=nline_bc+1
end do
!write(*,*) "nline_bc=",nline_bc 2017.07.25

!#[2]## Output line.msh
if (.false.) then
open(1,file="bcline.msh")
write(1,'(a)') "$MeshFormat"
write(1,'(a)')  "2.2 0 8"
write(1,'(a)')  "$EndMeshFormat"
write(1,'(a)')  "$Nodes"
write(1,*)  nodeg
do i=1,nodeg
write(1,*) i,xyzg(1:3,i)
end do
write(1,'(a)') "$EndNodes"
write(1,'(a)') "$Elements"
write(1,*)  nline_bc
icount=0
do i=1,nline
  if (line_bc(i)) then
    icount=icount+1
    write(1,*) icount, " 1 2 0", line_group(i), line(1,i), line(2,i)
  end if
end do
write(1,'(a)') "$EndElements"
close(1)
end if

return
end

!###############################################  SET_BC_3DJoint
! modified for 3DMT on 2021.09.14
! assumes the lines for dirichlet boundary are common for all the MT and CSEM problems
! Coded on Aug. 21, 2015
SUBROUTINE SET_BC_3DJoint(A,nline,nsr,rf,dirichlet_value,dirichlet,ip)
use iccg_var_takuto
implicit none
type(global_matrix), intent(inout) :: A
integer(4),          intent(in)    :: nline, ip, nsr             ! nsr= 2 for MT
complex(8),          intent(in)    :: dirichlet_value(nline,nsr) ! 2021.09.14
logical,             intent(in)    :: dirichlet(nline)           !whether dirichlet boundary
complex(8),          intent(inout) :: rf(nline,nsr)              ! 2021.09.16
integer(4)  :: i, j, ipos
complex(8)  ::  temp
!#[1]## Initialize
!rf(1:doftot)=(0.d0, 0.d0)
!doftot=nodtot
!dirichlet(:)=.false.
!dirichlet_value(:)=(0.d0,0.d0)
!do i=1,nline_bc
!  dirichlet(line_bc(i))=.true.
!  dirichlet_value(line_bc(i))=Avalue_bc(i)
!end do
!do i=1,doftot
!write(*,*) "dirichlet_value(i)=",dirichlet_value(i),"i=",i
!end do
!#[2]##  modify right-hand side vector
!     non-direchlet : {rf} = {rf} - [K]*{dirichlet_value}
!         direchlet : {rf} =            {dirichlet_value}
      do i=1,nline
        if (dirichlet(i)) then ! when line is for dirichlet boundary
	    do j=1,nsr                 ! j=1 for ex, j=2 for ey polarization
           rf(i,j) = dirichlet_value(i,j)
	    end do
        else
        do j=1,nsr  ! source type loop
          temp   =       A%D   (i   ) * dirichlet_value(i,j)
          do ipos=A%INL(i-1)+1,A%INL(i)
            temp = temp + A%AL(ipos) * dirichlet_value( A%IAL(ipos),j ) ! 2021.09.14
          enddo
          do ipos=A%INU(i-1)+1,A%INU(i)
            temp = temp + A%AU(ipos) * dirichlet_value( A%IAU(ipos),j ) ! 2021.09.14
          enddo
          rf(i,j) = rf(i,j) - temp ! 2017.07.11
        end do ! source loop end 2021.09.14
        endif
      enddo ! nline loop

!#[3]## modify coefficient matrix ## Note modification to A is common for all nsr souces
      do i=1,nline
        if (dirichlet(i)) then ! lines corresponding to boundary surfaces
              A%D   (i   ) = 1.d0
          do ipos=A%INL(i-1)+1,A%INL(i)
              A%AL(ipos) = 0.d0
          enddo
          do ipos=A%INU(i-1)+1,A%INU(i)
              A%AU(ipos) = 0.d0
          enddo
        else
          do ipos=A%INL(i-1)+1,A%INL(i)
            if (dirichlet(A%IAL(ipos))) then
              A%AL(ipos) = 0.d0
            endif
          enddo
          do ipos=A%INU(i-1)+1,A%INU(i)
            if (dirichlet(A%IAU(ipos))) then
              A%AU(ipos) = 0.d0
            endif
          enddo
        endif
      enddo
if (ip .eq. 0) write(*,*) "### SET_BC_3DMT END!! ###"
RETURN
END
!
!
