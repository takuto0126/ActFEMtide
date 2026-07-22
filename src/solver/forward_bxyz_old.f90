!# modified for multiple sources on Feb 22, 2018
!# Coded on March 4, 2016
subroutine forward_bxyz(A,h_mesh,l_line,doftot,nsr,bs,freq,sparam,g_param,g_cond,ip)
use mesh_type
use iccg_var_takuto
use line_type
use param
use constants ! added on 2016.10.17
implicit none
!include "../mesh_for_FEM/meshpara.f90" ! zmin, zmax, xout, yout! commented out on201610.17
!--------------- input and output variants ------------------
integer(4),intent(in) :: doftot,icomp,isource_for_B
integer(4),intent(inout) :: ip
real(8),intent(in)    :: freq
type(mesh),intent(in) :: h_mesh      ! see m_mesh_type.f90
type(line_info),intent(in) :: l_line ! see m_line_type.f90
type(param_forward),intent(in) :: g_param
type(param_source), intent(in) :: sparam
type(param_cond)  , intent(in) :: g_cond
complex(8),intent(in) :: bp(doftot)  ! [nT*km]
complex(8),intent(out):: bs(doftot)  ! [nT*km]
!--------------- internal variants
type(global_matrix),intent(inout) :: A ! see m_iccg_var_takuto.f90
complex(8),allocatable,dimension(:)  :: Avalue_bc
logical,   allocatable,dimension(:)  :: line_bc
complex(8),allocatable,dimension(:)  :: b_vec
integer(4) :: i,j,k,IB
real(8)    :: omega,zmin,zmax,xout,yout
allocate(Avalue_bc(doftot),line_bc(doftot),b_vec(doftot))

!#[1]## set
zmin = g_param%zbound(1)
zmax = g_param%zbound(4)
xout = g_param%xbound(4)
yout = g_param%ybound(4)

!#[3]## SET Coefficient matrix and Generate Matrix

!-------- start frequency loop ----------- start frequency loop--------
omega=2.d0*pi*freq
!write(*,*) "freq=",freq,"[Hz]"
!#[4]## Initialize the matrices and vectors
CALL INITIALIZE(A,b_vec,bs,doftot)

!#[4-2]## Generate Matrix for CRS format
CALL GENMAT(h_mesh,l_line,A,b_vec,omega,pi,dmu,bp,sparam,g_param,icomp,isource_for_B,g_cond)

!#[4-3]## Copy upper triangle to lower driangle
CALL COPY_UL_ICCG12(A,ip)

!#[5]## Generate Boundary Values
!#[5-1]## Dirichlet boundary at calculation boundaries
CALL GENBCCSEM(zmin,zmax,xout,yout,l_line,h_mesh,Avalue_bc,line_bc,IB)

!#[6]## Set Boundary Condition for dirichlet boundary
!#[6-1]## for CG
CALL SET_BC_ICCG(A, doftot, b_vec, Avalue_bc, line_bc, ip)

!#[7]## Solve
!call solveCG(doftot,A,b_vec,bs) ! b_vec, and bs are complex(8)
!call solveMUMPS(doftot,A,b_vec,bs,ip)  ! for MacbookPro 15inch
call solvePARDISO(doftot,A,b_vec,bs,ip) !　for EIC

write(*,*) "### forward_bz END !! ###"
return
end subroutine forward_bxyz
!########################################
subroutine initialize(A,b_vec,bs,doftot)
use iccg_var_takuto ! global_matrix
implicit none
type(global_matrix) :: A
integer(4),intent(in) :: doftot
complex(8),intent(out) :: b_vec(doftot),bs(doftot)

A%D(:)=(0.d0,0.d0)
A%AU(:)=(0.d0,0.d0)
A%AL(:)=(0.d0,0.d0)
bs(:)=(0.d0,0.d0)
b_vec(:)=(0.d0,0.d0)

return
end

!------------------------------------------------------------------- solveCG
subroutine solveCG(doftot,A,rhs,X)
use iccg_var_takuto ! use AU,AL,D,IAU,IAL,INU,INL,iau_tot,ial_tot
implicit none
type(global_matrix),intent(in) :: A
integer(4),intent(in)  :: doftot
complex(8),intent(in)  :: rhs(doftot)
complex(8),intent(out) :: X(doftot) ! solution vector
!##### iccg parameter #####
integer(4) :: iter, error, nset
character(20) :: PRECOND
real(8) :: resid, sigma_diag, sigma
X(:)=(0.d0,0.d0)
PRECOND="SSOR"
RESID=1.e-4       ! convergence tol.
SIGMA_DIAG=1.0     ! diag. mul.
SIGMA=0.0          ! non-diag. mul.
ERROR=0            ! error id (inout)
nset=1             ! precondition 0: no, 1: yes, 2: scaling
ITER=3000
CALL  COCG13_single (doftot,A%ial_tot,A%iau_tot,                       &
                       & A%D, A%AL, A%INL, A%IAL, A%AU, A%INU, A%IAU,  &
                       & rhs, X, precond, sigma_diag, sigma,&
                       & resid,  iter, error,  NSET) ! bs is b_secondary
return
end
!#################################################################### GENMAT8
! Coded on March 4, 2016
! replace integration by tetrai_table by analytical integration
subroutine GENMAT(h_mesh,l_line,A,b_vec,omega,pi,dmu,bp,sparam,g_param,icomp,isource_for_B,g_cond)
use  outerinnerproduct
use  iccg_var_takuto ! b_vec is not included, see m_iccg_var_takuto.f90
use  mesh_type       ! see m_mesh_type.f90
use  line_type       ! see m_line_type.f90
use  fem_util        ! for volume, intv, (see m_fem_utiil.f90 )
use  fem_edge_util   ! see fem_edge_util.f90
use  param
use  m_param_ana, only:cond,istructure ! see m_param_ana.f90  
implicit none
integer(4),intent(in) :: icomp ! 1 for solving by B, 2 for solving by E
integer(4),intent(in) :: isource_for_B ! added on May 12, 2016
! isource_for_B= 1 for B=bp+bs, 2 for mu*int[i_s cdot (rot w)]dv
!                3 for int[(rot w)cdot(rot w)]{Bl0}, where Bl0 is analytical solution
! on May 12, 2016, all of isource_for_B cannot create accurate solution...orz
type(mesh),intent(in)             :: h_mesh
type(line_info),intent(in)        :: l_line
type(param_forward),intent(in)    :: g_param
type(param_source),intent(in)     :: sparam
type(param_cond),  intent(in)     :: g_cond
type(global_matrix),intent(inout) :: A
real(8),   intent(in) :: omega, pi, dmu
complex(8),   intent(in) :: bp(l_line%nline) ! [nT*km]
complex(8),intent(out):: b_vec(l_line%nline)
real(8)    :: elm_xyz(3,4),xx(3,6), gn(3,4)
real(8)    :: w(6,3), S(6,6), v, sigma,sigma0, yy
complex(8) :: elm_k(6,6), S1(6,6),S2(6,6)
complex(8) :: iunit=(0.d0, 1.d0), rhs1
integer(4) :: table_dof_elm(6)
integer(4) :: iele, i, j, k, l, m, n, ii, jj, idirection(6),id_group
!---------------  scales ------------------------------------------------------
real(8), parameter  :: L0=1.d+3  ! [m]  scale length
real(8)    :: AA, a1(4),a2(4)
complex(8) :: BB, BB2
complex(8) :: b3(3,4),bl(6)
real(8)    :: x1(3),x2(3),x3(3)
logical :: itrue
real(8) :: x3s(3),x3e(3),PQ(3),unitPQ(3),x3p1(3),x3p2(3),localPQ(3),sigma_bell

!  sigma_bell = g_param%sigma_bell ! [S/m]
  x3s(1:3)=sparam%xs1(1:3)  ! [km]
  x3e(1:3)=sparam%xs2(1:3)  ! [km]
  PQ(1:3)=x3e(1:3)-x3s(1:3) ! [km]
  unitPQ(1:3)=PQ(1:3)/dsqrt(inner(PQ,PQ))

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
  BB2=iunit*omega*dmu*sigma *L0**2.d0 ! BB is complex in ww, dn_dx twice

  !# [4-2] ## assemble scheme No.2 ( analytical assembly)
  ! Since \nabla lambda_k =1/6/V*( x_ln \times x_lm )  is constant,
  ! int { w_i cdot w_j }dv can be calculated analytically
  S1(:,:)=(0.d0,0.d0) ; S2(:,:)=(0.d0,0.d0)
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
     S2(i,j)= yy*idirection(i)*idirection(j)*BB2 ! S1 [km*rad/s*S/m]
    end do
  end do

  !# [5] ## Construct elemnt matrix, elm_k
  elm_k(:,:)=S(:,:)+S1(:,:) ! elm_k (complex), S(real), S2(complex)

  !# [6] ## Set global matrix from elm_k
  do i=1,6
    table_dof_elm(i)=l_line%n6line(iele,i)*idirection(i) ! make n6line positive
  end do
  CALL sup_iccg(elm_k,table_dof_elm,6,A%D,A%INU,A%IAU,A%AU,l_line%nline,A%iau_tot)


  !# [7] ## set right hand side vector, b_vec  ############################### Right Hand Side
  !# [7-1]## method for E=Ep+Es or B=Bp+Bs
  !######################################################### icomp = 1 and isource_for_B = 1
  if ( icomp .eq. 1 .and. isource_for_B .eq. 1) then
   do i=1,6
    ii=l_line%n6line(iele,i)*idirection(i)
    do j=1,6
     jj=l_line%n6line(iele,j)*idirection(j)
     b_vec(ii)=b_vec(ii) - S2(i,j)*bp(jj) ! - is necessary, b_vec [nT*km]
    end do
   end do

  !####################################################  not (icomp = 1 and isource_for_B = 1)
  else
  !# [7-1] ## method for js expressed by Heaviseide function
  call checksourceelement(x3s,x3e,elm_xyz,itrue,x3p1,x3p2) ! check the penetration
  if (itrue) then ! when the cell includes the wire current
  !CALL checkcoeff(elm_xyz,xx,gn,elm_k,S,6,S1,6,6,v)
  call nodebasisfun(elm_xyz,x3p1,a1)
  call nodebasisfun(elm_xyz,x3p2,a2)
  localPQ(1:3)=x3p2(1:3) - x3p1(1:3) ! [km]

  !############################################################  icomp = 2 (solve by A)
  if ( icomp .eq. 2) then ! solve by A (E)
    do i=1,6
     k=kl(i,1)
     l=kl(i,2)
     ii=l_line%n6line(iele,i)*idirection(i)
     rhs1= dmu*sparam%I*( &
     &   (a1(k)+a2(k))/2.d0 * inner(gn(1:3,l),localPQ) &
     &  -(a1(l)+a2(l))/2.d0 * inner(gn(1:3,k),localPQ) &
     &   )*idirection(i)*1.d+6 ! [V/m] -> [mV/km]
     b_vec(ii)=b_vec(ii) + rhs1
    end do

!########################################################### icomp = 1 and isource_for_B = 2
    !#[int(mu*i_s cdot (rot w))dv]
    else if (icomp .eq. 1 .and. isource_for_B .eq. 2) then
     do i=1,6
      ii=l_line%n6line(iele,i)*idirection(i)
      rhs1=dmu*sparam%I*inner(xx(1:3,i)/3./v,localPQ) &
		& *idirection(i)/L0*1.d+9 ! [T] to [nT]
      b_vec(ii)=b_vec(ii) + rhs1
     end do

!############################################################  icomp = 1 and isource_for_B = 3
    else if (icomp .eq. 1 .and. isource_for_B .eq. 3 ) then ! solve by B  ! Bl [nT * km]
     do i=1,4
      x3(1:3)= elm_xyz(1:3,i)* 1.d3 ! [m]
      x1(1:3)=x3p1(1:3)*1.d3       ! [m]
      x2(1:3)=x3p2(1:3)*1.d3       ! [m]
      CALL magfromwirecurrent(x3,x1,x2,sparam%I,20,b3(1:3,i),0.d0) ! b3 [nT]
     end do
     do i=1,6 !##  prepare analytical B near electric wirecurrent
      k=kl(i,1) ; l=kl(i,2)
      bl(i)=inner(elm_xyz(1:3,l) - elm_xyz(1:3,k),&
	&      dreal((b3(1:3,l) +  b3(1:3,k))/2.d0))*idirection(i)! [nT*km]
     end do
     do i=1,6
     !#[int(mu*(rot w) cdot (rot w))dv{Bl}]
     ii=l_line%n6line(iele,i)*idirection(i)
      do j=1,6
       jj=l_line%n6line(iele,j)*idirection(j)
       b_vec(ii)=b_vec(ii) + S(i,j)*bl(j) ! int{(rot w) cdot (rot w)}dv{bl} [nT*km]
      end do
     end do
    !    write(*,*) "iele=",iele,"i=",i,"rhs1=",rhs1

!###########################################################
    end if ! icomp =2 and icomp =1 + isource_for_B = 2 or 3 end
  end if ! itrue is .true. end
  end if ! whether icomp =1 and isource_for_B = 1 or not ( RHS bransh end)

end do ! element loop end
write(*,*) "### GENMAT END !! ###"

!  do i=1,l_line%nline
!   if (b_vec(i) .ne. 0.d0) write(*,*) i,"b=",b_vec(i)
!  end do
!  stop

return
end subroutine GENMAT
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
!############################################### COPY_UL_ICCG12
! Copied from static_linear/copy_ul_iccg.f90 on June 10, 2015
! Converted aiccg_u, aiccg_l, diag to complex on Aug. 21, 2015
!      subroutine  copy_ul_iccg12( doftot,item_u_tot,item_l_tot,istack_u,item_u,aiccg_u,istack_l,item_l,aiccg_l,ip)
     subroutine  copy_ul_iccg12( A,ip)
!    copy upper triangle term value to lower for iccg solver
	use iccg_var_takuto
      implicit  none
	type(global_matrix),intent(inout) :: A
	integer(4),intent(in) :: ip
      integer(4 ) ::  i_gl, j_gl
      integer(4)  ::  num , num2, i_pos
      logical     ::  found
	logical,allocatable,dimension(:) :: assigned
	integer(4),allocatable,dimension(:) :: order
	complex(8) :: test0=0.d0, test1=0.d0
	A%AL(:)=(0.d0,0.d0)
	allocate (assigned(A%ial_tot),order(A%ial_tot))
	assigned(:)=.false.
!  * copy upper triangle to lower triangle  (i_gl,j_gl) -> (j_gl,i_gl)
!
      do i_gl=1,A%doftot
         do num=A%INU(i_gl-1)+1,A%INU(i_gl)
              j_gl=A%IAU(num)
              found = .false.
!                                                 .. search term pos.
!              do num2=istack_l(j_gl-1)+1,istack_l(j_gl)
               do num2=A%INL(j_gl-1)+1,A%INL(j_gl)
!                if ( i_gl == item_l(num2) ) then
                if ( i_gl == A%IAL(num2) ) then
                  i_pos =  num2
                  found = .true.
                  exit
                endif
              enddo
              if ( .not.found ) then
                write(*,*) ' ***** error on copy_ul_iccg: dof not found', i_gl, j_gl
!                call geofem_abort(34,'internal logic error')
                stop
              endif
!                                                 .. set the term
!              aiccg_l(i_pos) = aiccg_u(num)
               A%AL(i_pos) = A%AU(num)
		   assigned(i_pos)=.true.
		   order(num)=i_pos
!		   test0=test0+A%AU(num)*dconjg(A%AU(num))
!		   test1=test1+A%AL(i_pos)*dconjg(A%AL(i_pos))
               if (cdabs(A%AL(i_pos)) .ne. cdabs(A%AU(num)) ) then
                write(*,*) "GEGEGE"
		    stop
		   end if
        enddo
      enddo
      do i_gl=1,A%ial_tot
	 if ( .not. assigned(i_gl) ) then
	  write(*,*) "GEGEGE not assigned i_gl=",i_gl
	  stop
	 end if
!	 test0=test0+A%AU(i_gl)
!	 test1=test1+A%AL(i_gl)
      end do
!	write(*,*) "test0=",test0
!	write(*,*) "test1=",test1
      if (ip .eq. 0 ) write(*,*) "### COPY_UL_ICCG12 END!! ### ip=",ip
      return
      end
!########################################  GENBCCSEM
! Coded on October 19, 2015
! THis subroutine set the boundary condition for calculation boundaries
subroutine  GENBCCSEM(zmin,zmax,xout,yout,l_line,h_mesh,Avalue_bc,line_bc,IB)
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
integer(4),intent(out) :: IB ! # of boundary lines
allocate(line_group(l_line%nline))

!#[0]## set boudary coordinates
xout1=xout-1.d0; yout1=yout-1.d0; zmin1=zmin+1.d0;  zmax1=zmax-1.d0

!#[1]## Initialize Avalue_bc and line_bc
line_bc(:)=.false.
Avalue_bc(:)=(0.d0, 0.d0)
line_group(:)=0  ! 1:top, 2:north, 3:west, 4:south, 5:east, 6:bottom

!#[2]## Initialize Avalue_bc and line_bc
IB=0
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
	   IB=IB+1
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
write(*,*) "nline_bc=",nline_bc

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

!###############################################  SET_BC_ICCG
! This program is based on set_bc_iccg.f90 in GeoFEM
! Coded on Aug. 21, 2015
SUBROUTINE SET_BC_ICCG(A,doftot,rf,dirichlet_value,dirichlet,ip)
use iccg_var_takuto
implicit none
type(global_matrix),intent(inout) :: A
integer(4), intent(in) :: doftot, ip
complex(8), intent(in), dimension(doftot) :: dirichlet_value
logical,    intent(in), dimension(doftot) :: dirichlet
complex(8), intent(out) :: rf(doftot)
integer(4) :: i, ipos
complex(8)      ::  temp
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
!#[2]##  modify right-hand vector
!     non-direchlet : {rf} = {rf} - [K]*{dirichlet_value}
!         direchlet : {rf} =            {dirichlet_value}
      do i=1,doftot
        if (dirichlet(i)) then
          rf(i) = dirichlet_value(i)
        else
          temp   =       A%D   (i   ) * dirichlet_value(        i     )
          do ipos=A%INL(i-1)+1,A%INL(i)
            temp = temp + A%AL(ipos) * dirichlet_value( A%IAL(ipos) )
          enddo
          do ipos=A%INU(i-1)+1,A%INU(i)
            temp = temp + A%AU(ipos) * dirichlet_value( A%IAU(ipos) )
          enddo
            rf(i) = rf(i) - temp
        endif
      enddo

!#[3]## modify stiffness matrix
      do i=1,doftot
        if (dirichlet(i)) then
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
if (ip .eq. 0) write(*,*) "### SET_BC_ICCG END!! ###"
RETURN
END
!
!
!############################################# sup_iccg
! Coppied from static_linear/sup_iccg.f90 on June 10, 2015
! here elm_k, duag, aiccg_u are complex
subroutine  sup_iccg ( &
&      elm_k, table_dof_elm, edof,                                       &
&      diag , istack_u   , item_u , aiccg_u, doftot, item_u_tot     )
!    superpose symmetric element stiffness matrix to gloval stiffness
!    for iccg solver (lower element, upper whole stiffness)
implicit  none
integer(4),intent(in) :: edof, doftot, item_u_tot
complex   (8),   intent(in)   ::  elm_k(edof,edof) !  element stiffness matrix
integer(4), intent(in)   :: table_dof_elm(edof)
!    digree of freedom number table at a element table_dof_elm(i)
!      where i    : dof order of 'elm_k'
complex(8), intent(inout) ::  diag(doftot)                  ! diagonal value   (i-th dof)
integer(4), intent(in) ::  istack_u(0:item_u_tot) !  last term count at each freedom  (i-th dof)
integer(4), intent(in) ::  item_u(item_u_tot)       ! freedom number at each term      (i-th term)
complex(8), intent(inout) ::  aiccg_u(item_u_tot)      ! upper triangular k value         (i-th term)
integer(4)  ::  i_gl, j_gl, i_elm, j_elm
integer(4)  ::  num , i_pos
logical        ::  found
! * superpose element stiffness matrix
!      store diagonal and upper term for packed whole stiffness matrix
!       ( element stiffness matrix stored only lower triangle )
!
!    k(i_gl,j_gl) = k(i_gl,j_gl) + elm_k(i_elm,j_elm)
!
!   . make diagonal and upper triangle
!
      do i_elm=1,edof
         i_gl = table_dof_elm(i_elm)
         do j_elm=1,edof
            j_gl = table_dof_elm(j_elm)
            if      ( i_gl  == j_gl  ) then            ! diagonal
              diag(i_gl) = diag(i_gl) + elm_k(i_elm,j_elm)
            else if ( i_gl  <  j_gl  ) then            ! upper triangle
              found = .false.
!                                                 .. search term pos.
              do num=istack_u(i_gl-1)+1,istack_u(i_gl)
                if ( j_gl == item_u(num) ) then
                  i_pos =  num
                  found = .true.
                  exit
                endif
              enddo
              if ( .not.found ) then
                write(*,*) ' ***** error on sup_iccg : dof not found  ', i_gl, j_gl
!                call geofem_abort(34,'internal logic error')
                stop
              endif
!                                                 .. superpose the term
              if ( i_elm >  j_elm )  then              ! elm. lower
                aiccg_u(i_pos) = aiccg_u(i_pos) + elm_k(i_elm,j_elm)
              else                                     ! elm. upper(sym)
                aiccg_u(i_pos) = aiccg_u(i_pos) + elm_k(j_elm,i_elm)
              endif
            endif
         enddo
      enddo
      return
      end subroutine  sup_iccg
