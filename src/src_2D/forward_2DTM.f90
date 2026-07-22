!  <coordinate system>
!    z
!    ^
!    |
!    |
!  x -------> y   z=0 is sea level
!  (x,y,z) [km]
!
!      3
! (3) /|
!  / / | ^
! L /  | | (2)
!  /---|
! 1    2
!   ->
!  (1)
subroutine forward_2DTM(g_surface,freq,g_param,g_cond,ip,iface)
use mesh_type
use iccg_var_takuto
use line_type
use param
use  constants,   only:pi,dmu          ! see m_constants.f90, 2017.07.11
use surface_type ! 2021.05.21
implicit none
!include "../mesh_for_FEM/meshpara.f90" ! zmin, zmax, xout, yout! commented out on201610.17
!--------------- input and output variants ------------------
integer(4),         intent(in)       :: ip
real(8),            intent(in)       :: freq
type(param_forward),intent(in)       :: g_param
type(param_cond)  , intent(in)       :: g_cond
type(surface),      intent(inout)    :: g_surface ! see m_surface_type.f90 2021.05.31
integer(4),         intent(in)       :: iface ! 2022.01.04
!--------------- internal variants
type(global_matrix)                  :: A         ! see m_iccg_var_takuto.f90
complex(8),allocatable,dimension(:)  :: Avalue_bc
logical,   allocatable,dimension(:)  :: line_bc
complex(8),allocatable,dimension(:,:):: b_vec    ! 2017.07.11
complex(8),allocatable,dimension(:,:):: bs       ! 2017.07.11
real(8)                              :: omega
integer(4)                           :: nsr = 1 ,i ! for 2DTM problem
integer(4)                           :: doftot
character(1) :: num ! 2022.01.04
!#[1]## set
 doftot = g_surface%nline
 A = g_surface%A  ! 2021.06.08
 write(*,*) "nsr=",nsr,"ip",ip
 allocate(Avalue_bc(doftot),line_bc(doftot))
 allocate(b_vec(doftot,nsr)) ! 2017.07.11
 allocate(bs(   doftot,nsr)) ! 2021.06.01
 omega  = 2.d0*pi*freq

!-------- start frequency loop ----------- start frequency loop--------
!write(*,*) "freq=",freq,"[Hz]"

!#[4]## Initialize the matrices and vectors
 CALL INITIALIZE(A,b_vec,bs,doftot,nsr) ! see m_iccg_var_takuto.f90 2021.09.16

 if (.false.) then
 write(num,'(i1)') iface
 open(1,file="sigma2D"//num//".dat") ! ok
 do i=1,g_surface%ntri
  write(1,'(i6,g15.7)'),i,g_surface%cond(i)
 end do
 close(1)

 open(1,file="x1x2_2D"//num//".dat")! ok
 do i=1,g_surface%node
  write(1,'(i6,2g15.7)'),i,g_surface%x1x2_face(1:2,i)
 end do
 close(1)

 open(1,file="n3line_2D"//num//".dat")! ok
 do i=1,g_surface%ntri
  write(1,'(i6,3i5)'),i,g_surface%n3line(i,1:3)
 end do
 close(1)
end if

 !#[5]## Generate Matrix for CRS format
 CALL GENMAT2DTM(g_surface,A,omega) ! 2021.05.21 see below

 if ( .false.) then
 open(1,file="D2D"//num//".dat") ! ok
 write(1,*) A%doftot, A%intdoftot
 write(1,*) A%iau_tot,A%ial_tot
 write(1,'(i6,2g15.7)') (i,A%D(i),i=1,A%doftot)
 close(1)

 open(1,file="AU2D"//num//".dat") ! ok
 write(1,*) A%doftot, A%intdoftot
 write(1,*) A%iau_tot,A%ial_tot
 write(1,'(i6,2g15.7)') (i,A%AU(i),i=1,A%iau_tot)
 close(1)
 end if

!#[6]## Copy upper triangle to lower driangle
 CALL COPY_UL_ICCG12(A,ip) ! see m_iccg_var_takuto.f90 2021.09.16

 if (.false. ) then
 open(1,file="AL2D"//num//".dat") ! ok
 write(1,*) A%doftot, A%intdoftot
 write(1,*) A%iau_tot,A%ial_tot
 write(1,'(i6,2g15.7)') (i,A%AL(i),i=1,A%ial_tot)
 close(1)

open(1,file="INU2D"//num//".dat") ! ok
 write(1,'(3i7)') (i,A%INU(i),A%INL(i),i=0,A%doftot)
close(1)

open(1,file="IAU2D"//num//".dat") ! ok
 write(1,'(3i7)') (i,A%IAU(i),A%IAL(i),i=1,A%iau_tot)
close(1)
 end if

!#[7]## Dirichlet boundary at calculation boundaries
 CALL GENBC2D(g_surface)
 Avalue_bc = g_surface%Avalue_bc
 line_bc   = g_surface%line_bc
! do i=1,g_surface%nline
!  if (line_bc(i) ) write(*,12) i," Avalue",Avalue_bc(i)," line",g_surface%line(1:2,i)
! end do
12 format(i5,a,2g15.7,a,2i5)

if (.false.) then
write(num,'(i1)') iface
open(1,file="bc2D"//num//".dat") ! ok
do i=1,doftot
 if (line_bc(i) ) write(1,'(i6,2g20.12)'),i,Avalue_bc(i)
end do
close(1)

open(1,file="bvecbefore"//num//".dat") ! ok
write(1,'(i4,2g15.7)') (i,b_vec(i,1),i=1,doftot)
close(1)
end if

!#[8]## Set Boundary Condition for dirichlet boundary
 CALL SET_BC_ICCG(A, doftot, nsr, b_vec, Avalue_bc, line_bc, ip,iface) !2017.07.11

 if ( .false.) then
 open(1,file="bvecin"//num//".dat") !  problem!!
 write(1,'(i4,2g15.7)') (i,b_vec(i,1),i=1,doftot)
 close(1)

 open(1,file="D2Dafter"//num//".dat")
 write(1,*) A%doftot, A%intdoftot
 write(1,*) A%iau_tot,A%ial_tot
 write(1,'(i6,2g15.7)') (i,A%D(i),i=1,A%doftot)
 close(1)

 open(1,file="AU2Dafter"//num//".dat")
 write(1,*) A%doftot, A%intdoftot
 write(1,*) A%iau_tot,A%ial_tot
 write(1,'(i6,2g15.7)') (i,A%AU(i),i=1,A%iau_tot)
 close(1)
 end if

 !#[9]## Solve
 call solvePARDISO(doftot,nsr,A,b_vec,bs,ip) !　2017.07.11

 if (.false.) then
 open(1,file="2Din"//num//".dat")
 write(1,'(i4,2g15.7)') (i,bs(i,1),i=1,doftot)
 close(1)
 end if

 if (.not. allocated(g_surface%bs)) allocate(g_surface%bs(g_surface%nline))
 g_surface%bs(:) = bs(:,1)   ! store solution

write(*,*) "### forward_2DTM END !! ###"! 2021.10.14

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
!###################################################################
! 2021.06.08
subroutine PREPAOFSURFACE(g_surface,n,ip)
use surface_type
use iccg_var_takuto
implicit none
integer(4),   intent(in)    :: n,ip  ! # of surface
type(surface),intent(inout) :: g_surface(n)
integer(4) :: nline,dofn,ntri,j

dofn  = 1

do j=1,n
nline = g_surface(j)%nline
write(*,'(a,i3,a,i7)') "Surface #",j,"  nline =",nline ! 2021.10.14
ntri  = g_surface(j)%ntri

allocate(g_surface(j)%table_dof(nline,dofn))

CALL SET_TABLE_DOF(dofn,nline,g_surface(j)%table_dof,nline) ! see below
CALL set_iccg_var7_dofn(dofn,3,nline,nline,nline,ntri,&
&    g_surface(j)%n3line,g_surface(j)%table_dof,g_surface(j)%A,ip)

end do !surface loop

return
end

!############################################################ GENMAT2DTM
!# 2021.05.25
!# coded on 2020.10.28 for 2D TM mode calculation
!
!    [ int (rot w) dot (rot w) dS - i*omega*mu*sigma* int w dot w dS ]El = 0 <- El [V/m * m]
!    v
!    [ int (rot' w') dot (rot'w') dS' - i*omega*mu*sigma*L0*int w' dot w' dS']El' = 0 <- El' [mV/km *km]
!
! -> A = int (rot' w) dot (rot'w) dS'
! -> B = i*omega*mu*sigma*L0*int w' dot w' dS'
!
! w [1/m]       = 10^-3 * w' = \lambda_i \nabla' lambda_j = [1/km]
! int*dS [m^2]  = 10^6  * int dS' [km^2]
! rot w [1/m^2] = 10^-6 * rot'w = [1/km/km]
! El [V] = 10^-6 * El' [mV/km *km]
! L0 = 10^6
!
!# Takuto Minami
subroutine GENMAT2DTM(g_surface,A,omega)
use  outerinnerproduct
use  iccg_var_takuto ! b_vec is not included, see m_iccg_var_takuto.f90
use  mesh_type       ! src_mesh/m_mesh_type.f90
use  line_type       ! m_line_type.f90
use  fem_util        ! for volume, intv, (see solver/m_fem_utiil.f90 )
use  fem_edge_util   ! fem_edge_util.f90
use  surface_type    ! src_2D/m_surface_type.f90
use  param
use  constants,   only:pi,dmu          ! see m_constants.f90, 2017.07.11
implicit none
type(surface),      intent(inout)   :: g_surface
type(global_matrix),intent(inout)   :: A
real(8),            intent(in)      :: omega
real(8)                             :: y(3),z(3)
real(8)                             :: s,v, sigma,sigma0, yy
complex(8)                          :: iunit=(0.d0, 1.d0), rhs1
complex(8),         dimension(3,3)  :: elm_k
real(8),            dimension(3,3)  :: A1,B1 !2021.05.21
integer(4),         dimension(3)    :: table_dof_elm, idirection
integer(4) :: iele, i, j, k, l, m, n, ii, jj, id_group,nsource,iline
!---------------  scales ------------------------------------------------------
complex(8)                          :: BB
real(8), parameter                  :: L0=1.d+3  ! [m]  scale length 2021.09.14
integer(4)                          :: i_edge,j_edge
real(8)                             :: grcgr_i(3),grcgr_j(3)

!#[0]##
  
!#[1] ## left-hand side matrix
do iele=1, g_surface%ntri  ! start element loop 2021.05.25
!#
  !# [1] ## ! check the direction of edge, compared to the defined lines
  idirection(1:3)=1 ! 2020.05.25
  do j=1,3
    if ( g_surface%n3line(iele, j) .lt. 0 ) idirection(j)=-1 ! 2021.05.25
  end do

!if ( iele .eq. 30 .or. iele .eq. 31 .or.  iele .eq. 50 ) then
! write(*,*) "iele=",iele
!do j=1,3
!iline = abs(g_surface%n3line(iele, j))
!write(*,11) "j",j,"line",g_surface%n3line(iele, j),"node",g_surface%line(1:2,iline),"idirection",idirection(j)
!end do
!end if
!11 format(a,i5,a,i5,a,2i5,a,i3)

!# [2] ## Prepare the coordinates for 3 nodes of elements ! 2021.05.25
  do j=1,3 ! 2021.05.25
    y(j)=g_surface%x1x2_face(1,g_surface%n3(iele,j)) ! [km] 2021.05.25
    z(j)=g_surface%x1x2_face(2,g_surface%n3(iele,j)) ! [km] 2021.05.25
  end do

  !# [3] ## First term from the rot rot, S

  A1(:,:)=0.d0
  B1(:,:)=0.d0

  sigma = g_surface%cond(iele) ! 2021.05.21
!  sigma = 0.01d0
  BB    = iunit*omega*dmu*sigma *L0**2.d0 ! BB is complex in ww, dn_dx twice
  
  !# [4-2] ## assemble scheme No.2 ( analytical assembly)
  ! Since \nabla lambda_k =1/6/V*( x_ln \times x_lm )  is constant,
  ! int { w_i cdot w_j }dv can be calculated analytically
  call area_tri(y,z,s) ! 2021.05.25

  do i_edge=1,3                     ! 2020.10.28
    call grcgr(i_edge,y,z,grcgr_i)  ! 2021.05.25
!      write(*,*) "i_edge",i_edge,"grcgr_i",grcgr_i
      do j_edge=1,3                   ! 2020.10.28

      call grcgr(j_edge,y,z,grcgr_j)! 2021.05.25
!       write(*,*) "j_edge",j_edge,"grcgr_j",grcgr_j

      A1(i_edge,j_edge) = 4.*s*inner(grcgr_i,grcgr_j) ! 2021.05.25
      B1(i_edge,j_edge) = wiwjdS(i_edge,j_edge,y,z)   ! 2021.05.25
    end do
  end do

if (g_surface%cond(iele) .ge. 1.d-3  .and. .false.) then
   write(*,*) "g_surface%cond(iele)",g_surface%cond(iele)
   write(*,*) "s",s
   write(*,*) "1/s",1./s
   write(*,*) "BB",BB
   write(*,*) "A1",A1(1,1:3)
   write(*,*) "A1",A1(2,1:3)
   write(*,*) "A1",A1(3,1:3)
   write(*,*) "B1",B1(1,1:3)
   write(*,*) "B1",B1(2,1:3)
   write(*,*) "B1",B1(3,1:3)
  stop
  end if
  !# [5] ## Construct elemnt matrix, elm_k
  do i=1,3
   do j=1,3
    elm_k(i,j)=(A1(i,j)+BB*B1(i,j))*idirection(i)*idirection(j) ! elm_k (complex), S(real) ! 2021.05.25
   end do
  end do

if (g_surface%cond(iele) .ge. 1.d-3  .and. .false.) then
 write(*,*) "g_surface%cond(iele)",g_surface%cond(iele)
 write(*,*) "s",s
 write(*,*) "1/s",1./s
 write(*,*) "BB",BB
 write(*,*) "elm_k",elm_k(1,1:3)
 write(*,*) "elm_k",A1(2,1:3)
 write(*,*) "elm_k",A1(3,1:3)
 write(*,*) "elm_k",B1(1,1:3)
 write(*,*) "B1",B1(2,1:3)
 write(*,*) "B1",B1(3,1:3)
stop
end if

  !# [6] ## Set global matrix from elm_k
  do i=1,3
    table_dof_elm(i)=g_surface%n3line(iele,i)*idirection(i) ! make line positive
  end do
CALL sup_iccg(elm_k,table_dof_elm,3,A%D,A%INU,A%IAU,A%AU,g_surface%nline,A%iau_tot) !see m_iccg_var_takuto.f90 2021.09.16

  end do ! element loop end

write(*,*) "### GENMAT END !! ###"

!  do i=1,l_line%nline
!   if (b_vec(i) .ne. 0.d0) write(*,*) i,"b=",b_vec(i)
!  end do
!  stop

return
99 write(*,*) "GEGEGE z(1) .ne. z(2) "
   write(*,*) "line # ",i
   write(*,*) "y(1),z(1)=",y(1),z(1)
   write(*,*) "y(2),z(2)=",y(2),z(2)
stop
end subroutine GENMAT2DTM

!###############################################################################  GENBC2D
subroutine GENBC2D(g_surface)
use  surface_type    ! src_2D/m_surface_type.f90
implicit none
type(surface), intent(inout) :: g_surface
integer(4)                 :: i,j,nline
real(8)                    :: y(2),z(2)
integer(4)                 :: idirection
integer(4)                 :: iflag_bound ! 2021.06.07

!#[1]## set
nline = g_surface%nline
if ( .not. allocated(g_surface%Avalue_bc) ) then ! 2021.10.15
 allocate( g_surface%Avalue_bc(nline))
 allocate( g_surface%line_bc(nline))
end if                                           ! 2021.10.15
g_surface%line_bc = .false.
g_surface%Avalue_bc = 0.d0   ! 2022.01.04

!#[2]## set Avalue_bc, line_bc
!# Electric current is eigher in positive x or y direction
!#
  do i=1,nline
  iflag_bound = g_surface%iflag_bound(i)
  if ( iflag_bound .eq. 1 ) then ! meaning top boundary line (see m_surface_type.f90)
     do j=1,2 ! set coordinate 2021.05.25
      y(j)=g_surface%x1x2_face(1,g_surface%line(j,i)) ! [km] 2021.05.25
      z(j)=g_surface%x1x2_face(2,g_surface%line(j,i)) ! [km] 2021.05.25
     end do
     if ( z(1) .ne. z(2) ) goto 99
    ! note the source electric field is in positive x or y direction
    g_surface%line_bc(i)   = .true.
    ! edge is rightward, y1 < y2, then, Avalue is positive
    ! edge is leftward, y1 > y2, then, Avalue is negative
    g_surface%Avalue_bc(i) = (y(2) - y(1))*1.d0  ! El = Esrc * length [mV/km * km] Dhilichlet 2021.05.21
  else if ( iflag_bound .eq. 2 .or. iflag_bound .eq. 3 .or. iflag_bound .eq. 4 ) then ! left/right/bottom bound
    g_surface%line_bc(i)   = .true. ! 2021.06.07
    g_surface%Avalue_bc(i) = 0.d0   ! 2021.06.07
  end if

  end do ! line loop end


return
99 write(*,*) "GEGEGE"
write(*,*) "z(1)",z(1)
write(*,*) "z(2)",z(2)
stop
end


!###############################################  SET_BC_ICCG
! nsr for multiple sources are added on 2017.07.11
! This program is based on set_bc_iccg.f90 in GeoFEM
! Coded on Aug. 21, 2015
SUBROUTINE SET_BC_ICCG(A,doftot,nsr,rf,dirichlet_value,dirichlet,ip,iface)
use iccg_var_takuto
implicit none
type(global_matrix),          intent(inout) :: A
integer(4),                   intent(in)    :: doftot, ip, nsr ! 2017.07.11
complex(8),dimension(doftot), intent(in)    :: dirichlet_value
logical,   dimension(doftot), intent(in)    :: dirichlet
complex(8),                   intent(out)   :: rf(doftot,nsr)  ! 2017.07.11
integer(4) :: i, j, ipos
complex(8)      ::  temp
integer(4), intent(in) :: iface ! 2022.01.04
character(1) :: num
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
write(num,'(i1)') iface
!if ( ip.eq.0) open(1,file="set_bc2D"//num//".dat")
    do i=1,doftot
        if (dirichlet(i)) then
        do j=1,nsr                 ! 2017.07.11
           rf(i,j) = dirichlet_value(i)
        !   if ( ip.eq.0) write(1,*) i,rf(i,j)
        end do
        else
          temp   =       A%D   (i   ) * dirichlet_value(        i     )
          do ipos=A%INL(i-1)+1,A%INL(i)
            temp = temp + A%AL(ipos) * dirichlet_value( A%IAL(ipos) )
          enddo
          do ipos=A%INU(i-1)+1,A%INU(i)
            temp = temp + A%AU(ipos) * dirichlet_value( A%IAU(ipos) )
          enddo
          do j=1,nsr
          !  if ( ip.eq.0) write(1,*) i,temp
            rf(i,j) = rf(i,j) - temp ! 2017.07.11
          !  if ( ip.eq.0) write(1,*) i,rf(i,j)
          end do
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

