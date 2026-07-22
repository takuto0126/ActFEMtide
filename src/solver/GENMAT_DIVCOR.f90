!######################################## GENMAT_DIVCOR
!# Coded on Nov. 24, 2015
!#
subroutine GENMAT_DIVCOR(ntetg, n6line, xyz, n4, nodeg, nline)
use  outerinnerproduct
use  iccg_var_takuto ! b_vec is included
implicit none
integer(4),intent(in) :: ntetg, n6line(ntetg,6), n4(ntetg,5), nodeg, nline
real(8),   intent(in) :: xyz(3, nodeg)
real(8) :: elm_xyz(3,4)! nodal coordinates of element elm_xyz(i,j) i:x,y,z,j-th node
real(8) :: xx(3,6), v,  gn(3,4), S(4,4),S1(4,6),din,v1(3),v2(3),v3
complex(8) :: elm_k(4,4)
integer(4) :: table_dof_elm6(6), table_dof_elm4(4)
integer(4) :: iele, i, j, k, l, m, n, ii, jj, idirection(6), lm(6,2)
!---------------  scales ------------------------------------------------------
real(8), parameter  :: L0=1.d+3  ! [m]  scale length

!# [0] ## definition of line composition regarding a tetrahedron
  lm(1,1:2)=(/1,2/)
  lm(2,1:2)=(/2,3/)
  lm(3,1:2)=(/1,3/)
  lm(4,1:2)=(/1,4/)
  lm(5,1:2)=(/2,4/)
  lm(6,1:2)=(/3,4/)
  !write(*,*) "check0"

do iele=1, ntetg  ! start elemetn loop
!#
  !# [1] ## Prepare idirection
  !         check the direction of edge, compared to the defined lines
  idirection(1:6)=1
  do j=1,6
    if ( n6line(iele, j) .lt. 0 ) idirection(j)=-1
  end do
  !write(*,*) "check1"

  !# [2] ## Prepare the coordinates for 4 nodes of elements
  elm_xyz(1:3,1)=xyz(1:3,n4(iele,2)) ! [km]
  elm_xyz(1:3,2)=xyz(1:3,n4(iele,3))
  elm_xyz(1:3,3)=xyz(1:3,n4(iele,4))
  elm_xyz(1:3,4)=xyz(1:3,n4(iele,5))

  !# [3] ## Prepare volume of this tetrahedron
  xx(:,1)=elm_xyz(:, 4) - elm_xyz(:, 3) ! rot w
  xx(:,2)=elm_xyz(:, 4) - elm_xyz(:, 1)
  xx(:,3)=elm_xyz(:, 2) - elm_xyz(:, 4)
  xx(:,4)=elm_xyz(:, 3) - elm_xyz(:, 2)
  xx(:,5)=elm_xyz(:, 1) - elm_xyz(:, 3)
  xx(:,6)=elm_xyz(:, 2) - elm_xyz(:, 1)
  call volume(xx(:,1), xx(:,5), -xx(:,4), V) ! volume of this element
  if ( v .le. 0 ) goto 999

  !# [4] ## Gradient of node-base functions, gn(i,j) is i-th component of nablda N_j
  ! gradient of nodal shape function, where outer vector, gn(:,i), points to i-th node
  gn(:,1)=1.d0/6.d0/v*outer(elm_xyz(:,4)-elm_xyz(:,2), elm_xyz(:,3)-elm_xyz(:,2))
  gn(:,2)=1.d0/6.d0/v*outer(elm_xyz(:,3)-elm_xyz(:,1), elm_xyz(:,4)-elm_xyz(:,1))
  gn(:,3)=1.d0/6.d0/v*outer(elm_xyz(:,4)-elm_xyz(:,1), elm_xyz(:,2)-elm_xyz(:,1))
  gn(:,4)=1.d0/6.d0/v*outer(elm_xyz(:,2)-elm_xyz(:,1), elm_xyz(:,3)-elm_xyz(:,1))
  !write(*,*) "check4"

  !# [5] ## Left hand side matrix of Laprasian Phi, S(4,4) (real)
  ! int{ nablda(N) cdot nabla(N)^T }dv = nablda(N) cdot nabla(N)^T * v
  call mul_atb(gn(1:3,1:4), gn(1:3,1:4), S(1:4,1:4), 4, 3, 4)
  elm_k(:,:)=S(:,:)*v ! S (real) -> elm_k (complex)

  !# [6] ## Right hand side matrix of -int{ nabla {N} cdot {w}^T } dv {bsl}
  S1(:,:)=0.d0
  ! S1(i,j) = - int { nabla N cdot {w}^T }dv, where w is vector shape function
  !         = - int { gn(:,i) cdot [n_k*gn(:,l) - n_l*gn(:,k) ] } dv
  !         = -{  int ( n_k ) dv [ gn(:,i) cdot gn (:,l) ]      first
  !           - int ( n_l ) dv [ gn(:,i) cdot gn (:,k) ] }      second
  !         = -v/4 * gn(:,i) cdot { gn(:,l) - gn(:,k)}
  do i=1,4    ! node loop
    do j=1,6  ! line loop
      k=lm(j,1) ; l=lm(j,2)
      S1(i,j)= - v/4.d0*inner(gn(1:3,i),gn(1:3,l)-gn(1:3,k))*idirection(j)
    end do
  end do
  !CALL checkvalues(elm_xyz,xx,gn,elm_k,S,4,S1,4,6,v)

  !# [7] ## Set global matrix (LHS)
  do i=1,4
     table_dof_elm4(i)=n4(iele,i+1)
  end do
  CALL sup_iccg(elm_k,table_dof_elm4,4,D,INU,IAU,AU,nodeg,iau_tot)

  !# [8] ## Set global coefficient matrix nodeg*nline (RHS)
  do i=1,6
     table_dof_elm6(i)=n6line(iele,i)*idirection(i) ! make n6line positive
  end do
  CALL sup_iccg_asym(S1, table_dof_elm4, 4, table_dof_elm6, 6,&
 & istack, item, biccg, nodeg, nline, item_tot )

end do ! element loop end
write(*,*) "### GENMAT_DIVCOR END !! ###"
return
999 continue
write(*,*) "GEGEGE! volume is less than 0; iele=",iele
write(*,*)  "xx(:,1)=",xx(:,1)
write(*,*)  "xx(:,5)=", xx(:,5)
write(*,*)  "-xx(:,4)=",-xx(:,4)
stop
end subroutine GENMAT_DIVCOR
!####################################################  sup_iccg_asym
! Coppied from sup_iccg in n_ebfem.f90 on Nov. 24, 2015
! and modified for assymetric matrix
subroutine  sup_iccg_asym( elm_k, table_dof_elm1, edof1, table_dof_elm2, edof2,&
 & istack, item, biccg, doftot1, doftot2, item_tot )
!    superpose asymmetric element matrix (doftot1*doftot2 )to gloval stiffness
!    for iccg solver (lower element, upper whole stiffness)
implicit  none
integer(4),intent(in) :: edof1,edof2, doftot1, doftot2, item_tot
real(8),   intent(in)   ::  elm_k(edof1,edof2) !  element stiffness matrix
integer(4), intent(in)   :: table_dof_elm1(edof1),table_dof_elm2(edof2)
!    digree of freedom number table at a element table_dof_elm(i)
!      where i    : dof order of 'elm_k'
real(8),intent(inout) :: biccg(item_tot)   ! value   (i-th dof)
integer(4),intent(in) :: istack(0:doftot1) ! store # of terms at i-th row
integer(4),intent(in) :: item(item_tot)    ! storing colum
integer(4)  ::  i_gl, j_gl, i_elm, j_elm
integer(4)  ::  num , i_pos
logical     ::  found
! * superpose element stiffness matrix
!
!    k(i_gl,j_gl) = k(i_gl,j_gl) + elm_k(i_elm,j_elm)
!
!   . make diagonal and upper triangle
      do i_elm=1,edof1
         i_gl = table_dof_elm1(i_elm)
         do j_elm=1,edof2
            j_gl = table_dof_elm2(j_elm)
		  found = .false.
              do num=istack(i_gl-1)+1,istack(i_gl) ! .. search term pos.
                if ( j_gl == item(num) ) then
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
		  biccg(i_pos) = biccg(i_pos) + elm_k(i_elm,j_elm)
         enddo
      enddo
      return
      end subroutine  sup_iccg_asym
