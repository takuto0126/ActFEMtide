! modified for multiple source on 2017.07.13
! Coded by T.M. on 2017.05.15 originally from ../solver/solvePARDISO
! THis program solve normal forward problem as well as Au=P for Jacobian Matrix
INCLUDE 'mkl_pardiso.f90'
subroutine solvePARDISOinv(doftot,nsr,A,rhs,xout,PTR,ut,iflag_comp,ip,np) ! 2020.09.18
use mkl_pardiso
use iccg_var_takuto
use matrix         ! 2017.05.15
implicit none
!-- for forward ------
type(global_matrix),  intent(inout)  :: A           ! A is deallocated bfr pardiso
integer(4),           intent(in)     :: nsr         ! # of rhs vector, 2017.07.13
integer(4),           intent(in)     :: doftot
complex(8),           intent(in)     :: rhs(doftot,nsr) ! 2017.07.13
complex(8),           intent(out)    :: xout(doftot,nsr)! 2017.07.13
integer(4), optional, intent(in)     :: ip,np        ! 2020.09.18
!-- for inversion ---- added on 2017.05.15
type(real_crs_matrix),   intent(in)  :: PTR(5)       ! [nobs*ncomp,nline] ! 2018.10.04
integer(4),dimension(5), intent(in)  :: iflag_comp   ! 2018.10.04
type(complex_crs_matrix),intent(out) :: ut(5)        ! solution for Au = P,2018.10.04
type(complex_crs_matrix)             :: u,crscompout ! 2017.05.15
type(complex_ccs_matrix)             :: utccs        ! 2017.05.15
type(real_ccs_matrix)                :: ccsout
!==========================
!INTEGER, PARAMETER :: dp = KIND(1.0D0)
!.. Internal solver memory pointer 
TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)
!.. All other variables
INTEGER maxfct, mnum, mtype, phase, n, nrhs, error, msglvl, nnz
INTEGER error1
INTEGER,    ALLOCATABLE :: iparm( : )
INTEGER,    ALLOCATABLE :: ia( : )
INTEGER,    ALLOCATABLE :: ja( : )
COMPLEX(8), ALLOCATABLE :: amat( : )
COMPLEX(8), ALLOCATABLE :: b(:,:)  ! (ndof,nsr)  2018.10.04
COMPLEX(8), ALLOCATABLE :: b1(:,:) ! (ndof,nobs) 2018.10.04
COMPLEX(8), ALLOCATABLE :: x(:,:)
COMPLEX(8), ALLOCATABLE :: x1(:,:) ! 2018.10.05
INTEGER i, j,idum(1),j1,j2, icomp  ! 2018.10.04
COMPLEX(8) ddum(1)
integer(4)              :: nobs,ncolm ! 2017.05.15
type(real_ccs_matrix)   :: P
real(8) :: threshold = 1.d-10


!#[0]## set nobs and others
 do icomp=1,5                           ! 2018.10.04
  if (iflag_comp(icomp) .eq. 0 ) cycle  ! 2018.10.04
  nobs = PTR(icomp)%nrow
  if ( doftot .ne. PTR(icomp)%ncolm ) THEN
   write(*,*) "GEGEGE doftot=",doftot," .ne. PTR%nrow=",PTR(icomp)%nrow
   stop
  end if
 end do                      ! 2018.10.04
 if ( nobs .eq. 0 ) goto 998 ! 2018.10.04
 n      = doftot             ! number of equations
 nnz    = A%iau_tot + doftot ! upper triangle + diagonal
 nrhs   = nsr                ! number of right hand side vector, 2018.10.04
 maxfct = 1
 mnum   = 1

!---------------------------------------------------------
ALLOCATE( ia(n + 1),ja(nnz),amat(nnz))
ALLOCATE( b(n,nrhs), x(n,nrhs))

!### coefficient matrix
  IA(1)=1
  DO I = 2,n+1
   IA(I)=IA(I-1)+(A%INU(I-1)-A%INU(I-2))+1 !the first of Ith row
   JA(IA(I-1))=I-1 ! diagonal
   AMAT(IA(I-1))=A%D(I-1)
   DO J =1,(A%INU(I-1) - A%INU(I-2)) ! upper triangle of I-1 th row
   JA  (IA(I-1)+J) = A%IAU(A%INU(I-2)+J) !upper triangle
   AMAT(IA(I-1)+J) =  A%AU(A%INU(I-2)+J)
   END DO
  END DO

  b(:,:)     = 0.d0
  b(:,1:nsr) = rhs(:,1:nsr)      ! 2017.07.13
!

!write(*,*) "solvePARDISOinv start!!" ! commented out on 2017.12.22

!..
!.. Set up PARDISO control parameter
!..
ALLOCATE(iparm(64))

DO i = 1, 64
   iparm(i) = 0
END DO

iparm(1)  = 1 ! no solver default
iparm(2)  = 2 ! fill-in reordering from METIS
iparm(4)  = 0 ! no iterative-direct algorithm
iparm(5)  = 0 ! no user fill-in reducing permutation
iparm(6)  = 0 ! =0: solution is stored in x, while b is not changed
iparm(8)  = 2 ! numbers of iterative refinement steps
iparm(10) = 13 ! perturb the pivot elements with 1E-13
iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
iparm(13) = 1 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
! iparm(13) should be 1, if iparm(13)=0, the accuracy of solution dramatically fall down
! on May 13, 2016
iparm(14) = 0 ! Output: number of perturbed pivots
iparm(18) = -1 ! Output: number of nonzeros in the factor LU
iparm(19) = -1 ! Output: Mflops for LU factorization
iparm(20) = 0 ! Output: Numbers of CG Iterations

error     = 0 ! initialize error flag
msglvl    = 0 ! print statistical information
!msglvl = 1 ! print statistical information
mtype     = 6 ! complex and symmetric

!.. Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.

ALLOCATE (pt(64))
DO i = 1, 64
   pt(i)%DUMMY =  0 
END DO


!#[1]## only reordering and symbolic factorization
  !.. Reordering and Symbolic Factorization, This step also allocates
  ! all memory that is necessary for the factorization

  phase = 11 !

  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, AMAT, ia, ja, &
              idum, nrhs, iparm, msglvl, ddum, ddum, error)
    
  !WRITE(*,*) 'Reordering completed ... ' !commented out on 2017.12.22
  IF (error /= 0) THEN
    WRITE(*,*) 'The following ERROR was detected: ', error
    GOTO 1000
  END IF
  !WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
  !WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)

!#[2]## Factorization

  phase = 22 ! only factorization
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, AMAT, ia, ja, &
              idum, nrhs, iparm, msglvl, ddum, ddum, error)
  !WRITE(*,*) 'Factorization completed ... ' ! commented out on 2017.12.22
  IF (error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', error
   GOTO 1000
  ENDIF

!#[3]## solving: Back substitution and iterative refinement

  iparm(8) = 2 ! max numbers of iterative refinement steps
  phase = 33 ! only solving

  !#[3-1]## forward solving
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, AMAT, ia, ja, &
              idum, nrhs, iparm, msglvl, b, x, error)
  if ( present(ip) ) then ! 2020.09.18
  WRITE(*,'(a,i2)') ' Solve completed ... [Forward]  ip =',ip !," /",np ! 2020.09.18
  else   ! 2020.09.18
  WRITE(*,*) 'Solve completed ... [Forward]' ! 2018.10.04
  end if ! 2020.09.19
  xout(1:n,1:nsr)=x(1:n,1:nsr)  ! 2017.07.13

  IF (error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', error
   GOTO 1000
  ENDIF

  !#[3-2]## deallocate and allocate
  DEALLOCATE(b,         x          )   ! 2018.10.05
  ALLOCATE(  b1(n,nobs),x1(n,nobs) )   ! 2018.10.04

  !#[3-2]## 1 to 5 component solving for jacobian calculation
  do icomp = 1,5
  if ( iflag_comp(icomp) .eq. 0 ) cycle

  !#[3-2-1]## right hand side vector
  call trans_crs2ccs(PTR(icomp),P) ! see m_matrix.f90
  if (.false.) call CHECKP(P) ! see below 2018.10.05
  b1 = (0.d0,0.d0)
  do i=1,nobs                     ! row loop
   do j=P%stack(i-1)+1,P%stack(i) ! colum loop
    b1(P%item(j),i) = P%val(j)    ! 2017.07.18
   end do
  end do

 !#[3-2-2]## solve for Jacobian
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, AMAT, ia, ja, &
              idum, nobs, iparm, msglvl, b1, x1, error)

 !#[3-2-3]## u -> ut
 call convcomp_full2crs(x1(1:n,1:nobs),n,nobs,u,threshold)  ! m_matrix.f90,2017.07.13
 call transcomp_crs2ccs(u,utccs)                           ! m_matrix.f90
 call convcomp_ccs2crs(utccs,ut(icomp))                    ! m_matrix.f90

 end do
  if ( present(ip) ) then ! 2020.09.18
  WRITE(*,'(a,i2)') ' Solve completed ... [Jacobian] ip =',ip !" /",np ! 2020.09.18
  else   ! 2020.09.18
  WRITE(*,*) 'Solve completed ... [Jacobian]'      ! 2018.10.05
  end if ! 2020.09.18

!crscompout = ut
!write(*,*) "crs%nrow=",crscompout%nrow
!write(*,*) "crs%ncolm=",crscompout%ncolm
!write(*,*) "crs%ntot=",crscompout%ntot
!write(*,*) "crs%stack=",crscompout%stack
!write(*,*) "crs%item=",crscompout%item
!write(*,*) "crs%val=",crscompout%val
!DO i = 1, n
!   WRITE(*,*) ' x(',i,') = ', x(i)
!END DO


1000 CONTINUE
!.. Termination and release of memory
phase = -1 ! release internal memory
CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
              idum, nrhs, iparm, msglvl, ddum, ddum, error1)

IF (ALLOCATED(ia))      DEALLOCATE(ia)
IF (ALLOCATED(ja))      DEALLOCATE(ja)
IF (ALLOCATED(AMAT))    DEALLOCATE(AMAT)
!IF (ALLOCATED(b))       DEALLOCATE(b)
!IF (ALLOCATED(x))       DEALLOCATE(x)
IF (ALLOCATED(iparm))   DEALLOCATE(iparm)

IF (error1 /= 0) THEN
   WRITE(*,*) 'The following ERROR on release stage was detected: ', error1
   STOP 1
ENDIF

IF (error /= 0) STOP 1
!END PROGRAM pardiso_sym_f90



return
998 CONTINUE                 ! 2018.10.05
  write(*,*) "GEGEGE nobs=0" ! 2018.10.05
  stop                       ! 2018.10.05
end subroutine solvePARDISOinv
!######################################################## CHECKP
!# 2018.10.05
subroutine CHECKP(P)
use matrix
implicit none
type(real_crs_matrix),intent(in) :: P
type(real_crs_matrix) :: ccsout
integer(4) :: i,j,j1,j2
   ccsout = P
   write(*,*) "ccs%nrow=",ccsout%nrow
   write(*,*) "ccs%ncolm=",ccsout%ncolm
   write(*,*) "ccs%ntot=",ccsout%ntot
   do i=1,ccsout%ncolm
     ! if (ccsout%stack(i)-crsout%stack(i-1) .ne. 0) then
     write(*,*) i,"# of content",ccsout%stack(i)-ccsout%stack(i-1)
     j1=ccsout%stack(i-1)+1;j2=ccsout%stack(i)
     write(*,'(5g15.7)') (ccsout%item(j),j=j1,j2)
     write(*,'(5g15.7)') (ccsout%val(j),j=j1,j2)
     write(*,*) ""
     ! end if
   end do

return
end

