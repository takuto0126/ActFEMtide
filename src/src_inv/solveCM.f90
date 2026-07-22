! Coded on 2017.05.12 to gen CM, copied from ./solver/solvePARDISO.f90
!INCLUDE 'mkl_pardiso.f90'
subroutine solveCM(doftot,A,x)
use mkl_pardiso
use matrix
IMPLICIT NONE
!==========================
integer(4),intent(in)  :: doftot
type(real_crs_matrix),intent(in) :: A ! A is deallocated before pardiso
real(8),intent(out) :: x(doftot,doftot) ! solution will be stored
!==========================
!INTEGER, PARAMETER :: dp = KIND(1.0D0)
!.. Internal solver memory pointer 
TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)
!.. All other variables
real(8) :: rhs(doftot,doftot) ! right hand side matrix
INTEGER maxfct, mnum, mtype, phase, n, nrhs, error, msglvl, nnz
INTEGER error1
INTEGER, ALLOCATABLE :: iparm( : )
INTEGER, ALLOCATABLE :: ia( : )
INTEGER, ALLOCATABLE :: ja( : )
real(8), ALLOCATABLE :: amat( : )
real(8), ALLOCATABLE :: b(:,:)
!real(8), ALLOCATABLE :: x( : , :)
INTEGER i, j,idum(1)
real(8) ddum(1)
!write(*,*) "solvePARDISO start!!"
!.. Fill all arrays containing matrix data.
n = doftot     ! number of equations
nnz = A%ntot   ! all the components of matrix A
nrhs = n       ! number of right hand side vector
maxfct = 1 
mnum = 1
ALLOCATE(ia(n + 1))
ALLOCATE(ja(nnz))
ALLOCATE(amat(nnz))
ALLOCATE(b(n,n))
!ALLOCATE(x(n,n))
  IA(1)=1 ; b = 0.d0
  DO I = 2,n+1
   IA(I)=IA(I-1) + A%stack(I-1)-A%stack(I-2) !the first of Ith row
   b(I-1,I-1) = 1.d0
  END DO
  JA = A%item
  AMAT = A%val
!# example
!  AMAT= (/3,1,1,2,5,1,3,4,2,1,1,2,1,3,2,1,1,1,1/)*1.d0
!  IA  = (/1,5,9,13,18,20/)
!  JA  = (/1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,5,4,5/)
write(*,*) "solvePARDISO start!!"

!..
!.. Set up PARDISO control parameter
!..
ALLOCATE(iparm(64))

DO i = 1, 64
   iparm(i) = 0
END DO

iparm(1) = 1 ! no solver default
iparm(2) = 2 ! fill-in reordering from METIS
iparm(4) = 0 ! no iterative-direct algorithm
iparm(5) = 0 ! no user fill-in reducing permutation
iparm(6) = 0 ! =0: solution is stored in x, while b is not changed
iparm(8) = 2 ! numbers of iterative refinement steps
iparm(10) = 13 ! perturb the pivot elements with 1E-13
iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
iparm(13) = 1 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
! iparm(13) should be 1, if iparm(13)=0, the accuracy of solution dramatically fall down
! on May 13, 2016
iparm(14) = 0 ! Output: number of perturbed pivots
iparm(18) = -1 ! Output: number of nonzeros in the factor LU
iparm(19) = -1 ! Output: Mflops for LU factorization
iparm(20) = 0 ! Output: Numbers of CG Iterations

error  = 0 ! initialize error flag
msglvl = 0 ! print statistical information
!msglvl = 1 ! print statistical information
mtype  = 11 ! real and structually symmetric

!.. Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.

ALLOCATE (pt(64))
DO i = 1, 64
   pt(i)%DUMMY =  0 
END DO

!.. Reordering and Symbolic Factorization, This step also allocates
! all memory that is necessary for the factorization

phase = 13 ! factorization to solve

CALL pardiso (pt, maxfct, mnum, mtype, phase, n, AMAT, ia, ja, &
              idum, nrhs, iparm, msglvl, b, x, error)
WRITE(*,*) 'Solve completed ... '
IF (error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', error
   GOTO 1000
ENDIF
WRITE(*,*) 'The solution of the system is '

!xout(1:n)=x(1:n)
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
IF (ALLOCATED(b))       DEALLOCATE(b)
!IF (ALLOCATED(x))       DEALLOCATE(x)
IF (ALLOCATED(iparm))   DEALLOCATE(iparm)

IF (error1 /= 0) THEN
   WRITE(*,*) 'The following ERROR on release stage was detected: ', error1
   STOP 1
ENDIF

IF (error /= 0) STOP 1
!END PROGRAM pardiso_sym_f90

return
end subroutine solveCM
