! modified for multiple source on 2017.07.11
! Coded by T.M. on May 13, 2016
! based on pardiso_sym_f90.f90
!PROGRAM pardiso_sym_f90
INCLUDE 'mkl_pardiso.f90'
subroutine solvePARDISO(doftot,nsr,A,rhs,xout,ip) !2017.07.11
USE mkl_pardiso
use iccg_var_takuto
IMPLICIT NONE
!==========================
type(global_matrix),intent(inout) :: A ! A is deallocated before pardiso
integer(4),         intent(in)    :: doftot
integer(4),         intent(in)    :: nsr  ! # of rhs vector, 2017.07.11
complex(8),         intent(inout) :: rhs(doftot,nsr)  ! 2017.07.11
complex(8),         intent(out)   :: xout(doftot,nsr) ! 2017.07.11
integer(4),         intent(in)    :: ip
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
COMPLEX(8), ALLOCATABLE :: b( : , :) ! 2017.07.11
COMPLEX(8), ALLOCATABLE :: x( : , :) ! 2017.07.11
INTEGER i, j,idum(1)
COMPLEX(8) ddum(1)
!write(*,*) "solvePARDISO start!!"

!.. Fill all arrays containing matrix data.
n      = doftot               ! number of equations
nnz    = A%iau_tot + doftot  ! upper triangle + diagonal
nrhs   = nsr                 ! number of right hand side vector, 2017.07.11
maxfct = 1 
mnum   = 1
!################################################# for debug
!n = 8
!nnz = 18
!nrhs = 1
!maxfct = 1
!mnum = 1
!deallocate(A%AU,A%IAU,A%INU,A%D)
!allocate(A%AU(nnz-n),A%IAU(nnz),A%INU(0:n),A%D(n))
!A%D=(/7.d0,-4.d0,1.d0,7.d0,5.d0,-1.d0,11.d0,5.d0/)*(0.d0,1.d0)
!A%AU=(/1.d0,2.d0,7.d0,8.d0,2.d0,  5.d0,9.d0,1.d0,5.d0,5.d0/)*(0.d0,1.d0)
!A%INU=(/0,3,5,6,7,9,10,10,10/)
!A%IAU=(/3,6,7,3,5,8,7,6,7,8/)
!rhs(1:n)=(0.d0,0.d0)
!#################################################
ALLOCATE(ia(n + 1))
ALLOCATE(ja(nnz))
ALLOCATE(amat(nnz))
ALLOCATE(b(n,nrhs)) ! 2017.07.11
ALLOCATE(x(n,nrhs)) ! 2017.07.11
  IA(1)=1
  DO I = 2,n+1
   IA(I)=IA(I-1)+(A%INU(I-1)-A%INU(I-2))+1 !the first of Ith row
   JA(IA(I-1))=I-1 ! diagonal
   AMAT(IA(I-1))=A%D(I-1)
   b(I-1,:)=rhs(I-1,:) ! 2017.07.11
   DO J =1,(A%INU(I-1) - A%INU(I-2)) ! upper triangle of I-1 th row
   JA  (IA(I-1)+J) = A%IAU(A%INU(I-2)+J) !upper triangle
   AMAT(IA(I-1)+J) =  A%AU(A%INU(I-2)+J)
   END DO
  END DO
!do i=1,8
!write(*,'(a,3i7,2e15.7)') "I,INU,IAU,A",i,A%INU(i),A%IAU(i),A%AU(i)
!end do
!deallocate (A%IAU, A%IAL, A%INU, A%INL, A%AU, A%AL, A%D)
write(*,*) "solvePARDISO start!!"

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
iparm(14) = 0  ! Output: number of perturbed pivots
iparm(18) = -1 ! Output: number of nonzeros in the factor LU
iparm(19) = -1 ! Output: Mflops for LU factorization
iparm(20) = 0  ! Output: Numbers of CG Iterations

error  = 0  ! initialize error flag
msglvl = 0  ! print statistical information
!msglvl = 1 ! print statistical information
mtype  = 6  ! complex and symmetric

!.. Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.

ALLOCATE (pt(64))
DO i = 1, 64
   pt(i)%DUMMY =  0 
END DO

!.. Reordering and Symbolic Factorization, This step also allocates
! all memory that is necessary for the factorization

phase = 11 ! only reordering and symbolic factorization

CALL pardiso (pt, maxfct, mnum, mtype, phase, n, AMAT, ia, ja, &
              idum, nrhs, iparm, msglvl, ddum, ddum, error)
    
!WRITE(*,*) 'Reordering completed ... '
IF (error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', error
   GOTO 1000
END IF
!WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
!WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)

!.. Factorization.
phase = 22 ! only factorization
CALL pardiso (pt, maxfct, mnum, mtype, phase, n, AMAT, ia, ja, &
              idum, nrhs, iparm, msglvl, ddum, ddum, error)
!WRITE(*,*) 'Factorization completed ... '
IF (error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', error
   GOTO 1000
ENDIF

!.. Back substitution and iterative refinement
iparm(8) = 2 ! max numbers of iterative refinement steps
phase = 33 ! only solving
!############################## for debug
!DO i = 1, n
!   b(i) = 1.d0
!END DO
!##############################
CALL pardiso (pt, maxfct, mnum, mtype, phase, n, AMAT, ia, ja, &
              idum, nrhs, iparm, msglvl, b, x, error)
!WRITE(*,*) 'Solve completed ... '
IF (error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', error
   GOTO 1000
ENDIF
!WRITE(*,*) 'The solution of the system is '

xout(1:n,1:nrhs)=x(1:n,1:nrhs)
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
IF (ALLOCATED(x))       DEALLOCATE(x)
IF (ALLOCATED(iparm))   DEALLOCATE(iparm)

IF (error1 /= 0) THEN
   WRITE(*,*) 'The following ERROR on release stage was detected: ', error1
   STOP 1
ENDIF

IF (error /= 0) STOP 1
!END PROGRAM pardiso_sym_f90

return
end subroutine solvePARDISO
