!# coded on 2017.06.15
module solveCM_ap
use mkl_pardiso
use matrix

type PARDISO_PARAM
 TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)
 !.. All other variables
 INTEGER maxfct, mnum, mtype, n, nrhs, error, msglvl, nnz
 INTEGER error1
 INTEGER, ALLOCATABLE :: iparm( : )
 INTEGER, ALLOCATABLE :: ia( : )
 INTEGER, ALLOCATABLE :: ja( : )
 real(8), ALLOCATABLE :: amat( : )
 INTEGER(4) :: idum(1)
 REAL(8) ddum(1)
end type

contains
!######################################################## PARDISOphase1
subroutine PARDISOphase1(A,B)
IMPLICIT NONE
!==========================
type(real_crs_matrix),intent(in)    :: A ! A is deallocated before pardiso
type(PARDISO_PARAM),  intent(inout) :: B
integer(4) :: i,j,phase
!==========================
!INTEGER, PARAMETER :: dp = KIND(1.0D0)
!.. Internal solver memory pointer 
!.. Fill all arrays containing matrix data.
B%n    = A%nrow     ! number of equations
B%nnz  = A%ntot   ! all the components of matrix A
B%nrhs = 1       ! number of right hand side vector
B%maxfct = 1
B%mnum = 1
write(*,*) "B%n=",B%n
write(*,*) "B%nnz=",B%nnz

ALLOCATE(B%ia(B%n + 1))
ALLOCATE(B%ja(B%nnz))
ALLOCATE(B%amat(B%nnz))

  B%IA(1)=1
  DO I = 2,B%n+1
   B%IA(I)=B%IA(I-1) + A%stack(I-1)-A%stack(I-2) !the first of Ith row
  END DO
  B%JA   = A%item
  B%AMAT = A%val

write(*,*) "solvePARDISO start!!"

ALLOCATE(B%iparm(64))
DO i = 1, 64
   B%iparm(i) = 0
END DO

B%iparm(1) = 1 ! no solver default
B%iparm(2) = 2 ! fill-in reordering from METIS
B%iparm(4) = 0 ! no iterative-direct algorithm
B%iparm(5) = 0 ! no user fill-in reducing permutation
B%iparm(6) = 0 ! =0: solution is stored in x, while b is not changed
B%iparm(8) = 2 ! numbers of iterative refinement steps
B%iparm(10) = 13 ! perturb the pivot elements with 1E-13
B%iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
B%iparm(13) = 1 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
! iparm(13) should be 1, if iparm(13)=0, the accuracy of solution dramatically fall down
! on May 13, 2016
B%iparm(14) = 0 ! Output: number of perturbed pivots
B%iparm(18) = -1 ! Output: number of nonzeros in the factor LU
B%iparm(19) = -1 ! Output: Mflops for LU factorization
B%iparm(20) = 0 ! Output: Numbers of CG Iterations

B%error  = 0 ! initialize error flag
B%msglvl = 0 ! print statistical information
!msglvl = 1 ! print statistical information
B%mtype  = 11 ! real and structually symmetric

!.. Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.

ALLOCATE (B%pt(64))
DO i = 1, 64
   B%pt(i)%DUMMY =  0
END DO

!.. Reordering and Symbolic Factorization, This step also allocates
! all memory that is necessary for the factorization

phase = 11 ! only reordering and symbolic factorization

CALL pardiso (B%pt, B%maxfct, B%mnum, B%mtype, phase, B%n, B%AMAT, B%ia, B%ja, &
              B%idum, B%nrhs, B%iparm, B%msglvl, B%ddum, B%ddum, B%error)

WRITE(*,*) 'Reordering completed ... '
IF (B%error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', B%error
!   GOTO 1000
call PARDISOphase4(B)
END IF
WRITE(*,*) 'Number of nonzeros in factors = ',B%iparm(18)
WRITE(*,*) 'Number of factorization MFLOPS = ',B%iparm(19)

return
end subroutine PARDISOphase1

!########################################################### phase2
subroutine PARDISOphase2(B)
type(PARDISO_PARAM),intent(inout) :: B
integer(4) :: phase
write(*,*) "### PARDISOphase2 start!! ###"
!.. Factorization.
phase = 22 ! only factorization
CALL pardiso (B%pt, B%maxfct, B%mnum, B%mtype, phase, B%n, B%AMAT, B%ia, B%ja, &
              B%idum, B%nrhs, B%iparm, B%msglvl, B%ddum, B%ddum, B%error)
WRITE(*,*) 'Factorization completed ... '
IF (B%error /= 0) THEN
   WRITE(*,*) 'The following ERROR was detected: ', B%error
!   GOTO 1000
call PARDISOphase4(B)
ENDIF
write(*,*) "### PARDISOphase2 END!! ###"
return
end subroutine PARDISOphase2

!########################################################### phase3
subroutine PARDISOphase3(B,nmodel,CM)
type(PARDISO_PARAM),  intent(inout) :: B
integer(4),           intent(in)    :: nmodel
type(real_crs_matrix),intent(out)   :: CM
integer(4),           parameter     :: nrhsmax = 500
integer(4)                          :: phase,icount,nite,nres
real(8),allocatable,dimension(:,:)  :: bb,x
real(8),parameter                   :: threshold = 1.d-15
type(real_ccs_matrix),allocatable,dimension(:) :: xccs
type(real_ccs_matrix)                          :: CMCCS
integer(4) :: i,j

B%iparm(8) = 2 ! max numbers of iterative refinement steps
phase      = 33 ! only solving
B%msglvl   = 0 ! not print statistical info

nite = nmodel/nrhsmax
if ( nmodel .gt. nrhsmax*nite ) nite = nite + 1
write(*,*) "nite=",nite
allocate(bb(nmodel,nrhsmax),x(nmodel,nrhsmax))

icount=0
allocate(xccs(nite))

do i=1,nmodel/nrhsmax
  bb(:,:)=0.d0 ; x(:,:) = 0.d0
  do j=1,nrhsmax
   icount = icount + 1
   bb(icount,j)=1.d0
  end do
 write(*,*) "i=",i,"icount=",icount,"nmodel=",nmodel
 CALL pardiso (B%pt, B%maxfct, B%mnum, B%mtype, phase, B%n, B%AMAT, B%IA, B%JA, &
          &   B%idum, nrhsmax, B%iparm, B%msglvl, bb, x, B%error)
 write(*,*) 'i=',i,'Solve completed ... '
! if (icount .gt. 3000 ) write(*,*) "x(:,1)=",x(:,1)
 IF (B%error /= 0) goto 99
 call conv_full2ccs(x,nmodel,nrhsmax,xccs(i),threshold)
 write(*,*) "i=",i,"xccs(i) end"
end do

if ( nmodel .gt. int(nmodel/nrhsmax)*nrhsmax ) then
 nres = nmodel - int(nmodel/nrhsmax)*nrhsmax
 write(*,*) "nres=",nres
 bb(:,:)=0.d0; x(:,:) = 0.d0
 do i=1,nres
  icount = icount + 1
  bb(icount,i)=1.d0
 end do
 write(*,*) "icount=",icount,"nmodel=",nmodel
 CALL pardiso (B%pt, B%maxfct, B%mnum, B%mtype, phase, B%n, B%AMAT, B%ia, B%ja, &
              B%idum, nres, B%iparm, B%msglvl, bb(:,1:nres), x(:,1:nres), B%error)
 call conv_full2ccs(x(:,1:nres),nmodel,nres,xccs(nite),threshold)
end if
write(*,*) "## xccs are generated! ##"

!# assemble CM
CMCCS%nrow  = nmodel
CMCCS%ncolm = nmodel
CMCCS%ntot  = 0
do i=1,nite
 write(*,*) "i=",i,"xccs%ntot",xccs(i)%ntot
 CMCCS%ntot = CMCCS%ntot + xccs(i)%ntot
end do
write(*,*) "CMCCS%ntot =",CMCCS%ntot
allocate(CMCCS%stack(0:nmodel))
allocate(CMCCS%item(CMCCS%ntot))
allocate(CMCCS%val(CMCCS%ntot))
CMCCS%stack(0) = 0

do i=1,nite
 write(*,*) "i=",i
 call combine_real_ccs_mat(CMCCS,xccs(i),(i-1)*nrhsmax+1)
 write(*,*) "(i-1)*nrhsmax+1",(i-1)*nrhsmax+1,"CMCCS%stack(i)=",CMCCS%stack((i-1)*nrhsmax+1)
end do
write(*,*) "## CMCCS is generated! ##"
call conv_ccs2crs(CMCCS,CM)
write(*,*) "## CM is generated! ##"

write(*,*) "### PARDISOphase3 END!! ###"

return

99 continue
   WRITE(*,*) 'The following ERROR was detected: ', B%error
   CALL PARDISOphase4(B)
   stop

end subroutine PARDISOphase3

!########################################################### phase4
subroutine PARDISOphase4(B)
type(PARDISO_PARAM),intent(inout) :: B
integer(4) :: phase,error1
!.. Termination and release of memory
phase = -1 ! release internal memory
CALL pardiso (B%pt, B%maxfct, B%mnum, B%mtype, phase, B%n, B%ddum, B%idum, B%idum, &
              B%idum, B%nrhs, B%iparm, B%msglvl, B%ddum, B%ddum, error1)

IF (ALLOCATED(B%iparm))  DEALLOCATE(B%iparm)
IF (ALLOCATED(B%IA))     DEALLOCATE(B%IA)
IF (ALLOCATED(B%JA))     DEALLOCATE(B%JA)
IF (ALLOCATED(B%AMAT))   DEALLOCATE(B%AMAT)

IF (error1 /= 0) THEN
   WRITE(*,*) 'The following ERROR on release stage was detected: ', error1
   STOP 1
ENDIF

IF (B%error /= 0) STOP 1
!END PROGRAM pardiso_sym_f90

return
end subroutine PARDISOphase4

end module
