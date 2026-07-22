!# coded on 2017.12.12
module spectral
use matrix
use outerinnerproduct
use caltime ! 2017.12.22

contains
!################################################# alpha from spectral radius
!# 2017.12.25
!#
!#  alpha = rho(J^T[Cd^-1]J)/rho(BM)
!#
subroutine alphaspectralradius_v1(Cd,JJ,BM,alpha,kmax) ! 2017.12.22
implicit none
type(real_crs_matrix),intent(in)  :: Cd,JJ,BM
integer(4),           intent(in)  :: kmax     ! max lanczos procedure
real(8),              intent(out) :: alpha
real(8)                           :: denom,numer
type(real_crs_matrix)             :: CDI
integer(4)                        :: i
type(real_crs_matrix)             :: JTCDI,JTCDIJ
real(8)                           :: Emin,Emax,errmin,errmax
type(watch) :: t_watch,t_watch1 ! 2017.12.25

 call watchstart(t_watch) ! 2017.12.22

 !#[1]## numerator
 call watchstart(t_watch1) ! 2017.12.25
 CDI = CD
 do i=1,CD%nrow
  CDI%val(i) = 1.d0/CD%val(i)
 end do
 call mulreal_crs_crs_crs(JJ,CDI,JTCDI,"T") ! [JJ^T]CDI 2017.12.25
 call mulreal_crs_crs_crs(JTCDI,JJ,JTCDIJ)  ! 2017.12.25
! write(*,*) "Maximum Lanczos procedure, kmax=",kmax
 call LANCZOSMINMAX(JTCDIJ,Emin,Emax,errmin,errmax,kmax)
 numer = max(abs(Emax),abs(Emin)) ! 2017.12.19
 call watchstop(t_watch1)
 write(*,'(a,g15.7,a)') "### Lanzos for JTCDIJ END!! Time=",t_watch1%time,"[min] ###"!2017.12.25

 !#[2]## denominator
 call watchstart(t_watch1) ! 2017.12.25
 call LANCZOSMINMAX(BM,Emin,Emax,errmin,errmax,kmax)
 denom = max(abs(Emax),abs(Emin)) ! 2017.12.19
 call watchstop(t_watch1)
 write(*,'(a,g15.7,a)') "### Lanzos for BM     END!! Time=",t_watch1%time,"[min] ###"!2017.12.25

 !#[3]## alpha
 alpha = abs(numer/denom)
 write(*,*) "alpha=",alpha

call watchstop(t_watch) ! 2017.12.22
write(*,'(a,g15.7,a)')  "### Alphaspectralradius_v1 END!! ",t_watch%time,"[min] ###"!2017.12.22

return
end

!################################################# alpha from spectral radius
!# 2017.12.13
!#
!#  alpha = rho(JBMIJ^T)/rho(Cd)
!#
subroutine alphaspectralradius_v2(Cd,JJ,BMI,alpha,kmax) ! 2017.12.22
implicit none
type(real_crs_matrix),intent(in)  :: Cd,JJ,BMI
integer(4),           intent(in)  :: kmax     ! max lanczos procedure
real(8),              intent(out) :: alpha
real(8)                           :: denom,numer
integer(4)                        :: i
type(real_crs_matrix)             :: JCMJT
type(real_ccs_matrix)             :: CMJT,JT
real(8)                           :: Emin,Emax,errmin,errmax
type(watch) :: t_watch ! 2017.12.22

 call watchstart(t_watch) ! 2017.12.22
 !#[1]## numerator
 call trans_crs2ccs(JJ,JT)
 call mulreal_crs_ccs_ccs(BMI,JT,CMJT)
 call mulreal_crs_ccs_crs(JJ,CMJT,JCMJT)
 write(*,*) "Maximum Lanczos procedure, kmax=",kmax
 call LANCZOSMINMAX(JCMJT,Emin,Emax,errmin,errmax,kmax)
 numer = max(abs(Emax),abs(Emin)) ! 2017.12.19

 !#[2]## denominator
 denom = 0.d0
 do i=1,Cd%nrow
  denom = max(denom,Cd%val(i))
 end do

 !#[3]## alpha
 alpha = abs(numer/denom)
 write(*,*) "alpha=",alpha

call watchstop(t_watch) ! 2017.12.22
write(*,'(a,g15.7,a)') "### Alphaspectralradius_v2 END!!",t_watch%time,"[min] ###"!2017.12.22

return
end

!#####################################################  Lanczos algorighm
!# 2017.12.11
subroutine LANCZOSMINMAX(crsmat,Emin,Emax,errmin,errmax,kmax)
implicit none
type(real_crs_matrix),intent(in)   :: crsmat
integer(4),           intent(in)   :: kmax       ! max Lanczos procedure
real(8),              intent(out)  :: Emin,Emax
real(8),              intent(out)  :: errmin,errmax
integer(4)                         :: NA         ! size of NA
real(8),allocatable,dimension(:)   :: x          ! initial vector
real(8),allocatable,dimension(:,:) :: CQ         ! Capital Q=[q1,q2,..qkmax]
real(8),allocatable,dimension(:)   :: res        ! residual vector
real(8),allocatable,dimension(:)   :: Q,V        ! j-th orthnomal vector
real(8),allocatable,dimension(:)   :: alpha,beta ! for Lanczos
integer(4)                         :: i,j
!# for TQLI
real(8),allocatable,dimension(:,:) :: Z
real(8),allocatable,dimension(:)   :: D,E
integer(4)                         :: N,NP

!#[0]# set
!  kmax = crsmat%nrow
  NA   = crsmat%nrow
!  write(*,*) "crsmat%nrow",crsmat%nrow   ! commented out on 2017.12.22
!  write(*,*) "crsmat%ncolm",crsmat%ncolm ! commented out on 2017.12.22
!  write(*,*) "crsmat%ntot",crsmat%ntot   ! commented out on 2017.12.22
!  write(*,*) "NA",NA

!#[1]## Normal LANCZOS: Generate tridiagonal matrix
!# see algorithm 10.3 in Chapter 10 Arnoldi and Lanczos algorithms in ETH text
 allocate( CQ(NA,kmax),Q(NA),x(NA),res(NA),V(NA))
 allocate( alpha(kmax),beta(kmax) )
 q(:) = 1.d0 ! set initial q
 q(:) = q/sqrt(inner_n(q,q,NA))
 CQ(1:NA,1) = q(1:NA)
 call mul_matcrs_rv(crsmat,q,NA,res) ! res = Aq
 alpha(1)   = inner_n(q,res,NA)      ! alpha = (q,res)
 res        = res - alpha(1)*q
 beta(1)    = sqrt(inner_n(res,res,NA))    ! beta1 = |res|
 j=1
 do j = 2,kmax
  V = q
  q = res/beta(j-1)
  CQ(:,j)=q
  call mul_matcrs_rv(crsmat,q,NA,res) ! res = Aq
  res = res - beta(j-1)*V
  alpha(j) = inner_n(q,res,NA)
  res      = res - alpha(j)*q
  beta(j)  = sqrt(inner_n(res,res,NA))
  if ( beta(j) .eq. 0.d0 ) goto 100
!  write(*,*) "j=",j,"alpha",alpha(j),"beta",beta(j)
 end do

100 continue

 !#[2]## QL algorithm to calculate eigenvalues
 !# D(1:NP) : diagonal     (input)
 !# E(2:NP) : sub-diagonal (input)
 do j=2,kmax
  NP = j
  allocate(D(NP),E(NP),Z(NP,NP))
  E(1)   = 0.d0
  E(2:j) = beta(1:j-1)
  D(1:j) = alpha(1:j)
  N = J
  Z = RESHAPE( [(1,(0,i=1,NP), j = 1, n - 1),  1],  [NP, NP] )
  call TQLI(D,E,N,NP,Z) ! N : # of eigenvalues (output)
  call eigenminmax(D,NP,beta(j),Z(NP,:),Emin,Emax,errmin,errmax)
!  write(*,*) "j=",j,"alpha",alpha(j),"beta",beta(j)
!  write(*,*) "j=",j,"eigen values"
!  do i=1,j
!   write(*,*) "i",i,"D",D(i)
!  end do
!  write(*,*) "Emin",Emin,"errmin",errmin
!  write(*,*) "Emax",Emax,"errmax",errmax
  deallocate(D,E,Z)
 end do
 return

end subroutine LANCZOSMINMAX

!############################################
subroutine eigenminmax(D,NP,beta,S,Emin,Emax,errmin,errmax)
implicit none
real(8),   intent(out) :: Emin,Emax,errmin,errmax
integer(4),intent(in)  :: NP
real(8),   intent(in)  :: D(NP),beta
real(8),   intent(in)  :: S(NP)      ! NP-th components of eigenvectors
integer(4)             :: i
!#[3]## estimate largest and minimum eigenvalues and corresponding errors

  emin=D(1)
  errmin = abs(beta)*S(1)
  emax=D(1)
  errmax = abs(beta)*S(1)
  do i=2,NP
  if ( emin .gt. D(i) ) then
    emin = D(i)
    errmin = abs(beta*S(i))
  end if
  if ( emax .lt. D(i) ) then
    emax = D(i)
    errmax = abs(beta*S(i))
  end if
 end do

return
end
!##########
!# 2017.11.30
! From TQLI.for in numerical recipes
      SUBROUTINE TQLI(D,E,N,NP,Z)
!# C USES pythag
! QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of
! a real, symmetric, tridiagonal matrix, or of a real, symmetric matrix previously
! reduced by tred2 §11.2.
! -
! d is a vector of length np.
!   On input, its first n elements are the diagonal elements
!             of the tridiagonal matrix
!   On output, it returns the eigenvalues.
! -
! The vector e inputs the sub- diagonal elements of the tridiagonal
! matrix, with e(1) arbitrary.
! On output e is destroyed.
! -
! When finding only the eigenvalues, several lines may be omitted,
! as noted in the comments.
! If the eigenvectors of a tridiagonal matrix are desired,
! the matrix z (n by n matrix stored in np by np array) is input
! as the identity matrix.
! If the eigenvectors of a matrix that has been reduced by tred2 are required,
! then z is input as the matrix output by tred2.
! In either case, the kth column of z returns the normalized
! eigenvector corresponding to d(k).
      implicit real*8(a-h,o-z)  ! 2017.12.01
      integer(4),intent(in) :: NP
      integer(4)            :: N
      real(8) ::  D(NP),E(NP),Z(NP,NP) ! 2017.11.30
      IF (N.GT.1) THEN
        DO 11 I=2,N
          E(I-1)=E(I)
11      CONTINUE
        E(N)=0.
        DO 15 L=1,N
          ITER=0
1         DO 12 M=L,N-1
            DD=ABS(D(M))+ABS(D(M+1))
            IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12        CONTINUE
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30)PAUSE 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.*E(L))
            R=SQRT(G**2+1.)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.
            C=1.
            P=0.
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(ABS(F).GE.ABS(G))THEN
                C=G/F
                R=SQRT(C**2+1.)
                E(I+1)=F*R
                S=1./R
                C=C*S
              ELSE
                S=F/G
                R=SQRT(S**2+1.)
                E(I+1)=G*R
                C=1./R
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 13 K=1,N
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
13            CONTINUE
14          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.
            GO TO 1
          ENDIF
15      CONTINUE
      ENDIF
      RETURN
      END
end module

