! Copied from statick_linear/solver_CG.f90 and modified on Aug. 25, 2015
! modified for single core and complex coefficient matrix
!##
!     CG solves the linear system Ax = b using the
!     Conjugate Gradient iterative method with preconditioning.
      subroutine CG13_single_cpx  (NP,  NPL, NPU,                               &
     &                  D,  AL, INL, IAL, AU, INU, IAU,                 &
     &                  B,  X, PRECOND, SIGMA_DIAG, SIGMA,               &
     &                  RESID,  ITER, ERROR,  NSET)
!##
!      use  M_message_time
      implicit none
      integer(4)                 ::  n_log  = 0       !  convergence log file count
      character(len=20  ) ::  log_file_name ! log file name
	integer(4),  intent(in)       :: NP          ! doftot : total # of degree of freedom
      real   (8),     intent(inout) ::  RESID    ! residual
      real   (8),     intent(in)       ::  SIGMA_DIAG
      real   (8),     intent(in)      ::  SIGMA
      integer(4),  intent(inout) ::  ITER
      integer(4),  intent(inout) ::  ERROR
      integer(4),  intent(in)      :: NSET    ! solver parameter
      complex(8), dimension(NP )  , intent(inout)::  D
      complex(8), dimension(NP )  , intent(inout)::  B
      complex(8), dimension(NP )  , intent(inout)::  X
      complex(8), dimension(NPU)  , intent(inout)::  AU ! aiccg_u(item_u_tot)
      complex(8), dimension(NPU)  , intent(inout)::  AL ! aiccg_l(item_l_tot)
	integer(4), intent(in) :: NPU ! item_u_tot
	integer(4), intent(in) :: NPL ! item_l_tot
      integer(4 ), dimension(0:NP)  , intent(in) ::  INU ! istack_u(0:doftot)
      integer(4 ), dimension(  NPU) , intent(in) ::  IAU ! item_u(item_u_tot)
      integer(4 ), dimension(0:NP)  , intent(in) ::  INL ! istack_l(0:doftot)
      integer(4 ), dimension(  NPL) , intent(in) ::  IAL ! item_l(item_l_tot)
      character(len=20), intent(in) :: PRECOND
      complex   (8), dimension(:),   allocatable :: WS, WR
      complex   (8), dimension(:,:), allocatable :: WW
      complex   (8), dimension(:), allocatable, save :: DD
      real   (8), dimension(:), allocatable, save :: SCALE
      integer(4) :: P,Q,R,Z, MAXIT, IFLAG, id, ieL, ieU, isL, isU, i, j, k, kk, inod
      real   (8) :: TOL, DNRM20, DNRM2, BNRM20, BNRM2
	complex(8) :: W, WVAL, SS, SW, RHO, RHO0, RHO1, BETA, C1, C10, ALPHA
      data IFLAG /0/  ! add 0 to IFLAG

!C#[1]## initialization
      write(*,*) "## CG13_simple_cpx START! ###"
      ERROR= 0
      if (IFLAG.eq.0 .and. NSET.eq.0) then  ! For precondition, NSET =1 : yes, 2 : no
	 ERROR= 101
	 return
      endif
      allocate (WW(NP,3))  ! NP is the doftot, while N is the intdoftot
      allocate (WS(NP))
      allocate (WR(NP))
      !
      if (IFLAG.eq.0) then  !----------------------------------  usually go through
	 allocate (DD   (NP))   ! NP is the doftot, while N is the intdoftot
	 allocate (SCALE(NP))
	 IFLAG= 1                   ! After DD and SCALE are allocated, IFLAG becomes 1
	 SCALE(1:NP)= 1.d0 ! SCALE remain 1 if NSET .ne. 2, where NP is the doftot
      endif !------------------------------------------------------------------------
      ! way to use which colmun of WW
      R = 1  ! Residual column of WW
      Z = 2  !
      Q = 2
      P = 3
      MAXIT = ITER   ! max # of iteration
	TOL = RESID    ! tolerance for residual, TOL and RESID are real
      !if (my_rank .eq. 0 ) write(*,*) "MAXIT, TOL=", MAXIT, TOL

!C#[2]## SCALING
!     where A(i,j)*SCALE(i)*SCALE(j) and b(i)*SCALE(i) are used,
!     then X(i)/SCALE(i) will be obtained
      if (NSET.ge.1) then ! nset: 1 for precondition, 2 for precondition with scaling
	 write(*,*) "check0"
	 if (NSET.eq.2) then !------------------------  without  precondition
	  do i= 1, NP  ! N is intnodtot
	   SCALE(i)= 1.d0/dsqrt(cdabs(D(i)))! SCAlE(i)=1/sqrt(|D(i)|)
	  enddo
	  do j= 1, NP !  j is row id, where NP is intnodtot
	   D(j)= cdabs(D(j))/D(j)  !   original D(j) = current D(j) / SCALE (j)^2
!	   if (my_rank .eq. 0) write(*,*) "j=",j,"D(j)=",D(j)
	   isU= INU(j-1) + 1
	   ieU= INU(j  )
	   do i= isU, ieU
	    inod= IAU(i) ! IAU is colmun id of AU
	    AU(i)= AU(i)*SCALE(inod)*SCALE(j) ! SCALING OF UPPER TRIANGLE
	   enddo
	   isL= INL(j-1) + 1
	   ieL= INL(j  )
	   do i= isL, ieL
	    inod= IAL(i)
	    AL(i)= AL(i)*SCALE(inod)*SCALE(j) ! SCALING OF LOWER TRIANGLE
	   enddo
	  enddo
	 endif ! end if (NSET .eq. 2 )
	do i= 1, NP  ! N is the intdoftot
	 B(i)= B(i) * SCALE(i)    !  SCALE(1:NP)=1.d0 without scaling
	enddo


!C#[3]## Preconditioning
      write(*,*) "check1"
	DD(1:NP)= (0.d0,0.d0)
!C#[3-1]## -------IC or ILU----
      if (PRECOND(1:2).eq.'IC' .or. PRECOND(1:3).eq.'ILU') then
        do i= 1, NP
          isL= INL(i-1) + 1
          ieL= INL(i)
          W= D(i) * SIGMA_DIAG
          do k= isL, ieL
            SS=  AL(k)
            id= IAL(k)
            isU= INU(id-1) + 1
            ieU= INU(id)
            do kk= isU, ieU
              SS= SS + AU(kk) * SIGMA
            enddo
            W= W - AL(k)*SS*DD(id)
          enddo
          DD(i)= (1.d0, 0.d0) / W
        enddo
      endif
!C#[3-2]## --SSOR or DIAG----
      if (PRECOND(1:4).eq.'SSOR' .or. PRECOND(1:4).eq.'DIAG') then
        do i= 1, NP ! N is intnodtot
          DD(i)= (1.d0, 0.d0)/(SIGMA_DIAG*D(i))
        enddo
      endif
      endif  ! end if ( NSET .ge. 1 )

!C#[4]## calculate initial residual, r0
!C       +-----------------------+
!C       | {r0}= {b} - [A]{xini} |
!C       +-----------------------+

!C#[5]## BEGIN calculation
      do j= 1, NP    ! calculate WW = B - A*X
!	  write(*,*) "D(j)=",D(j),"B(j)=",B(j)
        WVAL= B(j) - D(j) * X(j)
        isU= INU(j-1) + 1  ! INU is istack_u
        ieU= INU(j  )
        do i= isU, ieU
          inod= IAU(i)                                !  IAU is item_u
          WVAL= WVAL - AU(i) * X(inod)  ! AU is aiccg_u
        enddo
        isL= INL(j-1) + 1
        ieL= INL(j  ) 
        do i= isL, ieL
          inod= IAL(i)
          WVAL= WVAL - AL(i) * X(inod)
        enddo
        WW(j,R)= WVAL     !  R = 1
        ! write(*,*) "WW(j,R)=",WW(j,R)
      enddo
      BNRM20 = 0.d0
      do i= 1, NP
        BNRM20 = BNRM20 + B(i)*conjg(B(i)) ! B is complex, BNRM20 is real
      enddo
	write(*,*) "BNRM20=", BNRM20
      BNRM2=BNRM20
!##  when the right hand side equals to 0
      if (BNRM2 .eq. 0.d0) then
	 X(:)=(0.d0,0.d0)        ! the solusion
	 RESID = 0.d0  ! residul
	 ITER=0           ! No iteration
	 goto 200
	 return
	end if
	write(*,*) "BNRM2=",BNRM2
      ITER = 0

      do iter= 1, MAXIT
!C********************************************* Conjugate Gradient Iteration
!C      +----------------+
!C#[1]##| {z}= [M^-1]{r} |
!C      +----------------+
      do i= 1, NP
        WW(i,Z)= WW(i,R) ! WW(:,R) is residual vector (complex)
      enddo
!      do i= 1+N, NP
!        WW(i,Z)= 0.d0
!      enddo

!C#[2] ## Preconditioning #################
!C#[2-1]##  incomplete CHOLESKY, "IC", "ILU", "SSOR"
      if (PRECOND(1:2).eq.'IC' .or. PRECOND(1:3).eq.'ILU' .or. PRECOND(1:4).eq.'SSOR') then
	 do i= 1, NP
        isL= INL(i-1) + 1
        ieL= INL(i  ) 
        WVAL= WW(i,R)
        do j= isL, ieL
	   inod = IAL(j)
	   WVAL=  WVAL -  AL(j) * WW(inod,Z)
        enddo
        WW(i,Z)= WVAL * DD(i)
       enddo
       do i= NP, 1, -1
        SW  = (0.d0, 0.d0)
        isU= INU(i-1) + 1 ! istack_u
        ieU= INU(i  ) 
        do j= isU, ieU
	   inod = IAU(j) ! IAU is item_u : colmun id
!	   write(*,*) "AU(j)=",AU(j),"WW(inod,Z)=",WW(inod,Z),"SW=",SW,"inod,j=",inod,j
	   SW= SW + AU(j) * WW(inod,Z)
        enddo
!	  write(*,*) "DD(i)=",DD(i),"SW=",SW,"i=",i
        WW(i,Z)= WW(i,Z) - DD(i) * SW
       enddo
      endif

!C#[2-2]##  if precond="DIAG"
      if (PRECOND(1:4).eq.'DIAG') then
       do i= 1, NP
        !write(*,*) "WW(i,R)=",WW(i,R),"DD(i)=",DD(i)
	  WW(i,Z)=  WW(i,R) * DD(i)
       enddo
      endif

!C#[3]## +---------------+
!C       | {RHO}= {r}{z} |
!C       +---------------+
      RHO0= (0.d0, 0.d0)
      do i= 1, NP
	 RHO0= RHO0 + WW(i,R)*WW(i,Z)
      enddo
      RHO=RHO0 ! RHO and RHO0 are complex

!C#[4]## +-----------------------------+
!C       | {p} = {z} if      ITER=1    |
!C       | BETA= RHO / RHO1  otherwise |
!C       +-----------------------------+
      if ( ITER.eq.1 ) then
	 do i= 1, NP
!	  write(*,*) "WW(i,Z)=",WW(i,Z)
	  WW(i,P)= WW(i,Z)  ! WW(NP,3) is allocated above
	 enddo
	else
	 BETA= RHO / RHO1 ! BETA, RHO, RHO1 are complex
	 do i= 1, NP   ! N is intdof tot
	  WW(i,P)= WW(i,Z) + BETA*WW(i,P)
	 enddo
      endif

!C#[5]## +-------------+
!C       | {q}= [A]{p} |
!C       +-------------+ BEGIN calculation
      do j= 1, NP                ! N is intdoftot, loop for j-th internal dof
	 WVAL= D(j) * WW(j,P)
	 isU= INU(j-1) + 1 ! INU is istack_u(0:doftot)
	 ieU= INU(j  )
	 do i= isU, ieU
	  inod= IAU(i)
	  WVAL= WVAL + AU(i) * WW(inod,P)
	 enddo
	 isL= INL(j-1) + 1
	 ieL= INL(j  )
	 do i= isL, ieL
	  inod= IAL(i)
	  WVAL= WVAL + AL(i) * WW(inod,P)
	 enddo
	 WW(j,Q)= WVAL
      enddo

!C#[6]## +---------------------+
!C       | ALPHA= RHO / {p}{q} |
!C       +---------------------+
      C10= (0.d0, 0.d0)
      do i= 1, NP
	 C10= C10 + WW(i,P)*WW(i,Q)
      enddo
      C1=C10
      ALPHA= RHO / C1

!C#[7]## +----------------------+
!C       | {x}= {x} + ALPHA*{p} |
!C       | {r}= {r} - ALPHA*{q} |
!C       +----------------------+
      do i= 1, NP
	 X(i)   = X (i)   + ALPHA * WW(i,P) ! X is solution vector       (complex)
	 WW(i,R)= WW(i,R) - ALPHA * WW(i,Q) ! WW(:,R) is residual vector (complex)
      enddo
      DNRM20 = 0.d0
      do i= 1, NP
	 DNRM20= DNRM20 + WW(i,R)*conjg(WW(i,R)) ! DNRM20 (real)
      enddo
	DNRM2=DNRM20
	RESID= dsqrt(DNRM2/BNRM2)
	write (12,'(i10,1pe16.7)')  ITER, RESID
	if ( mod(ITER,10) .eq. 0 ) write (*, *) "ITR=",ITER,"RESID=", RESID

!C#[8]## check if RESID is small enough or not
	if ( RESID .le. TOL   ) exit ! check whether to end iteration
	if ( ITER .eq.MAXIT ) ERROR= -100
	RHO1 = RHO

    enddo ! iteration loop end!

!C#[9]## recover X and B from scaling
      do i= 1, NP
        X(i)= X(i) * SCALE(i)
        B(i)= B(i) / SCALE(i)
      enddo

      deallocate (WW)
      deallocate (WR)
      deallocate (WS)
        close(12)
200 continue
      write(*,*) "Final RESID=", RESID, "FINAL ITER=",ITER
	write(*,*) "### CG13 END!! ###"
     return
100 continue
	stop
      end