!Coded on 2017.05.11
! test of m_matrix.f90
! complile this, for example, by "ifort m_matrix.f90 test.f90"
program test
use matrix ! 2017.09.04
implicit none
!      | 1 0 1 0 0|     | 1 0 0 0 1 |    | 2 0 1 0 1 |
! A1 = | 2 1 3 0 0|  A2=| 0 0 0 1 2 | A3=| 2 1 3 1 2 |
!      | 0 1 0 1 0|     | 1 0 1 0 0 |    | 1 1 1 1 0 |
!
!     | 1 3 1 0 0|      | 1 3 1 0 |
! B = | 1 0 0 2 0|  B2= | 1 0 0 2 |
!     | 2 0 3 1 0|      | 2 0 3 1 |
!     | 0 1 2 4 0|      | 0 1 2 4 |
!     | 0 0 0 0 0|      | 0 0 0 0 |
!
!     | 3 3  4 1 0|
! A*B=| 9 6 11 5 0|
!     | 1 1  2 6 0|
!
type(real_crs_matrix) :: crs, B,B2, C,crsout,A1,A2,A3,A4
type(real_ccs_matrix) :: ccs, ccs2,BCCS, ccsout,B1
type(complex_crs_matrix) :: AC,ACT,crscompout
type(complex_ccs_matrix) :: ACTCCS,ccscompout
complex(8) :: ACO(3,5)
real(8)    :: fullA2(3,5),fullA1(3,5)
real(8) :: threshold=1.d-10
integer(4) :: ntot,nrow,ncolm

! A1
fullA1(1,1:5)=(/1,0,1,0,0/)
fullA1(2,1:5)=(/2,1,3,0,0/)
fullA1(3,1:5)=(/0,1,0,1,0/)
! A2
fullA2(1,1:5)=(/1,0,0,0,1/)
fullA2(2,1:5)=(/0,0,0,1,2/)
fullA2(3,1:5)=(/1,0,1,0,0/)
call conv_full2crs(fullA1,3,5,A1,1.d-10)
call conv_full2crs(fullA2,3,5,A2,1.d-10)

call lateralcomb_real_crs_mat(A1,A2,crsout)
write(*,*) "crs%nrow=",crsout%nrow
write(*,*) "crs%ncolm=",crsout%ncolm
write(*,*) "crs%ntot=",crsout%ntot
write(*,*) "crs%stack=",crsout%stack
write(*,*) "crs%item=",crsout%item
write(*,*) "crs%val=",crsout%val

!stop

!# set B
nrow  = 5
ncolm = 5
ntot  = 11
B%nrow  = nrow
B%ncolm = ncolm
B%ntot  = ntot
allocate(B%item(ntot),B%val(ntot),B%stack(0:nrow))
B%stack(0:nrow) = (/ 0,3,5,8,11,11/)
B%item(1:ntot)  = (/ 1,2,3,1,4,1,3,4,2,3,4 /)
B%val(1:ntot)   = (/ 1,3,1,1,2,2,3,1,1,2,4 /)
B2=B
B2%nrow  = nrow
B2%ncolm = ncolm-1
B2%ntot  = ntot

call mulreal_crs_crs_crs(B2,B,crsout,"T")

write(*,*) "crs%nrow=",crsout%nrow
write(*,*) "crs%ncolm=",crsout%ncolm
write(*,*) "crs%ntot=",crsout%ntot
write(*,*) "crs%stack=",crsout%stack
write(*,*) "crs%item=",crsout%item
write(*,*) "crs%val=",crsout%val

stop

call conv_crs2ccs(B,ccsout)
write(*,*) "ccsout%ncrow=",ccsout%nrow
write(*,*) "ccsout%ncolm=",ccsout%ncolm
write(*,*) "ccsout%ntot=",ccsout%ntot
write(*,*) "ccsout%stack=",ccsout%stack
write(*,*) "ccsout%item=",ccsout%item
write(*,*) "ccsout%val=",ccsout%val
B1=ccsout

call mulreal_crs_ccs_ccs(A1,B1,ccsout)

write(*,*) "ccsout%ncrow=",ccsout%nrow
write(*,*) "ccsout%ncolm=",ccsout%ncolm
write(*,*) "ccsout%ntot=",ccsout%ntot
write(*,*) "ccsout%stack=",ccsout%stack
write(*,*) "ccsout%item=",ccsout%item
write(*,*) "ccsout%val=",ccsout%val

stop

write(*,*) "size(item)=",size(crsout%item)
write(*,*) "size(stack)=",size(crsout%stack)
write(*,*) "crs%nrow=",crsout%nrow
write(*,*) "crs%ncolm=",crsout%ncolm
write(*,*) "crs%ntot=",crsout%ntot
write(*,*) "crs%stack=",crsout%stack
write(*,*) "crs%item=",crsout%item
write(*,*) "crs%val=",crsout%val

stop



call add_crs_crs_crs(A1,A2,A3)

call conv_full2ccs(fullA1,3,5,ccs,1.d-10)
call conv_full2ccs(fullA2,3,5,ccs2,1.d-10)
allocate(ccsout%stack(0:10),ccsout%item(13),ccsout%val(13))
ccsout%nrow = 3 ; ccsout%ncolm=10 ; ccsout%ntot = 13
call combine_real_ccs_mat(ccsout,ccs,1)
call combine_real_ccs_mat(ccsout,ccs2,6)


crsout=A3


!# full matrix
ACO(1,1:5)=(/1,0,1,0,0/)*(0.,1.)
ACO(2,1:5)=(/2,1,3,0,0/)*(0.,1.)
ACO(3,1:5)=(/0,1,0,1,0/)*(0.,1.)
call convcomp_full2crs(ACO,3,5,AC,1.d-10)
call transcomp_crs2ccs(AC,ACTCCS)

ccscompout = ACTCCS
write(*,*) "size(ccsout%stack)=",size(ccscompout%stack)
write(*,*) "ccsout%ncrow=",ccscompout%nrow
write(*,*) "ccsout%ncolm=",ccscompout%ncolm
write(*,*) "ccsout%ntot=",ccscompout%ntot
write(*,*) "ccsout%stack=",ccscompout%stack
write(*,*) "ccsout%item=",ccscompout%item
write(*,*) "ccsout%val=",ccscompout%val
stop

call convcomp_ccs2crs(ACTCCS,ACT)


crscompout = ACT

write(*,*) "crs%nrow=",crscompout%nrow
write(*,*) "crs%ncolm=",crscompout%ncolm
write(*,*) "crs%ntot=",crscompout%ntot
write(*,*) "crs%stack=",crscompout%stack
write(*,*) "crs%item=",crscompout%item
write(*,*) "crs%val=",crscompout%val

stop





!# set A
nrow = 3
ntot = 7
crs%nrow =  nrow
crs%ncolm = 5
crs%ntot =  ntot
allocate(crs%item(ntot),crs%val(ntot),crs%stack(0:nrow))
crs%stack(0:nrow) = (/ 0,2,5,7/)
crs%item(1:ntot)  = (/ 1,3,1,2,3,2,4 /)
crs%val(1:ntot)   = (/ 1,1,2,1,3,1,1 /)

!call conv_crs2ccs(crs,ccsout)
!call conv_ccs2crs(ccsout,crsout)
call trans_crs2ccs(crs,ccsout)



call conv_ccs2crs(ccsout,crsout)

write(*,*) "crs%nrow=",crsout%nrow
write(*,*) "crs%ncolm=",crsout%ncolm
write(*,*) "crs%ntot=",crsout%ntot
write(*,*) "crs%stack=",crsout%stack
write(*,*) "crs%item=",crsout%item
write(*,*) "crs%val=",crsout%val
stop


call conv_crs2ccs(B,BCCS)
call conv_ccs2crs(BCCS,C)

call mulreal_crs_crs_crs(crs,C,crsout)



stop






end program
