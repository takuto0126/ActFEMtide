! 2017.09.08 pointer -> allocatable
! Coded on 2017.03.06
! to generate face info
module face_type
implicit none

type face_info
 integer(4) :: nface  ! # of face
 integer(4) :: node   ! # of node
 integer(4) :: ntet   ! # of tetrahedron
 integer(4) :: ntri   ! # of triangle
 integer(4),allocatable,dimension(:,:) :: face       ! face(1:3,i) is start and end node
 integer(4),allocatable,dimension(:)   :: face_stack ! line_stack(0:node)
 integer(4),allocatable,dimension(:,:) :: face_item  ! face_item(1:2,i) is 2 node
 integer(4),allocatable,dimension(:,:) :: n4face ! n4face(i,1:4) is face id for ith ele
 integer(4),allocatable,dimension(:,:) :: face2ele !face2ele(1:2,iface) face points 1->2
end type

contains
!######################################## READLINE added on April 21, 2016
subroutine READFACE(filename,g_face)
implicit none
character(50),intent(in) :: filename
type(face_info),intent(out) :: g_face
integer(4) :: nface,ntet,i,j
!
  open(10,file=filename)
   read(10,*) nface
   allocate(g_face%face(3,nface))
   read(10,'(3I10)') ((g_face%face(j,i),j=1,3),i=1,nface)
   read(10,*) ntet
   allocate(g_face%n4face(ntet,4))
   read(10,'(4I10)') ((g_face%n4face(i,j),j=1,4),i=1,ntet)
  close(10)
  g_face%nface=nface
  g_face%ntet =ntet
  write(*,*) "### READFACE     END!! ###" ! 2020.09.17
return
end subroutine READFACE

!######################################## MKFACE
!# Coded on May 6, 2016
subroutine MKFACE(l_face, node, ne, enod, n4)
implicit none
integer(4),            intent(in)     :: node, ne, enod, n4(ne,enod)
type(face_info),       intent(out)    :: l_face
integer(4),            parameter      :: maxnface=6000000
integer(4),allocatable,dimension(:,:) :: face_item_work
integer(4),allocatable,dimension(:)   :: face_stack
integer(4) :: i,j,k,l,m,n1,n2,n3,nface,icount,n1o,n2o,iele
integer(4),parameter :: maxnele=200
integer(4),allocatable,dimension(:,:) :: node2ele
integer(4),allocatable,dimension(:)   :: nele
allocate(nele(node),node2ele(maxnele,node))
allocate(face_item_work(2,maxnface),face_stack(0:node))

!#[0]## gen node-element connectivity data
nele(:)=0
do i=1,ne
 do j=1,enod
  n1=n4(i,j)
  if (maxnele .eq. nele(n1)) goto 998
  nele(n1)=nele(n1)+1     ! increase nele(n1) by 1
  node2ele(nele(n1),n1)=i ! set i-th element to element list of n1 node
 end do
end do
write(*,*) "node-element connectivity end!"

!#[1]## cal face_item_work, l_face%face_stack
nface=0 ; face_stack(:)=0
do n1=1,node ! element loop
  if (mod(n1,100000) .eq. 0) write(*,*) "n1=",n1,"node=",node
  face_stack(n1)=face_stack(n1-1)
  do j=1,nele(n1)     ! element loop for i-th node
    iele=node2ele(j,n1) ! element id
    !  write(*,*) "iele=",iele,"n4(iele,1:4)=",n4(iele,1:4)
    do k=1,enod ! node loop
     do m=1,enod
      n2=n4(iele,k) ; n3=n4(iele,m)
      if ( n1 .lt. n2 .and. n2 .lt. n3 ) then
       do l=face_stack(n1-1)+1, face_stack(n1)
	  ! check whether the face (n1, n2, n3) exist or not
        if ( n2 .eq. face_item_work(1,l) .and. &
	  &    n3 .eq. face_item_work(2,l) ) goto 100
       end do
       if ( nface .eq. maxnface ) goto 999
       face_stack(n1)=face_stack(n1)+1
	 face_item_work(1, face_stack(n1) ) = n2 !line_stack(n1) is latest line id
	 face_item_work(2, face_stack(n1) ) = n3 !line_stack(n2) is latest line id
	 nface=nface+1
      end if
	100 continue
    end do ! 4 node loop
    end do
   end do   ! element loop
!   write(*,*) "node=",n1,"# of lines",&
!       &      line_stack(n1)-line_stack(n1-1),"nele=",nele(n1)
!   write(*,*) "line_item_work=",line_item_work(line_stack(n1-1)+1:line_stack(n1))
!if (n1 .eq. 2) stop
end do      ! node loop
write(*,*) "nface=",nface

!do i=1,nodeg
!write(*,*) "i=",i,"line_stack(i)=",line_stack(i),"line_item(i)=",(line_item(j),j=line_stack(i-1)+1,line_stack(i))
!end do

!#[2]## set line_item and nline
     allocate( l_face%face(3,nface), l_face%face_item(2,nface))
     allocate( l_face%face_stack(0:node) )
     l_face%face_stack(0:node)=face_stack(0:node)
     l_face%face_item(1:2,1:nface)=face_item_work(1:2,1:nface)
     l_face%nface=nface
     l_face%node=node
     if ( enod .eq. 4 ) l_face%ntet=ne
     if ( enod .eq. 3 ) l_face%ntri=ne

!#[3]## set line
icount=0
do i=1,node
  do j=face_stack(i-1)+1, face_stack(i)
    icount=icount+1
    l_face%face(1:3,icount)=(/ i,l_face%face_item(1,j),l_face%face_item(2,j) /)
    !write(*,*) "icount=",icount, "line(1:2,icount)=",line(1:2,icount)
  end do
end do
if ( icount .ne. nface ) goto 99

write(*,*) "### MKFACE END ###"
return
99 continue
write(*,*) "icount=",icount,"is not equal to nface=",nface
stop
999 continue
write(*,*) "GEGEGE! # of face exceeds maxnface =", maxnface
if ( enod .eq. 4 ) write(*,*) "ntet=",ne,"i=",i
if ( enod .eq. 3 ) write(*,*) "ntri=",ne,"i=",i
stop
998 continue
write(*,*) "GEGEGE! # maxnele=",maxnele,"too small for node-element connectivity data"
stop
end subroutine MKFACE

!########################################################### MKN4FACE
! Coded on 2017.03.06
subroutine MKN4FACE(l_face,nodeg,ntetg,n4g) ! make
use fem_edge_util ! for kl (see m_fem_edge_util.f90)
implicit none
type(face_info),intent(inout) :: l_face
integer(4),intent(in) :: ntetg, n4g(ntetg, 4), nodeg
integer(4),allocatable,dimension(:,:) :: n4face
integer(4) :: nface
real(8)    :: r3(3)
integer(4) :: i,j,k, icount, faceid,idirection, n11,n12,n13,n3(3)
integer(4),allocatable,dimension(:)   :: face_stack
integer(4),allocatable,dimension(:,:) :: face_item
allocate(face_stack(0:l_face%node), face_item(2,1:l_face%nface))
allocate(n4face(ntetg,4))

!#[0]# set
nface=l_face%nface
face_stack(0:l_face%node) = l_face%face_stack(0:l_face%node)
face_item(1:2,1:nface)    = l_face%face_item(1:2,1:nface)

!#[1] cal n4face
do i=1,ntetg ! element loop
 do j=1,4    ! face loop
  !#[1-1] set 3 nodes
   n11=n4g(i,lmn(j,1)) ; n12=n4g(i,lmn(j,2)) ; n13=n4g(i,lmn(j,3))
  ! face is defined by MKFACE as { small node id -> large node id}
  if ( ( n11 .lt. n12 .and. n12 .lt. n13 ).or. &
  &    ( n12 .lt. n13 .and. n13 .lt. n11 ).or. &
  &    ( n13 .lt. n11 .and. n11 .lt. n12 ) ) then
   idirection = 1 ! face is outward
  else
   idirection = -1 ! face is inward
  end if
  ! sort
  n3(1:3)=(/n11, n12, n13/) ; r3(1:3)=n3(1:3)*1.d0
  call sort_index(3,n3,r3)
  icount=face_stack(n3(1)-1)
  !#[1-2] search for the node
  do k=face_stack(n3(1)-1)+1, face_stack(n3(1))
   icount=icount+1
   if ( face_item(1,k) .eq. n3(2) .and. &
   &    face_item(2,k) .eq. n3(3) ) goto 100
  end do
  goto 99
  100 continue
  faceid=icount
  n4face(i,j)=idirection*faceid
 end do
 !write(*,*) "i=",i,"ntetg=",ntetg,"n6line(i,1:6)=",n6line(i,1:6)
end do

!#[2]## set n6line
allocate(l_face%n4face(ntetg,4) )
l_face%n4face(1:ntetg,1:4)=n4face(1:ntetg,1:4)

write(*,*) "### MKN4FACE END!! ###"
return
99 write(*,*) "GEGEGE! i=",i,"/ntetg=",ntetg," face n11, n12=",n3(1),n3(2),N3(3), "was not found!"
stop
end subroutine MKN4FACE

!#################################################### FACE2ELEMENT
subroutine FACE2ELEMENT(l_face)
implicit none
type(face_info),intent(inout) :: l_face
integer(4) :: i,j,ntet,iface,nface

!#[1]## set
 ntet  = l_face%ntet
 nface = l_face%nface

!#[2]## cal face2ele
 allocate(l_face%face2ele(2,nface))
 l_face%face2ele = 0 ! 0 is default
 do i=1,ntet
  do j=1,4 ! face loop
   iface = l_face%n4face(i,j)
   if ( iface .gt. 0 ) then
     l_face%face2ele(1, iface) = i ! element where the face is outward
   else if ( iface .lt. 0 ) then
     l_face%face2ele(2,-iface) = i ! element where the face is inward
   else
     goto 100
   end if
  end do
 end do

 write(*,*) "### FACE2ELEMENT END!! ###"
return
100 continue
write(*,*) "GEGEGE"
stop
end subroutine

end module face_type
