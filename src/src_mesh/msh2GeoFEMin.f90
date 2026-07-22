! Coded by Takuto MINAMI on 2015.06.02
! This program converts the *.msh mesh file to GeoFEM input mesh file
! The format of GeoFEM mesh format is described in GeoFEMmanJ.pdf
! In this program, Recursive Coordinate Bisection (RCB) is adopted for partition
!--
! Modified on 2015.07.16 to partition ocean.msh using MeTiS (, which requires /usr/local/lib/libmetis.a)
! This program should be compliled as follows
! $ gfortran m_partition.f90 LOCAL_DATA.f90 msh2GeoFEMin.f90 /usr/local/lib/libmetis.a
! $ ./a.out
program msh2GeoFEMin
use partition
implicit none
character(50),dimension(2) :: meshfile, outmshfile, fileheader
integer(4) :: node,ntet,ntri,nlin,npoi,npower, NP, NP_RCB, NP_MeTiS, NELMGRP ! for ocean.msh
real(8),      allocatable, dimension(:,:) :: xyz
integer(4),   allocatable, dimension(:,:) :: n4,n3,n2
integer(4),   allocatable, dimension(:)   :: n1, IGROUP
integer(4),   allocatable, dimension(:)   :: ELMGRP_STACK, ELMGRP_ITEM
character(50),allocatable, dimension(:)   :: ELMGRPNAME
integer(4),   allocatable, dimension(:)   :: NPN, NPC, ISTACKN, ISTACKC
integer(4),   allocatable, dimension(:)   :: NEIBPETOT, INTNODTOT, NODTOT
integer(4),   allocatable, dimension(:,:) :: NEIBPE, STACK_IMPORT, STACK_EXPORT
character(50),allocatable, dimension(:)   :: WORKFIL, FILNAME
integer(4) :: iflag(2), i
!# RCB parameters idir_npower(1:NPOWER) = 1: xdir, 2: ydir, 3 : zdir
parameter ( NPOWER = 2 , NP_RCB=2**NPOWER) ! NP is number of local mesh after partitionig
integer(4),dimension(NPOWER) :: idir_npower=(/1, 2/)
!# MeTiS
parameter ( NP_MeTiS=4 )
integer(4),allocatable,dimension(:) :: eptr, eind, epart ! arrays for MeTiS
integer(4), pointer     :: vwgt=>null(), vsize=>null(), options=>null()
real(8), pointer     :: tpwgts=>null()
integer(4) :: n
!##
!##  Search which node/cell is related to ip-th local mesh
meshfile(1)="3dmesh.msh" ;             meshfile(2)="ocean.msh"
outmshfile(1)="RCBresult.msh" ; outmshfile(2)="MeTiSocean.msh"
fileheader(1)="GeoFEMmesh" ;    fileheader(2)="GeoFEMocean"
iflag(1)=2 ; iflag(2)=2  ! 1: RCB, 2 : MeTiS
!
do i=1,1
!#[1]##      read mshfile
write(*,*) "### Partition of ",meshfile(i), "start! ###"
CALL MSHCOUNT1(meshfile(i),node,ntet,ntri,nlin,npoi)
   allocate(xyz(node,3),n4(ntet,5),n3(ntri,4),n2(nlin,3),n1(npoi)) ! npoi = 0 makes no problems
CALL READMSH2(meshfile(i),xyz,node,n4,n3,n2,n1,ntet,ntri,nlin,npoi, nelmgrp)
   allocate(ELMGRP_STACK(0:nelmgrp), ELMGRP_ITEM(ntet), ELMGRPNAME(nelmgrp))
CALL GENGRPINFO(ntet, n4(1:ntet,1), nelmgrp, ELMGRP_STACK, ELMGRP_ITEM,ELMGRPNAME)

!#[2]## Partition the mesh by RCB or MeTiS (Recursive Coordinate Bisection)
     if ( iflag(i) .eq. 1) NP=NP_RCB
     if ( iflag(i) .eq. 2) NP=NP_MeTiS
     allocate(IGROUP(node))
     if ( iflag(i) .eq. 1) then !##########  RCB
       CALL RCB3(IGROUP, xyz, node, NPOWER, idir_npower)
     else if (iflag(i) .eq. 2) then !####### MeTiS
       allocate( eptr(ntet+1), eind(ntet*4), epart(ntet) )
       CALL SETADJ(ntet, eptr, eind, n4)
       eind(1:ntet*4)=eind(1:ntet*4)-1 ! node id should start with 0 in MeTiS
       CALL METIS_PartMeshNodal(ntet,node,eptr,eind,vwgt,vsize,NP,tpwgts,options,n,epart,IGROUP)
       write(*,*) "### MeTiS end ! ###"
       IGROUP(1:node)=IGROUP(1:node)+1 ! the group id should start with 1
       deallocate( eptr, eind, epart)
     else
     goto 100
     end if
     CALL OUTGRPNOD( outmshfile(i), xyz, node, n4, ntet, IGROUP, NP)

!#[3]##  create local data
!CALL EDGECUT <- not applied
     allocate( NPN(NP), NPC(NP), ISTACKN(0:NP), ISTACKC(0:NP) )
     allocate( NEIBPE(NP,NP), NEIBPETOT(NP) )
     allocate( STACK_IMPORT(0:NP,NP), STACK_EXPORT(0:NP,NP), INTNODTOT(NP), NODTOT(NP) )
     allocate(  WORKFIL(NP), FILNAME(NP) )
CALL CRE_LOCAL_DATA4(ISTACKN, ISTACKC, NPN,  NPC, NP, ntet, node, IGROUP, n4(:,2:5))
CALL CNT_OVERLAP5(ntet, IGROUP, node, n4(:,2:5))
CALL NEIB_PE6(node, NP, n4(:,2:5), ntet, IGROUP, ISTACKC, NEIBPE, NEIBPETOT)
CALL GENFILE(NP, WORKFIL, FILNAME, fileheader(i))
CALL INTERFACE_NODES(node, NP, STACK_IMPORT, STACK_EXPORT, NEIBPETOT, NEIBPE, WORKFIL,&
& ntet, n4(:,2:5), IGROUP, NPN, NPC, ISTACKN, ISTACKC, NODTOT, INTNODTOT)
!
CALL LOCAL_DATA(node, ntet, NELMGRP, NP, xyz,&
&                               WORKFIL, FILNAME, ELMGRPNAME, &
&                               ELMGRP_STACK, ELMGRP_ITEM, fileheader(i) )
!
!# deallocate
deallocate( xyz, n4, n3, n2, n1)
deallocate( ELMGRP_STACK, ELMGRP_ITEM, ELMGRPNAME)
deallocate( IGROUP)
deallocate( NPN, NPC, ISTACKN, ISTACKC)
deallocate( NEIBPE, NEIBPETOT )
deallocate( STACK_IMPORT, STACK_EXPORT, INTNODTOT, NODTOT )
deallocate(  WORKFIL, FILNAME )
deallocate(NPNID, NPCID)
end do
stop
!
100 continue
write(*,*) "iflag(i) is neither 1 nor 2 iflag(i)=",iflag(i)
end program msh2GeoFEMin
!########################################## SETADJ
subroutine setadj(ntet,eptr, eind, n4)
implicit none
integer(4),intent(in) :: ntet, n4(ntet,5)
integer(4),intent(out) :: eptr(ntet+1), eind(ntet*4)
integer(4) :: i
eptr(1)=0
do i=1,ntet
eptr(i+1)=eptr(i)+4
eind(eptr(i)+1:eptr(i+1))=n4(i,2:5)
end do
write(*,*) "### SETADJ END!! ###"
return
end subroutine setadj
!############################################  OUTGRPNOD
subroutine OUTGRPNOD(em3dfile,xyza,nodea,n4a,nteta, IGROUP, NP)
implicit none
integer(4),intent(in) :: nodea,nteta, NP
integer(4),intent(in) :: n4a(nteta,5), IGROUP(nodea)
real(8),dimension(nodea,3),intent(in) :: xyza
character(50),intent(in) :: em3dfile
integer(4) :: i, j, l, ip, icou, ishift
open(2, file=em3dfile )
write(*,*) "check2"
write(2,'(a11)') "$MeshFormat"
write(*,*) "check3"
write(2,'(a7)') "2.2 0 8"
write(2,'(a14)') "$EndMeshFormat"
write(2,'(a6)') "$Nodes"
write(2,*) nodea
do i=1,nodea
write(2,*) i, (xyza(i,j),j=1,3)
end do
write(2,'(a9)') "$EndNodes"
! out oceanfile ### tetrahedron elements
write(2,'(a9)') "$Elements"
write(2,*) nodea+nteta  !ntets+nteta
icou=0
do ip=1,NP
do i=1,nodea
  if ( IGROUP(i) .eq. ip ) then
  icou=icou+1
  write(2,*) icou,"  15   2   0 ",ip, i
  end if
end do
end do
if ( nodea .ne. icou ) goto 99
ishift=nodea
do i=1,nteta
write(2,*) i+ishift,"  4   2   0",(n4a(i,l),l=1,5) ! n4a(1)=1 for ocean, 2 for air, 3 for land
end do
write(2,'(a12)') "$EndElements"
close(2)
write(*,*) "### OUTGRPNOD END!! ###"
return
99 write(*,*) "## GEGEGE icou .ne. node in OUTGRPNOD ###"; stop
end subroutine OUTGRPNOD
!########################################## GENGRPINFO
subroutine GENGRPINFO(ntet, IELMGRP, nelmgrp, ELMGRP_STACK, ELMGRP_ITEM, ELMGRPNAME)
implicit real(selected_real_kind(8))(a-h,o-z)
integer(4),intent(in) :: ntet, nelmgrp, IELMGRP(ntet)
integer(4),intent(out) :: ELMGRP_STACK(0:nelmgrp), ELMGRP_ITEM(ntet)
character(50),intent(out) :: ELMGRPNAME(nelmgrp)
integer(4) :: IMASK(ntet) ! work array
icou=0 ; ELMGRP_STACK(0)=0
write(*,*) "check1"
icou=0
do ii=1,nelmgrp
  IMASK(:)=0
  do i=1,ntet
    if ( IELMGRP(i) .eq. ii ) IMASK(i)=1 ! check the GROUP ID = ii
  end do
  do i=1,ntet
    if ( IMASK(i) .eq. 1 ) then
    icou=icou+1
    ELMGRP_STACK(ii)=icou  ! replace ELEMGRP_STACK
    ELMGRP_ITEM(icou) = i   ! add node id to ELMGRP_ITEM
    end if
   end do
end do
ELMGRPNAME(1)="OCEAN"
if ( nelmgrp .ge. 2 ) ELMGRPNAME(2)="AIR"
if ( nelmgrp .ge. 3 ) ELMGRPNAME(3)="LAND"
do i=1, nelmgrp
write(*,*) "ELMGRP ID",i,ELMGRPNAME(i)(1:10)," # of elements is", ELMGRP_STACK(i) - ELMGRP_STACK(i-1)
end do
write(*,*) "### GENGRPINFO END!! ###"
return
end
!########################################## GENWORKFIL
subroutine GENFILE(NP, WORKFIL, FILNAME, fileheader)
implicit real(selected_real_kind(8))(a-h,o-z)
integer(4),intent(in) ::NP
character(50),intent(in) ::fileheader
character(50),intent(out) :: WORKFIL(NP), FILNAME(NP)
character(2) :: num
do ip=1,NP
write(num,'(i2.2)') ip
WORKFIL(ip)="work."//num(1:2)
FILNAME(ip)=fileheader(1:len_trim(fileheader))//"."//num(1:2)
end do
return
end
!########################################## Interface_nodes
! Coppied by T. MINAMI on 2015.06.03 from
!  geofem-v.6.0.1/src/util/partitioner/interface_nodes.f90 by RIST
subroutine interface_nodes(N, NP, STACK_IMPORT, STACK_EXPORT, NEIBPETOT, NEIBPE, WORKFIL,&
& IELMTOT, ICELNOD, IGROUP, NPN, NPC, ISTACKN, ISTACKC, NODTOT, INTNODTOT)
use partition ! share NPNID(N2n), NPCIC(N2c) global node/cell id <- ISTACKN(0:NP), ISTACKC(0:NP)
character(50),intent(in) :: WORKFIL(NP)
integer(4),intent(in)    :: N, NP, IELMTOT, ICELNOD(IELMTOT, 4), IGROUP(NP)
integer(4),intent(in)    :: NEIBPETOT(NP), NEIBPE(NP,NP), NPN(NP), NPC(NP)
integer(4),intent(in)    :: ISTACKN(0:NP), ISTACKC(0:NP) ! indicates stack count of NPNID, and NPCID
integer(4),intent(out)  :: STACK_IMPORT(0:NP, NP), STACK_EXPORT(0:NP,NP)
integer(4),intent(out)  :: NODTOT(NP), INTNODTOT(NP)
character*80 LINE
integer(4) :: WHOAMI(NP,NP), IWKM(NP), iwkbuf(100), INODLOCAL(N), NODELM(IELMTOT)
integer(4),allocatable,dimension(:) ::  NOD_IMPORT, IMASK
N2=max(N, IELMTOT)
allocate(IMASK(N2), NOD_IMPORT(N2))
       NODELM(:)=4
        STACK_IMPORT(0,1:NP)= 0
        WHOAMI(:,:)= 0
!#[1]## Generate WHOAMI
!C    pointer WHOAMI defines MY ID in my neighbor's neighboring PE array.
!C    WHOAMI(ip,j) : ip-th PE's ID in his/her j-th neighboring PE's  neighboring PE array
      do ip  = 1, NP
      do inei= 1, NEIBPETOT(ip)
        ip1= NEIBPE(ip,inei)
        do in1= 1, NEIBPETOT(ip1)
          if (NEIBPE(ip1,in1).eq.ip) WHOAMI(ip,inei)= in1
        enddo
      enddo
      enddo
!C +---------------------------------------------+
!C | create INITIAL FILE : LOCAL/IMPORT pointers |
!C +---------------------------------------------+
      do ip= 1, NP
!C-- define INTERIOR and EXTERIOR NODEs
        IMASK(1:N)=0
	  NODTOT(ip)= NPN(ip)       ! NODTOT(ip)  will be  # of internal + external nodes pf ip-th PE
        INTNODTOT(ip)= NPN(ip) ! INTNODTOT(ip) :  # of internal nodes of ip-th PE
!  add global node id of internal nodes for ip-th PE to INODLOCAL
        icou= 0
        do is= ISTACKN(ip-1)+1, ISTACKN(ip)
            in= NPNID(is)
          icou= icou + 1
          INODLOCAL(icou)= in   ! INODLOCAL(i) is global node id of i-th node in local PE
        enddo
! add global node id of external nodes of ip-th PE to INODELOCAL
        do is= ISTACKC(ip-1)+1, ISTACKC(ip)
          icel= NPCID(is)
        do k = 1, NODELM(icel)
          in= ICELNOD(icel,k)
          ig= IGROUP (in)
          if (ig.ne.ip ) then
!         originally subroutine MARK_INT add global node id of external nodes of ip-th PE to INODLOCAL
          if (IMASK(in).ne.1) then
            IMASK(in)= 1
            NODTOT(ip)= NODTOT(ip) + 1
            INODLOCAL( NODTOT(ip) )= in
          endif
	    end if
	  enddo
	  enddo
!   write the initial information on Workfile
          open (11,file=WORKFIL(ip), status='unknown',form='unformatted')
          write (11) ip
          write (11) NODTOT(ip)
          write (11) INTNODTOT(ip)
          write (11) NEIBPETOT(ip)
          write (11) (NEIBPE(ip,inei),inei= 1, NEIBPETOT(ip))
          do in= 1, NODTOT(ip)
            in1= INODLOCAL(in) ! global node id of internal and external nodes in ip-th PE
            write (11) in1
          enddo
!C-- ASSEMBLE EXTERIOR NODEs
!C            IMPORT pointers
        STACK_IMPORT(0,ip)= 0
        do inei= 1, NEIBPETOT(ip)
          icou= 0
            ig= NEIBPE(ip,inei)
          do it= INTNODTOT(ip)+1, NODTOT(ip) ! loop for external nodes of ip-th PE
            in= INODLOCAL(it)            ! in is the local node id / it is the global node id
            if (IGROUP(in).eq.ig) then
                    icou = icou + 1           !
              IMASK(icou)= it
            endif
          enddo

          STACK_IMPORT(inei,ip)= STACK_IMPORT(inei-1,ip) + icou
          do ic= 1, icou
            is= ic + STACK_IMPORT(inei-1,ip)
            NOD_IMPORT(is)= IMASK(ic)
          enddo
        enddo

          write (11) (STACK_IMPORT(inei,ip), inei= 1, NEIBPETOT(ip))


          do is= 1, STACK_IMPORT(NEIBPETOT(ip),ip)
            nodL= NOD_IMPORT(is)
            nodG= INODLOCAL(nodL)
            write (11) nodL, nodG, IGROUP(nodG)
          enddo

          close (11)
      enddo
!C===

!C
!C +-------------------------------+
!C | update FILE : EXPORT pointers |
!C +-------------------------------+
!C===

!C
!C-- INIT.

      do ip= 1, NP
        STACK_EXPORT(0,ip)= 0
      do inei= 1, NEIBPETOT(ip)
        ig  = NEIBPE(ip,inei)
        iobj= WHOAMI (ip,inei)
        STACK_EXPORT(inei,ip)= STACK_EXPORT(inei-1,ip) +                &
     &                         STACK_IMPORT(iobj,ig) -                  &
     &                               STACK_IMPORT(iobj-1,ig)
      enddo
      enddo

      do ip= 1, NP
        open   (11,file=WORKFIL(ip),status='unknown',form='unformatted')
        rewind (11)
          read (11) ID
          read (11) N1
          read (11) N2
          read (11) N3         
          read (11) (IWKM(i),i=1,N3)

          do in= 1, N1
            read (11) INODLOCAL(in)
          enddo

          read (11) (STACK_IMPORT(inei,ip), inei= 1, NEIBPETOT(ip))

          do is= 1, STACK_IMPORT(NEIBPETOT(ip),ip)
            read (11) NOD_IMPORT(is)
          enddo

        close (11)          

!C
!C-- "masking" with GLOBAL NODE ID

        open   (12,file=WORKFIL(ip),status='unknown',form='unformatted')
        rewind (12)
          write (12) ID
          write (12) N1
          write (12) N2
          write (12) N3         
          write (12) (IWKM(i),i=1,N3)

          do in= 1, N1
            in1= INODLOCAL(in)
            write (12) in1
          enddo

          write (12) (STACK_IMPORT(inei,ip), inei= 1, NEIBPETOT(ip))

          do is= 1, STACK_IMPORT(NEIBPETOT(ip),ip)
            nodL= NOD_IMPORT(is)
            nodG= INODLOCAL(nodL)
            write (12) nodL, nodG, IGROUP(nodG)
          enddo


        do i= 1, NODTOT(ip)
                in = INODLOCAL(i)
          IMASK(in)= i
        enddo
          
        write (12) (STACK_EXPORT(inei,ip), inei= 1, NEIBPETOT(ip))
       
        do inei= 1, NEIBPETOT(ip)
          ig= NEIBPE(ip,inei)
          open   (11,file=WORKFIL(ig),status='unknown',form='unformatted')
          rewind (11)

          iobj= WHOAMI (ip,inei)

          read (11) ID
          read (11) N1
          read (11) N2
          read (11) N3         
          read (11) (IWKM(i),i=1,N3)

          do i= 1, N1
            read (11) INODLOCAL(i)
          enddo

!C
!C-- process NEIGHBOR's IMPORT pointers

          STACK_IMPORT(0,ig)= 0
          read (11) (STACK_IMPORT(k,ig), k=1, N3)
          do is= 1, STACK_IMPORT(N3,ig)
            read (11) i1,NOD_IMPORT(is)          
          enddo

!C
!C-- TRANSFER INFO. into "MY" EXPORT pointers
          N2n= 2*max(N,IELMTOT)
          allocate (NOD_EXPORT(N2n))
  100     continue

          ISTART= STACK_EXPORT(inei-1,ip)
          do is= STACK_IMPORT(iobj-1,ig)+1, STACK_IMPORT(iobj,ig)
            icur= is - STACK_IMPORT(iobj-1,ig)
            inod= NOD_IMPORT(is)

            if (ISTART+icur.gt.N2n) then
              deallocate (NOD_EXPORT)
              N2n= N2n * 11/10 + 1
              allocate (NOD_EXPORT(N2n))
              goto 100
            endif
           
            NOD_EXPORT(ISTART+icur)= IMASK(is)
            write (12) IMASK(inod), inod, ig
          enddo
          deallocate (NOD_EXPORT)
        enddo

        write (12) NPC(ip)
        do is= ISTACKC(ip-1)+1, ISTACKC(ip)
           is0 = is - ISTACKC(ip-1)
           icel= NPCID(is)
           do kkk= 1, NODELM(icel)
             iwkbuf(kkk)= IMASK(ICELNOD(icel,kkk))
           enddo
           write (12) icel, (iwkbuf(k), k=1, NODELM (icel))
        enddo
       
        close (11)
        close (12)

      enddo
      return
      end
!##########################################
      subroutine MARK_INT (in,ip, NODTOT, INODLOCAL, G_NODE_MAX,IMASK)
      return
      end
!########################################## NEIB_PE6
! Coppied by T. MINAMI on 2015.06.03 from
!  geofem-v.6.0.1/src/util/partitioner/NEIB_PE.f90 by RIST
subroutine NEIB_PE6(N, NP, ICELNOD, IELMTOT, IGROUP, ISTACKC, NEIBPE, NEIBPETOT)
use partition ! share N2n, N2c, NPNID, NPCID
integer(4),intent(in) ::  N, NP, IELMTOT, ICELNOD(IELMTOT,4),IGROUP(N),ISTACKC(0:NP)
integer(4),intent(out) :: NEIBPETOT(NP), NEIBPE(NP, NP)
integer(4),dimension(IELMTOT) :: NODELM
NODELM(1:IELMTOT)=4
do ip= 1, NP
  NEIBPETOT(ip)= 0
  do is= ISTACKC(ip-1)+1, ISTACKC(ip)
	icel= NPCID(is)
	do  k= 1, NODELM(icel)
	   in= ICELNOD(icel,k)
	   ig= IGROUP (in)
	   if (ig.ne.ip) call FIND_NEIBPE (ig, ip,NEIBPE, NEIBPETOT, NP)
      enddo
  enddo
enddo
!
      write ( *,'(/," PE/NEIB-PE#    NEIB-PEs")')
      do ip= 1, NP
        write ( *,'(i3,i4,5x, 31i4)') ip, NEIBPETOT(ip), (NEIBPE(ip,k),k=1,NEIBPETOT(ip))
      enddo
	write(*,*) "### NEIB_PE6 END!! ###"
      return
end
!####################################################
      subroutine FIND_NEIBPE (ig, ip, NEIBPE, NEIBPETOT, NP)
      implicit real(selected_real_kind(8))(a-h,o-z)
	integer(4), intent(in) :: ig, ip
	integer(4), intent(inout) :: NEIBPE(NP, NP), NEIBPETOT(NP)
      do inei= 1, NEIBPETOT(ip)
        if (ig.eq.NEIBPE(ip,inei)) return
      enddo
      NEIBPETOT(ip) = NEIBPETOT(ip) + 1
      NEIBPE(ip,NEIBPETOT(ip))= ig
      return
      end
!########################################## CNT_OVERLAP5
! Coppied by T. MINAMI on 2015.06.03 from
!  geofem-v.6.0.1/src/util/partitioner/PROC_LOCAL.f90 by RIST
subroutine CNT_OVERLAP5(IELMTOT, IGROUP, N, ICELNOD)
implicit real(selected_real_kind(8))(a-h,o-z)
integer(4),intent(in) :: IELMTOT, N, IGROUP(N), ICELNOD(IELMTOT,4)
integer(4) :: NODELM(IELMTOT), ISTACK(IELMTOT)
      NODELM(:)=4
	ISTACK(:)= 0 ! work array
!
      do icel= 1, IELMTOT
        do k1= 1, NODELM(icel)
        do k2= 1, NODELM(icel)
          ig1= IGROUP(ICELNOD(icel,k1))
          ig2= IGROUP(ICELNOD(icel,k2))
          if (ig1.ne.ig2) ISTACK(icel)= 1
        enddo
        enddo
      enddo
!
      icou= 0
      do icel= 1, IELMTOT
        if (ISTACK(icel).eq.1) icou= icou + 1
      enddo
      write ( *,'(/,"OVERLAPPED ELEMENTS", i12)')  icou
      write(*,*) "### CNT_OVERLAP END!! ###"
return
end
!
!########################################## mshcount1
subroutine mshcount1(mshfile,node,ntet,ntri,nlin,npoi)
implicit none
integer(4),intent(out) :: node,ntet,ntri,nlin,npoi
integer(4) :: i,j,k,itype,ii,jj,kk,etype,nele,n44(4)
character(50),intent(in) :: mshfile
open(1,file=mshfile)
do i=1,4
read(1,*)
end do
read(1,*) node
!write(*,*) "node=",node
do i=1,node+2 ! include "$EndNodes", "$Elements" lines
read(1,*)
end do
read(1,*) nele ! total number of elements
npoi=0;nlin=0;ntri=0;ntet=0
do i=1,nele
read(1,*) j,itype,ii,jj,kk,(n44(k),k=1,etype(itype))
if ( itype .eq. 15) then   ! Point element
npoi=npoi+1
else if ( itype .eq. 1) then ! Line element
nlin=nlin+1
else if ( itype .eq. 2) then   ! Triangle element
ntri=ntri+1
else if ( itype .eq. 4) then ! Tetrahedral element
ntet=ntet+1
end if
end do
close(1)
!write(*,*) "# of Point elements is",npoi
!write(*,*) "# of Line elements is",nlin
!write(*,*) "# of Triangle elements is",ntri
!write(*,*) "# of Tetrahedron elements is",ntet
write(*,*) "### COUNT ", mshfile(1:len_trim(mshfile))," END!! ###"
return
end subroutine mshcount1
!
!##########################################  readmsh2
subroutine readmsh2(mshfile,xyz,node,n4,n3,n2,n1,ntet,ntri,nlin,npoi, IELMGRPTOT)
implicit none
integer(4),intent(in) :: node,ntet,ntri,nlin,npoi
real(8),dimension(node,3),intent(out) :: xyz
integer(4),intent(out) :: n4(ntet,5),n3(ntri,4),n2(nlin,3),n1(npoi), IELMGRPTOT
integer(4) :: itet,itri,ilin,ipoi,nele,inode, IELMGRP(100) ! 100 is the max of IELMGRPTOT
integer(4) :: i,j,k,ii,jj,kk,itype,etype
character(50) :: mshfile
integer(4),dimension(4) :: n44
open(1,file=mshfile)
!# [1] ### skip header
do i=1,4
read(1,*)
end do
!# [2] ### read node
read(1,*) inode
write(*,*) "# of nodes (node) =",inode
if ( inode .gt. node) goto 999
do i=1,node
read(1,*) j, (xyz(i,k),k=1,3)
end do
read(1,*)
!# [3] ### read elements
read(1,*) ! skip the starting line, "$Elements"
read(1,*) nele
write(*,*)"# of elements (nele)=",nele
ipoi=0;ilin=0;itri=0;itet=0
do i=1,nele
read(1,*) j,itype,ii,jj,kk,(n44(k),k=1,etype(itype))
!write(*,*) j,itype,ii,jj,kk,(n44(k),k=1,etype(itype))
if ( itype .eq. 15) then   ! Point element
ipoi=ipoi+1
n1(ipoi)=n44(1)
else if ( itype .eq. 1) then ! Line element
ilin=ilin+1
n2(ilin,1)=kk
n2(ilin,2:3)=n44(1:2)
else if ( itype .eq. 2 ) then ! Triangle element
itri=itri+1
n3(itri,1)=kk
n3(itri,2:4)=n44(1:3)
else if ( itype .eq. 4) then ! Tetrahedral element
itet=itet+1
n4(itet,1)=kk
n4(itet,2:5)=n44(1:4)
end if
end do
IELMGRPTOT=1; IELMGRP(1)=n4(1,1)
do itet=1,ntet
  do j=1, IELMGRPTOT !When the group id of itet element is already present, goto 100
    if (IELMGRP(j) .eq. n4(itet,1) ) goto 100
  end do
  IELMGRPTOT=IELMGRPTOT+1
  IELMGRP(IELMGRPTOT)=n4(itet,1)
100 continue
end do
write(*,*) "# of Point elements is",npoi
write(*,*) "# of Line elements is",nlin
write(*,*) "# of Triangle elements is",ntri
write(*,*) "# of Tetrahedron elements is",ntet
write(*,*) "# of Tetrahedron group is", IELMGRPTOT
close(1)
write(*,*) "### READ ", mshfile(1:len_trim(mshfile))," END!! ###"
return
999 write(*,*) "GEGEGE node .ne. inode, node=",node,"inode=",inode
stop
end subroutine readmsh2
!###########################################  function etype
function etype(itype)
implicit none
integer(4) :: itype,etype
etype=0
if (itype .eq. 15) etype=1! one point node
if (itype .eq. 1 ) etype=2 ! line
if (itype .eq. 2 ) etype=3 ! triangle
if (itype .eq. 4 ) etype=4 ! tetrahedron
return
end function
!##########################################
! This RCB subroutine is extracted by T. MINAMI from
! geofem-v.6.0.1/src/util/partitioner/geofem_tiger.f90 by RIST
subroutine RCB3(IGROUP, xyz, N, NPOWER, idir_npower)
implicit real(selected_real_kind(8))(a-h,o-z)
integer(4),intent(in) :: N, NPOWER, idir_npower(NPOWER)
real(8),intent(in) :: xyz(N,3)
integer(4),intent(out) :: IGROUP(N)
real(8) :: VAL(N)
integer(4) :: IS1(N), IS2(-N:+N)
!C
!C +-----+
!C | RCB |
!C +-----+
!C===
      do i= 1, N
        IGROUP(i)= 1
      enddo

      do iter= 1, NPOWER
	 idir = idir_npower(iter) ! set idir
        do ip0= 1, 2**(iter-1)
          icou= 0
          do i= 1, N
            if (IGROUP(i).eq.ip0) then
                icou= icou + 1
              IS1(icou)= i
              VAL(icou)= XYZ(i,idir)
            endif
          enddo

          call SORT (VAL, IS1, IS2, N, icou)

          do ic= 1, icou/2
            in= IS1(ic)
            IGROUP(in)= ip0 + 2**(iter-1)
          enddo

!          ic0= 0
!          do i= 1, N
!            if (IGROUP(i).eq.ip0) then
!              ic0= ic0 + 1
!              IGROUP(i)= ip0 + 2**(iter-1)
!              if (ic0.eq.icou/2) exit
!            endif
!          enddo

!C
!C-- KL optimization
         ip1= ip0 + 2**(iter-1)

!         call KL (ip0,ip1,icou)

!          NXP1=  5
!          NYP1= 11
!          do jj= NYP1,1,-1
!            do ii= 1, NXP1
!              LINE(ii)= CHAR(IGROUP((jj-1)*NXP1+ii))
!            enddo
!            write (*,'(5a2)')( LINE(k),k=1,NXP1)
!          enddo


        enddo

      enddo
      write(*,*) "### RCB3 END!! ###"
return
end subroutine RCB3
!#####################################################
!  Coppied by T. MINAMI on June 2, 2015 from
!  geofem-v.6.0.1/src/util/partitioner/geofem_tiger.f90 by RIST
!  for RCB partition of initial mesh
      subroutine SORT (STEM, INUM, ISTACK, NP, N)
      real   (kind=8), dimension(NP)      ::  STEM
      integer(kind=4), dimension(NP)      ::  INUM
      integer(kind=4), dimension(-NP:+NP) ::  ISTACK

      M     = 100
      NSTACK= NP

      jstack= 0
      l     = 1
      ir    = N

      ip= 0
 1    continue
      ip= ip + 1

      if (ir-l.lt.M) then
        do 12 j= l+1, ir
          ss= STEM(j)
          ii= INUM(j)

          do 11 i= j-1,1,-1
            if (STEM(i).le.ss) goto 2
            STEM(i+1)= STEM(i)
            INUM(i+1)= INUM(i)
 11       continue
          i= 0

 2        continue
            STEM(i+1)= ss
            INUM(i+1)= ii
 12     continue

        if (jstack.eq.0) return
        ir = ISTACK(jstack)
         l = ISTACK(jstack-1)
        jstack= jstack - 2
       else

        k= (l+ir) / 2
            temp = STEM(k)
        STEM(k)  = STEM(l+1)
        STEM(l+1)= temp

              it = INUM(k)
        INUM(k)  = INUM(l+1)     
        INUM(l+1)= it

        if (STEM(l+1).gt.STEM(ir)) then
              temp = STEM(l+1)
          STEM(l+1)= STEM(ir)
          STEM(ir )= temp
                it = INUM(l+1)
          INUM(l+1)= INUM(ir)
          INUM(ir )= it
        endif

        if (STEM(l).gt.STEM(ir)) then
             temp = STEM(l)
          STEM(l )= STEM(ir)
          STEM(ir)= temp
               it = INUM(l)
          INUM(l )= INUM(ir)
          INUM(ir)= it
        endif

        if (STEM(l+1).gt.STEM(l)) then
              temp = STEM(l+1)
          STEM(l+1)= STEM(l)
          STEM(l  )= temp
                it = INUM(l+1)
          INUM(l+1)= INUM(l)
          INUM(l  )= it
        endif

        i= l + 1
        j= ir

        ss= STEM(l)
        ii= INUM(l)

 3      continue
          i= i + 1
          if (STEM(i).lt.ss) goto 3

 4      continue
          j= j - 1
          if (STEM(j).gt.ss) goto 4     

        if (j.lt.i)        goto 5

        temp   = STEM(i)
        STEM(i)= STEM(j)
        STEM(j)= temp

        it     = INUM(i)
        INUM(i)= INUM(j)
        INUM(j)= it
 
        goto 3
      
 5      continue

        STEM(l)= STEM(j)
        STEM(j)= ss
        INUM(l)= INUM(j)
        INUM(j)= ii

        jstack= jstack + 2

        if (jstack.gt.NSTACK) then
          write (*,*) 'NSTACK overflow'
!          call  MPI_FINALIZE (errno)
          stop
        endif

        if (ir-i+1.ge.j-1) then
          ISTACK(jstack  )= ir
          ISTACK(jstack-1)= i
          ir= j-1
         else
          ISTACK(jstack  )= j-1
          ISTACK(jstack-1)= l
          l= i
        endif 

      endif     

      goto 1

      end subroutine SORT
!###################################################
! Copied by T. MINAMI on 2015.06.02 from
! !  geofem-v.6.0.1/src/util/partitioner/CRE_LOCAL_DATA.f90 by RIST
subroutine CRE_LOCAL_DATA4(ISTACKN, ISTACKC, NPN, NPC, NP, IELMTOT, N, IGROUP, ICELNOD)
use partition ! N2n, N2c, NPNID, NPCID are shared
integer(4),intent(in) :: NP, IELMTOT, N, IGROUP(N), ICELNOD(IELMTOT, 4)
integer(4),intent(out) :: ISTACKN(0:NP), ISTACKC(0:NP), NPN(NP), NPC(NP)
! MAXintN is the max number of local internal nodes
integer(4),allocatable,dimension(:) ::  ISTACK, NODELM
	N2n=2*max(N,IELMTOT)
	N2c=N2N
	allocate (NPNID(N2n), NPCID(N2c), ISTACK(N2n), NODELM(IELMTOT))
      NODELM(1:IELMTOT)=4 ! assume all the elemets are tetrahedrons
      ISTACKN(0)= 0
      ISTACKC(0)= 0
  100 continue
      do ip= 1, NP
	 ISTACK(1:IELMTOT)= 0
	 do icel= 1, IELMTOT
	 do    k= 1, NODELM(icel)
          in= ICELNOD(icel,k)
          ig= IGROUP (in)
          if (ig.eq.ip) ISTACK(icel)= 1
        enddo
        enddo
!
        icou= 0
        do icel= 1, IELMTOT
          if (ISTACK(icel).eq.1) then
            icou= icou + 1
              is= ISTACKC(ip-1) + icou
              if (is.gt.N2c) then
                deallocate (NPCID)
                N2c= N2c * 11/10 + 1
                  allocate (NPCID(N2c))
			goto 100
              endif
            ISTACKC(ip)= is
              NPC  (ip)= icou
              NPCID(is)= icel
          endif
        enddo
      enddo
!
      do ip= 1, NP
        NPN(ip)= 0
      enddo

      do i= 1, N
        ig= IGROUP(i)
        NPN(ig)= NPN(ig) + 1
      enddo
!
      do ip= 1, NP
        ISTACKN(ip)= ISTACKN(ip-1) + NPN(ip)
        NPN(ip)= 0
      enddo
!
      do i= 1, N
        ip = IGROUP(i)
        icou= NPN(ip) + 1
        is = ISTACKN(ip-1) + icou
              if (is.gt.N2n) then
                deallocate (NPNID)
                N2n= N2n * 11/10 + 1
                  allocate (NPNID(N2n))
                goto 100
              endif
          NPN(ip)= icou
        NPNID(is)= i
      enddo
!
      MAXN= NPN(1)
      MINN= NPN(1)
      MAXC= NPC(1)
      MINC= NPC(1)
      
      do ip= 2, NP
        MAXN= max (MAXN,NPN(ip))
        MINN= min (MINN,NPN(ip))
        MAXC= max (MAXC,NPC(ip))
        MINC= min (MINC,NPC(ip))
      enddo
      write ( *,'(/,"TOTAL NODE     #   ", i12)') N
      write ( *,'(  "TOTAL CELL     #   ", i12)') IELMTOT
      write ( *,'(/," PE    NODE#   CELL#")')
      do ip= 1, NP
        write ( *,'(i3,5i8)') ip, NPN(ip), NPC(ip)
      enddo
      write(*,*) "### CRE_LOCAL_DATA4 END!! ###"
      return
      end










