! coded on 2022.10.20
module freq_mpi_joint
use param
use param_mt
implicit none

type freq_info_joint
 integer(4) :: np
 integer(4) :: nfreq_tot    !  MT + ACT2022.10.20
 integer(4) :: nfreq_tot_ip ! 2022.10.20
 real(8),   allocatable,dimension(:)   :: freq_tot       ! (nfreq_tot)
 real(8),   allocatable,dimension(:,:) :: freq_tot_ip    ! (nfreq_tot_ip,0:np-1)
 integer(4),allocatable,dimension(:,:) :: ifreq_tot_ip   ! (nfreq_tot_ip,0:np-1)
 !# ACTIVE
 integer(4) :: nfreq_act    !  ACTIVE frequency 2022.10.20
 integer(4),allocatable,dimension(:) :: nfreq_act_ip ! (0:np-1)2022.10.20
 integer(4),allocatable,dimension(:,:) :: i_act_ip   ! (nfreq_tot_ip,0:np-1)
 integer(4),allocatable,dimension(:) :: ip_from_act  ! (nfreq_act)
 integer(4),allocatable,dimension(:) :: if_g2l_act   ! (nfreq_act)

 !# MT
 integer(4) :: nfreq_mt    !  MT frequency 2022.10.20
 integer(4),allocatable,dimension(:)   :: nfreq_mt_ip ! (0:np-1)2022.10.20
 integer(4),allocatable,dimension(:,:) :: i_mt_ip   ! (nfreq_tot_ip,0:np-1)
 integer(4),allocatable,dimension(:)   :: ip_from_mt  ! (nfreq_mt)
 integer(4),allocatable,dimension(:)   :: if_g2l_mt   ! (nfreq_mt)
end type

contains

!#################################################################################################
!# Coded on 2022.10.20
subroutine SETFREQIPJOINT(g_param,g_param_mt,ipp,np,g_freq_joint)
implicit none
integer(4),            intent(in)     :: np,ipp
type(param_forward),   intent(in)     :: g_param
type(param_forward_mt),intent(in)     :: g_param_mt
type(freq_info_joint), intent(inout)  :: g_freq_joint
!#
integer(4),allocatable,dimension(:)   :: if_act2tot,if_mt2tot,ip_from_tot
integer(4),allocatable,dimension(:,:) :: if_tot_ip
integer(4),allocatable,dimension(:,:) :: if_act_ip, if_mt_ip !(nf_tot_ip,0:np-1)
integer(4),allocatable,dimension(:,:) :: ind_f_tot, ind,ip_from_f_tot
integer(4),allocatable,dimension(:)   :: nf_act_ip,nf_mt_ip
integer(4),allocatable,dimension(:,:) ::  i_act_ip, i_mt_ip
integer(4),allocatable,dimension(:)   :: ip_from_act,ip_from_mt
integer(4),allocatable,dimension(:)   :: if_g2l_act, if_g2l_mt
real(8),   allocatable,dimension(:)   :: f_mt,f_act,f_tot
real(8),   allocatable,dimension(:)   :: freq
real(8),   allocatable,dimension(:,:) :: f_tot_ip !
integer(4)                            :: i,j,ii,ij,jj,ip
integer(4)                            :: nf_act, nf_mt       ! 2022.10.20
integer(4)                            :: nf_tot, nf_tot_ip   ! 2022.10.20
integer(4)                            :: nf_tot_n,nn,nn2,nf_tmp,ip_from
integer(4)                            :: icount_act,icount_mt

!#[1]##===================================================================  global frequency
!#[1-1]## set frequencies
  !# set ACTIVE frequency
  nf_act   = g_param%nfreq
  f_act    = g_param%freq
  !# set MT frequency
  nf_mt    = g_param_mt%nfreq
  f_mt     = g_param_mt%freq

!#[1-2]## set f_tot, nf_tot, ind_f_tot ## merge the same frequency between ACTIVE and MT)
  nn   = nf_act + nf_mt
  allocate( freq(nn),ind(2,nn))
  call mkfreqtot(nf_act,nf_mt,f_act,f_mt,nn,freq,nf_tot,ind,ipp) ! get nf_tot 2022.12.05
  ! nn total frequencies -> nf_tot frequencies by merging the same frequencies
  allocate( ind_f_tot(2,nf_tot), f_tot(nf_tot)) ! nf_tot is the final # of frequencies
  f_tot          = freq( 1:nf_tot)
  ind_f_tot(:,:) = ind(:,1:nf_tot) ! ind(1,i) : index for active freq, ind(2,i) : index for mt freq

!#[1-3]## set ip_from_f_tot(2,nfreq_tot) : i-th frequency -> ip 
   allocate( ip_from_tot(nf_tot))
   do i=1,nf_tot
    if (mod(i,np) .eq. 1 ) ip_from = -1
    ip_from = ip_from + 1 
    ip_from_tot(i) = ip_from
   end do
   call checkfreq_tot( nf_act,nf_mt,nn,nf_tot,f_tot,ind_f_tot,ip_from_tot,ipp) !2022.12.05

!#[1-4]## set if_act2tot, if_mt2tot
   allocate(if_act2tot(nf_act), if_mt2tot( nf_mt ))
   do i=1,nf_tot
    if ( ind_f_tot(1,i) .ne. 0 ) if_act2tot(ind_f_tot(1,i)) = i ! active freq index
    if ( ind_f_tot(2,i) .ne. 0 ) if_mt2tot( ind_f_tot(2,i)) = i ! mt freq index
   end do

!#[2]##=============================================== MPI node for global total frequency  
  !#[2-1]## nf_tot_ip     : # of frequencies to be caluclated at each node 
   nf_tot_ip = nf_tot/np  ! 2022.12.05 
   if ( nf_tot - nf_tot_ip*np .gt. 0 ) nf_tot_ip = nf_tot_ip + 1 ! nf_tot_ip is equal for all ip
   if ( ipp .eq. 0 ) write(*,'(a,i5)') " nf_tot_ip =",nf_tot_ip ! 2022.12.05

  !#[2-2]## f_tot_ip, if_tot_ip
   allocate(  f_tot_ip( nf_tot_ip,0:np-1) ) ! frequency for each node
   allocate( if_tot_ip( nf_tot_ip,0:np-1) ) ! freqency index (to global freq) for each node
   f_tot_ip( :,:)  = -1.d0 ! default
   if_tot_ip(:,:) = 0
   do ip = 0,np-1
     do i=1,nf_tot_ip ! # of freq for each node
       ii = ip+1 + (i-1)*np ! ii depends on ip
       if (ii .le. nf_tot)  f_tot_ip(i,ip)  = f_tot(ii) ! freqency               (depends on ip)
       if (ii .le. nf_tot)  if_tot_ip(i,ip) = ii        ! global frequency index (depends on ip)
     end do
   end do
   if (ipp .eq. 0) then    ! 2022.12.05
       do ip = 0,np-1
         write(*,'(a,i3)') " ip",ip
         do i=1,nf_tot_ip
           write(*,'(a,i3,a,f9.3,a,i3)') " i",i," freq",f_tot_ip(i,ip)," Hz       index",if_tot_ip(i,ip)
         end do
       end do
   end if

!#[3]##================================================== node for ACTIVE and MT frequency
  !#[3-1]## nf_act_ip, nf_mt_ip
   allocate(nf_act_ip(0:np-1),nf_mt_ip(0:np-1))
   nf_act_ip(:) = 0  ! : # of freq for ip node (ACTIVE)
   nf_mt_ip(:)  = 0  ! : # of freq for ip node (MT)
   do ip=0,np-1
     do i=1,nf_tot_ip   ! common over ip
       ii = if_tot_ip(i,ip) ! global frequency index (depends on ip)
       if ( ii .ne. 0 ) then
        if ( ind_f_tot(1,ii) > 0 ) nf_act_ip(ip) = nf_act_ip(ip) + 1
        if ( ind_f_tot(2,ii) > 0 ) nf_mt_ip(ip)  = nf_mt_ip(ip)  + 1
       end if
     end do
    end do
    if ( ipp .eq. 0) then ! 2022.12.05
     write(*,'(a)') "  ip  | nf_act_ip |  nf_mt_ip"
     do ip = 0,np-1
       write(*,'(i4,2i10)') ip,nf_act_ip(ip),nf_mt_ip(ip)
     end do
    end if

  !#[3-2]## i_act_ip, i_act_mt_ip
  !      frequency id in ip node -> ACTIVE/MT frequency id in ip node 
   allocate(i_act_ip(nf_tot_ip,0:np-1))  ! freq index in ip node for ACTIVE 
   allocate(i_mt_ip( nf_tot_ip,0:np-1))
   i_act_ip(:,:) = 0 
   i_mt_ip(:,:)  = 0 
   do ip=0,np-1
     icount_act = 0 ;  icount_mt  = 0
     do i=1,nf_tot_ip   ! common in ip
       ii = if_tot_ip(i,ip) ! ii: global freq index 
       if ( ii .ne. 0 ) then
       if ( ind_f_tot(1,ii) > 0 ) then ! ACTIVE
         icount_act     = icount_act+1
         i_act_ip(i,ip) = icount_act
       end if
       if ( ind_f_tot(2,ii) > 0 ) then
         icount_mt      = icount_mt + 1
         i_mt_ip(i,ip)  = icount_mt
       end if
       end if
     end do
    end do

  !#[3-3]## if_act_ip, if_mt_ip, to be used in [3-4] and [3-5]
   allocate(if_act_ip(0:nf_tot_ip,0:np-1)) ! nf_tot_ip >= nf_act_ip(:), nf_mt_ip(:)
   allocate(if_mt_ip( 0:nf_tot_ip,0:np-1 ))
   if_act_ip = 0 ; if_mt_ip = 0
   do ip=0,np-1
     do i=1,nf_tot_ip
       ii = if_tot_ip(i,ip) ! ii: grobal frequency index
       ij = i_act_ip(i,ip)  ! ij th freq for ACTIVE in ip node
       if (  ij > 0 )  if_act_ip(ij,ip) = ind_f_tot(1,ii) ! freq index for whole ACTIVE freq
       jj = i_mt_ip(i,ip)
       if (  jj > 0 )  if_mt_ip(jj,ip)  = ind_f_tot(2,ii) ! freq index for whole MT
     end do
   end do
    if ( ipp .eq. 0 ) then ! 2022.12.05
      write(*,'(a)') "  ip  i | i_act_ip (whole)  | i_mt_ip (whole) "
     do ip = 0,np-1
      do i=1,nf_tot_ip   ! common in ip
       ij = i_act_ip(i,ip)
       jj = i_mt_ip(i,ip)
       write(*,'(1x,2i3,2(6x,i5,a,i2,a))') ip,i,ij," (",if_act_ip(ij,ip),")",jj,'(', if_mt_ip(jj,ip) ,")"
      end do
     end do
    end if

  !#[3-4]## set ip_from_act, if_act_g2l_act
    allocate( ip_from_act(nf_act),if_g2l_act(nf_act))
    do i=1,nf_act
     ii = if_act2tot(i) ! global freq index
     ip = ip_from_tot(ii)
     ip_from_act(i) = ip
     do j=1,nf_act_ip(ip)
      if ( if_act_ip(j,ip) .eq. i ) if_g2l_act(i) = j
     end do
    end do
  
  !#[3-5]##
    allocate( ip_from_mt(nf_mt),if_g2l_mt(nf_mt))
    do i=1,nf_mt
     ii = if_mt2tot(i) ! global freq index
     ip = ip_from_tot(ii)
     ip_from_mt(i) = ip
     do j=1,nf_mt_ip(ip)
      if ( if_mt_ip(j,ip) .eq. i ) if_g2l_mt(i) = j
     end do
    end do


!#[3]## set output
  g_freq_joint%np           = np
  !# total frequency
  g_freq_joint%nfreq_tot    = nf_tot
  g_freq_joint%nfreq_tot_ip = nf_tot_ip
  g_freq_joint%freq_tot     = f_tot      ! (nf_tot)          allocate and fill 2022.10.20
  g_freq_joint%freq_tot_ip  = f_tot_ip   !(nf_tot_ip,0:np-1) allocate and fill 2022.10.20
  g_freq_joint%ifreq_tot_ip = if_tot_ip  !(nf_tot_ip,0:np-1) allocate and fill 2022.10.20

  !# active frequency
  g_freq_joint%nfreq_act    = nf_act
  g_freq_joint%nfreq_act_ip = nf_act_ip   ! allocate and fill
  g_freq_joint%i_act_ip     = i_act_ip    ! (nf_tot_ip) allocate and fill
  g_freq_joint%ip_from_act  = ip_from_act ! (nf_act)  allocate and fill
  g_freq_joint%if_g2l_act   = if_g2l_act  ! (nf_act)  allocate and fill

  !# mt frequency
  g_freq_joint%nfreq_mt    = nf_mt
  g_freq_joint%nfreq_mt_ip = nf_mt_ip   ! (0:np-1)    allocate and fill
  g_freq_joint%i_mt_ip     = i_mt_ip    ! (nf_tot_ip) allocate and fill
  g_freq_joint%ip_from_mt  = ip_from_mt ! (nf_mt)     allocate and fill
  g_freq_joint%if_g2l_mt   = if_g2l_mt  ! (nf_mt)     allocate and fill

!  write(*,*) "np",np
!  write(*,*) "nf_tot",nf_tot
!  write(*,*) "nf_tot_ip",nf_tot_ip
!  write(*,*) "f_tot",f_tot
!  write(*,*) "f_tot_ip",f_tot_ip
!  do ip=0,np-1
!    write(*,*) "ip",ip
!    write(*,*) "if_tot_ip(:,ip)",if_tot_ip(:,ip)
!  end do
!  do i=1,nf_act
!       write(*,*) "i",i,"ip_from",ip_from_act(i)
!       write(*,*) "i",i,"if_g2l_mt",if_g2l_act(i)
!  end do

return
end subroutine

!##################################################################################################
! 2022.12.05 ip is added in arguments
subroutine mkfreqtot(nfreq_act,nfreq_mt,freq_act,freq_mt,nfreq_tot,freq_tot,nfreq_tot_n,index_freq,ip)
    implicit none
    integer(4), intent(in)  :: ip ! 2022.12.05
    integer(4), intent(in)  :: nfreq_tot, nfreq_act, nfreq_mt
    real(8),    intent(in)  :: freq_act(nfreq_act),freq_mt(nfreq_mt)
    integer(4), intent(out) :: nfreq_tot_n
    real(8),    intent(out) :: freq_tot(nfreq_tot)
    integer(4), intent(out) :: index_freq(2,nfreq_tot) 
    
    real(8),   allocatable,dimension(:)   :: freq_work,freq_work2
    integer(4),allocatable,dimension(:,:) :: index,index_work,index2
    integer(4),allocatable,dimension(:)   :: order
    integer(4)                            :: i,ii,ij

    !#[1]##
    allocate(freq_work(nfreq_tot),freq_work2(nfreq_tot))
    allocate(index(2,  nfreq_tot),index2(2,  nfreq_tot))
    allocate(index_work(2,  nfreq_tot))
    allocate(order(nfreq_tot))
    index(:,:) = 0

    !#[2]## freq
    do i=1,nfreq_act
      freq_work(i) = freq_act(i)
      index(1,i) = i           ! ACT freq index
    end do
    do i=1,nfreq_mt
      freq_work(nfreq_act+i) = freq_mt(i)
      index(2,nfreq_act+i) = i !  MT freq index
    end do

    !#[3]# sort
    do i=1,nfreq_tot ! set order for nfreq_act + nfreq_mt
     order(i) = i
    end do
    call SORT_INDEX(nfreq_tot,order,freq_work)
    do i=1,nfreq_tot
        index_work(:,i)=index(:,order(i))
    end do
    if ( ip == 0 ) then ! output sorted frequencies 2022.12.05
      do i=1,nfreq_tot
        write(*,'(a,i5,a,f8.2,a,2i5)') " order",order(i)," freq",freq_work(i)," index",index_work(:,i)
      end do
    end if             ! 2022.12.05

    !#[4]## merge the same frequencies
       ii=1 ; ij =1
       freq_work2 = 0.d0
       do while ( ij < nfreq_tot )
!        write(*,*) "ii",ii,"ij",ij
        freq_work2(ii) = freq_work(ij)
        nfreq_tot_n = ii
        if ( abs(freq_work(ij) - freq_work(ij+1))/freq_work(ij) < 1.d-3 ) then ! difference is very small
         index2(1:2,ii) = index_work(1:2,ij) + index_work(1:2,ij+1)
         ij = ij + 2 ! skip the next
         ii = ii + 1
        else ! next is significanly larger 
         index2(1:2,ii) = index_work(1:2,ij)
         ij = ij + 1 ! set the next
         ii = ii + 1
        end if
       end do
       if ( ij .eq. nfreq_tot) then
!        write(*,*) "ii",ii,"ij",ij
        freq_work2(ii) = freq_work(ij)
        index2(1:2,ii) = index_work(1:2,ij)
        nfreq_tot_n = ii
       end if

       index_freq = index2
       freq_tot   = freq_work2

    return
    end
!##########################################################################################
!# Check sorted and merged frequencies
subroutine checkfreq_tot( nf_act,nf_mt,nf_tmp,nf_tot,freq_tot,index_freq_tot,ip_from_tot,ip)
implicit none
integer(4),intent(in) :: nf_act,nf_mt,nf_tmp,nf_tot,ip
integer(4),intent(in) :: index_freq_tot(2,nf_tot)
integer(4),intent(in) :: ip_from_tot(nf_tot)
real(8),   intent(in) :: freq_tot(nf_tot)
integer(4) :: i

  if ( ip .eq. 0 ) then
      write(*,*) "nfreq_act",nf_act,"nfreq_mt",nf_mt
      write(*,*) "nfreq_tot:",nf_tmp,"->",nf_tot
      write(*,'(a)') " Freq index |  Hz |  Index for  ACT / MT  (MPI node) --"
      do i=1,nf_tot
       write(*,'(a,i5,f8.2,a,2i5,i10)') "     ",i,freq_tot(i)," freq index",index_freq_tot(1:2,i),ip_from_tot(i)
      end do
  end if  
 return

end

end module freq_mpi_joint
