!=====================================================================
! AUTHORS:
!  Gary Egbert & Lana Erofeeva
!  College of Atmospheric and Oceanic Sciences
!  104 COAS Admin. Bldg.
!  Oregon State University
!  Corvallis, OR 97331-5503
!  
!  E-mail:  egbert@coas.oregonstate.edu                                      
!  Fax:     (541) 737-2064
!  Ph.:     (541) 737-2947                                        
!  https://www.tpxo.net/
!
! COPYRIGHT: OREGON STATE UNIVERSITY, 2010
! (see the file COPYRIGHT for lisence agreement)
!=====================================================================
      program extract_local_model
!
!cc   reads OTIS standard binary single constituent atlas files
!cc   (elevations, transports and grid),listed in  
!ccc  DATA/Model_tpxo9_atlas file, reads "limits" file
!cc   outputs multiple constituents "local" model
!cc   in standard OTIS binary format
! 
      implicit none
      include 'constit.h'
      complex, allocatable:: z1(:,:,:),u1(:,:,:),v1(:,:,:)
      real, allocatable:: depth(:,:),depth1(:,:)
      real th_lim(2),ph_lim(2),th_lim1(2),ph_lim1(2)
      real*8 dth,dph
      integer, allocatable:: mz(:,:),mz1(:,:)
!
      character*4 c_id(ncmx),c_id1(ncmx)
      character*80 modname,outname,hname,uname,gname,fname
      character*80 hname1,uname1,gname1
      character*80 rmCom,ctmp
      integer nc,n,m,nc1,n1,m1,n0,m0,i1,i2,j1,j2,ic
      integer i,ii,j,k,l,nca,funit(ncmx)
      complex rc
! read setup
      read(*,'(a80)')ctmp
      modname=rmCom(ctmp)
      read(*,'(a80)')ctmp
      outname=rmCom(ctmp)
      read(*,*)th_lim1
      read(*,*)ph_lim1
! read atlas file names
      call rd_mod_file(modname,hname,uname,gname,fname,nca,c_id)
! read local model file names
      call rd_mod_file(outname,hname1,uname1,gname1,fname,nc1,c_id1)
!
      call rd_mod_header(modname,'z',n,m,th_lim,ph_lim,nc,c_id,ctmp)
      write(*,*)'Atlas constituents:',c_id(1:nc)
!
      dth=(th_lim(2)-th_lim(1))/m
      dph=(ph_lim(2)-ph_lim(1))/n
      j1=(th_lim1(1)-th_lim(1))/dth
      j2=(th_lim1(2)-th_lim(1))/dth
      m1=j2-j1+1
      th_lim1(1)=th_lim(1)+dph*j1
      th_lim1(2)=th_lim1(1)+dph*m1
      j2=j1+m1-1
      if(ph_lim1(1).gt.0)then
       i1=(ph_lim1(1)-ph_lim(1))/dph
       i2=(ph_lim1(2)-ph_lim(1))/dph
       n1=i2-i1+1    
       ph_lim1(1)=ph_lim(1)+dph*i1
       ph_lim1(2)=ph_lim1(1)+dph*n1
       i2=i1+n1-1
      else
       i1=(ph_lim(2)+ph_lim1(1))/dph
       if(ph_lim1(2).lt.0)then
        i2=(ph_lim(2)+ph_lim1(2))/dph
       else
        i2=(ph_lim1(2)-ph_lim(1))/dph
       endif
       if(i1.lt.i2)then
        n1=i2-i1+1    
        ph_lim1(1)=ph_lim(1)+dph*i1-360.
        ph_lim1(2)=ph_lim1(1)+dph*n1-360.
        i2=i1+n1-1
       else ! pass through lon=0
        i2=i2+1
        n1=(n-i1+1)+i2
        ph_lim1(1)=ph_lim(2)-dph*(n1-i2)-360.
        ph_lim1(2)=ph_lim1(1)+dph*n1
       endif
      endif
      !write(*,*)i1,i2,n1,m1
      !write(*,*)ph_lim1,th_lim1
      !stop
!
      allocate(z1(nc,n1,m1),depth1(n1,m1),mz1(n1,m1))
      allocate(u1(nc,n1,m1),v1(nc,n1,m1))
!
       allocate(depth(n,m),mz(n,m))
       open(unit=1,file=trim(gname),status='old',&
            form='unformatted',err=4)
       read(1)n0,m0 ! ignore the rest of record
       if(n0.ne.n.or.m0.ne.m)then
         write(*,*)'Wrong grid in ',trim(gname)
         write(*,*)'Grid size and model size are different!'
         stop
       endif
       read(1) ! pass iobc
       read(1) depth
       read(1) mz
       close(1)
! write grid
       if(ph_lim1(1)*ph_lim1(2).gt.0)then
        depth1=depth(i1:i2,j1:j2)
        mz1=mz(i1:i2,j1:j2)
       else
        depth1(1:n-i1+1,:)=depth(i1:n,j1:j2)
        depth1(n-i1+2:n1,:)=depth(1:i2,j1:j2)
        mz1(1:n-i1+1,:)=mz(i1:n,j1:j2)
        mz1(n-i1+2:n1,:)=mz(1:i2,j1:j2)
       endif
!
       open(unit=1,file=trim(gname1),status='unknown',&
            form='unformatted',err=4)
       write(1)n1,m1,th_lim1,ph_lim1,12.,0
       write(1)0
       write(1)depth1
       write(1)mz1
       close(1)
       write(*,*)'Grid done:',trim(gname1)
       deallocate(depth1,mz1,depth,mz)
!
       if(i2.lt.i1)i2=i2+n
       fname=hname
       write(*,*)'Opening atlas h-files:'
       do ic=1,nc
        if(ic.gt.1)then
          k=index(fname,trim(c_id(ic-1)))
          l=len(trim(c_id(ic-1)))
          fname=fname(1:k-1)//trim(c_id(ic))//fname(k+l:80)
        endif
        write(*,*)trim(fname)
        funit(ic)=100+ic
        open(unit=funit(ic),file=trim(fname),status='old', &
             form='unformatted',recl=4,access='direct')
       enddo
       write(*,*)'done'
!     
      write(*,'(a,$)')' Reading atlas - z...'
      do ic=1,nc
       write(*,'(a,$)')c_id(ic)
       do i=i1,i2
        ii=i
        if(i.gt.n)ii=i-n
        if(i.lt.1)ii=i+n
        do j=j1,j2
         call rd_mod_value_da('z',funit(ic),ii,j,n0,m0,1,rc)
         z1(ic,i-i1+1,j-j1+1)=rc
        enddo
       enddo
       close(funit(ic))
      enddo
      write(*,*)'done'
!
       fname=uname
       write(*,*)'Opening atlas uv-files:'
       do ic=1,nc
        if(ic.gt.1)then
          k=index(fname,trim(c_id(ic-1)))
          l=len(trim(c_id(ic-1)))
          fname=fname(1:k-1)//trim(c_id(ic))//fname(k+l:80)
        endif
        write(*,*)trim(fname)
        funit(ic)=100+ic
        open(unit=funit(ic),file=trim(fname),status='old', &
             form='unformatted',recl=4,access='direct')
       enddo
       write(*,*)'done'

      write(*,'(a,$)')' Reading atlas - u, v ...'
      do ic=1,nc
       write(*,'(a,$)')c_id(ic)
       do i=i1,i2
        ii=i
        if(i.gt.n)ii=i-n
        if(i.lt.1)ii=i+n
        do j=j1,j2
         call rd_mod_value_da('u',funit(ic),ii,j,n0,m0,1,rc)
         u1(ic,i-i1+1,j-j1+1)=rc
         call rd_mod_value_da('v',funit(ic),ii,j,n0,m0,1,rc)
         v1(ic,i-i1+1,j-j1+1)=rc
        enddo
       enddo
       close(funit(ic))
      enddo
      write(*,*)'done'
!
      write(*,'(a,$)')' Writing local model z...'
      call write_z(hname1,1,z1,n1,m1,nc,th_lim1,ph_lim1,c_id)
      write(*,'(a,$)')'u,v...'
      call write_uv(uname1,1,u1,v1,n1,m1,nc,th_lim1,ph_lim1,c_id)
      write(*,*)'done'
      write(*,*)'OUTPUTS:'
      write(*,*)hname1
      write(*,*)uname1
      write(*,*)gname1
!
      deallocate(z1,u1,v1)
4     stop
      end
