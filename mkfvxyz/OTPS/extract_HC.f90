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
      program extract_HC_da
!cc   LANA, 2020 remake (direct access to binary model files)
!cc   needs min RAM and time to extract/predict by directly accessing
!cc   4 model nodes corresponding to given lat/lon cell 
!cc   reads OTIS standard binary complex model file
!cc   (elevations OR transports), reads a list of locations 
!cc   and outputs ASCII file with the complex amplitudes/amp,phases of
!cc   elevations/transports/currents interpolated at the locations
! 
      implicit none
      include 'constit.h'
      complex, allocatable:: z1(:),dtmp(:,:)
      complex, allocatable:: zl1(:)
      complex d1
      real, allocatable:: lat(:),lon(:),depth(:,:),x(:),y(:),lon0(:)
      real, allocatable:: phase(:)
      real latp,lonp,xt,yt
      real th_lim(2),ph_lim(2),dum,lth_lim(2),lph_lim(2)
      integer, allocatable:: cind(:),lcind(:),mz(:,:)
!
      character*4 c_id(ncmx),c_id_mod(ncmx),lc_id(ncmx)
      character*80 modname,lltname,outname,gname,ctmp,lname
      character*80 hname,uname,fname
      character*2000 fmt
      character*80 rmCom
      character*1 zuv,c1
      character*80 xy_ll_sub
      logical APRI,geo,interp_micon,ll_km
      integer ncon,nc,n,m,ndat,k,ierr,ierr1,ic,n0,m0
      integer ncl,nl,ml,ibl,funit(ncmx),funitl(ncmx),nca,l
!
      ll_km=.false.
      ibl=0
      funitl(:)=100
      funit(:)=101
      lname='DATA/load_file'
      call rd_inp(modname,lltname,zuv,c_id,ncon,APRI,geo, &
                  outname,interp_micon)
      call rd_mod_file(modname,hname,uname,gname,xy_ll_sub,nca,c_id_mod)
      write(*,*)
      write(*,*)'Lat/Lon file:',trim(lltname)
      if(ncon.gt.0)write(*,*)'Constituents to include: ',c_id(1:ncon)
      if(geo.and.zuv.eq.'z')then
         write(*,*)'Extract GEOCENTRIC tide HC'
      else
         write(*,*)'Extract OCEAN tide HC'
      endif
!
      open(unit=11,file=outname,status='unknown')
!
      write(*,*)
      call rd_mod_header(modname,zuv,n,m,th_lim,ph_lim,nc,c_id_mod,&
                          xy_ll_sub)
      write(*,*)'Model:        ',trim(modname(12:80))
      write(11,*)'Model:        ',trim(modname(12:80))
      if(trim(xy_ll_sub).eq.'')then
       write(*,*)'Lat limits:   ',th_lim
       write(*,*)'Lon limits:   ',ph_lim
      else
       ll_km=.true.
       if(trim(xy_ll_sub).ne.'xy_ll_N'.and.&
          trim(xy_ll_sub).ne.'xy_ll_S'.and.&
          trim(xy_ll_sub).ne.'xy_ll_CATs')then
        write(*,*)'No converting function ', trim(xy_ll_sub),&
                  'in the OTPS'
        stop 
       endif
       if(trim(xy_ll_sub).eq.'xy_ll_N')then
        call xy_ll_N(ph_lim(1),th_lim(1),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon lower left corner:',yt,xt
        call xy_ll_N(ph_lim(1),th_lim(2),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon upper left corner:',yt,xt
        call xy_ll_N(ph_lim(2),th_lim(1),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon lower right corner:',yt,xt
        call xy_ll_N(ph_lim(2),th_lim(2),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon upper right corner:',yt,xt
       elseif(trim(xy_ll_sub).eq.'xy_ll_S')then
        call xy_ll_S(ph_lim(1),th_lim(1),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon lower left corner:',yt,xt
        call xy_ll_S(ph_lim(1),th_lim(2),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon upper left corner:',yt,xt
        call xy_ll_S(ph_lim(2),th_lim(1),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon lower right corner:',yt,xt
        call xy_ll_S(ph_lim(2),th_lim(2),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon upper right corner:',yt,xt
       elseif(trim(xy_ll_sub).eq.'xy_ll_CATs')then
        call xy_ll_CATs(ph_lim(1),th_lim(1),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon lower left corner:',yt,xt
        call xy_ll_CATs(ph_lim(1),th_lim(2),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon upper left corner:',yt,xt
        call xy_ll_CATs(ph_lim(2),th_lim(1),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon lower right corner:',yt,xt
        call xy_ll_CATs(ph_lim(2),th_lim(2),xt,yt)
        if(xt.gt.180)xt=xt-360
        write(*,*)'Lat,Lon upper right corner:',yt,xt
       endif
      endif
      write(*,*)'Constituents: ',c_id_mod(1:nc)
      if(trim(xy_ll_sub).ne.'')then
           write(*,*)'Model is on uniform grid in km'
           write(*,*)'Function to convert x,y to lat,lon:',&
                      trim(xy_ll_sub)
      endif 
!
      if(zuv.eq.'z')write(*,*)'Output elevations (m)'
      if(zuv.eq.'U')write(*,*)'Output WE transport (m^2/s)'
      if(zuv.eq.'u')write(*,*)'Output WE velocity  (cm/s)'
      if(zuv.eq.'V')write(*,*)'Output SN transport (m^2/s)'
      if(zuv.eq.'v')write(*,*)'Output SN velocity  (cm/s)'
!
      if(zuv.eq.'z')write(11,*)'Elevations (m)'
      if(zuv.eq.'U')write(11,*)'WE transport (m^2/s)'
      if(zuv.eq.'u')write(11,*)'WE velocity  (cm/s)'
      if(zuv.eq.'V')write(11,*)'SN transport (m^2/s)'
      if(zuv.eq.'v')write(11,*)'SN velocity  (cm/s)'
!
      if(ncon.eq.0)then
       ibl=1
       ncon=nc
       c_id=c_id_mod
       write(*,*)'Constituents to include: ',c_id(1:ncon)
      endif
!
      allocate(cind(ncon))
      call def_con_ind(c_id,ncon,c_id_mod,nc,cind)
! 
      ndat=0
      open(unit=1,file=lltname,status='old',err=1)
3     read(1,*,end=2)dum,dum
      ndat=ndat+1
      go to 3
2     rewind(1)
      allocate(lat(ndat),lon(ndat),lon0(ndat))
      if(trim(xy_ll_sub).ne.'')allocate(x(ndat),y(ndat))
      do k=1,ndat
       read(1,*)lat(k),lon(k)
        if(trim(xy_ll_sub).eq.'xy_ll_N')then
         call ll_xy_N(lon(k),lat(k),x(k),y(k))
         !write(*,*)lat(k),lon(k),x(k),y(k)
        elseif(trim(xy_ll_sub).eq.'xy_ll_S')then
         call ll_xy_S(lon(k),lat(k),x(k),y(k))
        elseif(trim(xy_ll_sub).eq.'xy_ll_CATs')then
         call ll_xy_CATs(lon(k),lat(k),x(k),y(k))
        endif
        lon0(k)=lon(k)
        if(trim(xy_ll_sub).eq.'')then ! check on lon convention
         if(lon(k).gt.ph_lim(2))lon(k)=lon(k)-360
         if(lon(k).lt.ph_lim(1))lon(k)=lon(k)+360
        endif
      enddo ! k
      close(1)
!
      allocate(z1(ncon))
!
      if(zuv.eq.'z'.and.geo)then
       write(*,'(a,$)')'Reading load correction header...'
       call rd_mod_header1(lname,nl,ml,ncl,lth_lim,lph_lim,lc_id)
       allocate(lcind(ncon))
       call def_con_ind(c_id,ncon,lc_id,ncl,lcind)
       allocate(zl1(ncon))
       write(*,*)'done'
       open(unit=funitl(1),file=trim(lname),status='old', &
           form='unformatted',recl=4,access='direct')
      endif
      if(zuv.eq.'u'.or.zuv.eq.'v')then
       write(*,'(a,$)')'Reading grid file...'
       allocate(depth(n,m),dtmp(n,m),mz(n,m))
! currents case: need to read grid
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
       dtmp=depth
       deallocate(depth)
       write(*,*)'done'
      endif
      allocate(phase(ncon))
! output format
      fmt='(f9.3,x,f9.3,x'
      c1=','
      do ic=1,ncon
       if(ic.eq.ncon)c1=')'
       if(APRI)then
        fmt=trim(fmt)//'f8.3,x,f8.1,x'//c1
       else
        fmt=trim(fmt)//'f8.3,x,f8.3,x'//c1
       endif
      enddo
!
      if(APRI)then
        write(11,*)'  Lat     Lon       ',&
                    (trim(c_id(ic)),'_amp  ',&
                     trim(c_id(ic)),'_ph   ',ic=1,ncon)
      else
        write(11,*)'  Lat       Lon     ', &
                   (trim(c_id(ic)),'_Re   ', &
                    trim(c_id(ic)),'_Im   ',ic=1,ncon)
      endif
      c1=zuv ! since interp change zuv (U->u, V->v)
      latp=0.
      lonp=0.
!
      fname=hname
      if(zuv.ne.'z')fname=uname
      if(nca.eq.0)then
       open(unit=funit(1),file=trim(fname),status='old', &
            form='unformatted',recl=4,access='direct')
      else
       write(*,'(a,$)')'Opening atlas files:'
       do ic=1,ncon
        if(ic.gt.1)then
          k=index(fname,trim(c_id(ic-1)))
          l=len(trim(c_id(ic-1)))
          fname=fname(1:k-1)//trim(c_id(ic))//fname(k+l:80)
        endif
        write(*,'(a,$)')c_id(ic)
        funit(ic)=100+ic
        open(unit=funit(ic),file=trim(fname),status='old', &
             form='unformatted',recl=4,access='direct')
       enddo
       write(*,*)'done'
      endif  
!
      do k=1,ndat
       if(lat(k).ne.latp.or.lon(k).ne.lonp)then
        if(trim(xy_ll_sub).eq.'')then
         call interp_da(funit,n,m,th_lim,ph_lim, &
                        lat(k),lon(k),z1,ncon,cind,ierr,c1,nca)
        else
         z1(1)=-1
         call interp_da(funit,n,m,th_lim,ph_lim, &
                        y(k),x(k),z1,ncon,cind,ierr,c1,nca)
        endif
        if(ierr.eq.0)then
         if(zuv.eq.'u'.or.zuv.eq.'v')then
          if(trim(xy_ll_sub).eq.'')then
           call interp(dtmp,1,n,m,mz,th_lim,ph_lim, &
                       lat(k),lon(k),d1,ierr1,'z')
          else
           d1=-1
           call interp(dtmp,1,n,m,mz,th_lim,ph_lim, &
                       y(k),x(k),d1,ierr1,'z')
          endif
          z1=z1/real(d1)*100. ! currents cm/s
         elseif(zuv.eq.'z'.and.geo)then       
           call interp_da(funitl,nl,ml,lth_lim,lph_lim, &
                    lat(k),lon(k),zl1,ncon,lcind,ierr1,'z',0)
           z1=z1+zl1     ! apply load correction to get geocentric tide
         endif  
        endif
        if(ierr.eq.0)then
         if(APRI)then
          phase=atan2(-imag(z1),real(z1))*180/3.141593
          write(11,fmt)lat(k),lon0(k), &
                (abs(z1(ic)),phase(ic),ic=1,ncon)
         else
           write(11,fmt)lat(k),lon0(k), &
                    (real(z1(ic)),imag(z1(ic)),ic=1,ncon)
         endif
        else
          write(11,'(1x,f8.3,x,f8.3,a70)')lat(k),lon(k), &
       '************* Site is out of model grid OR land ***************'
        endif
        latp=lat(k)
        lonp=lon(k)
       endif  
      enddo
      deallocate(z1,cind,lat,lon,lon0)
      if(zuv.eq.'u'.or.zuv.eq.'v')deallocate(dtmp,mz)
      if(trim(xy_ll_sub).ne.'')deallocate(x,y)
      if(ibl.eq.1.and.geo)then
         ncon=0
         deallocate(zl1,lcind)
      endif
      close(11)
      close(funit(1))
      if(nca.ne.0.and.ncon.gt.1)then
       do ic=2,ncon
        close(funit(ic))
       enddo
      endif
      if(geo)close(funitl(1))
      write(*,*)'Results are in ',trim(outname)
      if(geo.and.ibl.eq.0)deallocate(zl1,lcind)
      stop
1     write(*,*)'Lat Lon file ',trim(lltname),' not found'
      write(*,*)'Check setup file, line 2.'
      stop
4     write(*,*)'Grid file ',trim(gname),' not found'
      write(*,*)'Check file ',trim(modname)
      stop
      end
