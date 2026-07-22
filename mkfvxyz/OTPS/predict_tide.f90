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
      program predict_tide_da
!cc 
!cc   LANA, 2020 remake (direct access to binary model files)
!cc   needs min RAM and time to extract/predict by directly accessing
!cc   4 model nodes corresponding to given lat/lon cell 
!
!cc   modified March 2011 to optimize for obtaining
!cc   time series at open boundaries
!cc
!cc   reads OTIS standard binary complex model file
!cc   (elevations OR transports), reads a list of locations,
!cc   reads list of times  
!cc   and outputs ASCII file with the tidal predictions of tidal 
!cc   elevations/transports/currents at the locations and times
!cc   
      implicit none
      include 'constit.h'
      complex, allocatable:: z1(:),dtmp(:,:)
      complex, allocatable:: zl1(:)
      complex, allocatable:: u1(:),v1(:)
      complex d1
      real, allocatable:: lat(:),lon(:),depth(:,:),x(:),y(:),&
                          lon0(:),zpred(:),upred(:),vpred(:)
      real xt,yt 
      real*8, allocatable:: time_mjd(:)
      real th_lim(2),ph_lim(2),dum,lth_lim(2),lph_lim(2)
      integer, allocatable:: cind(:),lcind(:),ccind(:),mz(:,:)
!
      character*4 c_id(ncmx),c_id_mod(ncmx),lc_id(ncmx),tcon(ncmx)
      character*80 modname,lltname,outname,ctmp,lname
      character*80 hname,uname,gname,fname
      character*2000 fmt
      character*80 rmCom,tfname
      character*1 zuv,c1,c2
      character*80 xy_ll_sub,arg
      character*10 cdate
      character*8 ctime
      character*10 deblank 
      logical APRI,geo,interp_micon,ll_km
      integer ncon,nc,n,m,ndat,i,j,k,k1,ierr,ierr1,ic,n0,m0,it
      integer ncl,nl,ml,nmod,imod,ibl,ntime,idum,mjd,julian
      integer yyyy1,mm1,dd1,iargc,narg
      integer funit(ncmx),funitl(ncmx),nca,l
      integer, allocatable:: yyyy(:),mm(:),dd(:),hh(:), &
                             mi(:),ss(:)
!
      ll_km=.false.
      narg=iargc()
      if(narg.gt.0)then
        call getarg(1,arg)
        read(arg(3:80),'(a80)')tfname
        ntime=0
        open(unit=1,file=tfname,status='old',err=16)
18      read(1,*,end=17)idum,idum,idum,idum,idum,idum
        ntime=ntime+1
        go to 18
17      rewind(1)
        ntime=ntime-1
        allocate(yyyy(ntime),mm(ntime),dd(ntime), &
               hh(ntime),mi(ntime),ss(ntime),time_mjd(ntime))
        call read_time(ntime,1,2,yyyy,mm,dd,hh,mi,ss,time_mjd)
        close(1)
      endif       
!
      nmod=1
      ibl=0
      funitl(:)=100
      funit(:)=101 
      lname='DATA/load_file'
      call rd_inp(modname,lltname,zuv,c_id,ncon,APRI,geo, &
                  outname,interp_micon)
      call rd_mod_file(modname,hname,uname,gname,xy_ll_sub,nca,c_id_mod)
      write(*,*)
      write(*,*)'Lat/Lon/Time file:',trim(lltname)
      if(ncon.gt.0)write(*,*)'Constituents to include: ',c_id(1:ncon)
      if(zuv.eq.'z')then
       if(geo)then
         write(*,*)'Predict GEOCENTRIC tide'
       else
         write(*,*)'Predict OCEAN tide'
       endif
      endif
      if(interp_micon)write(*,*)'Interpolate minor constituents'
!
      if(narg.eq.0)then ! default usage
! read times
       ntime=0
       open(unit=1,file=lltname,status='old',err=6)
8      read(1,*,end=7,err=7)dum,dum,idum,idum,idum,idum,idum,idum
       ntime=ntime+1
       go to 8
7      rewind(1)
       allocate(yyyy(ntime),mm(ntime),dd(ntime), &
                hh(ntime),mi(ntime),ss(ntime),time_mjd(ntime))
       call read_time(ntime,1,1,yyyy,mm,dd,hh,mi,ss,time_mjd)
       close(1)
      endif
!
      open(unit=11,file=outname,status='unknown')

      call rd_mod_header(modname,zuv,n,m,th_lim,ph_lim,nc,c_id_mod,&
                          xy_ll_sub)
      write(*,*)'Model:        ',trim(modname(12:80))
      write(11,'(60a1)')('-',i=1,60)
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
                  ' in the OTPS'
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
       else
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
      if(zuv.eq.'z')then
        write(*,*)'Predict elevations (m)'
      else
        write(*,*)'Predict transport (m^2/s) and currents (cm/s)'
      endif
!
      k1=1
      tcon=''  
      if(ncon.eq.0)then
       ibl=1
       ncon=nc
       c_id=c_id_mod
      else
! check if all required constituents are in the model
       do ic=1,ncon
        do k=1,nc
         if(c_id(ic).eq.c_id_mod(k))then
          tcon(k1)=c_id(ic)
          k1=k1+1
          go to 14
         endif
        enddo
        write(*,*)'Constituent ',c_id(ic), ' is NOT in the model'
14      continue
       enddo 
       ncon=k1-1
       c_id=tcon
      endif
 
      write(*,*)'Constituents to include: ',c_id(1:ncon)
      write(11,*)'Constituents included: ',c_id(1:ncon)
!
      allocate(cind(ncon),ccind(ncon))
      call def_con_ind(c_id,ncon,c_id_mod,nc,cind)
! find corresponding indices in constit.h
      call def_cid(ncon,c_id,ccind)
!
      if(narg.eq.0)then
       ndat=ntime
      else
       ndat=0
       open(unit=1,file=lltname,status='old',err=6)
 20    read(1,*,end=19)dum,dum
       ndat=ndat+1
       go to 20
 19    close(1)
      endif
!
      open(unit=1,file=lltname,status='old',err=6)
      allocate(lat(ndat),lon(ndat),lon0(ndat))
      if(zuv.eq.'z')then
       allocate(zpred(ndat))
      else
       allocate(upred(ndat),vpred(ndat))
      endif
      if(trim(xy_ll_sub).ne.'')allocate(x(ndat),y(ndat))
      do k=1,ndat
       read(1,*)lat(k),lon(k) ! ignore rest of record
       if(trim(xy_ll_sub).eq.'xy_ll_N')then
         call ll_xy_N(lon(k),lat(k),x(k),y(k))
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
      enddo
      close(1)
!
      if(zuv.eq.'z')then      
       allocate(z1(ncon))
       if(geo) then
        call rd_mod_header1(lname,nl,ml,ncl,lth_lim,lph_lim,lc_id)
        allocate(lcind(ncon))
        call def_con_ind(c_id,ncon,lc_id,ncl,lcind)
        allocate(zl1(ncon))
        open(unit=funitl(1),file=trim(lname),status='old', &
             form='unformatted',recl=4,access='direct')
       endif
      else
       allocate(u1(ncon),v1(ncon))
      endif
!
       allocate(depth(n,m),dtmp(n,m),mz(n,m))
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
!
      ctmp='    Lat       Lon        mm.dd.yyyy hh:mm:ss'
      write(11,*)'' 
      if(zuv.eq.'z')then
        write(11,*)trim(ctmp),'     z(m)   Depth(m)'
      else
        write(11,*)trim(ctmp),'   U(m^2/s)  V(m^2/s)',&
                              '   u(cm/s)   v(cm/s) Depth(m)'
      endif
      write(11,*)''
!
      c1=zuv ! since interp change zuv (U->u, V->v)
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
        xt=lon(k)
        yt=lat(k)
        if(trim(xy_ll_sub).ne.'')then
         xt=x(k)
         yt=y(k)
        endif 
        if(zuv.eq.'z')then
          if(ll_km)z1(1)=-1
          call interp_da(funit,n,m,th_lim,ph_lim, &
                         yt,xt,z1,ncon,cind,ierr,c1,nca)
         else
          if(ll_km)u1(1)=-1
          if(ll_km)v1(1)=-1
          c2='u'                           
          if(c1.eq.'U'.or.c1.eq.'V')c2='U'
          call interp_da(funit,n,m,th_lim,ph_lim, &
                        yt,xt,u1,ncon,cind,ierr,c2,nca)     
          c2='v'                           
          if(c1.eq.'U'.or.c1.eq.'V')c2='V'
          call interp_da(funit,n,m,th_lim,ph_lim, &
                        yt,xt,v1,ncon,cind,ierr,c2,nca)
        endif
        if(ierr.eq.0)then
          if(ll_km)d1=-1
          call interp(dtmp,1,n,m,mz,th_lim,ph_lim, &
                       yt,xt,d1,ierr1,'z')
         if(zuv.eq.'z'.and.geo)then
          call interp_da(funitl,nl,ml,lth_lim,lph_lim, &
                lat(k),lon(k),zl1,ncon,lcind,ierr1,'z',0)
          z1=z1+zl1    ! apply load correction to get geocentric tide
         endif  
! predict tide
         if(narg.eq.0)then ! default usage
          if(zuv.eq.'z')then
           call ptide(z1,c_id,ncon,ccind,lat(k),time_mjd(k),1,&
                     interp_micon,zpred(k))
          else
           call ptide(u1,c_id,ncon,ccind,lat(k),time_mjd(k),1,&
                     interp_micon,upred(k))
           call ptide(v1,c_id,ncon,ccind,lat(k),time_mjd(k),1,&
                     interp_micon,vpred(k))
          endif
!
          write(cdate,'(i2,a1,i2,a1,i2)')hh(k),':',mi(k),':',ss(k)
          ctime=deblank(cdate)
          write(cdate,'(i2,a1,i2,a1,i4)')mm(k),'.',dd(k),'.',yyyy(k)
          cdate=deblank(cdate)
          if(zuv.eq.'z')then
          write(11,'(1x,f10.4,f10.4,5x,a10,1x,a8,f10.3,f10.3)')&
            lat(k),lon0(k),cdate,ctime,zpred(k),real(d1)
          else
          write(11,'(1x,f10.4,f10.4,5x,a10,1x,a8,5(f10.3))')&
           lat(k),lon0(k),cdate,ctime,upred(k),vpred(k),&
           upred(k)/real(d1)*100,vpred(k)/real(d1)*100,real(d1)
          endif
! OB usage March 2011
         else
          do it=1,ntime
           if(zuv.eq.'z')then
            call ptide(z1,c_id,ncon,ccind,lat(k),time_mjd(it),1,&
                      interp_micon,zpred(k))
           else
            call ptide(u1,c_id,ncon,ccind,lat(k),time_mjd(it),1,&
                      interp_micon,upred(k))
            call ptide(v1,c_id,ncon,ccind,lat(k),time_mjd(it),1,&
                     interp_micon,vpred(k))
           endif
!
           write(cdate,'(i2,a1,i2,a1,i2)')hh(it),':',mi(it),':',ss(it)
           ctime=deblank(cdate)
           write(cdate,'(i2,a1,i2,a1,i4)')mm(it),'.',dd(it),'.',yyyy(it)
           cdate=deblank(cdate)
           if(it.eq.1)then
            write(11,'(1x,f10.4,f10.4)')lat(k),lon0(k)
           endif
           if(zuv.eq.'z')then
            write(11,'(26x,a10,1x,a8,f10.3,f10.3)')&
                  cdate,ctime,zpred(k),real(d1)
           else
            write(11,'(26x,a10,1x,a8,5(f10.3))')&
                  cdate,ctime,upred(k),vpred(k),&
             upred(k)/real(d1)*100,vpred(k)/real(d1)*100,real(d1)
           endif
          enddo
         endif
        else
          write(11,'(1x,f10.4,f10.4,a)')lat(k),lon(k),&
       '***** Site is out of model grid OR land *****'
        endif  
      enddo
      deallocate(cind,ccind,lat,lon,lon0,dtmp,mz)
      if(zuv.eq.'z')then
       deallocate(z1,zpred)
      else
       deallocate(u1,v1,upred,vpred)
      endif
      if(trim(xy_ll_sub).ne.'')deallocate(x,y)
      if(zuv.eq.'z'.and.ibl.eq.1.and.geo)then
         ncon=0
         deallocate(zl1,lcind)
      endif
      close(11)
      close(12)
      close(funit(1))
      if(nca.ne.0.and.ncon.gt.1)then
       do ic=2,ncon
        close(funit(ic))
       enddo
      endif
      if(geo)close(funitl(1))
      write(*,*)'Results are in ',trim(outname)
      stop
1     write(*,*)'Lat lon file ',trim(lltname),' not found'
      write(*,*)'Check setup file, line 2.'
      stop
4     write(*,*)'Grid file ',trim(gname),' not found'
      write(*,*)'Check file ',trim(modname),', line 3'
      stop
6     write(*,*)'File ',trim(lltname),' not found'
      write(*,*)'Check setup file, line 7.'
      stop
16    write(*,*)'File ',trim(tfname),' not found'
      write(*,*)'Check spelling in command line'
      call usage()
      stop
11    write(*,*)'File ''model.list'' was NOT found...'
      write(*,*)'TO CREATE please do:'
      write(*,*)'ls -1 DATA/Model_*>model.list'
      stop
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine usage()
      write(*,*)'Usage:'
      write(*,*)'predict_tide [-t<time_file>]<setup.inp'
      write(*,*)'Default: lat_lon_time file is used for input'
      write(*,*)'          if option -t is given, then'
      write(*,*)'          lats/lons are read from lat_lon file'
      write(*,*)'          set in setup.inp, snd times are read'
      write(*,*)'          from <time_file>.'
      write(*,*)'          Use, when need output for time series'
      write(*,*)'          for open boundaries, i.e. times are'
      write(*,*)'          the same in all nodes'
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_time(ntime,iunit,iopt,yyyy,mm,dd,hh,&
                           mi,ss,time_mjd)
      implicit none
      integer k,iunit,iopt,julian,mjd,mm1,dd1,yyyy1
      integer ntime,yyyy(ntime),mm(ntime),dd(ntime)
      integer hh(ntime),mi(ntime),ss(ntime)
      real dum
      real*8 time_mjd(ntime)
      character*10 cdate,deblank
!
      do k=1,ntime
       if(iopt.eq.1)then
        read(iunit,*)dum,dum,yyyy(k),mm(k),dd(k),hh(k),mi(k),ss(k)
       else
        read(iunit,*)yyyy(k),mm(k),dd(k),hh(k),mi(k),ss(k)
       endif
! convert to mjd
       call date_mjd(mm(k),dd(k),yyyy(k),mjd)
! check if exists such a date
        julian=mjd+2400001        
        call CALDAT (JULIAN,MM1,DD1,YYYY1)
        if(mm(k).ne.mm1.or.dd(k).ne.dd1.or.yyyy(k).ne.yyyy1)then
         write(cdate,'(i2,a1,i2,a1,i4)')mm(k),'.',dd(k),'.',yyyy(k)
         cdate=deblank(cdate)
         write(*,*)'Wrong date in (lat_lon)_time file:',cdate
         stop
        endif
        time_mjd(k)=dble(mjd)+dble(hh(k))/24.D0+ &
                    dble(mi(k))/(24.D0*60.D0)  + &
                    dble(ss(k))/(24.D0*60.D0*60.D0)
      enddo
      return
      end
