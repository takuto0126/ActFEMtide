!#[3]## Prepare primary magnetic field
!  CALL PREPXYS(g_param,sparam)         ! commented out on 2017.05.14
!  CALL PREPZOBS(g_param,g_mesh,sparam) ! calculate z of obs and source
  nline   = g_line%nline
  ntet    = g_line%ntet
  nsr     = sparam%nsource     ! 2018.02.22
  ixyflag = g_param%ixyflag    ! 2018.02.22
  nfreq   = g_param%nfreq      ! 2018.02.22

!#[4]## allocate respdata and open output files for each observatory

!#[5]## allocate global matrix
  dofn=1; ip=0
  ! nline is the number of triangular mesh 2020.10.29
  allocate( table_dof(nline,dofn))
  CALL SET_TABLE_DOF(dofn,nline,table_dof,nline)
  CALL set_iccg_var7_dofn(dofn,6,nline,nline,nline,ntet,g_line%n6line,table_dof,A,ip)

!#[6]## prepare coefficients
  CALL PREPOBSCOEFF(g_param,g_mesh,g_line,coeffobs) ! for x,y,z component

!#[7]## set resp5
  allocate( resp5(5,nsr,nfreq)) ! 2018.02.22
  CALL ALLOCATERESP(g_param%nobs,nsr,resp5,ip,nfreq) ! 2018.02.22, see below

!#[3]## Prepare coefficients for values at observatories
  if ( ixyflag .eq. 1 .and. ip .lt. nfreq ) then
   CALL PREPOBSCOEFF_XY(g_param,g_mesh,h_mesh,g_line,obs_xy) ! 2018.02.22
   allocate( resp_xy(5,nsr,nfreq) )                          ! 2018.02.22
   CALL ALLOCATERESP(obs_xy%nobs,nsr,resp_xy,ip,nfreq)       ! 2018.02.22
  end if
  if ( ixyflag .eq. 2) then ! 2017.10.11
   write(*,*) "GEGEGE not suported in the current version"   ! 2017.10.12
   stop
  end if

 allocate( fp(nline,nsr),fs(nline,nsr) ) ! 2018.02.22 for multiple source

!==================================================================  freq loop
do i=1,nfreq

  freq  = g_param%freq(i)   ! 2018.02.22
  omega = 2.d0*(4.d0*datan(1.d0))*freq
  write(*,*) "frequency =",g_param%freq(i),"[Hz]"

 !#[7]## conduct forward calculation with model
 CALL forward_bxyz(A,g_mesh,g_line,nline,nsr,fs,g_param%freq(i),sparam,g_param,g_cond,ip)

 !#[8]## calculate response at every observation point
    !# calculate bx,by,bz
    CALL CALOBSEBCOMP(fp,fs,nline,nsr,omega,coeffobs,resp5(:,:,i)) !2017.07.11 see below
    ! E field at xy plane
    if ( ixyflag .eq. 1 ) then ! 2017.10.12
     CALL CALOBSEBCOMP(fp,fs,nline,nsr,omega,obs_xy%coeff,resp_xy(:,:,i)) !2018.02.22
     CALL OUTFREQFILES2(freq,nsr,resp_xy(:,:,i),g_param,obs_xy)!m_outresp.f90
    end if

!#[9]## output resp to frequency file
    j=1
    if (.false.) CALL OUTFREQ(freq,g_param,resp5(:,j,i)) !2017.07.11, see below

end do ! freq loop end

  !#[8]## output resp to obs file
   CALL OUTOBSFILESFWD(g_param,sparam,nsr,resp5,nfreq) !2017.07.11 m_outresp.f90

