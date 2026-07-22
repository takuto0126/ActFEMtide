!cc   This file contains the standard parameters which define the
!cc   amplitudes, frequencies, etc. for the primary tidal constituents
      integer, parameter :: ncmx = 48
      integer, parameter :: ngmx = 4 ! max # constituents/group
      character*4 constid(ncmx),cid8(8)
      data constid                              &
                  /'m2  ','s2  ','k1  ','o1  ', &
                   'n2  ','p1  ','k2  ','q1  ', &
                   '2n2 ','mu2 ','nu2 ','l2  ', &
                   't2  ','j1  ','m1  ','oo1 ', &
                   'rho1','mf  ','mm  ','ssa ', &
                   'm4  ','ms4 ','mn4 ','m6  ', &
                   'm8  ','mk3 ','s6  ','msf ', &
            '2mk3','m2_1','m2_2','o1_1','k1_1', &
            's2_1','sa  ','sta ','msm ','mst ', &
  'mt  ','msq ','mq ','msp  ','2q1 ','node','m2a ','m2b ','s1  ','m3  '/
! these 8 constituents used to infer minor constituents
       data cid8/'m2  ','s2  ','n2  ','k2  ', &
                 'k1  ','o1  ','p1  ','q1  '/
!    FOR EACH POSSIBLE CONSTIUENT, these parameters are given:
!    alpha = correction factor for first order load tides
!    amp = amplitude of equilibrium tide in m
!    ph = Currently set to zero ...   phases for
!             each constituent are referred to the time
!             when the phase of the forcing for that
!             constituent is zero on the Greenich meridian.)
!
!    omega = angular frequency of constituent, in radians
      real*8 alpha_d(ncmx),ph_d(ncmx),amp_d(ncmx),omega_d(ncmx) &
            ,phase_mkB(ncmx),beta_SE(ncmx)
      integer ispec_d(ncmx)

!     Tidal parameters taken from Rodney's constituent.h, 2/23/96:
!     (except for ispec).
      data ispec_d/&
          2,2,1,1, &
          2,1,2,1, &
          2,2,2,2, &
          2,1,1,1, &
          1,0,0,0, &
          4,4,4,6, &
          8,3,6,0, &
         3,2,2,1,1, &
         2,0,0,0,0, &
       0,0,0,0,1,0,2,2,1,3/
!cc     note: for now I am just leaving ispec for M4 set to 0 (ispec
!cc     is only used to define forcing in atgf, and this is always  0
!cc     for M4)

      data alpha_d/&
          0.693,0.693,0.736,0.695, &
          0.693,0.706,0.693,0.695, &
          0.693,0.693,0.693,0.693, &
          0.693,0.695,0.695,0.695, &
          0.695,0.693,0.693,0.693, &
          0.693,0.693,0.693,0.693, &
          0.693,0.693,0.693,0.693, &
    0.693,0.693,0.693,0.695,0.736, &
    0.693,0.693,0.693,0.693,0.693, &
    0.693,0.693,0.693,0.693,0.693,0.693,0.693,0.693,0.693,0.802/

!
!cpy=omega_m2/(2*pi)*3600*24*365;omp1=(cpy+1)/3600/24/365*2*pi=1.407181e-04
!omm1=(cpy-1)/3600/24/365*2*pi=1.403197e-04
!omega_m2=1.405189e-04 (initial value)
!
      data omega_d/&
          1.405189e-04,1.454441e-04,7.292117e-05,6.759774e-05, &
          1.378797e-04,7.252295e-05,1.458423e-04,6.495854e-05, &
          1.352405e-04,1.355937e-04,1.382329e-04,1.431581e-04, &
          1.452450e-04,7.556036e-05,7.028195e-05,7.824458e-05, &
          6.531174e-05,0.053234e-04,0.2639204e-05,0.003982e-04, &
          2.810377e-04,2.859630e-04,2.783984e-04,4.215566e-04, &
          5.620755e-04,2.134402e-04,4.363323e-04,4.925200e-06, &
         2.081166e-04,1.405189e-04,1.405189e-04,6.759774e-05,7.292117e-05, &
     1.454441e-04,0.1990970E-06,0.5973098E-06,0.2285998E-05,0.7609364E-05, &
          0.7962619E-05,0.1024862E-04,0.1060182E-04,0.1288782E-04, &
          0.6231934E-04,0.1069693E-07,1.403197e-04,1.407181e-04,7.2722e-05, &
          2.107783523e-04/

      data ph_d/48*0.0/

      data amp_d/ &
          0.244102,0.113568,0.142435,0.101270, &
          0.046735,0.047129,0.030879,0.019387, &
          0.006184,0.007408,0.008811,0.006931, &
          0.006608,0.007965,0.007915,0.004361, &
          0.003661,0.042041,0.022193,0.019547, &
!cc       amplitude for M4 etc. is zero
          0.,0.,0.,0., &
          0.,0.,0.,0.003681, &
          0.,0.242334,0.242334,0.101270,0.142435, &
          0.113568,0.003104,0.001141,0.004244,0.001528, &
          0.008044,0.001285,0.001064,0.000310,0.002565,0.017617,0.244102,0.244102,7.6464e-04,0.003192/
 
! Astronomical arguments, obtained with Richard Ray's
! "arguments" and "astrol", for Jan 1, 1992, 00:00 Greenwich time
! Corrected July 12, 2000
       data phase_mkB/ &
          1.731572000,0.000000000,0.173006000,1.558566000, &
          6.050735000,6.110179000,3.487605000,5.877729000, &
          4.086713000,3.463115091,5.427136701,0.553986502, &
          0.052841931,2.137028000,2.436575100,1.929039000, &
          5.254133027,1.756033000,1.964022000,3.487605000, &
          3.463115091,1.731557546,1.499093481,5.194672637, &
          6.926230184,1.904561220,0.000000000,4.551613000, &
          3.809122439,1.731557546,1.731557546,1.558566000,0.173006000, &
	  0.000000000,4.885395,3.437166,2.587591,4.343624, &
  3.720055,0.024461,5.684077,1.988483,3.913707,4.541314,1.7843995692782,1.6787156184449,0.000000000,5.738991/
! I am putting 0 for ms2,mn4 etc. for now: correct later
! Now this correction is done using the SAL file (h_TPXO3_90-90.load)
! I replace beta_SE with units for now (on case we decide to switch back
! to old version) and comment the old numbers - this way I do NOT change
! anything in subroutines  
! This was in weights.h before - placed here not to mix with w!
! to remove solid Earth tide multily by beta:
       data beta_SE/ &
          0.9540,0.9540,0.9400,0.9400, &
          0.9540,0.9400,0.9540,0.9400, &
          0.9540,0.9540,0.9540,0.9540, &
          0.9540,0.9400,0.9400,0.9400, &
          0.9400,0.9400,0.9400,0.9400, &
!ccc      for M4 just using value for semi-diurnals (no good reason!)
          0.9540,0.9540,0.9540,0.954, &
          0.9540,0.9540,0.9540,0.954, &
	  0.9540,0.9540,0.9540,0.954,0.954, &
          0.9540,0.9540,0.9540,0.9540,0.954, &
  0.9540,0.9540,0.9540,0.954,0.954,0.954,0.954,0.954,0.94,0.954/
!       data beta_SE/48*1./
