## lines starting with "##" work as lines for comments ! 2020.09.28
## this file is inversion control file for tide-only inversion (ijoint = 1)
##-----10!-------20! inv201408_ap.ctl
face info file     !../../mesh_joint/faceinfo.dat
1:cond,2:model     !1
ref  cond file     !../structure/cond_homo.msh
1:cond,2:model     !1
init cond file     !../structure/cond_homo.msh
output folder      |./result_inv/
Roughness type     |1
1:L,2:Cl,3:Mi,4:Gr !2
alpha init         !100.
factor(10^factor)  !-0.25
iflag_replace 0,1  |0
## iboundflag can set upper and lower limit of the conductivity value in inversion
## iboundflag = 0 : off : no boundary for conductivity value
## iboundflag = 1 : simple upper and lower limit will be specified by cutting
## iboundflag = 2 : transformed model variable that automatically satisfies upper
##                  and lower limit are used. (see Grayver et al. 2013)
## Note that if you use Smooth constraint, 0 is recommended.
 iboundflag =      !0
## ACTIVE data #########################################################
# of srces for inv !1
## select components if one or more data will be used for the component
Bx,By,Bz,Ex,Ey(5i2)!1 1 1 0 0
Act err floor [0-1]!0.03
index of source(S1)!1
# of obesrvatories !14
P5  amp data       !   1 ../fwd_tide_H/data_inv/P5_tide_amp.dat
P6  amp data       !   2 ../fwd_tide_H/data_inv/P6_tide_amp.dat
P7  amp data       !   3 ../fwd_tide_H/data_inv/P7_tide_amp.dat
P8  amp data       !   4 ../fwd_tide_H/data_inv/P8_tide_amp.dat
P9  amp data       !   5 ../fwd_tide_H/data_inv/P9_tide_amp.dat
P10 amp data       !   6 ../fwd_tide_H/data_inv/P10_tide_amp.dat
P11 amp data       !   7 ../fwd_tide_H/data_inv/P11_tide_amp.dat
P12 amp data       !   8 ../fwd_tide_H/data_inv/P12_tide_amp.dat
P13 amp data       !   9 ../fwd_tide_H/data_inv/P13_tide_amp.dat
P14 amp data       !  10 ../fwd_tide_H/data_inv/P14_tide_amp.dat
P15 amp data       !  11 ../fwd_tide_H/data_inv/P15_tide_amp.dat
P16 amp data       !  12 ../fwd_tide_H/data_inv/P16_tide_amp.dat
P17 amp data       !  13 ../fwd_tide_H/data_inv/P17_tide_amp.dat
P18 amp data       !  14 ../fwd_tide_H/data_inv/P18_tide_amp.dat
## phase
P5  pha data       !   1 ../fwd_tide_H/data_inv/P5_tide_pha.dat
P6  pha data       !   2 ../fwd_tide_H/data_inv/P6_tide_pha.dat
P7  pha data       !   3 ../fwd_tide_H/data_inv/P7_tide_pha.dat
P8  pha data       !   4 ../fwd_tide_H/data_inv/P8_tide_pha.dat
P9  pha data       !   5 ../fwd_tide_H/data_inv/P9_tide_pha.dat
P10 pha data       !   6 ../fwd_tide_H/data_inv/P10_tide_pha.dat
P11 pha data       !   7 ../fwd_tide_H/data_inv/P11_tide_pha.dat
P12 pha data       !   8 ../fwd_tide_H/data_inv/P12_tide_pha.dat
P13 pha data       !   9 ../fwd_tide_H/data_inv/P13_tide_pha.dat
P14 pha data       !  10 ../fwd_tide_H/data_inv/P14_tide_pha.dat
P15 pha data       !  11 ../fwd_tide_H/data_inv/P15_tide_pha.dat
P16 pha data       !  12 ../fwd_tide_H/data_inv/P16_tide_pha.dat
P17 pha data       !  13 ../fwd_tide_H/data_inv/P17_tide_pha.dat
P18 pha data       !  14 ../fwd_tide_H/data_inv/P18_tide_pha.dat
## MT data section is not read when ijoint = 1 (tide only); do not write it here
########################################################################
## icombine = 0: normal
## icombine = 1: integrate the outside blocks to one and assign one model parameter
## icombine = 2: integrate the outside blocks to one and fix the modelparameter with given cond
## icombine = -1: inner outer mode
icombine:0,1,2:fix !1
27
-250.0
-150.0
-100.0
-80.0
-70.0
-60.0
-50.0
-42.0
-35.0
-28.0
-21.0
-14.0
-7.0
0.0
7.0
14.0
21.0
28.0
35.0
42.0
50.0
60.0
70.0
80.0
100.0
150.0
250.0
27
-250.0
-150.0
-100.0
-80.0
-70.0
-60.0
-50.0
-42.0
-35.0
-28.0
-21.0
-14.0
-7.0
0.0
7.0
14.0
21.0
28.0
35.0
42.0
50.0
60.0
70.0
80.0
100.0
150.0
250.0
27
-350.0
-300.0
-250.0
-200.0
-170.0
-140.0
-110.0
-80.0
-60.0
-50.0
-40.0
-35.0
-30.0
-25.0
-20.0
-17.0
-14.0
-12.0
-11.0
-10.0
-9.0
-8.0
-6.0
-4.0
-2.0
-1.0
0.0
ioutlevel:0,1:Jacob|0
final rms          |1.0

