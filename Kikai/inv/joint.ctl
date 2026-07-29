## lines starting with "##" work as lines for comments ! 2020.09.28
## this file is joint inversion control file
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
Act err floor [0-1]!0.01
index of source(S1)!1
# of obesrvatories !9
P10 amp data       !   1 ../fwd/data_inv/P10_tide_amp.dat
P11 amp data       !   2 ../fwd/data_inv/P11_tide_amp.dat
P12 amp data       !   3 ../fwd/data_inv/P12_tide_amp.dat
P13 amp data       !   4 ../fwd/data_inv/P13_tide_amp.dat
P14 amp data       !   5 ../fwd/data_inv/P14_tide_amp.dat
P15 amp data       !   6 ../fwd/data_inv/P15_tide_amp.dat
P16 amp data       !   7 ../fwd/data_inv/P16_tide_amp.dat
P17 amp data       !   8 ../fwd/data_inv/P17_tide_amp.dat
P18 amp data       !   9 ../fwd/data_inv/P18_tide_amp.dat
## phase
P10 pha data       !   1 ../fwd/data_inv/P10_tide_pha.dat
P11 pha data       !   2 ../fwd/data_inv/P11_tide_pha.dat
P12 pha data       !   3 ../fwd/data_inv/P12_tide_pha.dat
P13 pha data       !   4 ../fwd/data_inv/P13_tide_pha.dat
P14 pha data       !   5 ../fwd/data_inv/P14_tide_pha.dat
P15 pha data       !   6 ../fwd/data_inv/P15_tide_pha.dat
P16 pha data       !   7 ../fwd/data_inv/P16_tide_pha.dat
P17 pha data       !   8 ../fwd/data_inv/P17_tide_pha.dat
P18 pha data       !   9 ../fwd/data_inv/P18_tide_pha.dat
## MT data #############################################################
## initial version support only impedance data 2021.12.13
imp: 0, 1:amp,pha  !0
unit:1Ohm,2mV/km/nT!2
MT errfloor SSQ*   !0.01
# of observatories |5
P10 impedance      !   1 ../fwd_MT/result/P10_MT_imp.dat
P11 impedance      !   2 ../fwd_MT/result/P11_MT_imp.dat
P12 impedance      !   3 ../fwd_MT/result/P12_MT_imp.dat
P13 impedance      !   4 ../fwd_MT/result/P13_MT_imp.dat
P14 impedance      !   5 ../fwd_MT/result/P14_MT_imp.dat
P15 impedance      !   6 ../fwd_MT/result/P15_MT_imp.dat
P16 impedance      !   7 ../fwd_MT/result/P16_MT_imp.dat
P17 impedance      !   8 ../fwd_MT/result/P17_MT_imp.dat
P18 impedance      !   9 ../fwd_MT/result/P18_MT_imp.dat
P10 impedance err  !   1 ../fwd_MT/data_inv/P10_MT_imp_err.dat
P11 impedance err  !   2 ../fwd_MT/data_inv/P11_MT_imp_err.dat
P12 impedance err  !   3 ../fwd_MT/data_inv/P12_MT_imp_err.dat
P13 impedance err  !   4 ../fwd_MT/data_inv/P13_MT_imp_err.dat
P14 impedance err  !   5 ../fwd_MT/data_inv/P14_MT_imp_err.dat
P15 impedance err  !   6 ../fwd_MT/data_inv/P15_MT_imp_err.dat
P16 impedance err  !   7 ../fwd_MT/data_inv/P16_MT_imp_err.dat
P17 impedance err  !   8 ../fwd_MT/data_inv/P17_MT_imp_err.dat
P18 impedance err  !   9 ../fwd_MT/data_inv/P18_MT_imp_err.dat
########################################################################
iflag_tipper       !0
## icombine = 0: normal
## icombine = 1: integrate the outside blocks to one and assign one model parameter
## icombine = 2: integrate the outside blocks to one and fix the modelparameter with given cond
## icombine = -1: inner outer mode
icombine:0,1,2:fix !1
19
-1.5
-1.2
-0.9
-0.75
-0.6
-0.45
-0.3
-0.2
-0.1
0.0
0.1
0.2
0.3
0.45
0.6
0.75
0.9
1.2
1.5
19
-1.5
-1.2
-0.9
-0.75
-0.6
-0.45
-0.3
-0.2
-0.1
0.0
0.1
0.2
0.3
0.45
0.6
0.75
0.9
1.2
1.5
21
0.0
0.25
0.5
0.6
0.7
0.8
0.85
0.9
0.95
1.0
1.025
1.05
1.075
1.1
1.125
1.15
1.175
1.2
1.225
1.25
1.50
ioutlevel:0,1:Jacob|0
final rms          |1.0

