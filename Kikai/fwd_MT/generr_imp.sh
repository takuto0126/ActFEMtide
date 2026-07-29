#!/bin/bash

#ifort -g generr_mt_imp.f90  -o generr.exe
gfortran generr_mt_imp.f90  -o ./generr.exe

cat <<EOF > err.ctl
        10     20->|header
# of obs files     |9
# of frequencies   |5
input folder       |./result/
output folder      |./mt_err/
in file1 ***.dat   |MT1_MT_imp.dat
errfile1           |MT1_MT_imp_err.dat
in file2 ***.dat   |MT2_MT_imp.dat
errfile2           |MT2_MT_imp_err.dat
in file3 ***.dat   |MT3_MT_imp.dat
errfile3           |MT3_MT_imp_err.dat
in file4 ***.dat   |MT4_MT_imp.dat
errfile4           |MT4_MT_imp_err.dat
in file5 ***.dat   |MT5_MT_imp.dat
errfile5           |MT5_MT_imp_err.dat
in file6 ***.dat   |MT6_MT_imp.dat
errfile6           |MT6_MT_imp_err.dat
in file7 ***.dat   |MT7_MT_imp.dat
errfile7           |MT7_MT_imp_err.dat
in file8 ***.dat   |MT8_MT_imp.dat
errfile8           |MT8_MT_imp_err.dat
in file9 ***.dat   |MT9_MT_imp.dat
errfile9           |MT9_MT_imp_err.dat
ratio              !0.01
EOF

./generr.exe < err.ctl

rm err.ctl

