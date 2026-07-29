# generated on 2026.03.01
#!/bin/bash

gfortran generr_mt_tip.f90  -o generr_tip.exe

cat <<EOF > err_tip.ctl
        10     20->|header
# of obs files     |9
# of frequencies   |5
input folder       |./result/
output folder      |./mt_err/
in file1 ***.dat   |MT1_TIP.dat
errfile1           |MT1_TIP_err.dat
in file2 ***.dat   |MT2_TIP.dat
errfile2           |MT2_TIP_err.dat
in file3 ***.dat   |MT3_TIP.dat
errfile3           |MT3_TIP_err.dat
in file4 ***.dat   |MT4_TIP.dat
errfile4           |MT4_TIP_err.dat
in file5 ***.dat   |MT5_TIP.dat
errfile5           |MT5_TIP_err.dat
in file6 ***.dat   |MT6_TIP.dat
errfile6           |MT6_TIP_err.dat
in file7 ***.dat   |MT7_TIP.dat
errfile7           |MT7_TIP_err.dat
in file8 ***.dat   |MT8_TIP.dat
errfile8           |MT8_TIP_err.dat
in file9 ***.dat   |MT9_TIP.dat
errfile9           |MT9_TIP_err.dat
ratio              !0.005
EOF

./generr_tip.exe < err_tip.ctl

rm err_tip.ctl

