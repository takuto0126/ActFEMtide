#!/bin/bash

#ifort -g generr_mt_imp.f90  -o generr.exe
gfortran src/generr_mt_imp.f90  -o ./src/generr.exe

cat <<EOF > err.ctl
        10     20->|header
# of obs files     |14
# of frequencies   |5
input folder       |./result/
output folder      |./data_inv/
infile1 ***.dat    |P5_MT_imp.dat
errfile1           |P5_MT_imp_err.dat
infile2 ***.dat    |P6_MT_imp.dat
errfile2           |P6_MT_imp_err.dat
infile3 ***.dat    |P7_MT_imp.dat
errfile3           |P7_MT_imp_err.dat
infile4 ***.dat    |P8_MT_imp.dat
errfile4           |P8_MT_imp_err.dat
infile5 ***.dat    |P9_MT_imp.dat
errfile5           |P9_MT_imp_err.dat
infile6 ***.dat    |P10_MT_imp.dat
errfile6           |P10_MT_imp_err.dat
infile7 ***.dat    |P11_MT_imp.dat
errfile7           |P11_MT_imp_err.dat
infile8 ***.dat    |P12_MT_imp.dat
errfile8           |P12_MT_imp_err.dat
infile9 ***.dat    |P13_MT_imp.dat
errfile9           |P13_MT_imp_err.dat
infile10 ***.dat   |P14_MT_imp.dat
errfile10          |P14_MT_imp_err.dat
infile11 ***.dat   |P15_MT_imp.dat
errfile11          |P15_MT_imp_err.dat
infile12 ***.dat   |P16_MT_imp.dat
errfile12          |P16_MT_imp_err.dat
infile13 ***.dat   |P17_MT_imp.dat
errfile13          |P17_MT_imp_err.dat
infile14 ***.dat   |P18_MT_imp.dat
errfile14          |P18_MT_imp_err.dat
ratio              !0.03
EOF

./src/generr.exe < err.ctl

rm err.ctl

