#!/bin/bash

gfortran src/generr_tide.f90  -o ./src/generr_tide.exe

cat <<EOF > err.ctl
        10     20->|header
# of obs files     |14
# of frequencies   |1
input folder       |./result/
output folder      |./data_inv/
in file1 ***.dat   |P5_S1.dat
data errfile head  |P5_tide
in file2 ***.dat   |P6_S1.dat
data errfile head  |P6_tide
in file3 ***.dat   |P7_S1.dat
data errfile head  |P7_tide
in file4 ***.dat   |P8_S1.dat
data errfile head  |P8_tide
in file5 ***.dat   |P9_S1.dat
data errfile head  |P9_tide
in file6 ***.dat   |P10_S1.dat
data errfile head  |P10_tide
in file7 ***.dat   |P11_S1.dat
data errfile head  |P11_tide
in file8 ***.dat   |P12_S1.dat
data errfile head  |P12_tide
in file9 ***.dat   |P13_S1.dat
data errfile head  |P13_tide
in file10 ***.dat  |P14_S1.dat
data errfile head  |P14_tide
in file11 ***.dat  |P15_S1.dat
data errfile head  |P15_tide
in file12 ***.dat  |P16_S1.dat
data errfile head  |P16_tide
in file13 ***.dat  |P17_S1.dat
data errfile head  |P17_tide
in file14 ***.dat  |P18_S1.dat
data errfile head  |P18_tide
ratio              !0.1
EOF

./src/generr_tide.exe < err.ctl

rm err.ctl
