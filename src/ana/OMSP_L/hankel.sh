# Coded on 2016.10.21
#
#!/bin/bash
make clean
make
./hankel_takuto_L.exe < ana1.ctl
./hankel_takuto_L.exe < ana2.ctl
./hankel_takuto_L.exe < ana3.ctl
./hankel_takuto_L.exe < ana4.ctl
./hankel_takuto_L.exe < ana5.ctl