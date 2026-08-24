#!/bin/bash

ite=10
b=bz
echo $b
F=result_inv

for S in S1 S2
do
gmt begin ${b}_${S}_fit pdf
## Amplitude
# Amp comp
gmt gmtset FONT_ANNOT_PRIMARY = "10p,Helvetica,black"
Bxa=xa1f3g1
Bya=ya1f3g1
Bxp=xa1f3g1+l"Frequency[Hz]" # x axis for phase plot
Byp=ya45g15 # y axis for phase plot
RA=0.8/110/0.005/0.2
RP=0.8/110/90/190
E="-Sc0.3 -W1.5,red -Ey+p0.5+w0.5"
FF="-JX6/8 -R0/6/0/8 -F+f12p,Helvetica,black+jLM -Gwhite -W1,black"
#A01
gmt basemap -JX6l/8l -R$RA -B${Bxa} -B${Bya}+l"Amplitude[nT/A]" -X3 -Y13
cat ${F}/A01_${S}_${b}ampIN.dat | awk '{print($1,$3,$4)}' | gmt plot $E
cat ${F}/A01_${S}_${b}amp01.dat | awk '{print($1,$3)}'    | gmt plot -W2,green 
cat ${F}/A01_${S}_${b}amp${ite}.dat | awk '{print($1,$3)}'| gmt plot -W2,blue
echo "0.5 0.5 A01 $b $S" | gmt text $FF

#A02
gmt basemap -JX6l/8l -R$RA -B${Bxa} -B${Bya} -X7
cat ${F}/A02_${S}_${b}ampIN.dat | awk '{print($1,$3,$4)}' | gmt plot $E
cat ${F}/A02_${S}_${b}amp01.dat | awk '{print($1,$3)}'    | gmt plot -W2,green 
cat ${F}/A02_${S}_${b}amp${ite}.dat | awk '{print($1,$3)}'| gmt plot -W2,blue 
echo "0.5 0.5 A02 $b $S" | gmt text $FF

#A03
gmt basemap -JX6l/8l -R$RA -B${Bxa} -B${Bya} -X7
cat ${F}/A03_${S}_${b}ampIN.dat | awk '{print($1,$3,$4)}' | gmt plot $E
cat ${F}/A03_${S}_${b}amp01.dat | awk '{print($1,$3)}'    | gmt plot -W2,green 
cat ${F}/A03_${S}_${b}amp${ite}.dat | awk '{print($1,$3)}'| gmt plot -W2,blue 
echo "0.5 0.5 A03 $b $S" | gmt text $FF

#A04
gmt basemap -JX6l/8l -R$RA -B${Bxa} -B${Bya} -X7
cat ${F}/A04_${S}_${b}ampIN.dat | awk '{print($1,$3,$4)}' | gmt plot $E
cat ${F}/A04_${S}_${b}amp01.dat | awk '{print($1,$3)}'    | gmt plot -W2,green 
cat ${F}/A04_${S}_${b}amp${ite}.dat | awk '{print($1,$3)}'| gmt plot -W2,blue 
echo "0.5 0.5 A04 $b $S" | gmt text $FF

## Phase
gmt basemap -JX6l/8 -R$RP -B${Bxp} -B${Byp}+l"Phase[deg]" -X-21 -Y-9
#A01
cat ${F}/A01_${S}_${b}phaIN.dat | awk '{print($1,$3,$4)}' | gmt plot $E
cat ${F}/A01_${S}_${b}pha01.dat | awk '{print($1,$3)}'    | gmt plot -W1,green 
cat ${F}/A01_${S}_${b}pha${ite}.dat | awk '{print($1,$3)}'| gmt plot -W1,blue 
echo "0.5 0.5 A01 $b $S" | gmt text $FF

#A02
gmt basemap -JX6l/8 -R$RP -B${Bxp} -B${Byp} -X7
cat ${F}/A02_${S}_${b}phaIN.dat | awk '{print($1,$3,$4)}' | gmt plot $E
cat ${F}/A02_${S}_${b}pha01.dat | awk '{print($1,$3)}'    | gmt plot -W1,green 
cat ${F}/A02_${S}_${b}pha${ite}.dat | awk '{print($1,$3)}'    | gmt plot -W1,blue 
echo "0.5 0.5 A02 $b $S" | gmt text $FF
#A03
gmt basemap -JX6l/8 -R$RP -B${Bxp} -B${Byp} -X7
cat ${F}/A03_${S}_${b}phaIN.dat | awk '{print($1,$3,$4)}' | gmt plot $E
cat ${F}/A03_${S}_${b}pha01.dat | awk '{print($1,$3)}'    | gmt plot -W1,green 
cat ${F}/A03_${S}_${b}pha${ite}.dat | awk '{print($1,$3)}'    | gmt plot -W1,blue 
echo "0.5 0.5 A03 $b $S" | gmt text $FF
#A04
gmt basemap -JX6l/8 -R$RP -B${Bxp} -B${Byp} -X7
cat ${F}/A04_${S}_${b}phaIN.dat | awk '{print($1,$3,$4)}' | gmt plot $E
cat ${F}/A04_${S}_${b}pha01.dat | awk '{print($1,$3)}'    | gmt plot -W1,green 
cat ${F}/A04_${S}_${b}pha${ite}.dat | awk '{print($1,$3)}'    | gmt plot -W1,blue 
echo "0.5 0.5 A04 $b $S" | gmt text $FF

gmt end


done


open ${b}_S1_fit.pdf &
open ${b}_S2_fit.pdf &
