region=-700/700/-700/700
gmt surface checkvxyz -Gvxyz.grid -I5 -R$region
#gmt surface result_By -GBx.grid -I5 -R$region
#region=-400/401/-400/401 # map region east/west/south/north                            # map projection and scale
ticks=a5f1g5                            # boundary tick info
frame=WSne+tcolor_fill                     # boundary frameinfo
climit=-100/100/1
#climit=0/10/1  # color table min/max/interval
grdfile=vxyz.grid # input bathymetry grid file
cptfile=haxby_grad.cpt                  # color table
psfile=vxyz.ps                 # output postscript file nameproj=x1

gmt makecpt -Chaxby -T$climit -Z > $cptfile
#gmt makecpt -Chaxby -Qi > $cptfile

gmt grdimage $grdfile -R$region -C$cptfile -Jx0.01 -K -V > $psfile
gmt psscale -D12/3/9/0.5 -C$cptfile -Bf1a1 -K -O >> $psfile
#gmt pscoast -R$region -J$proj -Di -Ggray -Wthin,black -K -V -O >> $psfile
#gmt psbasemap -R$region -J$proj -B$ticks -B$framm -O  -V >> $psfile
gmt psbasemap -R$region -Jx0.01 -Bafg -O  -V >> $psfile


