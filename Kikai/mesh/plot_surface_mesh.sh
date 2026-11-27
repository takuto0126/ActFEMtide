#!/bin/bash
# Plot sea-surface and ground-surface triangles of em3d.msh with GMT       2026.11.27 CC
#  em3d.msh has only tetrahedra (type 4), so the surface triangles are      2026.11.27 CC
#  extracted as faces shared by two tetrahedra of different regions:        2026.11.27 CC
#   sea surface    = air (region 1) / ocean (region 2)                      2026.11.27 CC
#   ground surface = air (region 1) / earth (region 3)                      2026.11.27 CC
#  The region is the 2nd tag of each element (same as n4flag(iele,2)).      2026.11.27 CC
# Usage: ./plot_surface_mesh.sh [half width of zoom area in km (default 200)] 2026.11.27 CC
# Output: surface_mesh.png                                                  2026.11.27 CC
set -e                                                     # 2026.11.27 CC

msh=em3d.msh                                               # 2026.11.27 CC
zoom=${1:-200}                                             # 2026.11.27 CC
out=surface_mesh                                           # 2026.11.27 CC
sea=surface_sea.xy                                         # 2026.11.27 CC
land=surface_land.xy                                       # 2026.11.27 CC
nonconf=surface_nonconforming.xy                           # 2026.11.27 CC
obs=obs_xyz.dat                                            # written by tmtgem mesh programs 2026.11.27 CC

[ -s $msh ] || { echo "$msh not found. Run tetmeshgen.sh first."; exit 1; } # 2026.11.27 CC

#[1]## extract surface triangles as GMT multi-segment polygons  2026.11.27 CC
python3 - $msh $sea $land $nonconf <<'EOF'
import sys
msh, fsea, fland, fnon = sys.argv[1:5]
xyz = {}
faces = {}   # sorted node triple -> set of regions
with open(msh) as f:
    mode = None
    for line in f:
        if line.startswith('$Nodes'):
            mode = 'n'; next(f); continue
        if line.startswith('$EndNodes'):
            mode = None; continue
        if line.startswith('$Elements'):
            mode = 'e'; next(f); continue
        if line.startswith('$EndElements'):
            break
        a = line.split()
        if mode == 'n':
            xyz[int(a[0])] = (float(a[1]), float(a[2]), float(a[3]))
        elif mode == 'e' and a[1] == '4':
            nt = int(a[2]); reg = int(a[4]) if nt >= 2 else int(a[3])
            n = [int(v) for v in a[3 + nt:7 + nt]]
            for tri in ((n[0], n[1], n[2]), (n[0], n[1], n[3]), (n[0], n[2], n[3]), (n[1], n[2], n[3])):
                k = tuple(sorted(tri))
                s = faces.get(k)
                if s is None:
                    faces[k] = reg
                elif s != reg:
                    faces[k] = (s, reg)
                else:
                    faces[k] = -1   # interior face of one region
nsea = nland = 0
with open(fsea, 'w') as fs, open(fland, 'w') as fl:
    for k, s in faces.items():
        if not isinstance(s, tuple):
            continue
        pair = tuple(sorted(s))
        if pair == (1, 2):
            fo = fs; nsea += 1
        elif pair == (1, 3):
            fo = fl; nland += 1
        else:
            continue
        zc = sum(xyz[i][2] for i in k) / 3.0
        fo.write('> -Z%.6f\n' % zc)
        for i in k + (k[0],):
            fo.write('%.4f %.4f\n' % (xyz[i][0], xyz[i][1]))
# 2026.11.27 CC: sea-surface triangles of the ocean mesh (z = 0) that are not shared by
# any air tetrahedron = non-conforming interface (gmsh refined the air-side surface)
nnon = 0
with open(fnon, 'w') as fn:
    for k, s in faces.items():
        if s == 2 and all(abs(xyz[i][2]) < 1e-6 for i in k):
            nnon += 1
            fn.write('>\n')
            for i in k + (k[0],):
                fn.write('%.4f %.4f\n' % (xyz[i][0], xyz[i][1]))
xs = [p[0] for p in xyz.values()]; ys = [p[1] for p in xyz.values()]
print('sea-surface triangles   :', nsea)
print('ground-surface triangles:', nland)
print('non-conforming sea-surface triangles of ocean mesh:', nnon)
print('REGION %.4f %.4f %.4f %.4f' % (min(xs), max(xs), min(ys), max(ys)))
EOF

#[2]## plot with GMT  2026.11.27 CC
read xmin xmax ymin ymax < <(python3 -c "
import sys
xs=[];ys=[]
for f in ('$sea','$land'):
    for l in open(f):
        if l[0]!='>':
            a=l.split(); xs.append(float(a[0])); ys.append(float(a[1]))
print(min(xs),max(xs),min(ys),max(ys))")
Rall=$xmin/$xmax/$ymin/$ymax
Rzoom=-$zoom/$zoom/-$zoom/$zoom
[ -s $obs ] && awk '$1!="name" && NF>=2{print $1,$2}' $obs > surface_obs.xy || : > surface_obs.xy

gmt begin $out png E200
gmt set FONT_TITLE 12p FONT_LABEL 10p FONT_ANNOT_PRIMARY 9p
gmt makecpt -Cgeo -T-1/3/0.1 -H > surface_land.cpt   # land elevation [km]
gmt subplot begin 1x2 -Fs12c/12c -M1c/0.5c -T"Surface triangles of em3d.msh (blue: sea surface, colored: ground, pink: non-conforming sea surface)"
#--- whole domain
gmt subplot set 0
gmt plot $land -R$Rall -JX12c -Csurface_land.cpt -L -W0.1p,gray30 -Bxaf+l"Easting [km]" -Byaf+l"Northing [km]" -BWSne+t"Whole domain"
gmt plot $sea -W0.1p,dodgerblue3 -L
[ -s $nonconf ] && gmt plot $nonconf -Gpink -W0.1p,magenta -L   # 2026.11.27 CC
gmt plot -W0.6p,red -A <<EOT
-$zoom -$zoom
$zoom -$zoom
$zoom $zoom
-$zoom $zoom
-$zoom -$zoom
EOT
#--- zoom around the focus area
gmt subplot set 1
gmt plot $land -R$Rzoom -JX12c -Csurface_land.cpt -L -W0.2p,gray30 -Bxaf+l"Easting [km]" -Byaf+l"Northing [km]" -BWSne+t"Zoom (+-$zoom km)"
gmt plot $sea -W0.2p,dodgerblue3 -L
[ -s $nonconf ] && gmt plot $nonconf -Gpink -W0.2p,magenta -L   # 2026.11.27 CC
[ -s surface_obs.xy ] && gmt plot surface_obs.xy -St0.25c -Gred -W0.3p
gmt subplot end
gmt colorbar -Csurface_land.cpt -DJBC+w10c/0.3c+o0c/1.2c+h -Bxa1+l"Ground elevation [km]"
gmt end
rm -f surface_land.cpt
echo "output: $out.png"
