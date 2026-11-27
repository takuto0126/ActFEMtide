#!/usr/bin/env python3
# Coded on 2026.11.27 by Claude Code (CC) for ActFEMtide
# Interpolate HAMTIDE11a tidal currents (HAMcurrent11a_<const>.nc) to the nodes in
# lat_lon_mesh and write u_xyz_mesh / v_xyz_mesh in the same format as the output of
# OTPS extract_HC (lat, lon, amplitude [cm/s], Greenwich phase [deg]), so that
# mkvxyzfile.exe can be used without changes.
#  - u: eastward (UAMP, UPHA), v: northward (VAMP, VPHA)
#  - Interpolation is bilinear on the real and imaginary parts (A cos G, A sin G),
#    not on amplitude and phase, to avoid errors where the phase wraps at 0/360 deg.
#  - Points whose surrounding HAMTIDE cells include land (no data) are written as
#    "Site is out of model grid OR land", like OTPS. They are filled afterwards by
#    fill_nearest_uv.py with the values of the nearest valid point.
#  - One line is written for each run of consecutive identical (lat, lon) in
#    lat_lon_mesh, as OTPS does (mkvxyzfile.exe relies on this).
# Requires GMT (gmt grdmath, gmt grdtrack).
# Usage: hamtide2uv.py <lat_lon_mesh> <HAMcurrent11a_const.nc> <const> <u_out> <v_out>
import math
import os
import subprocess
import sys
import tempfile

llfile, ncfile, const, uout, vout = sys.argv[1:6]

#[1]## read lat_lon_mesh and keep one point per run of identical (lat, lon)
pts = []
prev = None
with open(llfile) as f:
    for line in f:
        a = line.split()
        if len(a) < 2:
            continue
        lat, lon = float(a[0]), float(a[1])
        if (lat, lon) != prev:
            pts.append((lat, lon))
            prev = (lat, lon)
if not pts:
    sys.exit('GEGEGE no points in ' + llfile)

#[2]## real/imaginary grids in the region of the points (HAMTIDE longitude is 0-360)
lons = [lon % 360.0 for _, lon in pts]
lats = [lat for lat, _ in pts]
w = max(0.0, math.floor(min(lons)) - 1); e = min(360.0, math.ceil(max(lons)) + 1)
s = max(-90.0, math.floor(min(lats)) - 1); n = min(90.0, math.ceil(max(lats)) + 1)
region = '-R%g/%g/%g/%g' % (w, e, s, n)

tmp = tempfile.mkdtemp(prefix='hamtide2uv_', dir='.')
ptsfile = os.path.join(tmp, 'pts.txt')
with open(ptsfile, 'w') as f:
    for lat, lon in pts:
        f.write('%.10f %.10f\n' % (lon % 360.0, lat))

def gmt(*args):
    r = subprocess.run(['gmt'] + list(args), capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit('GEGEGE gmt %s failed:\n%s' % (args[0], r.stderr))
    return r.stdout

vals = {}
for comp, amp, pha in (('u', 'UAMP', 'UPHA'), ('v', 'VAMP', 'VPHA')):
    for part, op in (('re', 'COSD'), ('im', 'SIND')):
        grd = os.path.join(tmp, '%s_%s.nc' % (comp, part))
        gmt('grdmath', region, '%s?%s' % (ncfile, amp), '%s?%s' % (ncfile, pha), op, 'MUL', '=', grd)
        out = gmt('grdtrack', ptsfile, '-G' + grd, '-nl')   # bilinear; NaN if a cell has no data
        vals[(comp, part)] = [float(l.split()[2]) for l in out.splitlines() if l.strip()]
        if len(vals[(comp, part)]) != len(pts):
            sys.exit('GEGEGE grdtrack returned %d values for %d points' % (len(vals[(comp, part)]), len(pts)))

#[3]## write OTPS-like files
for comp, fname, label in (('u', uout, 'WE velocity  (cm/s)'), ('v', vout, 'SN velocity  (cm/s)')):
    nmiss = 0
    with open(fname, 'w') as f:
        f.write(' Model:        HAMTIDE11a (%s)\n' % os.path.basename(ncfile))
        f.write(' %s\n' % label)
        f.write('   Lat     Lon       %s_amp  %s_ph   \n' % (const, const))
        for (lat, lon), re_, im_ in zip(pts, vals[(comp, 're')], vals[(comp, 'im')]):
            if math.isnan(re_) or math.isnan(im_):
                nmiss += 1
                f.write('%9.3f %9.3f       ************* Site is out of model grid OR land ***************\n' % (lat, lon))
            else:
                amp = math.hypot(re_, im_)
                ph = math.degrees(math.atan2(im_, re_)) % 360.0
                f.write('%9.3f %9.3f %8.3f %8.1f\n' % (lat, lon, amp, ph))
    print('%s: %d points, %d out of model grid OR land' % (fname, len(pts), nmiss))

for fn in os.listdir(tmp):
    os.remove(os.path.join(tmp, fn))
os.rmdir(tmp)
