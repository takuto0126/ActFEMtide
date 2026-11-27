#!/usr/bin/env python3
# Coded on 2026.11.27 by Claude Code (CC) for ActFEMtide
# Fill the points marked "Site is out of model grid OR land" in u_xyz_mesh / v_xyz_mesh
# (output of OTPS extract_HC or hamtide2uv.py) with the amplitude and phase of the
# nearest valid point in the same file (great-circle distance).
# Such points are ocean nodes of the mesh that fall on land or outside the tidal model
# grid (mainly near coastlines). Without filling, mkvxyzfile.exe sets their velocity to 0.
# The original files are kept as <file>.raw.
# Usage: fill_nearest_uv.py <u_xyz_mesh> [<v_xyz_mesh> ...]
import math
import shutil
import sys

def read(fname):
    head, rows = [], []
    with open(fname) as f:
        lines = f.read().splitlines()
    head = lines[:3]
    for l in lines[3:]:
        if not l.strip():
            continue
        a = l.split()
        lat, lon = float(a[0]), float(a[1])
        try:
            amp, ph = float(a[2]), float(a[3])
            ok = True
        except (IndexError, ValueError):
            amp, ph, ok = 0.0, 0.0, False
        rows.append([lat, lon, amp, ph, ok])
    return head, rows

def unit(lat, lon):
    la, lo = math.radians(lat), math.radians(lon)
    return (math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo), math.sin(la))

for fname in sys.argv[1:]:
    head, rows = read(fname)
    valid = [r for r in rows if r[4]]
    miss = [r for r in rows if not r[4]]
    if not valid:
        sys.exit('GEGEGE no valid points in ' + fname)
    shutil.copyfile(fname, fname + '.raw')
    vu = [unit(r[0], r[1]) for r in valid]
    dmax = 0.0
    for r in miss:
        p = unit(r[0], r[1])
        best, bi = -2.0, 0
        for i, q in enumerate(vu):       # largest dot product = smallest distance
            d = p[0] * q[0] + p[1] * q[1] + p[2] * q[2]
            if d > best:
                best, bi = d, i
        r[2], r[3] = valid[bi][2], valid[bi][3]
        dmax = max(dmax, 6371.0 * math.acos(max(-1.0, min(1.0, best))))
    with open(fname, 'w') as f:
        for h in head:
            f.write(h + '\n')
        for r in rows:
            f.write('%9.3f %9.3f %8.3f %8.1f\n' % (r[0], r[1], r[2], r[3]))
    print('%s: %d points, %d filled with the nearest valid point (max distance %.1f km)'
          % (fname, len(rows), len(miss), dmax))
