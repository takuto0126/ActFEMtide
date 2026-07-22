src="../mkfvxyz/src"  # 2024.03.07 Tatkuo Minami
OTPS="../mkfvxyz/OTPS"
igrf="../mkfvxyz/igrf"

cd $src  # 2024.03.07 Tatkuo Minami
make clean
make
cd -
cd $OTPS # 2024.03.07 Tatkuo Minami
make clean
make
cd -

fldr=/home/nonoyama/tide/Kikai/mesh_light # 2024.03.07 Tatkuo Minami
meshctl=${fldr}/mesh_kikai.ctl            # 2024.03.07 Tatkuo Minami
oceanmesh=${fldr}/ocean.msh               # 2024.03.07 Tatkuo Minami
em3dmesh=${fldr}/em3d.msh                 # 2024.03.07 Tatkuo Minami

#[1]## generate lat_lon_mesh
${src}/mkmshlonlat.exe <<EOF
"${meshctl}"
"${oceanmesh}"
"${em3dmesh}"
EOF

#[2]## generate fxyz_mesh
${src}/meshfxyz.exe # input is lat_lon_mesh

#[3]## generate vxyz_mesh
#[3-1]## generate tide model list
cat <<EOF > Model_tpxo9_atlas_v3
${OTPS}/DATA/TPXO9_atlas_v3_test/h_*_tpxo9_atlas_30_v3
${OTPS}/DATA/TPXO9_atlas_v3_test/u_*_tpxo9_atlas_30_v3
${OTPS}/DATA/TPXO9_atlas_v3_test/grid_tpxo9_atlas_30_v3
EOF

#[3-2]## ge
${OTPS}/extract_HC << finish
Model_tpxo9_atlas_v3
lat_lon_mesh
u
m2
AP
oce
1
u_xyz_mesh
finish

${OTPS}/extract_HC << finish
Model_tpxo9_atlas_v3
lat_lon_mesh
v
m2
AP
oce
1
v_xyz_mesh
finish

${src}/mkvxyzfile.exe # input is u_vxyz_mesh and v_vxyz_mesh

exit

gfortran m_readfvxyz.f90
./a.out



