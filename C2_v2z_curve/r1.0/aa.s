for file in i*
do
#echo $file
ex $file <<_end
1,\$s/.*orbital_symmetries(1:norb)/xx/
wq
_end
done

# 1,\$s/d2h                               point_group/dih                               point_group/
#.s;.*;\&hf_det irreps=1,2,3,4 irrep_occs_up=3,1,1,1 irrep_occs_dn=3,1,1,1 /;
