2726516565436524 1346563466354361 irand_seed
hci                               run_type
1e-3  1e-7      1.e-4   1         eps_var, eps_pt, pt_error, n_states
f                                 dump_wf_var
'chem'  1                         hamiltonian_type,ipr
8        4                        nelec, nup
d2h                               point_group
.true.                            time_sym
-1                                z
26                                norb
1,5,3,2,1,7,6,5,1,2,3,5,1,6,7,5,3,2,1,4,8,5,1,7,6,5, orbital_symmetries(1:norb)
1                                 spatial_symmetry_wf
1                                 diagonalize_ham
&selected_ci eps_var_sched=2*2e-3 eps_pt_big=0.e-5 n_mc=0 /
!&hf_det lz=1, u=T/
&hf_det irreps=1,2,3,5 irrep_occs_up=2,1,0,1 irrep_occs_dn=2,0,1,1 /