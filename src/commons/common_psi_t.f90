module common_psi_t
  use types, only: rk, ik, ik_vec
  implicit none ; save
  integer :: ndeg ! number of degenerate HF ground states (for hubbard) 
  integer :: ndet_psi_t
  integer :: ndet_psi_t_in
  logical :: use_psit_in,use_psit_out,print_psit_wo_sqmc
  logical :: use_psit_con_in,use_psit_con_out
  logical :: use_elems_in,use_elems_out
  character(len=100) :: psit_in_file="psit.in",psit_out_file="psit.out"
  character(len=100) :: psit_con_in_file="psit_connections.in",psit_con_out_file="psit_connections.out"
  character(len=100) :: dtm_elems_in_file="dtm_projector.in",dtm_elems_out_file="dtm_projector.out"
  real(rk) :: dtm_energy ! deterministic space energy (only for printing into dtm_elems_out_file)
  real(rk) :: e_trial
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), allocatable :: dets_up_psi_t(:), dets_dn_psi_t(:)
  type(ik_vec), allocatable :: dets_up_psi_t_in(:), dets_dn_psi_t_in(:)
  type(ik_vec),allocatable :: psi_t_connected_dets_up(:),psi_t_connected_dets_dn(:)
#else
  integer(ik), allocatable :: dets_up_psi_t(:), dets_dn_psi_t(:)
  integer(ik), allocatable :: dets_up_psi_t_in(:), dets_dn_psi_t_in(:)
  integer(ik),allocatable :: psi_t_connected_dets_up(:),psi_t_connected_dets_dn(:)
#endif
  integer, allocatable :: iwdet_psi_t(:)
  real(rk), allocatable :: cdet_psi_t(:), psi_g(:),real_wt(:),cdet_psi_t_in(:)
  real(rk) :: psi_g_energy,psi_g_epsilon,psi_g_epsilon_inv ! parameters used in guiding wave function
  integer :: trial_wf_iters ! number of times find connected dets and truncation are performed
  integer,allocatable :: norb_trial_wf(:),n_initiators_trial_wf(:), n_truncate_trial_wf(:) ! only used in construction of psi t

  integer :: ndet_psi_t_connected
  real(rk),allocatable :: psi_t_connected_e_loc_num(:),psi_t_connected_e_loc_den(:)
  logical :: hf_to_psit ! if true, use psi trial instead of HF for the 1st state
end module common_psi_t
