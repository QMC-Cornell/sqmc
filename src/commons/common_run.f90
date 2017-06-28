module common_run
  use types, only: rk, ik, ik_vec
  use, intrinsic :: iso_c_binding, only: c_ptr
  implicit none ; save
  integer nstep, nblk, nblk_eq, ipr, nwalk, importance_sampling!, n_excite
  real(rk) tau, tau_multiplier, population_control_exponent, reweight_factor_inv_max_multiplier, reweight_factor_inv_max, initiator_rescale_power,partial_node_eps !, e_excite
  character*16 proposal_method
  character*16 run_type
  integer :: max_connected_dets
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),allocatable :: connected_dets_up(:),connected_dets_dn(:)
#else
  integer(ik),allocatable :: connected_dets_up(:),connected_dets_dn(:)
#endif
  integer :: n_connected_dets_hf
  real(rk),allocatable :: connected_matrix_elements(:),connected_matrix_elements_fn(:)
  logical :: first_time = .true.
  integer :: debug_counter = 0
  integer :: ndet_outside_ct=0 ! number of determinants outside C(T) at the end of an iteration (these are already sorted by label, so no need to re-sort them)
  ! Following line is for interface in to C++ MPS code
  type(c_ptr) :: h_psi_pointer,psi_pointer,qp_pointer,ps_pointer
  integer,allocatable :: occ(:) ! For the occupation number representation used by MPS code
  logical :: use_efficient_heatbath

  ! For quickly calculating diagonal elements in O(N) time rather than naive O(N^2) time
  type diag_elem_info
    real(rk) :: old_diag_elem
    integer :: p,q,r,s
  end type diag_elem_info

  type(diag_elem_info),allocatable :: connected_diag_elems_info(:)
  integer, parameter :: input_copy_unit = 111

end module common_run
