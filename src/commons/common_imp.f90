module common_imp
  use types, only: rk, ik, ik_vec, i8b
  implicit none ; save
  logical :: semistochastic ! true if semistochastic FCI-QMC
  integer :: imp_iters
  integer,allocatable :: norb_imp(:), n_imp_initiators(:), n_imp_truncate(:)
  integer :: n_imp
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),allocatable :: imp_up(:),imp_dn(:)
#else
  integer(ik),allocatable :: imp_up(:),imp_dn(:)
#endif
  integer(i8b),allocatable :: minus_tau_H_indices(:),minus_tau_H_nonzero_elements(:)
  real(rk),allocatable :: minus_tau_H_values(:) ! Sparse pieces of projector matrix
  ! These are for generating the deterministic space and psi trial simultaneously:
  logical :: diff_from_psi_t ! if false, generate deterministic space same as psi trial
  integer :: size_deterministic ! if above false, truncate to this many CSFs on the final step
end module common_imp
