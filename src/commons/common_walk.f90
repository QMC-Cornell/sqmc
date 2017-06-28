module common_walk
! Note MWALK is 8-byte but nwalk is 4-byte because we calculate MWALK and stop if it is larger than can be stored in 4-byte integer = 2^31-1
  use types, only: rk, ik, ik_vec, i1b, i8b
  implicit none ; save
  integer(i8b) :: MWALK
  integer nwalk, n_connected_dets, m_connected_dets
  integer, allocatable :: walk_det(:), connected_dets(:) ! These are for "toy" problems, i.e., where the projection matrix can be stored
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), allocatable :: walk_dets_up(:), walk_dets_dn(:)!, connected_dets_up(:), connected_dets_dn(:)
#else
  integer(ik), allocatable :: walk_dets_up(:), walk_dets_dn(:)!, connected_dets_up(:), connected_dets_dn(:)
#endif
  integer(i1b), allocatable :: imp_distance(:)
  integer, allocatable ::  connected_dets_sign(:)
! integer, allocatable ::  walk_wt(:), connected_dets_sign(:)
  real(rk), allocatable ::  walk_wt(:),walk_wt_fn(:)
  real(rk), allocatable ::  matrix_elements(:)
  real(rk), allocatable ::  e_num_walker(:),e_den_walker(:)
end module common_walk
