module common_ham
  use types, only: rk
  implicit none ; save
  integer nelec, nup, ndn, norb, n_core_orb, ndet, diagonalize_ham
  real(rk) hf_energy, max_energy, energy_exact, diagonal_ham_lowest, diagonal_ham_highest
  character*16 hamiltonian_type
end module common_ham
