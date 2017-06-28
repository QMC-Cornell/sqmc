module semistoch

use tools,  only : merge_sort2_up_dn, merge2_up_dn, merge_original_with_spawned3, alloc ! Note: this merge_original_with_spawned3 is very different from the routine in do_walk.f90!
use more_tools, only : print_int_matrix,print_real_matrix, matrix_lanczos
use hubbard, only : hamiltonian_hubbard_k,hamiltonian_hubbard_k_space_sym,k_vectors,ktot,l_x,l_y,nsites,k_energies,k_hf_deg_up,k_hf_deg_dn,ndeg,energy_pieces_hubbard_k,check_momentum_conservation,generate_sparse_ham_hubbardk,space_sym,symmetry_reduce_and_replace_hubbardk
use chemistry, only : hamiltonian_chem,hamiltonian_chem_time_sym,norb,energy_pieces_chem,filter_dets,n_core_orb,spatial_symmetry_wf,generate_sparse_ham_chem
use heg, only : hamiltonian_heg, generate_sparse_ham_heg
use types, only : ik, ik_vec, i8b, rk
#ifdef NUM_ORBITALS_GT_127
use overload
#endif
use common_ham, only : nelec,nup, ndn,diagonalize_ham,hamiltonian_type,norb
use common_run, only: run_type,tau,connected_dets_up,connected_dets_dn,connected_matrix_elements,max_connected_dets,diag_elem_info,connected_diag_elems_info
use common_psi_t, only : dets_up_psi_t,dets_dn_psi_t,ndet_psi_t,e_trial,cdet_psi_t
use common_walk, only : MWALK ! For stopping after generating final trial wavefn if MWALK < 0
use common_imp, only : n_imp, norb_imp, imp_up, imp_dn, diff_from_psi_t, size_deterministic

implicit none
save
private
public :: generate_psi_t_connected_e_loc, generate_space_iterate, perform_selected_ci, perform_truncated_lanczos, hamiltonian, find_doubly_excited, estimate_n_connections

contains

!------------------------------------------------------------------------------------------------------------------------------

  subroutine generate_psi_t_connected_e_loc(ndet_psi_t_connected,psi_t_connected_dets_up,psi_t_connected_dets_dn,psi_t_connected_e_loc_num,psi_t_connected_e_loc_den)

   ! Saves the local energies for all states connected to psi trial by the Hamiltonian.

      use common_run, only: ipr,importance_sampling
      use common_psi_t, only : cdet_psi_t,use_psit_con_out,psit_con_out_file,print_psit_wo_sqmc,use_elems_out
      use types, only : num_words
      use common_ham, only : n_core_orb
      use mpi_routines, only : master_core

      implicit none

      integer, intent(out) :: ndet_psi_t_connected
#ifdef NUM_ORBITALS_GT_127
      type(ik_vec),allocatable,intent(out) :: psi_t_connected_dets_up(:),psi_t_connected_dets_dn(:)
      type(ik_vec) :: tmp_det
#else
      integer(ik),allocatable,intent(out) :: psi_t_connected_dets_up(:),psi_t_connected_dets_dn(:)
      integer(ik) :: tmp_det
#endif
      real(rk),allocatable,intent(out) :: psi_t_connected_e_loc_num(:),psi_t_connected_e_loc_den(:)
      integer :: norb_include
      integer :: i,j
      character(len=8) :: fmt
      integer :: orb_up(nup),orb_dn(ndn)
      integer :: ndet_connections_nonzero

      call my_second(1,'generate energy pieces') ; write(6,*); call flush(6)

      if (hamiltonian_type.eq.'chem')  norb_include=norb
      if (hamiltonian_type.eq.'hubbardk')  norb_include=nsites
      if (hamiltonian_type.eq.'heg')  norb_include=norb

      write (6,*) "generating local energies; importance_sampling=",importance_sampling; call flush(6)
      call find_doubly_excited(ndet_psi_t_connected,psi_t_connected_dets_up,psi_t_connected_dets_dn,dets_up_psi_t(1:ndet_psi_t),dets_dn_psi_t(1:ndet_psi_t),norb_include,n_core_orb,ref_coeffs=cdet_psi_t(1:ndet_psi_t),e_mix_num=psi_t_connected_e_loc_num,e_mix_den=psi_t_connected_e_loc_den,importance_sampling=importance_sampling)
      ! importance_sampling is a global variable stored in common_run

      write (6,'(''Number of saved local energies, ndet_psi_t_connected='',i10,'' ='',es11.4)') ndet_psi_t_connected, real(ndet_psi_t_connected)
      call flush(6)

      write (6,*) "First 1000 (sorted by determinant label):"
      write (6,*) "  Det up     Det dn      E_Num        E_Den           E_L"
      write (fmt, '(i2)') 2*num_words
      do i=1,min(size(psi_t_connected_e_loc_num),1000)
        write(6,'(' // trim(fmt) // 'i10 f15.5 f15.5 f15.5)') psi_t_connected_dets_up(i),psi_t_connected_dets_dn(i),psi_t_connected_e_loc_num(i),psi_t_connected_e_loc_den(i),psi_t_connected_e_loc_num(i)/psi_t_connected_e_loc_den(i)
      enddo

      call my_second(2,'generate energy pieces') ; write(6,*)
      call flush(6)

      ! Write out connections to Psi_T
      ! If core orbitals are used, the labels of the filled orbitals start with the first non-core orbital as number 1
      if (use_psit_con_out) then ! if print Psi_T connections out into an output file
        if (master_core) then
          write (6,*) "Dumping the local energies into file ",psit_con_out_file
          call my_second(1,'dump local energies'); call flush(6)
          open(8, file=psit_con_out_file,status='new')
          ndet_connections_nonzero = ndet_psi_t_connected
          do i=1,ndet_psi_t_connected ! don't bother printing out dets with numerators < 1e-10
            if (.not. (abs(psi_t_connected_e_loc_num(i))>1.e-10_rk)) then
              ndet_connections_nonzero=ndet_connections_nonzero-1
            endif
          enddo
          if (hamiltonian_type.eq.'chem') then
            write(8,'(i8,i12,2i4,i6,i3,f15.8,'' ndet_psi_t, ndet_connections_nonzero, nup-n_core_orb, ndn-n_core_orb, norb, 0, E_T'')') &
         &  ndet_psi_t, ndet_connections_nonzero, nup-n_core_orb, ndn-n_core_orb, norb, 0, psi_t_connected_e_loc_num(1)/psi_t_connected_e_loc_den(1)
          elseif (hamiltonian_type.eq.'hubbardk') then
            write(8,'(i8,i12,2i4,i6,i3,f15.8,'' ndet_psi_t, ndet_connections_nonzero, nup-n_core_orb, ndn-n_core_orb, nsites, 0, E_T'')') &
         &  ndet_psi_t, ndet_connections_nonzero, nup-n_core_orb, ndn-n_core_orb, nsites, 0, psi_t_connected_e_loc_num(1)/psi_t_connected_e_loc_den(1)
          endif
          write(8,'(''orb_up          orb_dn           e_loc_num            e_loc_den'')')
          write (fmt, '(i5)') nup+ndn-2*n_core_orb-1
          do i = 1, ndet_psi_t_connected
            if (abs(psi_t_connected_e_loc_num(i))>1.e-10_rk) then
#ifdef NUM_ORBITALS_GT_127
              tmp_det = psi_t_connected_dets_up(i)-maskr_vec(n_core_orb)
#else
              tmp_det = psi_t_connected_dets_up(i)-maskr(n_core_orb,ik)
#endif

              do j=n_core_orb+1,nup
                orb_up(j) = trailz(tmp_det)+1
                tmp_det = ibclr(tmp_det,orb_up(j)-1)
              enddo
#ifdef NUM_ORBITALS_GT_127
              tmp_det = psi_t_connected_dets_dn(i)-maskr_vec(n_core_orb)
#else
              tmp_det = psi_t_connected_dets_dn(i)-maskr(n_core_orb,ik)
#endif

              do j=n_core_orb+1,ndn
                orb_dn(j) = trailz(tmp_det)+1
                tmp_det = ibclr(tmp_det,orb_dn(j)-1)
              enddo
              do j=n_core_orb+1,nup
                orb_up(j) = orb_up(j) - n_core_orb
              enddo
              do j=n_core_orb+1,ndn
                orb_dn(j) = orb_dn(j) - n_core_orb
              enddo
              write(8,'(i3,' // trim(fmt) // 'i4,f22.15,f19.15)') orb_up(n_core_orb+1:nup), orb_dn(n_core_orb+1:ndn), psi_t_connected_e_loc_num(i), psi_t_connected_e_loc_den(i)
            endif
          enddo
          call my_second(2,'dump local energies'); call flush(6)
        endif
      endif
      close(8)

      if (print_psit_wo_sqmc.and..not.use_psit_con_out.and..not.use_elems_out) then
        write (6,*) "Printed trial wave function connections in ",psit_con_out_file
        write (6,*) "Terminating SQMC"
        call flush(6)
        stop "Printed trial wave function connections. Terminating SQMC."
      endif

  end subroutine generate_psi_t_connected_e_loc
!-----------------------------------------------------------------------------------------------------------------------------------

  subroutine generate_space_iterate(iters,norb_include,n_initiators,n_truncate,n_det,dets_up,dets_dn,values,H_indices,H_nonzero_elements,e_trial,init_up,init_dn,init_wts,dtm_energy)

    ! Generates space as follows. Starting with HF, find connected dets out of norb_include(1) orbitals, then truncates to the n_sym_uniq_det
    ! highest weight csfs out of those connections, then repeats, using that truncated list as the reference, etc. iters is the number of times
    ! this process is repeated, and norb_include and n_sym_uniq_det are vectors of length iters corresponding to the norb_include and
    ! n_sym_uniq_det values at each iteration. Finally computes the sparse Hamiltonian in this space.
    ! Rediagonalize in the final space.
    ! Used for both important space and trial wave function generation.
    ! A Holmes, 30 Apr 2012
    ! Modified by A Holmes, 7 Nov 2012. Generate deterministic space H using guiding wave function (for importance sampling).
    ! Modified by A Holmes, 14 Feb 2013. Use previous iteration's (truncated) wts as starting vector for Lanczos.
    ! Modified by A Holmes, 23 Jul 2013. Optional inputs init_up,init_dn,init_wts allow you to use a different
    !                                    vector than HF for the first iteration. Also, if iters=0, then just
    !                                    diagonalize in this space.
    ! Modified by A Holmes, 3 Feb 2015. Use perturbation theory to reduce number of determinants at each step before rediagonalizing. (Default: reduce by a factor of 10)

      use common_run, only : ipr,tau_multiplier,tau,importance_sampling
      use chemistry, only  : diagonalize_sparse_hamiltonian_chem,time_sym,is_connected_chem,norb,find_connected_dets_chem,generate_sparse_ham_chem_upper_triangular
      use hubbard, only    : diagonalize_sparse_hamiltonian_hubbard,space_sym,is_connected_hubbard_fast,nsites,find_connected_dets_hubbard_k,c_sym_psi_t,generate_sparse_ham_hubbardk_upper_triangular
      use heg, only        : diagonalize_sparse_hamiltonian_heg,find_connected_dets_heg
      use more_tools, only : real_symmetric_diagonalize,binary_search,linear_search_list
      use generic_sort, only : sort
      use common_imp, only : minus_tau_H_indices,minus_tau_H_nonzero_elements,minus_tau_H_values
      use common_psi_t, only : hf_to_psit,psit_out_file,print_psit_wo_sqmc,use_psit_out,use_psit_con_in,use_psit_con_out,use_elems_out
      use types, only : num_words,i8b

      implicit none

      integer,intent(in) :: iters
      integer,intent(in) :: norb_include(:)
      integer,intent(inout) :: n_initiators(:) ! find connections to this many dets each iteration
      integer,intent(in) :: n_truncate(:) ! keep this many dets each iteration
      integer,intent(out) :: n_det
      real(rk),allocatable,intent(out) :: values(:) ! if the following inputs are present, this is nonzero elements of sparse Hamiltonian; else, it is lowest eigenvector of sparse Hamiltonian.
      integer(i8b),allocatable,optional,intent(out) :: H_indices(:),H_nonzero_elements(:)
      real(rk),optional,intent(out) :: e_trial
      real(rk),optional,intent(in) :: init_wts(:)
#ifdef NUM_ORBITALS_GT_127
      type(ik_vec),allocatable,intent(out) :: dets_up(:),dets_dn(:)
      type(ik_vec),optional,intent(in) :: init_up(:),init_dn(:)
      type(ik_vec),allocatable :: tmp_up(:),tmp_dn(:)
      type(ik_vec),allocatable :: tmp_up2(:),tmp_dn2(:)
      type(ik_vec), allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
      type(ik_vec) :: tmp_det
#else
      integer(ik),allocatable,intent(out) :: dets_up(:),dets_dn(:)
      integer(ik),optional,intent(in) :: init_up(:),init_dn(:)
      integer(ik),allocatable :: tmp_up(:),tmp_dn(:)
      integer(ik),allocatable :: tmp_up2(:),tmp_dn2(:)
      integer(ik), allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
      integer(ik) :: tmp_det
#endif
      real(rk),optional,intent(out) :: dtm_energy
      integer :: i,j,k,n_csfs
      real(rk) :: lowest_eigenvalue
      real(rk),allocatable :: lowest_eigenvector(:)
      integer :: iter,n_det_tmp
      logical :: rediagonalize
      real(rk), parameter :: epsilon = 1.e-10_rk ! Although the wavefn coefs. may not be converged to this accuracy, the agreement between symmetry related coefs. exists to higher precision
      real(rk) :: previous_value
      real(rk) :: tau_out_for_diag_ham ! Just for the output of diagonalize_hamiltonian routine for the iterations other than the final one. Not used otherwise.
      integer,allocatable :: iorder(:),temp_i_2(:)
      integer             :: nnzero
      real(rk)            :: normalization
      integer             :: n_csfs_dtm
      real(rk)            :: sign_condition_number_numerator,sign_condition_number_denominator
      real(rk)            :: H_psi,H_psi_abs
      real(rk)            :: element
      logical             :: connected
      real(rk),allocatable :: lanczos_initial_vector(:)
      character(len=2)    :: fmt
      integer :: orb_up(nup),orb_dn(ndn)
      real(rk),allocatable :: new_wts(:)
      integer :: n_keep
      real(rk) :: h_ii
      integer,allocatable :: indices(:)

      rediagonalize=.true.

      if (.not.present(init_up)) then
        ! Use HF (or linear combination for hubbard if ground state is degenerate)
        n_det_tmp=1
        if (hamiltonian_type.eq.'hubbardk')  n_det_tmp = ndeg

        allocate(tmp_up(n_det_tmp))
        allocate(tmp_dn(n_det_tmp))
        if (hamiltonian_type.eq.'hubbardk') then
          tmp_up=k_hf_deg_up
          tmp_dn=k_hf_deg_dn
        else
          tmp_up(1)=2_ik**nup-1_ik
          tmp_dn(1)=2_ik**ndn-1_ik
        endif

        allocate(lowest_eigenvector(n_det_tmp))
        lowest_eigenvector(:) = 0._rk
        lowest_eigenvector(1) = 1._rk
        if (hamiltonian_type.eq.'hubbardk')  lowest_eigenvector=c_sym_psi_t
      endif

      if (iters>0) then

        do iter=1,iters ! iters = number of iterations of finding connections, diagonalizing, and truncating
          call my_second(1,'find_doubly_excited')
          call flush(6)
          ! Single application of H:
          if (iter==1.and.present(init_up)) then
            call find_doubly_excited(n_det_tmp,dets_up,dets_dn,init_up,init_dn,norb_include(iter),n_core_orb,init_wts,ninitiator=n_initiators(iter))
          else
            call find_doubly_excited(n_det_tmp,dets_up,dets_dn,tmp_up,tmp_dn,norb_include(iter),n_core_orb,lowest_eigenvector,ninitiator=n_initiators(iter),e_mix_num=new_wts)
          endif
          call my_second(2,'find_doubly_excited')
          call flush(6)

          ! Use perturbation theory to estimate the values of the new determinants and only use the 10% with largest expected absolute weight.
          if (.not.(iter==1.and.present(init_up))) then
            n_keep = size(tmp_up) + (n_det_tmp-size(tmp_up))/10 ! size = old + 10% of the change from old to new
            n_keep = max(n_keep,nint(1.5*n_truncate(iter)))
            if (iter==iters.and..not.diff_from_psi_t)  n_keep = max(n_keep,nint(1.5*size_deterministic))
            n_keep = min(n_keep,n_det_tmp)
            if (n_keep==n_det_tmp) then
              write (6,*) "Perturbation theory not needed to reduce connected determinants list"; call flush(6)
            else
              write (6,*) "First-order perturbation theory used to reduce",n_det_tmp,"dets to",n_keep; call flush(6)
              allocate(iorder(size(tmp_up)))
              allocate(temp_i16_up((size(tmp_up)+1)/2))
              allocate(temp_i16_dn((size(tmp_up)+1)/2))
              allocate(temp_i_2((size(tmp_up)+1)/2))
              do j=1,size(tmp_up)
                iorder(j)=j
              enddo
              call merge_sort2_up_dn(tmp_up(1:size(tmp_up)),tmp_dn(1:size(tmp_up)), iorder, size(tmp_up), temp_i16_up, temp_i16_dn, temp_i_2)
              lowest_eigenvector(1:size(tmp_up)) = lowest_eigenvector(iorder(1:size(tmp_up)))
              deallocate(iorder)
              deallocate(temp_i16_up)
              deallocate(temp_i16_dn)
              deallocate(temp_i_2)
              allocate(indices(size(dets_up)))
              call linear_search_list(dets_up,dets_dn,tmp_up,tmp_dn,indices)
              do i=1,n_det_tmp
               if (indices(i).ne.0) then ! new det in old det list
                 new_wts(i) = lowest_eigenvector(indices(i))
                !new_wts(i) = 10 ! This line forces all dets from last iteration to be in current iteration
               else
                 call hamiltonian(dets_up(i),dets_dn(i),dets_up(i),dets_dn(i),h_ii)
                 new_wts(i) = new_wts(i)/(lowest_eigenvalue-h_ii)
               endif
              enddo
              deallocate(indices)
              call sort(new_wts, dets_up(1:n_det_tmp), dets_dn(1:n_det_tmp))
              allocate(tmp_up2(n_keep))
              allocate(tmp_dn2(n_keep))
              tmp_up2(1:n_keep) = dets_up(1:n_keep)
              tmp_dn2(1:n_keep) = dets_dn(1:n_keep)
              deallocate(dets_up,dets_dn)
              allocate(dets_up(n_keep))
              allocate(dets_dn(n_keep))
              dets_up(1:n_keep) = tmp_up2(1:n_keep)
              dets_dn(1:n_keep) = tmp_dn2(1:n_keep)
              allocate(iorder(n_keep))
              allocate(temp_i16_up((n_keep+1)/2))
              allocate(temp_i16_dn((n_keep+1)/2))
              allocate(temp_i_2((n_keep+1)/2))
              do j=1,n_keep
                iorder(j)=j
              enddo
              call merge_sort2_up_dn(dets_up(1:n_keep),dets_dn(1:n_keep), iorder, n_keep, temp_i16_up, temp_i16_dn, temp_i_2)
              new_wts(1:n_keep) = new_wts(iorder(1:n_keep))
              deallocate(iorder)
              deallocate(temp_i16_up)
              deallocate(temp_i16_dn)
              deallocate(temp_i_2)
              deallocate(tmp_up2,tmp_dn2)
              n_det_tmp = n_keep
              deallocate(new_wts) ! Can make starting Lanczos vector later
            endif
          endif

          ! diagonalization
          call my_second(1,'diagonalize_sparse_hamiltonian in generate_space_iterate')
          call flush(6)

          if (.not.((hamiltonian_type.eq.'chem'.or.hamiltonian_type.eq.'hubbardk').and.iter>1)) then
            if (.not.present(init_up)) then
              deallocate(tmp_up,tmp_dn)
              deallocate(lowest_eigenvector)
            endif
            allocate(lowest_eigenvector(n_det_tmp))
          endif

          if (hamiltonian_type.eq.'hubbardk') then ! Need to call different routines because within the diagonalize routine they call different Hamiltonians.  Could be combined.
            if (iter>1) then
              ! Seed with last iteration's values.
              ! Note: lowest_eigenvector,tmp_up,tmp_dn (from last iteration) are currently sorted by abs wt,
              ! whereas dets_up,dets_dn (from this iteration) are sorted by label!
              allocate(lanczos_initial_vector(n_det_tmp))
              lanczos_initial_vector(:) = 0._rk
              do i=1,size(tmp_up)
                call binary_search(tmp_up(i),tmp_dn(i),dets_up(1:n_det_tmp),dets_dn(1:n_det_tmp),k)
                if (k>0)  lanczos_initial_vector(k)=lowest_eigenvector(i)
              enddo
              deallocate(lowest_eigenvector)
              allocate(lowest_eigenvector(n_det_tmp))
              call diagonalize_sparse_hamiltonian_hubbard(dets_up(1:n_det_tmp),dets_dn(1:n_det_tmp),lowest_eigenvector,lowest_eigenvalue,tau_out_for_diag_ham,initial_vector=lanczos_initial_vector)
              deallocate(lanczos_initial_vector)
              deallocate(tmp_up,tmp_dn)
            else
              call diagonalize_sparse_hamiltonian_hubbard(dets_up(1:n_det_tmp),dets_dn(1:n_det_tmp),lowest_eigenvector,lowest_eigenvalue,tau_out_for_diag_ham)
            endif
           !call diagonalize_sparse_hamiltonian_hubbard(dets_up(1:n_det_tmp),dets_dn(1:n_det_tmp),lowest_eigenvector,lowest_eigenvalue,tau_out_for_diag_ham)
          elseif (hamiltonian_type.eq.'chem') then
            if (iter>1) then
              ! Seed with last iteration's values.
              ! Note: lowest_eigenvector,tmp_up,tmp_dn (from last iteration) are currently sorted by abs wt,
              ! whereas dets_up,dets_dn (from this iteration) are sorted by label!
              allocate(lanczos_initial_vector(n_det_tmp))
              lanczos_initial_vector(:) = 0._rk
              do i=1,size(tmp_up)
                call binary_search(tmp_up(i),tmp_dn(i),dets_up(1:n_det_tmp),dets_dn(1:n_det_tmp),k)
                if (k>0)  lanczos_initial_vector(k)=lowest_eigenvector(i)
              enddo
              deallocate(lowest_eigenvector)
              allocate(lowest_eigenvector(n_det_tmp))
              call diagonalize_sparse_hamiltonian_chem(dets_up(1:n_det_tmp),dets_dn(1:n_det_tmp),n_det_tmp,lowest_eigenvector,lowest_eigenvalue,tau_out_for_diag_ham,.true.,initial_vector=lanczos_initial_vector)
              deallocate(lanczos_initial_vector)
              deallocate(tmp_up,tmp_dn)
            else
              call diagonalize_sparse_hamiltonian_chem(dets_up(1:n_det_tmp),dets_dn(1:n_det_tmp),n_det_tmp,lowest_eigenvector,lowest_eigenvalue,tau_out_for_diag_ham,.true.)
            endif
          elseif (hamiltonian_type.eq.'heg') then
            call diagonalize_sparse_hamiltonian_heg(dets_up(1:n_det_tmp),dets_dn(1:n_det_tmp),n_det_tmp,lowest_eigenvector,lowest_eigenvalue)
          else
            stop 'This hamiltonian type has not been coded (only hubbardk, chem, and heg)'
          endif
          ! lowest_eigenvector is now sorted according to abs wt
          call my_second(2,'diagonalize_sparse_hamiltonian in generate_space_iterate')
          call flush(6)

          !pick the number of desired unique (by symmetry) determinants
          !Actually for molecules (unlike the heg) this does not guarantee it because the coefficients of the determinants in a CSF need not be equal in absolute magnitude.

          n_det = 0
          n_csfs = 0
          previous_value = 0._rk

          do i = 1, n_det_tmp
            if ( abs(abs(previous_value) - abs(lowest_eigenvector(i))) > epsilon ) then
              n_csfs=n_csfs+1
              previous_value = lowest_eigenvector(i)
              if (i>n_truncate(iter)) then
                n_csfs=n_csfs-1
                n_det = i-1
                exit
              endif
            endif
          enddo
          if (n_det.eq.0)  n_det=n_det_tmp

          if (iter==iters.and..not.diff_from_psi_t) then
            ! On last iteration truncate differently for deterministic space than for psi trial
            n_imp = 0
            n_csfs_dtm = 0
            previous_value = 0.d0
            do i = 1, n_det_tmp
              if ( abs(abs(previous_value) - abs(lowest_eigenvector(i))) > epsilon ) then
                n_csfs_dtm=n_csfs_dtm+1
                previous_value = lowest_eigenvector(i)
                if (i>size_deterministic) then
                  n_csfs_dtm=n_csfs_dtm-1
                  n_imp = i-1
                  exit
                endif
              endif
            enddo
            if (n_imp.eq.0)  n_imp=n_det_tmp
          endif

          normalization = 1 / sqrt(dot_product(lowest_eigenvector(1:n_det),lowest_eigenvector(1:n_det)))
          lowest_eigenvector(:) = lowest_eigenvector(:)*normalization

          if (diff_from_psi_t.or.iter<iters) then
            if (present(H_indices)) then
              write(6,'(/,''Deterministic space after iteration '',i2,'' has:'',i8,'' dets., with energy (before rediag)'',f12.6)') iter, n_det_tmp, lowest_eigenvalue
              write(6,'(/,''Truncated Deterministic space after iteration '',i2,'' has: '',i8,'' CSFs,'',i8,'' dets.'')') iter, n_csfs, n_det
            else
              write(6,'(/,''Trial wavefunction after iteration '',i2,'' has:'',i8,'' dets., with energy (before rediag)'',f12.6)') iter, n_det_tmp, lowest_eigenvalue
              write(6,'(/,''Truncated Trial wavefunction after iteration '',i2,'' has: '',i8,'' CSFs,'',i8,'' dets.'')') iter, n_csfs, n_det
            endif
            write(6,'(''Largest'',i6,'' renormalized abs coefs. are:'')') min(n_det,100)
            write(6,*) "        up-det         dn-det"
            write (fmt, '(i2)') 2*num_words
            do i = 1, min(n_det,100)
              write(6,'(' // trim(fmt) // 'i22,f18.14)') dets_up(i), dets_dn(i), lowest_eigenvector(i)
            enddo
          else
            ! generating simultaneously
            write(6,'(/,''Deterministic space after iteration '',i2,'' has:'',i8,'' dets., with energy (before rediag)'',f12.6)') iter, n_det_tmp, lowest_eigenvalue
            write(6,'(/,''Truncated Deterministic space after iteration '',i2,'' has: '',i8,'' CSFs,'',i8,'' dets.'')') iter, n_csfs_dtm, n_imp
            write(6,'(/,''Trial wavefunction after iteration '',i2,'' has:'',i8,'' dets., with energy (before rediag)'',f12.6)') iter, n_det_tmp, lowest_eigenvalue
            write(6,'(/,''Truncated Trial wavefunction after iteration '',i2,'' has: '',i8,'' CSFs,'',i8,'' dets.'')') iter, n_csfs, n_det
            write(6,'(''Largest'',i6,'' renormalized abs coefs. are:'')') min(max(n_det,n_imp),100)
            write(6,*) "        up-det         dn-det"
            write (fmt, '(i2)') 2*num_words
            do i = 1, min(max(n_det,n_imp),100)
              write(6,'(' // trim(fmt) // 'i22,f18.14)') dets_up(i), dets_dn(i), lowest_eigenvector(i)
            enddo
          endif
          call flush(6)

          ! Calculate sign condition number, sum_i(abs(sum_j(H_ij*c_j)))/sum_i(sum_j(abs(H_ij*c_j)))
          sign_condition_number_numerator = 0._rk
          sign_condition_number_denominator = 0._rk
          do i=1,n_det
            H_psi = 0._rk
            H_psi_abs = 0._rk
            do j=1,n_det
              call hamiltonian(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j),element,connected)
              if (connected) then
                H_psi = H_psi + element*lowest_eigenvector(j)
                H_psi_abs = H_psi_abs + abs(element*lowest_eigenvector(j))
              endif
            enddo
            sign_condition_number_numerator = sign_condition_number_numerator + abs(H_psi)
            sign_condition_number_denominator = sign_condition_number_denominator + H_psi_abs
          enddo

          write(6,'(/,''Sign condition number after iteration '',i2,'' is: '',f12.6)') iter, sign_condition_number_numerator/sign_condition_number_denominator

          if (iter<iters) then
            allocate(tmp_up(n_det))
            allocate(tmp_dn(n_det))
            tmp_up = dets_up(1:n_det)
            tmp_dn = dets_dn(1:n_det)
            deallocate(dets_up,dets_dn)
          endif

          ! Following code forces all connections to psi trial to be in deterministic space, so that walkers on psi trial never have to spawn (commented out because we've been having trouble with it)
         !if (iter==iters.and.present(H_indices)) then ! important space on final iteration
         !  allocate(tmp_up(n_imp))
         !  allocate(tmp_dn(n_imp))
         !  tmp_up = dets_up(1:n_imp)
         !  tmp_dn = dets_dn(1:n_imp)
         !  deallocate(dets_up,dets_dn)
         !  if (hamiltonian_type.eq.'hubbardk') then
         !    call find_doubly_excited(n_det=n_det,dets_up=dets_up,dets_dn=dets_dn,ref_up=tmp_up,ref_dn=tmp_dn,norb=nsites,n_core_orb=n_core_orb,ninitiator=1)
         !  else
         !    call find_doubly_excited(n_det=n_det,dets_up=dets_up,dets_dn=dets_dn,ref_up=tmp_up,ref_dn=tmp_dn,norb=norb,n_core_orb=n_core_orb,ninitiator=1)
         !  endif
         !if (allocated(tmp_up)) deallocate(tmp_up,tmp_dn)
         !endif
          ! End

          if (iter==iters)  deallocate(lowest_eigenvector)
         !if ((rediagonalize.or.iter<iters.or.present(H_indices)))  deallocate(lowest_eigenvector)

        enddo ! iter
      else
        ! copy over tmp_up,dn into dets_up,dn
        n_det = size(init_up)
        allocate(dets_up(n_det))
        allocate(dets_dn(n_det))
        dets_up(1:n_det)=init_up(1:n_det)
        dets_dn(1:n_det)=init_dn(1:n_det)
      endif

     !if (.not.diff_from_psi_t.and.iters>0) then
     !  ! Following code forces all connections to psi trial to be in deterministic space, so that dets on psi trial never have to spawn
     !  allocate(tmp_up(n_imp))
     !  allocate(tmp_dn(n_imp))
     !  tmp_up = dets_up(1:n_imp)
     !  tmp_dn = dets_dn(1:n_imp)
     !  ! up to now, dets_up is the same for dtm space and psit
     !  if (hamiltonian_type.eq.'hubbardk') then
     !    call find_doubly_excited(n_det=n_imp,dets_up=imp_up,dets_dn=imp_dn,ref_up=tmp_up,ref_dn=tmp_dn,norb=nsites,n_core_orb=n_core_orb,ninitiator=1)
     !  else
     !    call find_doubly_excited(n_det=n_imp,dets_up=imp_up,dets_dn=imp_dn,ref_up=tmp_up,ref_dn=tmp_dn,norb=norb,n_core_orb=n_core_orb,ninitiator=1)
     !  endif
     !  deallocate(tmp_up,tmp_dn)
     ! !n_det = n_imp
     !  ! End
     !elseif (.not.diff_from_psi_t) then
      if (.not.diff_from_psi_t) then
        !deallocate(imp_up,imp_dn)
        allocate(imp_up(n_imp))
        allocate(imp_dn(n_imp))
        imp_up(1:n_imp)=dets_up(1:n_imp)
        imp_dn(1:n_imp)=dets_dn(1:n_imp)
      endif
      if (size(dets_up).ne.n_det) then
        if (allocated(tmp_up))  deallocate(tmp_up,tmp_dn)
        allocate(tmp_up(n_det))
        allocate(tmp_dn(n_det))
        tmp_up(1:n_det)=dets_up(1:n_det)
        tmp_dn(1:n_det)=dets_dn(1:n_det)
        deallocate(dets_up,dets_dn)
        allocate(dets_up(n_det))
        allocate(dets_dn(n_det))
        dets_up(1:n_det)=tmp_up(1:n_det)
        dets_dn(1:n_det)=tmp_dn(1:n_det)
        deallocate(tmp_up,tmp_dn)
      endif

      ! Sort again, since the dets are currently in order of absolute magnitude of their coefficients.
      allocate(iorder(n_det))
      allocate(temp_i16_up((n_det+1)/2))
      allocate(temp_i16_dn((n_det+1)/2))
      allocate(temp_i_2((n_det+1)/2))
      do j=1,n_det
        iorder(j)=j
      enddo

      call merge_sort2_up_dn(dets_up,dets_dn, iorder, n_det, temp_i16_up, temp_i16_dn, temp_i_2)

      deallocate(iorder)
      deallocate(temp_i16_up)
      deallocate(temp_i16_dn)
      deallocate(temp_i_2)

      if (.not.diff_from_psi_t) then

        allocate(iorder(n_imp))
        allocate(temp_i16_up((n_imp+1)/2))
        allocate(temp_i16_dn((n_imp+1)/2))
        allocate(temp_i_2((n_imp+1)/2))
        do j=1,n_imp
          iorder(j)=j
        enddo

        call merge_sort2_up_dn(imp_up,imp_dn, iorder, n_imp, temp_i16_up, temp_i16_dn, temp_i_2)

        deallocate(iorder)
        deallocate(temp_i16_up)
        deallocate(temp_i16_dn)
        deallocate(temp_i_2)

      endif

      ! Finally, rediagonalize

      call my_second(1,'rediagonalize')
      call flush(6)

      if (present(H_indices)) then ! important space
        if (hamiltonian_type.eq.'hubbardk') then
          call generate_sparse_ham_hubbardk_upper_triangular(dets_up(1:n_det),dets_dn(1:n_det),H_indices,H_nonzero_elements,values,hf_to_psit,importance_sampling=(importance_sampling==1))
         !call generate_sparse_ham_hubbardk_upper_triangular(dets_up(1:n_det),dets_dn(1:n_det),H_indices,H_nonzero_elements,values,importance_sampling=(importance_sampling==1))
          ! Diagonalize just to get the energy in the deterministic space
          if (n_det>1) then
            allocate(lowest_eigenvector(n_det))
            call diagonalize_sparse_hamiltonian_hubbard(dets_up(1:n_det),dets_dn(1:n_det),lowest_eigenvector,lowest_eigenvalue,tau_out_for_diag_ham)
            n_csfs = count_csfs(n_det,lowest_eigenvector(1:n_det))
            deallocate(lowest_eigenvector)
          else
            if (space_sym) then
                call hamiltonian_hubbard_k_space_sym(dets_up(1),dets_dn(1),dets_up(1),dets_dn(1),lowest_eigenvalue,nnzero)
            else
                call hamiltonian_hubbard_k(dets_up(1),dets_dn(1),dets_up(1),dets_dn(1),lowest_eigenvalue)
            endif
          endif
        elseif (hamiltonian_type.eq.'chem') then
          call generate_sparse_ham_chem_upper_triangular(dets_up(1:n_det),dets_dn(1:n_det),H_indices,H_nonzero_elements,values,.true.,hf_to_psit)
         !call generate_sparse_ham_chem_upper_triangular(dets_up(1:n_det),dets_dn(1:n_det),H_indices,H_nonzero_elements,values,.true.)
          if (n_det>1) then
            allocate(lowest_eigenvector(n_det))
            call diagonalize_sparse_hamiltonian_chem(dets_up(1:n_det),dets_dn(1:n_det),n_det,lowest_eigenvector,lowest_eigenvalue,tau_out_for_diag_ham,.true.)
            n_csfs = count_csfs(n_det,lowest_eigenvector(1:n_det))
            !call diagonalize_sparse_hamiltonian(dets_up(1:n_det),dets_dn(1:n_det),n_det,lowest_eigenvector,lowest_eigenvalue,tau_out_for_diag_ham)
            ! Calculate sign condition number, sum_i(abs(sum_j(H_ij*c_j)))/sum_i(sum_j(abs(H_ij*c_j)))
            sign_condition_number_numerator = 0._rk
            sign_condition_number_denominator = 0._rk
            do i=1,n_imp
              H_psi = 0._rk
              H_psi_abs = 0._rk
              do j=1,n_imp
                call hamiltonian(imp_up(i),imp_dn(i),imp_up(j),imp_dn(j),element,connected)
                if (connected) then
                  H_psi = H_psi + element*lowest_eigenvector(j)
                  H_psi_abs = H_psi_abs + abs(element*lowest_eigenvector(j))
                endif
              enddo
              sign_condition_number_numerator = sign_condition_number_numerator + abs(H_psi)
              sign_condition_number_denominator = sign_condition_number_denominator + H_psi_abs
            enddo
            write(6,'(/,''Sign condition of final deterministic space is: '',f12.6)') sign_condition_number_numerator/sign_condition_number_denominator
            deallocate(lowest_eigenvector)
          else
            if (time_sym) then
              call hamiltonian_chem_time_sym(dets_up(1),dets_dn(1),dets_up(1),dets_dn(1),lowest_eigenvalue)
            else
              call hamiltonian_chem(dets_up(1),dets_dn(1),dets_up(1),dets_dn(1),0,lowest_eigenvalue)
            endif
          endif
        elseif (hamiltonian_type.eq.'heg') then
          write (6,*) "Generating sparse ham heg"
          call flush(6)
          call generate_sparse_ham_heg(dets_up(1:n_det),dets_dn(1:n_det),H_indices,H_nonzero_elements,values)
          write (6,*) "Finished generating sparse ham heg"
          call flush(6)
          ! Diagonalize just to get the energy in the deterministic space
          if (n_det>1) then
            allocate(lowest_eigenvector(n_det))
            call diagonalize_sparse_hamiltonian_heg(dets_up(1:n_det),dets_dn(1:n_det),n_det,lowest_eigenvector,lowest_eigenvalue)
            n_csfs = count_csfs(n_det,lowest_eigenvector(1:n_det))
            deallocate(lowest_eigenvector)
          else
            call hamiltonian_heg(dets_up(1),dets_dn(1),dets_up(1),dets_dn(1),lowest_eigenvalue)
          endif
        else
          stop 'This hamiltonian type has not been coded (only hubbardk, chem, and heg)'
        endif
        write(6,'(/,''Final deterministic space has: '',i7,'' CSFs'',i8,'' dets., ground state energy='',f12.6)') n_csfs, n_det, lowest_eigenvalue
        if (present(dtm_energy))  dtm_energy=lowest_eigenvalue
        values=-tau*values

      else ! psi trial
        allocate(values(n_det))
        if(rediagonalize) then
          if (iters>=0) then
            if (hamiltonian_type.eq.'hubbardk') then
              call diagonalize_sparse_hamiltonian_hubbard(dets_up(1:n_det),dets_dn(1:n_det),values,lowest_eigenvalue,tau_out_for_diag_ham)
            elseif (hamiltonian_type.eq.'chem') then
              call diagonalize_sparse_hamiltonian_chem(dets_up(1:n_det),dets_dn(1:n_det),n_det,values,lowest_eigenvalue,tau_out_for_diag_ham,.true.)
            elseif (hamiltonian_type.eq.'heg') then
              call diagonalize_sparse_hamiltonian_heg(dets_up(1:n_det),dets_dn(1:n_det),n_det,values,lowest_eigenvalue)
            else
              stop 'This hamiltonian type has not been coded (only hubbardk, chem, and heg)'
            endif
          else
            dets_up(1:n_det) = init_up(1:n_det)
            dets_dn(1:n_det) = init_dn(1:n_det)
            values(1:n_det) = init_wts(1:n_det)
          endif

          n_csfs = count_csfs(n_det,values(1:n_det))
          if (run_type.eq.'selected_ci') then
            write(6,'(/,''Deterministically selected determinants: '',i8,'' CSFs'',i8,'' dets., ground state energy='',f12.6)') n_csfs, n_det, lowest_eigenvalue
          else
            write(6,'(/,''Final trial wavefunction has: '',i8,'' CSFs'',i8,'' dets., ground state energy='',f12.6)') n_csfs, n_det, lowest_eigenvalue
          endif
          if (use_psit_out) then
            open(8, file=psit_out_file,status='new')
            if (hamiltonian_type.eq.'chem') then
              write(8,'(i6,2i4,i6,i3)')  n_det,nup-n_core_orb,ndn-n_core_orb,norb,0
            elseif (hamiltonian_type.eq.'hubbardk') then
              write(8,'(i6,2i4,i6,i3)')  n_det,nup-n_core_orb,ndn-n_core_orb,nsites,0
            endif
            write(8,'(''occ'')')
            write (fmt, '(i2)') nup+ndn-2*n_core_orb
            do i = 1, n_det !min(n_det,100)
#ifdef NUM_ORBITALS_GT_127
               tmp_det = dets_up(i)-maskr_vec(n_core_orb)
#else
               tmp_det = dets_up(i)-maskr(n_core_orb,ik)
#endif

               do j=n_core_orb+1,nup
                 orb_up(j) = trailz(tmp_det)+1
                 tmp_det = ibclr(tmp_det,orb_up(j)-1)
               enddo
#ifdef NUM_ORBITALS_GT_127
               tmp_det = dets_dn(i)-maskr_vec(n_core_orb)
#else
               tmp_det = dets_dn(i)-maskr(n_core_orb,ik)
#endif

               do j=n_core_orb+1,ndn
                 orb_dn(j) = trailz(tmp_det)+1
                 tmp_det = ibclr(tmp_det,orb_dn(j)-1)
               enddo
               do j=n_core_orb+1,nup
                 orb_up(j) = orb_up(j) - n_core_orb
               enddo
               do j=n_core_orb+1,ndn
                 orb_dn(j) = orb_dn(j) - n_core_orb
               enddo
               write(8,'(' // trim(fmt) // 'i4,f16.10)') orb_up(n_core_orb+1:nup), orb_dn(n_core_orb+1:ndn), values(i)
            enddo
          endif
          if (print_psit_wo_sqmc.and..not.use_psit_con_out.and..not.use_elems_out) then
            write (6,*) "Printed trial wave function in ",psit_out_file
            write (6,*) "Terminating SQMC"
            call flush(6)
            stop "Printed trial wave function. Terminating SQMC."
          else
            if (.not.diff_from_psi_t.and.iters<=0) then
              n_imp = 0
              n_csfs_dtm = 0
              previous_value = 0.d0
             !write (6,*) "About to find CSFs for deterministic space. n_det_tmp =",n_det_tmp
              do i = 1, n_det
                if ( abs(abs(previous_value) - abs(values(i))) > epsilon ) then
                  n_csfs_dtm=n_csfs_dtm+1
                  previous_value = values(i)
                  if (i>size_deterministic) then
                 !if (n_csfs_dtm>size_deterministic) then
                    n_csfs_dtm=n_csfs_dtm-1
                    n_imp = i-1
                    exit
                  endif
                endif
              enddo
              if (n_imp.eq.0)  n_imp=n_det

              allocate(imp_up(n_imp))
              allocate(imp_dn(n_imp))
              imp_up(1:n_imp)=dets_up(1:n_imp)
              imp_dn(1:n_imp)=dets_dn(1:n_imp)
            endif
            if (.not.diff_from_psi_t) then

              allocate(iorder(n_imp))
              allocate(temp_i16_up((n_imp+1)/2))
              allocate(temp_i16_dn((n_imp+1)/2))
              allocate(temp_i_2((n_imp+1)/2))
              do j=1,n_imp
                iorder(j)=j
              enddo

              call merge_sort2_up_dn(imp_up,imp_dn, iorder, n_imp, temp_i16_up, temp_i16_dn, temp_i_2)

              deallocate(iorder)
              deallocate(temp_i16_up)
              deallocate(temp_i16_dn)
              deallocate(temp_i_2)
            endif

          endif
          write(6,'(''Largest'',i6,'' abs coefs. are:'')') min(n_det,100)
          write(6,*) "      dets_up             dets_dn        dets_up      dets_dn  excit_level coefficient"
          write (fmt, '(i2)') 2*num_words
          do i = 1, min(n_det,100)
            write(6,'(' // trim(fmt) // 'b65,' // trim(fmt) // 'i15,i3,f13.9)') dets_up(i), dets_dn(i), dets_up(i), dets_dn(i), popcnt(iand(dets_up(i),not(dets_up(1))))+popcnt(iand(dets_dn(i),not(dets_dn(1)))), values(i)
          enddo
        else
          values=lowest_eigenvector(1:n_det)
          deallocate(lowest_eigenvector)
        endif

        ! Calculate sign condition number, sum_i(abs(sum_j(H_ij*c_j)))/sum_i(sum_j(abs(H_ij*c_j)))
        sign_condition_number_numerator = 0._rk
        sign_condition_number_denominator = 0._rk
        do i=1,n_det
          H_psi = 0._rk
          H_psi_abs = 0._rk
          do j=1,n_det
            call hamiltonian(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j),element,connected)
            if (connected) then
              H_psi = H_psi + element*values(j)
              H_psi_abs = H_psi_abs + abs(element*values(j))
            endif
          enddo
          sign_condition_number_numerator = sign_condition_number_numerator + abs(H_psi)
          sign_condition_number_denominator = sign_condition_number_denominator + H_psi_abs
        enddo

        write(6,'(/,''Sign condition number of final trial wave function is: '',f12.6)') sign_condition_number_numerator/sign_condition_number_denominator

        if(MWALK.lt.0) then
          write(6,'(''stop because MWALK < 0'')')
          stop 'stop because MWALK < 0'
        endif

        if (use_psit_con_in.or..not.diff_from_psi_t) then
          if (hamiltonian_type.eq.'hubbardk') then
            call generate_sparse_ham_hubbardk_upper_triangular(imp_up(1:n_imp),imp_dn(1:n_imp),minus_tau_H_indices,minus_tau_H_nonzero_elements,minus_tau_H_values,hf_to_psit,importance_sampling=(importance_sampling==1))
           !call generate_sparse_ham_hubbardk_upper_triangular(imp_up(1:n_imp),imp_dn(1:n_imp),minus_tau_H_indices,minus_tau_H_nonzero_elements,minus_tau_H_values,importance_sampling=(importance_sampling==1))
            ! Diagonalize just to get the energy in the deterministic space
            if (n_imp>1) then
              allocate(lowest_eigenvector(n_imp))
              call diagonalize_sparse_hamiltonian_hubbard(imp_up(1:n_imp),imp_dn(1:n_imp),lowest_eigenvector,lowest_eigenvalue,tau_out_for_diag_ham,sort_or_not=.true.)
              n_csfs_dtm = count_csfs(n_imp,lowest_eigenvector(1:n_imp))
              deallocate(lowest_eigenvector)
            else
              if (space_sym) then
                  call hamiltonian_hubbard_k_space_sym(imp_up(1),imp_dn(1),imp_up(1),imp_dn(1),lowest_eigenvalue,nnzero)
              else
                  call hamiltonian_hubbard_k(imp_up(1),imp_dn(1),imp_up(1),imp_dn(1),lowest_eigenvalue)
              endif
            endif
          elseif (hamiltonian_type.eq.'chem') then
            call generate_sparse_ham_chem_upper_triangular(imp_up(1:n_imp),imp_dn(1:n_imp),minus_tau_H_indices,minus_tau_H_nonzero_elements,minus_tau_H_values,.true.,hf_to_psit=hf_to_psit)
           !call generate_sparse_ham_chem_upper_triangular(imp_up(1:n_imp),imp_dn(1:n_imp),minus_tau_H_indices,minus_tau_H_nonzero_elements,minus_tau_H_values,.true.)
            if (n_imp>1) then
              allocate(lowest_eigenvector(n_imp))
              call diagonalize_sparse_hamiltonian_chem(imp_up(1:n_imp),imp_dn(1:n_imp),n_imp,lowest_eigenvector,lowest_eigenvalue,tau_out_for_diag_ham,.true.,.false.)
              ! Calculate sign condition number, sum_i(abs(sum_j(H_ij*c_j)))/sum_i(sum_j(abs(H_ij*c_j)))
              sign_condition_number_numerator = 0._rk
              sign_condition_number_denominator = 0._rk
              do i=1,n_imp
                H_psi = 0._rk
                H_psi_abs = 0._rk
                do j=1,n_imp
                  call hamiltonian(imp_up(i),imp_dn(i),imp_up(j),imp_dn(j),element,connected)
                  if (connected) then
                    H_psi = H_psi + element*lowest_eigenvector(j)
                    H_psi_abs = H_psi_abs + abs(element*lowest_eigenvector(j))
                  endif
                enddo
                sign_condition_number_numerator = sign_condition_number_numerator + abs(H_psi)
                sign_condition_number_denominator = sign_condition_number_denominator + H_psi_abs
              enddo
              write(6,'(/,''Sign condition of final deterministic space is: '',f12.6)') sign_condition_number_numerator/sign_condition_number_denominator
              n_csfs_dtm = count_csfs(n_imp,lowest_eigenvector(1:n_imp))
              deallocate(lowest_eigenvector)
            else
              if (time_sym) then
                  call hamiltonian_chem_time_sym(imp_up(1),imp_dn(1),imp_up(1),imp_dn(1),lowest_eigenvalue)
              else
                  call hamiltonian_chem(imp_up(1),imp_dn(1),imp_up(1),imp_dn(1),0,lowest_eigenvalue)
              endif
            endif
          elseif (hamiltonian_type.eq.'heg') then
            call generate_sparse_ham_heg(imp_up(1:n_imp),imp_dn(1:n_imp),minus_tau_H_indices,minus_tau_H_nonzero_elements,minus_tau_H_values)
            ! Diagonalize just to get the energy in the deterministic space
            if (n_imp>1) then
              allocate(lowest_eigenvector(n_imp))
              call diagonalize_sparse_hamiltonian_heg(imp_up(1:n_imp),imp_dn(1:n_imp),n_imp,lowest_eigenvector,lowest_eigenvalue)
              n_csfs_dtm = count_csfs(n_imp,lowest_eigenvector(1:n_imp))
              deallocate(lowest_eigenvector)
            else
              call hamiltonian_heg(imp_up(1),imp_dn(1),imp_up(1),imp_dn(1),lowest_eigenvalue)
            endif
          else
            stop 'This hamiltonian type has not been coded (only hubbardk, chem, and heg)'
          endif
          write(6,'(/,''Final deterministic space has: '',i7,'' CSFs'',i8,'' dets., ground state energy='',f12.6)') n_csfs_dtm, n_imp, lowest_eigenvalue
          if (present(dtm_energy))  dtm_energy=lowest_eigenvalue

         !write(6,'(''Largest'',i6,'' abs coefs. are:'')') min(n_imp,100)
         !write(6,*) "      imp_up             imp_dn        imp_up      imp_dn  excit_level coefficient"
         !write (fmt, '(i2)') 2*num_words
         !do i = 1, min(n_imp,100)
         !   write(6,'(' // trim(fmt) // 'b65,' // trim(fmt) // 'i15,i3,f13.9)') imp_up(i), imp_dn(i), imp_up(i), imp_dn(i), popcnt(iand(imp_up(i),not(imp_up(1))))+popcnt(iand(imp_dn(i),not(imp_dn(1)))), values(i)
         !enddo

          minus_tau_H_values=-tau*minus_tau_H_values

        endif
      endif

      if (present(e_trial))  e_trial=lowest_eigenvalue

      call my_second(2,'generate_space_iterate')
      call flush(6)

  end subroutine generate_space_iterate
!-----------------------------------------------------------------------------------------------------------------------------------

  subroutine perform_selected_ci(iters,max_num_of_initiators,n_sym_uniq_det)

    ! Performs a Selected CI, where the determinants are selected in a way inspired by initiator-based FCIQMC:
    ! Starting with HF, each iteration, find all connections to a list of "initiators", i.e., a set of the
    ! largest-magnitude determinants from the previous iteration. Then, apply a projector that connects
    ! the previous iteration's coefficients to the new list of connections. Finally, sort and truncate the
    ! list of connections, and this becomes next iteration's starting vector. This process is repeated as
    ! specified by the inputs, and then the Hamiltonian is diagonalized in the set of determinants left on
    ! the final iteration.
    ! A Holmes, 30 Jul 2012

      use common_run, only : ipr,tau_multiplier,tau
      use chemistry, only  : diagonalize_sparse_hamiltonian_chem,time_sym,is_connected_chem,norb,find_connected_dets_chem
      use hubbard, only    : diagonalize_sparse_hamiltonian_hubbard,space_sym,is_connected_hubbard_fast,nsites,find_connected_dets_hubbard_k
      use heg, only        : diagonalize_sparse_hamiltonian_heg,find_connected_dets_heg
      use more_tools, only : real_symmetric_diagonalize,binary_search
      use generic_sort, only : sort
      use common_ham, only : max_energy
      use common_selected_ci, only : max_nonzero_elements,too_big_to_store,log_num_nonzero_elements
      use types, only : num_words

      implicit none

      integer,intent(in) :: iters,n_sym_uniq_det(:)
      integer,intent(in) :: max_num_of_initiators(:)

      integer :: i,j,k,n_csfs,n_det
      real(rk) :: lowest_eigenvalue
      real(rk),allocatable :: lowest_eigenvector(:)
#ifdef NUM_ORBITALS_GT_127
      type(ik_vec),allocatable :: tmp_up(:),tmp_dn(:)
      type(ik_vec), allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
      type(ik_vec),allocatable :: dets_up(:),dets_dn(:)
      type(ik_vec) :: det_tmp,max_up,max_dn ! for finding max energy
#else
      integer(ik),allocatable :: tmp_up(:),tmp_dn(:)
      integer(ik), allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
      integer(ik),allocatable :: dets_up(:),dets_dn(:)
      integer(ik) :: det_tmp,max_up,max_dn ! for finding max energy
#endif
      integer :: iter,n_det_tmp
      logical :: rediagonalize,connected
      real(rk), parameter :: epsilon = 1.e-8_rk
      real(rk) :: previous_value
      real(rk) :: tau_out_for_diag_ham ! Just for the output of diagonalize_hamiltonian routine for the iterations other than the final one. Not used otherwise.
      integer,allocatable :: iorder(:),temp_i_2(:)
      integer             :: ninitiator,n_det_old
      real(rk),allocatable :: tmp_wt(:)
      real(rk) :: matrix_element,normalization
     !real(rk) :: tau_deterministic
     !real(rk) :: max_diag,min_diag
      real(rk),allocatable :: ham(:,:),eigenvalues(:),eigenvectors(:,:)
      integer :: n_connected_dets
      real(rk) :: mixed_energy,element
      real(rk) :: initiator
      integer :: n_initiators
      integer :: highest_occ
      real(rk),allocatable :: lanczos_initial_vector(:)
      character(len=2) :: fmt

      write (6,*) "max_num_of_initiators=",max_num_of_initiators
      write (6,*) "truncate to=",n_sym_uniq_det
      call flush(6)

      rediagonalize=.true.

      n_det_tmp=1
      if (hamiltonian_type.eq.'hubbardk')  n_det_tmp = ndeg

      allocate(tmp_up(n_det_tmp))
      allocate(tmp_dn(n_det_tmp))
      if (hamiltonian_type.eq.'hubbardk') then
        tmp_up=k_hf_deg_up
        tmp_dn=k_hf_deg_dn
      else
        tmp_up(1)=2_ik**nup-1_ik
        tmp_dn(1)=2_ik**ndn-1_ik
      endif

      initiator = 0._rk
      allocate(lowest_eigenvector(n_det_tmp))
      lowest_eigenvector(:) = 0._rk
      lowest_eigenvector(1:n_det_tmp) = initiator + 1._rk ! so it definitely passes initiator test
      if (hamiltonian_type.eq.'hubbardk')  norb=nsites

      do iter=1,iters
        call my_second(1,'find_doubly_excited')
        call flush(6)
        ! Single application of H:
        n_det_old = n_det_tmp
        n_initiators = min(n_det_tmp,max_num_of_initiators(iter))
        write (6,*) "n_det_tmp",n_det_tmp
        write (6,*) "norb=",norb
        write (6,*) "ninitiators",n_initiators;call flush(6)
        call find_doubly_excited(n_det=n_det_tmp,dets_up=dets_up,dets_dn=dets_dn,ref_up=tmp_up,ref_dn=tmp_dn,norb=norb,n_core_orb=n_core_orb,ref_coeffs=lowest_eigenvector,ninitiator=n_initiators)
        write (6,*) "n_det_tmp",n_det_tmp
        call flush(6)
        ! diagonalization

       !call hamiltonian(tmp_up(1),tmp_dn(1),tmp_up(1),tmp_dn(1),matrix_element,connected)
       !max_diag = matrix_element
       !min_diag = matrix_element
       !if (n_det_old.ge.2) then
       !  do i=2,n_det_old
       !    call hamiltonian(tmp_up(i),tmp_dn(i),tmp_up(i),tmp_dn(i),matrix_element,connected)
       !    if (matrix_element>max_diag)  max_diag = matrix_element
       !    if (matrix_element<min_diag)  min_diag = matrix_element
       !  enddo
       !endif
       !write (6,*) "taum,max,min",tau_multiplier,max_diag,min_diag
       !if (iter==1) then
       !  tau_deterministic = tau
       !else
       !  tau_deterministic = tau_multiplier/(max_diag-min_diag)
       !endif
       !write (6,*) "tau_deterministic",tau_deterministic
       !write (6,*) "maxval(lowest_eigenvector) before proj",maxval(lowest_eigenvector)
        call my_second(1,'Find contributions from currently occupied dets')
        call flush(6)
        allocate(tmp_wt(n_det_tmp))

        ! For each det in dets_up,dn(1:n_det_tmp), calculate its weight by summing over the contributions from dets in tmp_up,tmp_dn(1:n_det_old).
        if (iter==1) then
          ! diagonalize
          allocate(ham(n_det_tmp,n_det_tmp))
          ham(:,:) = 0._rk
          do i=1,n_det_tmp
            do j=i,n_det_tmp
              call hamiltonian(dets_up(i),dets_dn(i),dets_up(j),dets_dn(j),matrix_element,connected)
              if (connected) then
                ham(i,j) = matrix_element
                ham(j,i) = matrix_element
              endif
            enddo
          enddo
          if (n_det_tmp<1000) then
            allocate(eigenvectors(n_det_tmp,n_det_tmp))
            allocate(eigenvalues(n_det_tmp))
            call real_symmetric_diagonalize(n_det_tmp,ham,eigenvectors,eigenvalues)
            tmp_wt(:) = eigenvectors(:,1)
            lowest_eigenvalue = eigenvalues(1)
            deallocate(ham,eigenvectors,eigenvalues)
          else
            call matrix_lanczos(n_det_tmp,tmp_wt,lowest_eigenvalue,ham)
          endif
        else

        ! For each old, generate connections; binary search new for overlap (Scales as N_truncate*N_connected*log(N_initiator*N_connected))

          ! sort dets_up,dn(1:n_det_old) by label
          allocate(iorder(n_det_tmp))
          allocate(temp_i16_up((n_det_tmp+1)/2))
          allocate(temp_i16_dn((n_det_tmp+1)/2))
          allocate(temp_i_2((n_det_tmp+1)/2))
          do j=1,n_det_tmp
            iorder(j)=j
          enddo

          call merge_sort2_up_dn(dets_up(1:n_det_tmp),dets_dn(1:n_det_tmp), iorder, n_det_tmp, temp_i16_up, temp_i16_dn, temp_i_2)

          deallocate(iorder)
          deallocate(temp_i16_up)
          deallocate(temp_i16_dn)
          deallocate(temp_i_2)

          ! Apply projector
          if (n_det_tmp>10) then
            det_tmp = max(maxval(tmp_up(1:n_det_old)),maxval(tmp_dn(1:n_det_old)))
            do i=1,100
              highest_occ = trailz(det_tmp) + 1
              det_tmp = ibclr(det_tmp,highest_occ-1)
              if (det_tmp==0_ik)  exit
            enddo
            if (i>90)  write (6,*) "infinite loop!"
            write (6,*) "highest occupied orbital=",highest_occ
            if (hamiltonian_type.eq.'hubbardk') then
              max_energy = 0._rk
              do i=1,nup
                max_energy = max_energy + k_energies(highest_occ-i+1)
              enddo
              do i=1,ndn
                max_energy = max_energy + k_energies(highest_occ-i+1)
              enddo
            elseif (hamiltonian_type.eq.'chem') then
              max_up = 0_ik
              do i=1,nup
                max_up = ibset(max_up,highest_occ-i)
              enddo
              max_dn = 0_ik
              do i=1,ndn
                max_dn = ibset(max_dn,highest_occ-i)
              enddo
              call hamiltonian(max_up,max_dn,max_up,max_dn,max_energy)
            endif
            write (6,*) "Max energy=",max_energy; call flush(6)
          endif

          tmp_wt(:) = 0._rk
          do i=1,n_det_old
            ! generate connections = con_up,dn
            if (hamiltonian_type.eq.'hubbardk') then
              call find_connected_dets_hubbard_k(tmp_up(i),tmp_dn(i),n_connected_dets,connected_dets_up,connected_dets_dn,nsites)
            else
              call find_connected_dets_chem(tmp_up(i), tmp_dn(i), n_connected_dets, connected_dets_up, connected_dets_dn,norb)
            endif
            ! loop over connections: binary search the tmp list for each connection.
            ! if exists, tmp_wt(i) = tmp_wt(i) + projection element between connection j and dets_up(i)
            do j=1,n_connected_dets
              call binary_search(connected_dets_up(j),connected_dets_dn(j),dets_up(1:n_det_tmp),dets_dn(1:n_det_tmp),k)
              if (k<0) then
                write (6,*) "Infinite loop!"
                call flush(6)
                stop
              endif
              if (k>0) then
                ! New projector: H - max(H_ii)
                call hamiltonian(tmp_up(i),tmp_dn(i),dets_up(k),dets_dn(k),matrix_element)
                tmp_wt(k) = tmp_wt(k) - matrix_element * lowest_eigenvector(i)
                if (tmp_up(i)==dets_up(k).and.tmp_dn(i)==dets_dn(k)) then
                  tmp_wt(k) = tmp_wt(k) + max_energy * lowest_eigenvector(i)
                endif
              endif
            enddo
          enddo

        endif

        deallocate(tmp_up,tmp_dn)

        call my_second(2,'Find contributions from currently occupied dets')
        call flush(6)

        !pick the number of desired unique (by symmetry) determinants
        !Actually for molecules (unlike the heg) this does not guarantee it because the coefficients of the determinants in a CSF need not be equal in absolute magnitude.

        ! Estimate energy
        ! Currently commented out because this is approximately 5 times as expensive as applying the projector (in the tests done so far)

       !if (iter>1) then
       !  lowest_eigenvalue = 0._rk
       !  do i=1,n_det_tmp
       !    ! generate connections = connected_dets_up,dn
       !    if (hamiltonian_type.eq.'hubbardk') then
       !      call find_connected_dets_hubbard_k(dets_up(i),dets_dn(i),n_connected_dets,connected_dets_up,connected_dets_dn,nsites)
       !    else
       !      call find_connected_dets_chem(dets_up(i), dets_dn(i), n_connected_dets, connected_dets_up, connected_dets_dn,norb)
       !    endif
       !    ! loop over connections: binary search the tmp list for each connection.
       !    ! if exists, lowest_eigenvector(i) = lowest_eigenvector(i) + projection element between connection j and dets_up(i)
       !    do j=1,n_connected_dets
       !      call binary_search(connected_dets_up(j),connected_dets_dn(j),dets_up(1:n_det_tmp),dets_dn(1:n_det_tmp),k)
       !      if (k>0) then
       !        call hamiltonian(dets_up(i),dets_dn(i),dets_up(k),dets_dn(k),matrix_element,connected)
       !        if (connected)  lowest_eigenvalue = lowest_eigenvalue + matrix_element * tmp_wt(k) * tmp_wt(i)
       !      endif
       !    enddo
       !  enddo
       !  lowest_eigenvalue = lowest_eigenvalue/(dot_product(tmp_wt(1:n_det_tmp),tmp_wt(1:n_det_tmp)))
       !endif
       !write (6,*) "Pure energy=",lowest_eigenvalue
       !call flush(6)

        call sort(tmp_wt, dets_up(1:n_det_tmp), dets_dn(1:n_det_tmp))

        ! Calculate mixed energy for current wave function
        ! Significantly cheaper than the variational energy
        mixed_energy = 0._rk
        do i=1,n_det_tmp
          call hamiltonian(dets_up(i),dets_dn(i),dets_up(1),dets_dn(1),element,connected)
          if (connected)  mixed_energy=mixed_energy+element*tmp_wt(i)
        enddo
        mixed_energy = mixed_energy / tmp_wt(1)
        write (6,*) "mixed energy at iteration number",iter,"is",mixed_energy
        call flush(6)

        ! for selected CI, we should truncate by number of dets instead of number of CSFs, since applying the projector only once will cause us to incorrectly count CSFs.

        ! scan weights, starting with n_sym_uniq_det(iter)'th one, to include all dets with wt equal to that one
        if (n_sym_uniq_det(iter)<n_det_tmp) then
          previous_value = tmp_wt(n_sym_uniq_det(iter))
          do i=n_sym_uniq_det(iter),n_det_tmp
            if ( abs(abs(previous_value) - abs(tmp_wt(i))) > epsilon ) then
              n_det = i-1
              exit
            endif
          enddo
        else
          n_det = n_det_tmp
        endif

        n_csfs=-1 ! because it's not a useful quantity here.
        n_det_tmp=n_det

       !deallocate(tmp_up,tmp_dn)
        deallocate(lowest_eigenvector)
        allocate(lowest_eigenvector(n_det))
        lowest_eigenvector = tmp_wt(1:n_det)
        deallocate(tmp_wt)

        normalization = 1 / sqrt(dot_product(lowest_eigenvector(1:n_det),lowest_eigenvector(1:n_det)))
       !normalization = 1 / sqrt(dot_product(lowest_eigenvector(1:n_det_tmp),lowest_eigenvector(1:n_det_tmp)))
        lowest_eigenvector(:) = lowest_eigenvector(:)*normalization

        ! reset initiator cutoff value
        initiator = abs(lowest_eigenvector(min(max_num_of_initiators(iter),n_det)))
       !initiator = abs(lowest_eigenvector(min(max_num_of_initiators(iter),n_det_tmp)))

        ! count number of initiators for this iteration.
        ninitiator=0
        do i=1,n_det
       !do i=1,n_det_tmp
          if (abs(lowest_eigenvector(i)).ge.initiator-epsilon) then
            ninitiator = ninitiator + 1
          endif
        enddo

        write(6,*) "initiator=", initiator
        write(6,'(/,''Selected CI wavefunction after iteration '',i2,'' has:'',i6,'' initiators'',i8,'' dets., with mixed energy'',f12.6)') iter, ninitiator, n_det_tmp, mixed_energy
        write(6,'(/,''Truncated Selected CI wavefunction after iteration '',i2,'' has: '',i8,'' initiators (dets),'',i8,'' CSFs,'',i8,'' dets.'')') iter, ninitiator, n_csfs, n_det
        write(6,'(/,''Smallest abs coefficient'',2es12.4)') lowest_eigenvector(n_det), minval(abs(lowest_eigenvector(1:n_det)))
        write(6,'(''Largest'',i6,'' abs coefs. are:'')') min(n_det,100)
        write(6,*) "        up-det         dn-det"
        write (fmt, '(i2)') 2*num_words
        do i = 1, min(n_det,100)
          write(6,'(' // trim(fmt) // 'i15,f13.9)') dets_up(i), dets_dn(i), lowest_eigenvector(i)
        enddo
        call flush(6)

        if (iter<iters) then
          allocate(tmp_up(n_det))
          allocate(tmp_dn(n_det))
          tmp_up = dets_up(1:n_det)
          tmp_dn = dets_dn(1:n_det)
          deallocate(dets_up,dets_dn)
        endif

      enddo

      if (n_det.ge.2000) then
        allocate(lanczos_initial_vector(n_det))
        lanczos_initial_vector(:)=lowest_eigenvector(:)
      endif
      deallocate(lowest_eigenvector)

      if (size(dets_up).ne.n_det) then
        allocate(tmp_up(n_det))
        allocate(tmp_dn(n_det))
        tmp_up(1:n_det)=dets_up(1:n_det)
        tmp_dn(1:n_det)=dets_dn(1:n_det)
        deallocate(dets_up,dets_dn)
        allocate(dets_up(n_det))
        allocate(dets_dn(n_det))
        dets_up(1:n_det)=tmp_up(1:n_det)
        dets_dn(1:n_det)=tmp_dn(1:n_det)
        deallocate(tmp_up,tmp_dn)
      endif

      if(rediagonalize) then

        ! Sort again, since the dets are currently in order of absolute magnitude of their coefficients.
        allocate(iorder(n_det))
        allocate(temp_i16_up((n_det+1)/2))
        allocate(temp_i16_dn((n_det+1)/2))
        allocate(temp_i_2((n_det+1)/2))
        do j=1,n_det
          iorder(j)=j
        enddo

        call merge_sort2_up_dn(dets_up,dets_dn, iorder, n_det, temp_i16_up, temp_i16_dn, temp_i_2)

        if (n_det.ge.2000) lanczos_initial_vector(1:n_det)=lanczos_initial_vector(iorder(1:n_det))

        deallocate(iorder)
        deallocate(temp_i16_up)
        deallocate(temp_i16_dn)
        deallocate(temp_i_2)

        ! Finally, rediagonalize

        call my_second(1,'rediagonalize')
        call flush(6)

        allocate(lowest_eigenvector(n_det))
        if (hamiltonian_type.eq.'hubbardk') then
          call diagonalize_sparse_hamiltonian_hubbard(dets_up(1:n_det),dets_dn(1:n_det),lowest_eigenvector,lowest_eigenvalue,tau_out_for_diag_ham,initial_vector=lanczos_initial_vector)
        elseif (hamiltonian_type .eq. 'chem') then
          call diagonalize_sparse_hamiltonian_chem(dets_up(1:n_det),dets_dn(1:n_det),n_det,lowest_eigenvector,lowest_eigenvalue,tau_out_for_diag_ham,.true.,initial_vector=lanczos_initial_vector)
        elseif (hamiltonian_type .eq. 'heg') then
          call diagonalize_sparse_hamiltonian_heg(dets_up(1:n_det),dets_dn(1:n_det),n_det,lowest_eigenvector,lowest_eigenvalue)
        endif
        write(6,'(/,''Final selected ci ground state wavefunction has: '',i8,'' CSFs'',i8,'' dets., ground state energy='',f12.6)') n_csfs, n_det, lowest_eigenvalue
        write(6,'(''Largest'',i6,'' abs coefs. are:'')') min(n_det,100)
        write(6,*) "      dets_up             dets_dn        dets_up      dets_dn  excit_level coefficient"
        write (fmt, '(i2)') 2*num_words
        do i = 1, min(n_det,100)
           write(6,'(' // trim(fmt) // 'b65,' // trim(fmt) // 'i15,i3,f13.9)') dets_up(i), dets_dn(i), dets_up(i), dets_dn(i), popcnt(iand(dets_up(i),not(dets_up(1))))+popcnt(iand(dets_dn(i),not(dets_dn(1)))), lowest_eigenvector(i)
        enddo
      endif

      write(6,'(/,''Deterministically selected determinants: '',i8,'' CSFs'',i8,'' dets., ground state energy='',f12.6)') n_csfs, n_det, lowest_eigenvalue
      call my_second(2,'selected_ci')
      call flush(6)

  end subroutine perform_selected_ci
!-------------------------------------------------------------------------------

  subroutine perform_truncated_lanczos(iters,n_initiators,n_truncate)
    ! Perform a Lanczos diagonalization, truncating the current vector at every step to n_sym_uniq_det,
    ! and only applying the Hamiltonian to the highest-weight max_num_of_initiators dets.

    ! A Holmes, 6 Aug 2012

      use common_run, only : ipr,tau_multiplier,tau
      use chemistry, only  : diagonalize_sparse_hamiltonian_chem,time_sym,is_connected_chem,norb,find_connected_dets_chem
      use hubbard, only    : diagonalize_sparse_hamiltonian_hubbard,space_sym,is_connected_hubbard_fast,nsites,find_connected_dets_hubbard_k,c_sym_psi_t
      use more_tools, only : real_symmetric_diagonalize,binary_search
      use generic_sort, only : sort
      use tools, only : merge_sort2_up_dn
      use common_ham, only : n_core_orb

      implicit none

      integer,intent(in) :: iters,n_initiators,n_truncate

      integer :: i,j,k,m
      integer :: iter,n_det_tmp
      logical :: connected
      real(rk), parameter :: epsilon = 1.e-8_rk
      real(rk) :: normalization
      real(rk),allocatable :: eigenvalues(:),eigenvectors(:,:)
      real(rk) :: element
#ifdef NUM_ORBITALS_GT_127
      type(ik_vec),allocatable :: dets_up(:),dets_dn(:)
      type(ik_vec),allocatable :: up(:,:),dn(:,:),up_iter(:),dn_iter(:)
      type(ik_vec),allocatable :: temp_i16_up(:),temp_i16_dn(:)
#else
      integer(ik),allocatable :: dets_up(:),dets_dn(:)
      integer(ik),allocatable :: up(:,:),dn(:,:),up_iter(:),dn_iter(:)
      integer(ik),allocatable :: temp_i16_up(:),temp_i16_dn(:)
#endif
      integer :: n_connected_dets

      real(rk),allocatable :: v(:,:)
      integer,allocatable :: len_v(:)
      real(rk),allocatable :: H(:,:)
      real(rk),allocatable :: dotproducts(:),v_tmp(:),krylov_ground_states(:)
      integer :: actual_n_initiators
      integer,allocatable :: iorder(:)
      integer,allocatable :: temp_i_2(:)
      real(rk) :: tmp1
      real(rk),allocatable :: A(:,:),x(:),b(:)
      integer,allocatable :: ipiv(:)
      real(rk) :: initiator

      allocate(dotproducts(iters))
      allocate(krylov_ground_states(iters))

      write (6,*) "Performing truncated Lanczos for",iters,"iterations, truncating to",n_truncate,"dets, and applying H to",n_initiators,"initiators."
      call flush(6)

      ! Allocate set of all vectors used in truncated Lanczos
      allocate(v(n_truncate,iters))
      allocate(len_v(iters))
      allocate(up(n_truncate,iters))
      allocate(dn(n_truncate,iters))
      allocate(up_iter(n_truncate))
      allocate(dn_iter(n_truncate))
      v(:,:)=0._rk
      len_v(:)=0
      up(:,:)=0_ik
      dn(:,:)=0_ik

      ! Assign the first vector to be HF
      if (hamiltonian_type.eq.'hubbardk') then
        len_v(1) = ndeg
        v(1:len_v(1),1) = c_sym_psi_t
        write (6,*) "c_sym_psi_t",c_sym_psi_t
        up(1:len_v(1),1) = k_hf_deg_up
        dn(1:len_v(1),1) = k_hf_deg_dn
      else
        len_v(1) = 1
        v(1,1) = 1._rk
        up(1,1) = 2_ik**nup-1_ik
        dn(1,1) = 2_ik**ndn-1_ik
      endif
      if (hamiltonian_type.eq.'hubbardk')  norb=nsites

      ! Allocate the Hamiltonian in the Krylov space
      allocate(H(iters,iters))
      H(:,:)=0._rk
      do i=1,len_v(1)
        do j=1,len_v(1)
          call hamiltonian(up(i,1),dn(i,1),up(j,1),dn(j,1),element,connected)
          if (connected)  H(1,1)=H(1,1)+element*v(i,1)*v(j,1)
        enddo
      enddo

      write (6,*) "Hamiltonian in Krylov space after iteration",1,":"
      call print_real_matrix(iter,iter,H(1:iter,1:iter))

      ! 1 Apply H to the first n_initiators dets (and sort by label)
      ! 2 Orthogonalize current vector relative to all previous
      ! 3 Sort by abs wt, then truncate to n_truncate, then sort by label
      ! 4 Re-orthogonalize truncated vector relative to all previous (shouldn't do much)
      ! 5 Normalize truncated vector
      ! 6 Compute matrix elements of H in truncated Krylov space
      ! 7 Perform exact diagonalization on Krylov H
      ! 8 Repeat until either convergence or max number of iterations

      do iter=2,iters

        ! adjust number of initiators to include all that have equal wt to one another.
        actual_n_initiators = min(len_v(iter-1),n_initiators)
        if (actual_n_initiators.ne.len_v(iter-1)) then
          initiator = abs(v(min(actual_n_initiators,len_v(iter-1)),iter-1))
          actual_n_initiators = 0
          do i=1,len_v(iter-1)
            if (abs(v(i,iter-1))-initiator>-epsilon) then
              actual_n_initiators = actual_n_initiators + 1
            else
              exit
            endif
          enddo
        endif
        write (6,*) "Number of initiators:",actual_n_initiators
        call flush(6)

        ! 1 Apply H to the first n_initiators dets
        call find_doubly_excited(n_det=n_det_tmp,dets_up=dets_up,dets_dn=dets_dn,ref_up=up(1:len_v(iter-1),iter-1),ref_dn=dn(1:len_v(iter-1),iter-1),norb=norb,n_core_orb=n_core_orb,ref_coeffs=v(1:len_v(iter-1),iter-1),ninitiator=actual_n_initiators,e_mix_num=v_tmp)

        ! Apply H to all noninitiators
        if (len_v(iter-1)>actual_n_initiators) then
          do i=actual_n_initiators+1,len_v(iter-1)
            ! generate connections = con_up,dn
            if (hamiltonian_type.eq.'hubbardk') then
              call find_connected_dets_hubbard_k(up(i,iter-1),dn(i,iter-1),n_connected_dets,connected_dets_up,connected_dets_dn,nsites)
            else
              call find_connected_dets_chem(up(i,iter-1),dn(i,iter-1), n_connected_dets, connected_dets_up, connected_dets_dn,norb)
            endif
            ! loop over connections: binary search the tmp list for each connection.
            ! if exists, tmp_wt(i) = tmp_wt(i) + projection element between connection j and dets_up(i)
            do j=1,n_connected_dets
              call binary_search(connected_dets_up(j),connected_dets_dn(j),dets_up(1:n_det_tmp),dets_dn(1:n_det_tmp),k)
              if (k>0) then
                call hamiltonian(up(i,iter-1),dn(i,iter-1),dets_up(k),dets_dn(k),element)
                v_tmp(k) = v_tmp(k) - element * v(i,iter-1)
              endif
            enddo
          enddo
        endif

       !if (iter==2) then
       !  if (hamiltonian_type.eq.'hubbardk') then ! Need to call different routines because within the diagonalize routine they call different Hamiltonians.  Could be combined.
       !    call diagonalize_sparse_hamiltonian_hubbard(dets_up(1:n_det_tmp),dets_dn(1:n_det_tmp),v_tmp,lowest_eigenvalue,tau_out_for_diag_ham)
       !  else
       !    call diagonalize_sparse_hamiltonian_chem(dets_up(1:n_det_tmp),dets_dn(1:n_det_tmp),n_det_tmp,v_tmp,lowest_eigenvalue,tau_out_for_diag_ham,.true.)
       !  endif
       !  write (6,*) "CISD energy=",lowest_eigenvalue
       !endif

        ! Orthogonalize current vector relative to all previous
        ! Current vector: dets_up,dn(1:n_det_tmp), all previous: v(1:len_v(i),i), for all i=1,...,iter-1

        allocate(A(iter-1,iter-1))
        allocate(x(iter-1))
        allocate(b(iter-1))
        allocate(ipiv(iter-1))

        v_tmp(1:n_det_tmp) = v_tmp(1:n_det_tmp) / sqrt(dot_product(v_tmp(1:n_det_tmp),v_tmp(1:n_det_tmp)))

        call orthogonalize(v_tmp,dets_up,dets_dn,v(:,1:iter-1),up(:,1:iter-1),dn(:,1:iter-1),len_v(1:iter-1),A,x,b,ipiv)

        do i=1,iter-1
          tmp1=dotproduct(v(1:len_v(i),i),up(1:len_v(i),i),dn(1:len_v(i),i),v_tmp,dets_up,dets_dn)
          if (abs(tmp1)>1.e-12)  write (6,*) "Not orthogonal after first orthogonalization!",i,tmp1
        enddo

        ! Truncate to n_truncate
        call sort(v_tmp,dets_up,dets_dn)
        len_v(iter) = min(n_truncate,n_det_tmp)

        ! adjust number of truncated dets to include all that have equal wt to one another.
        if (actual_n_initiators.ne.len_v(iter-1)) then
          initiator = abs(v_tmp(n_truncate))
          do i=len_v(iter)-1,1,-1
            if (abs(v_tmp(i))-initiator>epsilon) then
              len_v(iter) = i
              exit
            endif
          enddo
        endif

        v(1:len_v(iter),iter) = v_tmp(1:len_v(iter))
        up(1:len_v(iter),iter) = dets_up(1:len_v(iter))
        dn(1:len_v(iter),iter) = dets_dn(1:len_v(iter))

        ! Re-orthogonalize truncated vector relative to all previous
        ! First, re-sort by label
        allocate(iorder(len_v(iter)))
        allocate(temp_i16_up((len_v(iter)+1)/2))
        allocate(temp_i16_dn((len_v(iter)+1)/2))
        allocate(temp_i_2((len_v(iter)+1)/2))
        do j=1,len_v(iter)
          iorder(j)=j
        enddo

        call merge_sort2_up_dn(up(1:len_v(iter),iter),dn(1:len_v(iter),iter), iorder, len_v(iter), temp_i16_up, temp_i16_dn, temp_i_2)
        v(1:len_v(iter),iter)=v(iorder(1:len_v(iter)),iter)

        deallocate(iorder)
        deallocate(temp_i16_up)
        deallocate(temp_i16_dn)
        deallocate(temp_i_2)

        v(1:len_v(iter),iter) = v(1:len_v(iter),iter) / sqrt(dot_product(v(1:len_v(iter),iter),v(1:len_v(iter),iter)))
        call orthogonalize(v(1:len_v(iter),iter),up(1:len_v(iter),iter),dn(1:len_v(iter),iter),v(:,1:iter-1),up(:,1:iter-1),dn(:,1:iter-1),len_v(1:iter-1),A,x,b,ipiv)
        deallocate(A,x,b,ipiv)

        do i=1,iter-1
          tmp1=dotproduct(v(1:len_v(i),i),up(1:len_v(i),i),dn(1:len_v(i),i),v(1:len_v(iter),iter),up(1:len_v(iter),iter),dn(1:len_v(iter),iter))
          if (abs(tmp1)>1.e-12)  write (6,*) "Not orthogonal after second orthogonalization!",i,tmp1
        enddo

        ! Normalize truncated vector
        normalization = 1/sqrt(dot_product(v(1:len_v(iter),iter),v(1:len_v(iter),iter)))
        v(1:len_v(iter),iter) = normalization * v(1:len_v(iter),iter)

        ! Compute matrix elements of H in truncated Krylov space
        do i=1,iter
          up_iter(1:len_v(iter)) = up(1:len_v(iter),iter)
          dn_iter(1:len_v(iter)) = dn(1:len_v(iter),iter)
          ! compute H between i and iter
          do m=1,len_v(i)
            if (hamiltonian_type.eq.'hubbardk') then
              call find_connected_dets_hubbard_k(up(m,i),dn(m,i),n_connected_dets,connected_dets_up,connected_dets_dn,nsites)
            else
              call find_connected_dets_chem(up(m,i), dn(m,i), n_connected_dets, connected_dets_up, connected_dets_dn,norb)
            endif
            do j=1,n_connected_dets
              call binary_search(connected_dets_up(j),connected_dets_dn(j),up_iter(1:len_v(iter)),dn_iter(1:len_v(iter)),k)
              if (k>0) then
                call hamiltonian(up(m,i),dn(m,i),up(k,iter),dn(k,iter),element,connected)
                if (connected)  H(i,iter)=H(i,iter)+element*v(m,i)*v(k,iter)
              endif
            enddo
          enddo
          if (i.ne.iter)  H(iter,i)=H(i,iter)
        enddo

        ! Perform exact diagonalization on Krylov H
        allocate(eigenvectors(iter,iter))
        allocate(eigenvalues(iter))
        call real_symmetric_diagonalize(iter,H(1:iter,1:iter),eigenvectors,eigenvalues)
        krylov_ground_states(iter) = eigenvalues(1)
        if (iter<=10) then
          write (6,*) "Hamiltonian in Krylov space after iteration",iter,":"
          call print_real_matrix(iter,iter,H(1:iter,1:iter))
        endif
        write (6,*) "Ground state after iteration",iter,"is",krylov_ground_states(iter)
        call flush(6)
        deallocate(eigenvectors,eigenvalues)
        ! Repeat until either convergence or max number of iterations
        !! STILL NEED TO DO THIS !!

        deallocate(v_tmp)
        call my_second(2,"truncated lanczos iteration")

      enddo

      write (6,*) "checking orthogonality"
      do i=1,iters
        do j=1,iters
          if (i.ne.j) then
            tmp1=dotproduct(v(1:len_v(i),i),up(1:len_v(i),i),dn(1:len_v(i),i),v(1:len_v(j),j),up(1:len_v(j),j),dn(1:len_v(j),j))
            write (6,*) i,j,tmp1
          endif
        enddo
      enddo
      write (6,*) ""

      write (6,*) "Truncated Lanczos energy:",krylov_ground_states(iters)
      call flush(6)

  end subroutine perform_truncated_lanczos

!-------------------------------------------------------------------------------
  subroutine find_doubly_excited(n_det, dets_up, dets_dn, ref_up, ref_dn, norb, n_core_orb, ref_coeffs, ninitiator, e_mix_num, e_mix_den, term1_big, term2_big, importance_sampling, pt, eps_var_pt, eps_var_pt_big, ref_diag_elems, new_diag_elems, n_mc, w_over_p, min_H_already_done,core_up,core_dn,virt_up,virt_dn,active_only)
    !---------------------------------------------------------------------------
    ! Description : Generates determinants connected to reference. This was initially written to handle arbitrary excitation level.
    !               Returns dets_up,dn sorted by label.
    !               Written by A. Holmes to handle general case of excitation level. Rewritten by F. Petruzielo, 20 Feb 2012
    ! Modified    : A Holmes, 5 Jul 2012. Optional input ninitiator only finds dets connected to that many of the reference states
    ! Modified    : A Holmes, 19 Nov 2012. If optional input importance_sampling is 1, generate guiding wave function here, and use it to modify e_mix_num and e_mix_den.
    ! Modified    : A Holmes, 22 Jul 2013. Added input n_core_orb
    ! Modified    : A Holmes, 8 Feb 2016. Added optional pt, which if included, returns H_ij*C_j/(E-H_ii), rather than H_ij*C_j, for e_mix_num
    ! Modified    : A Holmes, 24 Feb 2016. Added optional input eps_var_pt, which if included, causes find_important_connected_dets_chem to be called with that eps_var_pt
    !             :                        eps_var_pt in the calling routine can be either eps_var or eps_pt
    ! Modified    : A Holmes, 9 Mar 2016. Added optional inputs ref_diag_elems, new_diag_elems
    !               If included, use the ref_diag_elems to generate the new_diag_elems in O(N) time,
    !               instead of computing them from scratch in O(N^2) time
    ! Modified    : A Holmes, 17 May 2017. Added optional arguments core_up/dn,
    !               virt_up/dn, active_only for active space variational calculation
    !---------------------------------------------------------------------------
    use common_run, only: tau_multiplier, ipr, n_connected_dets_hf !,importance_sampling
    use hubbard, only : find_connected_dets_hubbard_k
    use chemistry, only : find_connected_dets_chem,find_important_connected_dets_chem,get_new_diag_elem
    use heg, only : find_connected_dets_heg, find_important_connected_dets_heg, get_new_diag_elem_heg
    use more_tools, only : binary_search
    use common_psi_t, only : ndet_psi_t,dets_up_psi_t,dets_dn_psi_t,cdet_psi_t,psi_g,psi_g_energy,psi_g_epsilon,psi_g_epsilon_inv
    use types, only : i8b
    use mpi_routines, only : mpi_merge_sort2,whoami,mpi_barr,ncores

    implicit none

    integer,intent(out) :: n_det
    integer,intent(in) :: norb
    integer,intent(in) :: n_core_orb
    real(rk),optional,intent(in) :: ref_coeffs(:)
    integer,optional,intent(in) :: ninitiator !requested number of initiators
    real(rk),allocatable,optional,intent(out) :: e_mix_num(:),e_mix_den(:)
    real(rk),allocatable,optional,intent(out) :: term1_big(:),term2_big(:)
    integer,optional,intent(in) :: importance_sampling
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), allocatable,intent(out) :: dets_up(:), dets_dn(:)
    type(ik_vec),intent(in) :: ref_up(:),ref_dn(:)
    type(ik_vec),optional,intent(in) :: core_up,core_dn,virt_up,virt_dn
    type(ik_vec),allocatable :: old_dets_up(:),old_dets_dn(:)
    type(ik_vec), allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#else
    integer(ik), allocatable,intent(out) :: dets_up(:), dets_dn(:)
    integer(ik),intent(in) :: ref_up(:),ref_dn(:)
    integer(ik),optional,intent(in) :: core_up,core_dn,virt_up,virt_dn
    integer(ik),allocatable :: old_dets_up(:),old_dets_dn(:)
    integer(ik), allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#endif
    real(rk),optional,intent(in) :: pt ! If present, the energy from the last iteration in HCI
    real(rk),optional,intent(in) :: eps_var_pt
    real(rk),optional,intent(in) :: eps_var_pt_big
    real(rk),optional,intent(in) :: ref_diag_elems(:)
    real(rk),allocatable,optional,intent(out) :: new_diag_elems(:)
    integer,optional,intent(in) :: n_mc ! If present, then numerators become term1, denominators become term2
    real(rk),optional,intent(in) :: w_over_p(:) ! If present, then numerators become term1, denominators become term2
    real(rk),optional,intent(in) :: min_H_already_done(:)
    logical,optional,intent(in) :: active_only

    !local variables
    integer ::  i,j,num_mpi_sort
    integer :: n_connected_dets
    integer :: n_ref, modnum
    integer(i8b) :: n_allocate,n_allocate_new
    integer,allocatable :: iorder(:),temp_i_2(:)
    type(diag_elem_info),allocatable :: diag_elems_info(:)
    real(rk) :: initiator
    real(rk),parameter :: epsilon = 1e-10_rk
    logical :: imp,core_filled
    integer :: n_det_core_filled
    real(rk) :: H_ii, avail, avail2, avail3
    logical :: is_connected
    real(rk),allocatable :: connected_matrix_elements_big(:)
    integer is_done, is_not_done, n_samples, ind

    if (present(eps_var_pt_big))  allocate(connected_matrix_elements_big(size(connected_matrix_elements)))

    imp = .false.
    if (present(importance_sampling))  imp=(importance_sampling==1)

    n_ref=size(ref_up,1)
    if (ipr.ge.1)  write (6,'(''n_ref='',i8)') n_ref ; call flush(6)
   !if (present(ninitiator))  n_ref=min(ninitiator,n_ref)
    if (ipr.ge.1)  write (6,'(''Find_doubly_excited called with norb='',i5)') norb ; call flush(6)

    if (present(ninitiator)) then
      if (n_ref.ne.ninitiator) then
        if (present(ref_coeffs)) then
          initiator = abs(ref_coeffs(min(ninitiator,n_ref)))
          n_ref = 0
          do i=1,size(ref_up)
            if (abs(ref_coeffs(i))-initiator>-epsilon) then
              n_ref = n_ref + 1
            else
              exit
            endif
          enddo
        else
          n_ref = ninitiator
        endif
      endif
      if (ipr.ge.1)  write (6,*) "In find_doubly_excited, Number of initiators:",n_ref
      call flush(6)
    endif

! Sample upto 100 variational determinants to estimate the number of connected determinants
    if (present(eps_var_pt)) then
      if (n_ref<=1) then
        if (n_ref>0) then
          n_allocate=int(n_connected_dets_hf,i8b)+size(ref_up,1,i8b)
        else
          n_allocate=1
        endif
      else
        if (present(core_up)) then
          n_allocate = estimate_n_connections(n_ref,ref_up,ref_dn,ref_coeffs,eps_var_pt,core_up,core_dn,virt_up,virt_dn,active_only)
        else
          n_allocate = estimate_n_connections(n_ref,ref_up,ref_dn,ref_coeffs,eps_var_pt)
        endif
       !n_allocate=nint(0.5*real(n_connected_dets)*n_ref*(1.0_rk-ref_coeffs(1)**2),i8b)+size(ref_up,1) ! This is to estimate the effective number of reference dets
      endif
    else
      n_allocate=int(n_connected_dets_hf,i8b)*int(n_ref,i8b)+size(ref_up,1,i8b)
    endif
    if (ipr.ge.1)  write(6,'(''allocating dets_up etc. arrays of size'',i10,'' ='',es11.4)') n_allocate, real(n_allocate)
    write(6,'(''find_doubly_excited: allocating dets_up etc. arrays of size'',i10,'' ='',es11.4)') n_allocate, real(n_allocate)
    call flush(6)
    allocate(dets_up(n_allocate))
    allocate(dets_dn(n_allocate))
    if (n_core_orb > 0) then
      allocate(old_dets_up(n_allocate))
      allocate(old_dets_dn(n_allocate))
    endif

    if (present(e_mix_num)) then
      allocate(e_mix_num(n_allocate))
    endif
    if (present(e_mix_den)) then
      allocate(e_mix_den(n_allocate))
    endif
    if (present(ref_diag_elems)) then
      allocate(diag_elems_info(n_allocate))
    endif
    if (present(term1_big)) then
      allocate(term1_big(n_allocate))
    endif
    if (present(term2_big)) then
      allocate(term2_big(n_allocate))
    endif
    if (n_ref.ne.size(ref_up,1)) then
      ! add in the previous iteration's dets
      n_det = size(ref_up,1) - n_ref
      dets_up(1:n_det) = ref_up(n_ref+1:n_ref+n_det)
      dets_dn(1:n_det) = ref_dn(n_ref+1:n_ref+n_det)
    endif

    write(6,'(''Initially allocating (for merging) iorder etc. arrays of size'',i10,'' ='',es11.4,i4)') n_allocate, real(n_allocate)
    call flush(6)
    allocate(iorder(n_allocate))
    allocate(temp_i16_up((n_allocate+1)/2))
    allocate(temp_i16_dn((n_allocate+1)/2))
    allocate(temp_i_2((n_allocate+1)/2))

    n_det = 0
    modnum = n_ref/ncores !MJO we want every core to do the same number of mpi_sort calls
    num_mpi_sort = 0
    if (modnum<1) then
      modnum = 1
    endif

! Find the dets connected to all the reference dets
    do i=1,n_ref
       if(ipr.ge.1) write(6,'(''In find_doubly_excited, n_ref,i='',9i9)') n_ref,i
       if(ipr.ge.1) call system("free -m | grep Mem")
       !call system("free -m | grep Mem >> mem_output")
       call flush(6)

       !find determinants connected to this determinant by the hamiltonian
       if (present(e_mix_num)) then
         if (hamiltonian_type.eq.'hubbardk') then
           call find_connected_dets_hubbard_k(ref_up(i), ref_dn(i), n_connected_dets, connected_dets_up, connected_dets_dn, norb, connected_matrix_elements)
         elseif (hamiltonian_type.eq.'chem') then
           if (present(eps_var_pt)) then
             if (abs(ref_coeffs(i)).ne.0._rk) then
               if (present(core_up)) then
                 if (present(ref_diag_elems)) then
                   if (present(eps_var_pt_big)) then
                     call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements, ref_diag_elems(i), connected_diag_elems_info, eps_var_pt_big/abs(ref_coeffs(i)), connected_matrix_elements_big,core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn,active_only=active_only)
                   else
                     call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements, ref_diag_elems(i), connected_diag_elems_info,core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn,active_only=active_only)
                   endif
                 else
                   if (present(eps_var_pt_big)) then
                     call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements, eps_big=eps_var_pt_big/abs(ref_coeffs(i)), matrix_elements_big=connected_matrix_elements_big,core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn,active_only=active_only)
                   else
                     call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements,core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn,active_only=active_only)
                   endif
                 endif
               else
                 if (present(ref_diag_elems)) then
                   if (present(eps_var_pt_big)) then
                     call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements, ref_diag_elems(i), connected_diag_elems_info, eps_var_pt_big/abs(ref_coeffs(i)), connected_matrix_elements_big)
                   else
                     call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements, ref_diag_elems(i), connected_diag_elems_info)
                   endif
                 else
                   if (present(eps_var_pt_big)) then
                     call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements, eps_big=eps_var_pt_big/abs(ref_coeffs(i)), matrix_elements_big=connected_matrix_elements_big)
                   else
                     call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements)
                   endif
                 endif
               endif
             endif
           else
             call find_connected_dets_chem(ref_up(i), ref_dn(i), n_connected_dets, connected_dets_up, connected_dets_dn,norb,connected_matrix_elements)
           endif
         elseif (hamiltonian_type.eq.'heg') then
           if (present(eps_var_pt)) then
             if (abs(ref_coeffs(i)).ne.0._rk) then
               if (present(ref_diag_elems)) then
                 if (present(eps_var_pt_big)) then
                   call find_important_connected_dets_heg(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements, ref_diag_elems(i), connected_diag_elems_info, eps_var_pt_big/abs(ref_coeffs(i)), connected_matrix_elements_big)
                 else
                   call find_important_connected_dets_heg(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements, ref_diag_elems(i), connected_diag_elems_info)
                 endif
               else
                 if (present(eps_var_pt_big)) then
                   call find_important_connected_dets_heg(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements, eps_big=eps_var_pt_big/abs(ref_coeffs(i)), matrix_elements_big=connected_matrix_elements_big)
                 else
                   call find_important_connected_dets_heg(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements)
                 endif
               endif
             endif
           else
             call find_connected_dets_heg(ref_up(i), ref_dn(i), n_connected_dets, connected_dets_up, connected_dets_dn, connected_matrix_elements)
           endif
         else
           stop 'This hamiltonian type has not been coded (only hubbardk, chem, and heg)'
         endif
       else ! e_mix_num not present
         if (hamiltonian_type.eq.'hubbardk') then
           call find_connected_dets_hubbard_k(ref_up(i),ref_dn(i),n_connected_dets,connected_dets_up,connected_dets_dn,norb)
         elseif (hamiltonian_type.eq.'chem') then
           if (present(eps_var_pt)) then
             if (present(min_H_already_done)) then
               if (abs(ref_coeffs(i))*min_H_already_done(i)>eps_var_pt) then ! If the absolute value of the variational coef has not gone up, it will not generate new connections
                is_not_done=is_not_done+1
                !write (6,*) "Generating connections! i,ref_coeffs(i),eps_var_pt,ref_coeffs(i)/eps_var_pt,min_H_already_done=",i,ref_coeffs(i),eps_var_pt,abs(ref_coeffs(i))/eps_var_pt,min_H_already_done(i)
                 if (present(core_up)) then
                   if (present(ref_diag_elems)) then
!                    call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, ref_diag_elem=ref_diag_elems(i), diag_elems_info=connected_diag_elems_info)
                     call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, ref_diag_elem=ref_diag_elems(i), diag_elems_info=connected_diag_elems_info, min_H_already_done_elem=min_H_already_done(i),core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn,active_only=active_only)
                   else
!                    call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn)
                     call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, min_H_already_done_elem=min_H_already_done(i),core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn,active_only=active_only)
                   endif
                 else
                   if (present(ref_diag_elems)) then
!                    call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, ref_diag_elem=ref_diag_elems(i), diag_elems_info=connected_diag_elems_info)
                     call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, ref_diag_elem=ref_diag_elems(i), diag_elems_info=connected_diag_elems_info, min_H_already_done_elem=min_H_already_done(i))
                   else
!                    call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn)
                     call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, min_H_already_done_elem=min_H_already_done(i))
                   endif
                 endif
                !min_H_already_done(i) = eps_var_pt/abs(ref_coeffs(i))
               else
                 is_done=is_done+1
                 n_connected_dets = 1
                 connected_dets_up(1) = ref_up(i)
                 connected_dets_dn(1) = ref_dn(i)
              !  write (6,*) "Skipping! i,ref_coeffs(i),eps_var_pt,ref_coeffs(i)/eps_var_pt,min_H_already_done=",i,ref_coeffs(i),eps_var_pt,abs(ref_coeffs(i))/eps_var_pt,min_H_already_done(i)
               endif
             else ! min_H_already_done not present
               if (abs(ref_coeffs(i)).ne.0._rk) then
                 if (present(core_up)) then
                   if (present(ref_diag_elems)) then
                     call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, ref_diag_elem=ref_diag_elems(i), diag_elems_info=connected_diag_elems_info,core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn,active_only=active_only)
                   else
                     call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn,core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn,active_only=active_only)
                   endif
                 else
                   if (present(ref_diag_elems)) then
                     call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, ref_diag_elem=ref_diag_elems(i), diag_elems_info=connected_diag_elems_info)
                   else
                     call find_important_connected_dets_chem(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn)
                   endif
                 endif
               endif
             endif
           else
             call find_connected_dets_chem(ref_up(i), ref_dn(i), n_connected_dets, connected_dets_up, connected_dets_dn, norb)
           endif
         elseif (hamiltonian_type.eq.'heg') then
           if (present(eps_var_pt)) then
             if (present(min_H_already_done)) then
               if (abs(ref_coeffs(i))*min_H_already_done(i)>eps_var_pt) then
                !write (6,*) "Generating connections! i,ref_coeffs(i),eps_var_pt,ref_coeffs(i)/eps_var_pt,min_H_already_done=",i,ref_coeffs(i),eps_var_pt,abs(ref_coeffs(i))/eps_var_pt,min_H_already_done(i)
                 if (present(ref_diag_elems)) then
                   call find_important_connected_dets_heg(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, ref_diag_elem=ref_diag_elems(i), diag_elems_info=connected_diag_elems_info)
                 else
                   call find_important_connected_dets_heg(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn)
                 endif
                !min_H_already_done(i) = eps_var_pt/abs(ref_coeffs(i))
               else
                 n_connected_dets = 1
                 connected_dets_up(1) = ref_up(i)
                 connected_dets_dn(1) = ref_dn(i)
              !  write (6,*) "Skipping! i,ref_coeffs(i),eps_var_pt,ref_coeffs(i)/eps_var_pt,min_H_already_done=",i,ref_coeffs(i),eps_var_pt,abs(ref_coeffs(i))/eps_var_pt,min_H_already_done(i)
               endif
             else ! min_H_already_done not present
               if (abs(ref_coeffs(i)).ne.0._rk) then
                 if (present(ref_diag_elems)) then
                   call find_important_connected_dets_heg(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn, ref_diag_elem=ref_diag_elems(i), diag_elems_info=connected_diag_elems_info)
                 else
                   call find_important_connected_dets_heg(ref_up(i), ref_dn(i), eps_var_pt/abs(ref_coeffs(i)), n_connected_dets, connected_dets_up, connected_dets_dn)
                 endif
               endif
             endif
           else
             call find_connected_dets_heg(ref_up(i), ref_dn(i), n_connected_dets, connected_dets_up, connected_dets_dn)
           endif
         else
           stop 'This hamiltonian type has not been coded (only hubbardk, chem, and heg)'
         endif
       endif
       if(ipr.ge.1) write(6,'(''n_det,n_connected_dets,n_allocate,size(connected_dets_up)='',9i11)') n_det,n_connected_dets,n_allocate,size(connected_dets_up)
       call flush(6)
       !MJO We want to do a parallel sort / merge every 1000 determinants to save memory

       if (present(e_mix_num).and.present(e_mix_den).and.present(ref_diag_elems).and..not.present(term1_big)) then
         if (mod(i,modnum).eq.0.and.ncores>1) then !MJO TEMP PARALLEL
           num_mpi_sort = num_mpi_sort + 1
           if (num_mpi_sort.le.ncores) then
             write(6,'(''Performing'',i4,'' of'',i4,'' mpi sort/merge of size'',i11,'' ='',es11.4)') num_mpi_sort, ncores, n_det, real(n_det) ; call flush(6)
             !We only want to do ncores. If the mod didn't work, exit here
             do j=1,n_det
               iorder(j)=j
             enddo
             !MJO Local sort and merge
             ! sort the determinants
             call merge_sort2_up_dn(dets_up(1:n_det),dets_dn(1:n_det), iorder, n_det, temp_i16_up, temp_i16_dn, temp_i_2)

             ! merge the determinants.  The e_mix_num, e_mix_den, term1_big, term2_big contributions get added.
             e_mix_num(1:n_det) = e_mix_num(iorder(1:n_det))
             e_mix_den(1:n_det) = e_mix_den(iorder(1:n_det))
             diag_elems_info(1:n_det) = diag_elems_info(iorder(1:n_det))
             call merge_original_with_spawned3(n_det, dets_up, dets_dn, e_mix_num, e_mix_den, diag_elems_info)
             !MJO Remote sort and merge
             call mpi_merge_sort2(n_det,dets_up,dets_dn,n_allocate,e_mix_num,e_mix_den,diag_elems_info,iorder,temp_i16_up,temp_i16_dn,temp_i_2)
           endif
         endif
       endif


       if (n_det+n_connected_dets>n_allocate) then ! Reduce the number by sorting and merging

         if (ipr.ge.1)  write(6,'(''allocating (for merging) iorder etc. arrays of size'',i10,'' ='',es11.4)') n_det, real(n_det)
         write(6,'(''find_doubly_excited: Performing sort and merge of size '',i10,'' ='',es11.4)') n_det, real(n_det); call flush(6)
         call flush(6)
         do j=1,n_det
           iorder(j)=j
         enddo
         !MJO Local sort and merge

         ! sort the determinants
         call merge_sort2_up_dn(dets_up(1:n_det),dets_dn(1:n_det), iorder, n_det, temp_i16_up, temp_i16_dn, temp_i_2)

         ! merge the determinants.  The e_mix_num, e_mix_den, term1_big, term2_big contributions get added.
         if (present(e_mix_num)) then
           e_mix_num(1:n_det) = e_mix_num(iorder(1:n_det))
           if (present(e_mix_den)) then
             e_mix_den(1:n_det) = e_mix_den(iorder(1:n_det))
             if (present(ref_diag_elems)) then
               diag_elems_info(1:n_det) = diag_elems_info(iorder(1:n_det))
               if (present(term1_big)) then
                 term1_big(1:n_det) = term1_big(iorder(1:n_det))
                 term2_big(1:n_det) = term2_big(iorder(1:n_det))
                 call merge_original_with_spawned3(n_det, dets_up, dets_dn, e_mix_num, e_mix_den, diag_elems_info, term1_big, term2_big)
               else
                 call merge_original_with_spawned3(n_det, dets_up, dets_dn, e_mix_num, e_mix_den, diag_elems_info)
               endif
             else
               call merge_original_with_spawned3(n_det, dets_up, dets_dn, e_mix_num, e_mix_den)
             endif
           else
             call merge_original_with_spawned3(n_det, dets_up, dets_dn, e_mix_num)
           endif
         else
           call merge_original_with_spawned3(n_det, dets_up, dets_dn)
         endif

         if (n_det+n_connected_dets>n_allocate) then ! If there are too many dets even after merging, then reallocate

           if (present(eps_var_pt)) then
             n_allocate_new = n_allocate+n_det+n_connected_dets
           else
             n_allocate_new = n_allocate+int(n_connected_dets,i8b)*(n_ref-i+1)
           endif
          !if (ipr.ge.1)  write(6,'(''reallocating dets_up etc. arrays of size'',i10,'' ='',es11.4)') n_allocate_new, real(n_allocate_new)
           write(6,'(''find_doubly_excited: reallocating dets_up etc. arrays of size'',i10,'' ='',es11.4)') n_allocate_new, real(n_allocate_new)
           call flush(6)

           call alloc('iorder',iorder,n_allocate_new)
           call alloc('temp_i16_up',temp_i16_up,(n_allocate_new+1)/2)
           call alloc('temp_i16_dn',temp_i16_dn,(n_allocate_new+1)/2)
           call alloc('temp_i_2',temp_i_2,(n_allocate_new+1)/2)

           call alloc('dets_up',dets_up,n_allocate_new)
           call alloc('dets_dn',dets_dn,n_allocate_new)

           if (present(e_mix_num)) then
             call alloc('e_mix_num',e_mix_num,n_allocate_new)
           endif

           if (present(e_mix_den)) then
             call alloc('e_mix_den',e_mix_den,n_allocate_new)
           endif

           if (present(ref_diag_elems)) then
             call alloc('diag_elems_info',diag_elems_info,n_allocate_new)
           endif

           if (present(term1_big)) then
             call alloc('term1_big',term1_big,n_allocate_new)
           endif

           if (present(term2_big)) then
             call alloc('term2_big',term2_big,n_allocate_new)
           endif

           n_allocate = n_allocate_new

           if (n_core_orb > 0) then
             call alloc('old_dets_up',old_dets_up,n_allocate)
             call alloc('old_dets_dn',old_dets_dn,n_allocate)
           endif

           if (ipr.ge.1)  write(6,'(''done reallocating dets_up etc. arrays'')')
           call flush(6)

         endif

       endif


       if (run_type.eq.'selected_ci') then
         if (abs(ref_coeffs(i)).ge.initiator) then
           dets_up(n_det+1:n_det+n_connected_dets) = connected_dets_up(1:n_connected_dets)
           dets_dn(n_det+1:n_det+n_connected_dets) = connected_dets_dn(1:n_connected_dets)
         else
           dets_up(n_det+1) = ref_up(i)
           dets_dn(n_det+1) = ref_dn(i)
           n_connected_dets = 1
         endif
       else
        !dets_up(n_det+1) = ref_up(i)
        !dets_dn(n_det+1) = ref_dn(i)
         dets_up(n_det+1:n_det+n_connected_dets) = connected_dets_up(1:n_connected_dets)
         dets_dn(n_det+1:n_det+n_connected_dets) = connected_dets_dn(1:n_connected_dets)
       endif


       if (present(e_mix_num)) then
         if (present(w_over_p)) then
      ! term1(k) = sum_i^{N_MC_diff} H_{ki} c_i w_i / p_i
           e_mix_num(n_det + 1:n_det + n_connected_dets) = connected_matrix_elements(1:n_connected_dets)*ref_coeffs(i)*w_over_p(i)
           if (present(term1_big)) then
             term1_big(n_det + 1:n_det + n_connected_dets) = connected_matrix_elements_big(1:n_connected_dets)*ref_coeffs(i)*w_over_p(i)
           endif ! present(term1_big)
         else
           e_mix_num(n_det + 1:n_det + n_connected_dets) = connected_matrix_elements(1:n_connected_dets)*ref_coeffs(i)
         endif
       endif
       if (present(e_mix_den)) then
         if (present(w_over_p)) then
           ! term2(k) = sum_i^{N_MC_diff} (H_{ki} c_i)**2 * ( (n_mc-1) * w_i / p_i - w_i**2 / p_i**2 )
           e_mix_den(n_det + 1:n_det + n_connected_dets) = (connected_matrix_elements(1:n_connected_dets)*ref_coeffs(i))**2 * ( (n_mc-1)*w_over_p(i) - w_over_p(i)**2 )
           if (present(term2_big)) then
             term2_big(n_det + 1:n_det + n_connected_dets) = (connected_matrix_elements_big(1:n_connected_dets)*ref_coeffs(i))**2 * ( (n_mc-1)*w_over_p(i) - w_over_p(i)**2 )
           endif ! present(term2_big)
         else
           e_mix_den(n_det + 1:n_det + n_connected_dets) = 0._rk
           e_mix_den(n_det + 1) = ref_coeffs(i) ! Diagonal "connected det"
         endif
       endif
       if (present(ref_diag_elems)) then
         diag_elems_info(n_det + 1:n_det + n_connected_dets) = connected_diag_elems_info(1:n_connected_dets)
       endif

       n_det = n_det + n_connected_dets

    enddo ! i=1,n_ref
    if (.not.present(e_mix_num)) write(6,'(''is_done, is_not_done'',9i10)') is_done, is_not_done, is_done+is_not_done

    !MJO Local sort and merge
    if (n_ref > 0) then
      write(6,'(''Performing final sort/merge of size         '',i11,'' ='',es11.4)') n_det, real(n_det)
      call flush(6)
      if (present(e_mix_num)) then
        do j=1,n_det
          iorder(j)=j
        enddo
      endif
    endif

    call merge_sort2_up_dn(dets_up(1:n_det),dets_dn(1:n_det), iorder, n_det, temp_i16_up, temp_i16_dn, temp_i_2)

    ! merge the determinants.  The e_mix_num, e_mix_den, term1_big, term2_big contributions get added.
    if (present(e_mix_num)) then
      e_mix_num(1:n_det) = e_mix_num(iorder(1:n_det))
      if (present(e_mix_den)) then
        e_mix_den(1:n_det) = e_mix_den(iorder(1:n_det))
        if (present(ref_diag_elems)) then
          diag_elems_info(1:n_det) = diag_elems_info(iorder(1:n_det))
          if (present(term1_big)) then
            term1_big(1:n_det) = term1_big(iorder(1:n_det))
            term2_big(1:n_det) = term2_big(iorder(1:n_det))
            if (n_det > 0 ) then
               call merge_original_with_spawned3(n_det, dets_up, dets_dn, e_mix_num, e_mix_den, diag_elems_info, term1_big, term2_big)
            endif
            !MJO To use global populations, we do an mpi_merge_sort here.
            !    If we wanted to do local, independent populations, we would not
            !    merge_sort here
            if (ncores>1) then !MJO TEMP PARALLEL
               call mpi_merge_sort2(n_det,dets_up,dets_dn,n_allocate,e_mix_num,e_mix_den,diag_elems_info,term1_big,term2_big,iorder,temp_i16_up,temp_i16_dn,temp_i_2)
            endif
          else
            if (n_det > 0) then
               call merge_original_with_spawned3(n_det, dets_up, dets_dn, e_mix_num, e_mix_den, diag_elems_info)
            endif
            !MJO Remote sort and merge
            if (ncores>1) then !MJO TEMP PARALLEL
               call mpi_merge_sort2(n_det,dets_up,dets_dn,n_allocate,e_mix_num,e_mix_den,diag_elems_info,iorder,temp_i16_up,temp_i16_dn,temp_i_2)
            endif
          endif
        else
          call merge_original_with_spawned3(n_det, dets_up, dets_dn, e_mix_num, e_mix_den)
        endif
      else
!        !MJO FIXME: Add remote sort and merge
        call merge_original_with_spawned3(n_det, dets_up, dets_dn, e_mix_num)
      endif
    else
       if (n_det > 0) then
          call merge_original_with_spawned3(n_det, dets_up, dets_dn)
       endif
       !MJO Remote sort and merge
      if (ncores>1) then
         call mpi_merge_sort2(n_det,dets_up,dets_dn,n_allocate,iorder,temp_i16_up,temp_i16_dn,temp_i_2)
      endif
    endif

    ! Go through list and weed out dets without the core filled.
    if (n_core_orb>0) then
      n_det_core_filled = 0
      do i=1,n_det
        core_filled=.true.
        do j=1,n_core_orb
          if (btest(dets_up(i),j-1).eqv..false..or.btest(dets_dn(i),j-1).eqv..false.)  core_filled = .false.
        enddo
        if (core_filled) then
          n_det_core_filled = n_det_core_filled + 1
          old_dets_up(n_det_core_filled) = dets_up(i)
          old_dets_dn(n_det_core_filled) = dets_dn(i)
        endif
      enddo
      deallocate(dets_up)
      deallocate(dets_dn)
      allocate(dets_up(n_det_core_filled))
      allocate(dets_dn(n_det_core_filled))
      do i=1,n_det_core_filled
        dets_up(i) = old_dets_up(i)
        dets_dn(i) = old_dets_dn(i)
      enddo
      n_det = n_det_core_filled
    endif

    if (imp) then
      ! store psi_g to make things easier
      ! e_num -> e_num/psi_g
      ! e_den -> e_den/psi_g
      psi_g_epsilon_inv = 1._rk/psi_g_epsilon
      allocate(psi_g(n_det))
      do i=1,n_det
        psi_g(i) = e_mix_num(i)/psi_g_energy
        if (abs(e_mix_num(i))<1.e-10)  write (6,'(''i,emix_num='',i6,es12.4)') i,e_mix_num(i)
        e_mix_num(i) = psi_g_energy
        e_mix_den(i) = e_mix_den(i) / psi_g(i) ! Probably not necessary since this should be zero, but since we take sign(+/- 0), weird things can happen.
      enddo
      do i=1,ndet_psi_t
        call binary_search(dets_up_psi_t(i),dets_dn_psi_t(i),dets_up(1:n_det),dets_dn(1:n_det),j)
        if (j>0) then
          e_mix_num(j) = e_mix_num(j) * psi_g(j) / cdet_psi_t(i)
          e_mix_den(j) = e_mix_den(j) * psi_g(j) / cdet_psi_t(i)
          psi_g(j) = cdet_psi_t(i)
        endif
      enddo
    endif

    ! Sort and merge has already happened at this point

    if (present(ref_diag_elems)) then
      allocate(new_diag_elems(n_det))
      do i=1,n_det
        if (diag_elems_info(i)%old_diag_elem>1.e50_rk) then ! single excitation; compute the easy way for now
          call hamiltonian(dets_up(i),dets_dn(i),dets_up(i),dets_dn(i),new_diag_elems(i))
        else
          if (hamiltonian_type .eq. 'chem') then
            call get_new_diag_elem(diag_elems_info(i),dets_up(i),dets_dn(i),new_diag_elems(i))
          elseif (hamiltonian_type .eq. 'heg') then
            ! In HEG, this cost is negligible comparing to finding important connected dets.
            call hamiltonian(dets_up(i), dets_dn(i), dets_up(i), dets_dn(i), new_diag_elems(i))
            ! call get_new_diag_elem_heg(diag_elems_info(i),dets_up(i),dets_dn(i),new_diag_elems(i))
          endif
        endif
      enddo
      if (present(pt)) then
        ! Divide above weight by (lowest_energy - H_ii)
        if (n_ref==1) then
          e_mix_num(1) = 1._rk
        else
          e_mix_num(1) = e_mix_num(1) / (pt - new_diag_elems(1)) + e_mix_den(1)
        endif
        do i=2,n_det
          e_mix_num(i) = e_mix_num(i) / (pt - new_diag_elems(i)) ! 1st order PT in wavefunction coefficients
          e_mix_num(i) = e_mix_num(i) + e_mix_den(i) ! Add the PT correction to the previous iteration's coeff
        enddo
      endif

    else ! not (present(ref_diag_elems))
      if (present(pt)) then
        ! Divide above weight by (lowest_energy - H_ii)
        if (n_ref==1) then
          e_mix_num(1) = 1._rk
        else
          call hamiltonian(dets_up(1),dets_dn(1),dets_up(1),dets_dn(1),H_ii)
          e_mix_num(1) = e_mix_num(1) / (pt - H_ii) + e_mix_den(1)
        endif
        do i=2,n_det
          call hamiltonian(dets_up(i),dets_dn(i),dets_up(i),dets_dn(i),H_ii)
          e_mix_num(i) = e_mix_num(i) / (pt - H_ii) ! 1st order PT in wavefunction coefficients
          e_mix_num(i) = e_mix_num(i) + e_mix_den(i) ! Add the PT correction to the previous iteration's coeff
        enddo
      endif
    endif

    if (n_core_orb > 0) then
      deallocate(old_dets_up, old_dets_dn)
    endif

  end subroutine find_doubly_excited

!-------------------------------------------------------------------------------
  subroutine hamiltonian(up1,dn1,up2,dn2,element,connected)
  ! Finds whether states are connected.
  ! If so, calculates matrix element.
  ! A Holmes, 19 Jul 2012.

  use chemistry, only  : hamiltonian_chem,hamiltonian_chem_time_sym,time_sym,is_connected_chem,excitation_level
  use hubbard, only    : hamiltonian_hubbard_k,hamiltonian_hubbard_k_space_sym,space_sym,is_connected_hubbard_fast
  use heg, only        : hamiltonian_heg,is_connected_heg

  implicit none

  real(rk),intent(out) :: element
  logical,optional,intent(out) :: connected
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: up1,dn1,up2,dn2
  type(ik_vec) :: up1_tmp,up2_tmp,dn1_tmp,dn2_tmp
#else
  integer(ik),intent(in) :: up1,dn1,up2,dn2
  integer(ik) :: up1_tmp,up2_tmp,dn1_tmp,dn2_tmp
#endif
  integer :: excite_level,nnzero

    up1_tmp=up1
    up2_tmp=up2
    dn1_tmp=dn1
    dn2_tmp=dn2

    if (present(connected))  connected=.false.
 
    if (hamiltonian_type.eq.'hubbardk') then
      if (space_sym) then
        call hamiltonian_hubbard_k_space_sym(up1_tmp,dn1_tmp,up2_tmp,dn2_tmp,element,nnzero)
        if (present(connected)) then
          if (abs(element) .gt. 1.0e-10) connected=.true.
        endif
      else
        if (is_connected_hubbard_fast(up1_tmp,dn1_tmp,up2_tmp,dn2_tmp)) then
          call hamiltonian_hubbard_k(up1_tmp,dn1_tmp,up2_tmp,dn2_tmp,element,.true.)
          if (present(connected))  connected=.true.
        endif
      endif
    elseif (hamiltonian_type.eq.'chem') then
      element = 0._rk
      if (present(connected)) then
        if (time_sym) then
          call is_connected_chem(up1_tmp,dn1_tmp,up2_tmp,dn2_tmp,connected,excite_level)
          if (.not.connected)  call is_connected_chem(up1_tmp,dn1_tmp,dn2_tmp,up2_tmp,connected,excite_level)
          if (connected)  call hamiltonian_chem_time_sym(up1_tmp, dn1_tmp, up2_tmp, dn2_tmp, element)
        else
          call is_connected_chem(up1_tmp,dn1_tmp,up2_tmp,dn2_tmp,connected,excite_level)
          if (excite_level >= 0) call hamiltonian_chem(up1_tmp, dn1_tmp, up2_tmp, dn2_tmp, excite_level, element)
        endif
      else
        element = 0._rk
        if (time_sym) then
          call hamiltonian_chem_time_sym(up1_tmp, dn1_tmp, up2_tmp, dn2_tmp, element)
        else
          call excitation_level(up1_tmp, dn1_tmp, up2_tmp, dn2_tmp, excite_level)
          if (excite_level >= 0) call hamiltonian_chem(up1_tmp, dn1_tmp, up2_tmp, dn2_tmp, excite_level, element)
        endif
      endif
    elseif (hamiltonian_type.eq.'heg') then
      if (is_connected_heg(up1_tmp,dn1_tmp,up2_tmp,dn2_tmp)) then
        call hamiltonian_heg(up1_tmp,dn1_tmp,up2_tmp,dn2_tmp,element)
        if (present(connected))  connected=.true.
      endif
    else
      stop 'This hamiltonian type has not been coded (only hubbardk, chem, and heg)'
    endif

  end subroutine hamiltonian

  function dotproduct(v1,up1,dn1,v2,up2,dn2)
  ! Dot product of v1 and v2
  ! Assumes that at least the longer of the two vectors is
  ! already sorted by label.
  ! v1, v2 are vectors of coefficients
  ! up/dn1, up/dn2 are vectors of determinant labels
  ! A Holmes, 7 Aug 2012

    use more_tools, only : binary_search

    implicit none

    real(rk),intent(in) :: v1(:),v2(:)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: up1(:),dn1(:),up2(:),dn2(:)
#else
    integer(ik),intent(in) :: up1(:),dn1(:),up2(:),dn2(:)
#endif
    real(rk) :: dotproduct

    integer :: i,j

    dotproduct = 0._rk

    if (size(v1)<size(v2)) then
      do i=1,size(v1)
        call binary_search(up1(i),dn1(i),up2,dn2,j)
        if (j>0) then
          dotproduct=dotproduct+v1(i)*v2(j)
        endif
      enddo
    else
      do j=1,size(v2)
        call binary_search(up2(j),dn2(j),up1,dn1,i)
        if (i>0) then
          dotproduct=dotproduct+v1(i)*v2(j)
        endif
      enddo
    endif

  end function dotproduct
!------------------------------------------------------------------------------------------------------------------------------

  subroutine addvectors(v1,up1,dn1,v2,up2,dn2)
  ! Add the elements of v2 in common with v1 to v1.
  ! Assumes that v1 is sorted by label.
  ! v1, v2 are vectors of coefficients
  ! up/dn1, up/dn2 are vectors of determinant labels
  ! Not presently being used.
  ! A Holmes, 7 Aug 2012

    use more_tools, only : binary_search

    implicit none

    real(rk),intent(inout) :: v1(:)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(inout) :: up1(:),dn1(:)
    type(ik_vec),intent(in) :: up2(:),dn2(:)
#else
    integer(ik),intent(inout) :: up1(:),dn1(:)
    integer(ik),intent(in) :: up2(:),dn2(:)
#endif
    real(rk),intent(in) :: v2(:)

    integer :: i,j

    do j=1,size(v2)
      call binary_search(up2(j),dn2(j),up1,dn1,i)
      if (i>0)  v1(i)=v1(i)+v2(j)
    enddo

  end subroutine addvectors
!------------------------------------------------------------------------------------------------------------------------------

  subroutine orthogonalize(vnew,upnew,dnnew,v,up,dn,len_v,A,x,b,ipiv)
  ! Orthogonalize vnew with respect to all columns of v.
  ! This is a modified version of Gram-Schmidt procedure,
  ! where the number of nonzero elements in vnew stays the same.
  ! Assumes that all columns of v are sorted by label.
  ! A Holmes, 21 Aug 2012

  use more_tools, only : binary_search
  use tools, only : matrix_inversion

  implicit none

  real(rk),intent(inout) :: vnew(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: upnew(:),dnnew(:)
  type(ik_vec),intent(in) :: up(:,:),dn(:,:)
#else
  integer(ik),intent(in) :: upnew(:),dnnew(:)
  integer(ik),intent(in) :: up(:,:),dn(:,:)
#endif
  real(rk),intent(in) :: v(:,:)
  integer,intent(in) :: len_v(:)
  real(rk),intent(inout) :: A(:,:),x(:),b(:)
  integer,intent(inout) :: ipiv(:)

  integer :: i,j,k,m,n,ctr
  real(rk),allocatable :: A_tmp(:,:),work(:)
  real(rk),parameter :: epsilon=1e-17_rk
  integer,allocatable :: ind(:)

  n = size(v,2)
  allocate(ind(n))

  allocate(A_tmp(n,n))

  A(:,:) = 0._rk
  b(:) = 0._rk

  ! compute A,b
  ctr = 0
  do i=1,n
    do k=1,len_v(i)
      call binary_search(up(k,i),dn(k,i),upnew(:),dnnew(:),m)
      if (m>0) then
        b(i) = b(i) + v(k,i)*vnew(m)
        do j=i,n
          call binary_search(up(k,i),dn(k,i),up(1:len_v(j),j),dn(1:len_v(j),j),m)
          if (m>0) then
            if (i.eq.j) then
              A(i,j)=A(i,j)+v(k,i)*v(m,j)
            else
              A(i,j)=A(i,j)+v(k,i)*v(m,j)
              A(j,i)=A(j,i)+v(k,i)*v(m,j)
            endif
          endif
        enddo
      endif
    enddo
  enddo
  do i=1,n
    A(:,i) = A(:,i) * b(i)
    if (abs(b(i))<epsilon) then
      ctr = ctr+1
      ind(ctr) = i ! indices of vectors already orthogonal to current vector
    endif
  enddo
 !write (6,*) "b=",b
 !write (6,*) "A=",A

  if (ctr==n) then
    write (6,*) "Still orthogonal after truncation"
    x(:)=0._rk
    return
  endif
  if (ctr>0) then
    ! one or more old vectors are already orthogonal to current vector; remove these
    do i=1,ctr
      do j=ind(ctr)+1,n
        b(j-1) = b(j)
        A(:,j-1) = A(:,j)
      enddo
      do j=ind(ctr)+1,n
        A(j-1,:) = A(j,:)
      enddo
    enddo
  endif
  n = n - ctr
 !write (6,*) "b=",b(1:n)
 !write (6,*) "A=",A(1:n,1:n)

  ! compute x
  A_tmp(1:n,1:n)=A(1:n,1:n)
  allocate(work(n**2))
  call matrix_inversion(n,A(1:n,1:n),x(1:n),b(1:n),A_tmp(1:n,1:n),ipiv(1:n),work)
  deallocate(work)
  write (6,*) "x=",x

  ! Now put the 0's back in x
  n = n + ctr
  if (ctr>0) then
    ! one or more old vectors are already orthogonal to current vector; remove these
    do i=ctr,1,-1
      do j=n,ind(i)-i+2,-1
        x(j) = x(j-1)
        b(j) = b(j-1)
      enddo
      x(ind(i)-i+1) = 0._rk
      b(ind(i)-i+1) = 0._rk
    enddo
  endif

  ! orthogonalize
  do i=1,n
    write (6,*) "b(i),x(i)=",b(i),x(i)
    do k=1,len_v(i)
      call binary_search(up(k,i),dn(k,i),upnew(:),dnnew(:),m)
      if (m>0)  vnew(m)=vnew(m)-b(i)*x(i)*v(k,i)
    enddo
  enddo

  end subroutine orthogonalize

  integer function count_csfs(size_vector,vector)
  ! Counts the number of CSFs, where a CSF is assumed to be the set of all determinants
  ! whose coefficients are equal in magnitude (up to epsilon defined below)
  ! A Holmes, 23 May 2014

    integer,intent(in) :: size_vector
    real(rk),intent(in) :: vector(size_vector)
    integer :: i
    real(rk) :: previous_value
    real(rk), parameter :: epsilon = 1.e-10_rk ! Although the wavefn coefs. may not be converged to this accuracy, the agreement between symmetry related coefs. exists to higher precision

!   write(6,'(''size_vector='',i9)') size_vector
!   write(6,'(''vector='',100f14.10)') vector(1:size_vector)
    count_csfs = 0
    previous_value = 0._rk
    do i = 1,size_vector
      if ( abs(abs(previous_value) - abs(vector(i))) > epsilon ) then
        count_csfs=count_csfs+1
        previous_value = vector(i)
      endif
    enddo

  end function count_csfs
!-------------------------------------------------------------------------------------

  integer(i8b) function estimate_n_connections(ndets,dets_up,dets_dn,coeffs,eps_pt,core_up,core_dn,virt_up,virt_dn,active_only)
  ! Estimate the number of connections that exceed eps2 for the reference wavefunction
  ! A Holmes, 24 June 2016
  ! Modified : A Holmes, 23 May 2017. Improve estimate by taking more samples
  !            from quickly-varying regions

    use generic_sort, only : sort
    use chemistry, only : find_important_connected_dets_chem,time_sym
    use heg, only: find_important_connected_dets_heg
    use tools, only : sort_and_merge

    integer,intent(in) :: ndets
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
    type(ik_vec),optional,intent(in) :: core_up,core_dn,virt_up,virt_dn
    type(ik_vec),allocatable :: dets_up_sorted(:),dets_dn_sorted(:)
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
    integer(ik),optional,intent(in) :: core_up,core_dn,virt_up,virt_dn
    integer(ik),allocatable :: dets_up_sorted(:),dets_dn_sorted(:)
#endif
    real(rk),intent(in) :: coeffs(:),eps_pt
    logical,optional,intent(in) :: active_only
    real(rk),allocatable :: coeffs_sorted(:)

    integer :: n_connected_dets_left,n_connected_dets_right

      ! Can be zero variational dets if the variational dets are all on other
      ! cores
      if (ndets==0) then
        estimate_n_connections = 0_i8b
        return
      endif

      ! Sort by abs coeff first because connections can vary a lot and a random sample is very noisy!
      allocate(dets_up_sorted(ndets))
      allocate(dets_dn_sorted(ndets))
      allocate(coeffs_sorted(ndets))
      dets_up_sorted(1:ndets) = dets_up(1:ndets)
      dets_dn_sorted(1:ndets) = dets_dn(1:ndets)
      coeffs_sorted(1:ndets) = coeffs(1:ndets)
      call sort(coeffs_sorted, dets_up_sorted, dets_dn_sorted) ! Sort in order of decreasing absolute coef

      ! Now call recursive function on sorted lists
      if (hamiltonian_type .eq. 'chem') then
        if (present(core_up)) then
          call find_important_connected_dets_chem(dets_up_sorted(1), dets_dn_sorted(1), eps_pt/abs(coeffs_sorted(1)), n_connected_dets_left, connected_dets_up, connected_dets_dn, core_up=core_up, core_dn=core_dn, virt_up=virt_up, virt_dn=virt_dn, active_only=active_only)
        else
          call find_important_connected_dets_chem(dets_up_sorted(1), dets_dn_sorted(1), eps_pt/abs(coeffs_sorted(1)), n_connected_dets_left, connected_dets_up, connected_dets_dn)
        endif
        if (time_sym) then
          ! Sort and merge because there could be repeats
          call sort_and_merge(n_connected_dets_left,connected_dets_up,connected_dets_dn)
        endif
        if (present(core_up)) then
          call find_important_connected_dets_chem(dets_up_sorted(ndets), dets_dn_sorted(ndets), eps_pt/abs(coeffs_sorted(ndets)), n_connected_dets_right, connected_dets_up, connected_dets_dn, core_up=core_up, core_dn=core_dn, virt_up=virt_up, virt_dn=virt_dn, active_only=active_only)
        else
          call find_important_connected_dets_chem(dets_up_sorted(ndets), dets_dn_sorted(ndets), eps_pt/abs(coeffs_sorted(ndets)), n_connected_dets_right, connected_dets_up, connected_dets_dn)
        endif
        if (time_sym) then
          ! Sort and merge because there could be repeats
          call sort_and_merge(n_connected_dets_right,connected_dets_up,connected_dets_dn)
        endif
      elseif (hamiltonian_type .eq. 'heg') then
        call find_important_connected_dets_heg(dets_up_sorted(1), dets_dn_sorted(1), eps_pt/abs(coeffs_sorted(1)), n_connected_dets_left, connected_dets_up, connected_dets_dn)
        call find_important_connected_dets_heg(dets_up_sorted(ndets), dets_dn_sorted(ndets), eps_pt/abs(coeffs_sorted(ndets)), n_connected_dets_right, connected_dets_up, connected_dets_dn)
      endif
      write (6,*) "j,n_connected_dets=",1,n_connected_dets_left
      write (6,*) "j,n_connected_dets=",ndets,n_connected_dets_right
      if (present(core_up)) then
        estimate_n_connections = estimate_connections(ndets,dets_up_sorted,dets_dn_sorted,coeffs_sorted,eps_pt,n_connected_dets_left,n_connected_dets_right,0, core_up=core_up, core_dn=core_dn, virt_up=virt_up, virt_dn=virt_dn, active_only=active_only)
      else
        estimate_n_connections = estimate_connections(ndets,dets_up_sorted,dets_dn_sorted,coeffs_sorted,eps_pt,n_connected_dets_left,n_connected_dets_right,0) - ndets
      endif

  end function estimate_n_connections
!-------------------------------------------------------------------------------------

  recursive function estimate_connections(ndets,dets_up,dets_dn,coeffs,eps_pt,n_connected_dets_left,n_connected_dets_right,offset,core_up,core_dn,virt_up,virt_dn,active_only) result (answer)
  ! Estimate number of connections using a binary sampling to take more samples
  ! of quickly-varying regions
  ! Assumes input lists are sorted in decreasing order by abs coeff
  ! Also assumes that numbers of connections also decrease roughly monotonically
  ! with decreasing abs coeff
  ! A Holmes, 24 May 2017

    use chemistry, only : find_important_connected_dets_chem,time_sym
    use heg, only: find_important_connected_dets_heg
    use common_run, only : connected_dets_up,connected_dets_dn
    use tools, only : sort_and_merge

    integer,intent(in) :: ndets
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
    type(ik_vec),optional,intent(in) :: core_up,core_dn,virt_up,virt_dn
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
    integer(ik),optional,intent(in) :: core_up,core_dn,virt_up,virt_dn
#endif
    real(rk),intent(in) :: coeffs(:),eps_pt
    integer,intent(in) :: n_connected_dets_left,n_connected_dets_right
    integer,intent(in) :: offset
    logical,optional,intent(in) :: active_only

    integer :: left,right,mid
    integer :: n_connected_dets_mid
    real(rk),parameter :: max_variation = 2 ! Maximum variation before more samples are taken
    integer(i8b) :: answer

      ! If too few elements:
      if (ndets==1) then
        answer = n_connected_dets_left
        return
      elseif (ndets==2) then
        answer = n_connected_dets_left + n_connected_dets_right
        return
      endif

      ! Get n_connections at beginning, middle, and end of interval
      left = 1
      right = ndets
      mid = ndets/2
      if (hamiltonian_type .eq. 'chem') then
        if (present(core_up)) then
          call find_important_connected_dets_chem(dets_up(mid), dets_dn(mid), eps_pt/abs(coeffs(mid)), n_connected_dets_mid, connected_dets_up, connected_dets_dn, core_up=core_up, core_dn=core_dn, virt_up=virt_up, virt_dn=virt_dn, active_only=active_only)
        else
          call find_important_connected_dets_chem(dets_up(mid), dets_dn(mid), eps_pt/abs(coeffs(mid)), n_connected_dets_mid, connected_dets_up, connected_dets_dn)
        endif
        if (time_sym) then
          ! Sort and merge because there could be repeats
          call sort_and_merge(n_connected_dets_mid,connected_dets_up,connected_dets_dn)
        endif
      elseif (hamiltonian_type .eq. 'heg') then
        call find_important_connected_dets_heg(dets_up(mid), dets_dn(mid), eps_pt/abs(coeffs(mid)), n_connected_dets_mid, connected_dets_up, connected_dets_dn)
      endif
      write (6,*) "j,n_connected_dets=",mid+offset,n_connected_dets_mid

      ! Get better estimate of part of interval if it varies a lot and has large magnitudes
      ! Left side
      if ((real(n_connected_dets_left)/n_connected_dets_mid>max_variation.or.real(n_connected_dets_mid)/n_connected_dets_left>max_variation)) then !.and.(n_connected_dets_left+n_connected_dets_mid)/2>(mid-left)) then
        if (present(core_up)) then
          answer = estimate_connections(mid,dets_up(left:mid),dets_dn(left:mid),coeffs(left:mid),eps_pt,n_connected_dets_left,n_connected_dets_mid,offset,core_up=core_up, core_dn=core_dn, virt_up=virt_up, virt_dn=virt_dn, active_only=active_only)
        else
          answer = estimate_connections(mid,dets_up(left:mid),dets_dn(left:mid),coeffs(left:mid),eps_pt,n_connected_dets_left,n_connected_dets_mid,offset)
        endif
      else
        answer = int(int(n_connected_dets_left+n_connected_dets_mid,i8b)/2.0*mid,i8b)
      endif
      ! Right side
      if ((real(n_connected_dets_mid)/n_connected_dets_right>max_variation.or.real(n_connected_dets_right)/n_connected_dets_mid>max_variation)) then !.and.(n_connected_dets_right+n_connected_dets_mid)/2>(right-mid)) then
        if (present(core_up)) then
          answer = answer + estimate_connections(right-mid,dets_up(mid:right),dets_dn(mid:right),coeffs(mid:right),eps_pt,n_connected_dets_mid,n_connected_dets_right,offset+mid-1,core_up=core_up, core_dn=core_dn, virt_up=virt_up, virt_dn=virt_dn, active_only=active_only)
        else
          answer = answer + estimate_connections(right-mid,dets_up(mid:right),dets_dn(mid:right),coeffs(mid:right),eps_pt,n_connected_dets_mid,n_connected_dets_right,offset+mid-1)
        endif
      else
        answer = answer + int(int(n_connected_dets_mid+n_connected_dets_right,i8b)/2.0*(right-mid))
      endif

  end function estimate_connections


end module semistoch
