module hci
! Heat-bath Configuration Interaction
! Input two parameters: eps_var,eps_pt
! Starting with HF, each iteration,
! 1. Find lowest eigenvector in the current set of determinants.
! 2. Add all determinants i for which H_{ij}*c_j>eps_var for at least one
!    determinant j in current set.
! Repeat until det list increases by less than 0.001% on an iteration (used to be 1%) or energy changes by < 1e-5.
! Once the converged wavefunction is generated, perform second-order perturbation theory to
! improve the energy estimate, using eps_pt to avoid wasting time on small contributions
! to the sum.
! A Holmes, Started on 7 Feb 2016

use types, only : rk, ik, ik_vec, i8b, num_words
use tools, only : alloc, round_r, round_i
use common_ham, only : nelec, nup, ndn, norb, n_core_orb, hamiltonian_type
use common_selected_ci, only: n_max_connections
#ifdef NUM_ORBITALS_GT_127
!use overload, only : maskr_vec
use overload
#endif

implicit none
save
private
public :: perform_hci

  type neighbor_table
    integer :: n_hash_functions,n_bins_per_hash
    logical,allocatable :: is_empty(:,:)
    integer,allocatable :: bin_counts(:,:)
    integer,allocatable :: start_index(:,:),end_index(:,:)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),allocatable :: dets_up(:,:),dets_dn(:,:)
#else
    integer(ik),allocatable :: dets_up(:,:),dets_dn(:,:)
#endif
    real(rk),allocatable :: coeffs(:,:)
  end type neighbor_table

  integer :: RandomHash1(0:255) !,RandomHash2(0:255)
  integer,allocatable :: n_per_bin_singles(:)
  integer,allocatable :: n_per_bin_doubles(:)
 !integer,allocatable :: n_per_bin(:)
  integer,allocatable :: hash_ind_singles(:),hash_info_singles(:)
  integer,allocatable :: hash_ind_doubles(:),hash_info_doubles(:)
 !integer,allocatable :: hash_ind(:),hash_info(:)
  integer,allocatable :: start_index_singles(:),end_index_singles(:)
  integer,allocatable :: start_index_doubles(:),end_index_doubles(:)
 !integer,allocatable :: start_index(:),end_index(:)

!!integer :: n_bins = 1000000007 ! prime number of bins in hash tables. can change, but this seems reasonable
!!integer :: n_bins = 100000007 ! prime number of bins in hash tables. can change, but this seems reasonable
! integer :: n_bins = 1000003 ! prime number of bins in hash tables. can change, but this seems reasonable
  ! Wolfram alpha can be used to find primes near a number (search "primes near x")
  ! In python, import sympy and then nextprime(x)
  logical :: use_mpi
  integer :: max_sum_terms
  logical :: wf_file_exist
  real(rk), pointer, dimension(:),contiguous :: diag_elems_global
  integer :: ndets_global_old
  real(rk),parameter :: inv_sqrt2=1._rk/sqrt(2._rk)

contains

  subroutine perform_hci(eps_pt, eps_pt_big, target_error, n_states, hf_up, hf_dn)

    use chemistry, only : orbital_energies,time_sym,z
    use types, only : i8b, int_vec
    use semistoch, only : hamiltonian, estimate_n_connections
    use tools, only : merge_sort2_up_dn
    use heg, only : energy_hf, energy_madelung
    use common_psi_t, only : trial_wf_iters
    use common_run, only : n_connected_dets_hf
    use common_selected_ci, only : eps_var_sched, dump_wf_var, sparse_ham, n_mc, get_natorbs, get_greens_function, n_w, w_min, w_max, use_pt, eps_pt_big_energy, n_var_orbs, var_orbs, n_var_e_up, n_var_e_dn, n_energy_batch, use_hash_generation
    use mpi_routines, only : master_core_node, master_core, init_hash_owners, get_det_owner, whoami, &
         &det_map, mpi_barr, det_map_l, mpi_barr_in_node, shmem_allocate, shmem_deallocate, sharedComm, mpierr, &
         &nnodes, ncores, ncores_on_node, mpi_allred, mpi_bsend, mpi_distribute_remote_det_map, mpi_bsend_between_nodes, &
         get_det_owner, mpi_stop, mpi_gath
    use fhash_module__ik_int_list
#ifdef MPI
    INCLUDE 'mpif.h'
#endif

    real(rk),intent(in)    :: eps_pt, target_error
    real(rk)               :: eps_var
    real(rk),intent(inout) :: eps_pt_big

    integer :: iter, iter_done=0, max_iters=50
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),optional, intent(in) :: hf_up, hf_dn
    type(ik_vec),allocatable :: old_dets_up(:), old_dets_dn(:), new_dets_up(:), new_dets_dn(:),tmp_dets_up(:),tmp_dets_dn(:)
    type(ik_vec),allocatable :: temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
    type(ik_vec) :: core_up,core_dn,virt_up,virt_dn

#else
    integer(ik),optional, intent(in) :: hf_up, hf_dn
    integer(ik),allocatable :: old_dets_up(:), old_dets_dn(:), new_dets_up(:), new_dets_dn(:),tmp_dets_up(:),tmp_dets_dn(:)
    integer(ik),allocatable :: temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
    integer(ik) :: core_up,core_dn,virt_up,virt_dn
#endif
    integer,intent(in) :: n_states ! Number of low lying states to compute (1 for just ground state, 2 for ground and 1st excited, etc)
    real(rk),allocatable :: old_wts(:,:),starting_wts(:,:),tmp_wts(:,:)
    integer :: ndets_old,ndets_new,ndets_tmp
    real(rk),allocatable :: energy(:),old_energy(:)
    integer :: i,j,k
    integer,allocatable :: iorder(:),temp_i_2(:)
!   integer :: tdets=9999999 ! max number of dets.  Warning: this should be redone in a more robust way.
    real(rk),allocatable :: diag_elems(:),ave_connections_remote(:)
    integer :: tdets=9999 ! Initial number of dets.
    integer(i8b) :: n_conn_estimate
    integer :: n_batches, ierr=0
    real(rk) :: rdm(norb,norb)
    real(rk) :: delta_e_2pt
    real(rk) :: pt_big, pt_diff, pt_diff_std_dev, eps_pt_big_active, eps_pt_big_sav, avail, avail2, avail3
    real(rk) :: active_pt_big, active_pt_diff, active_pt_diff_std_dev
    integer :: core_id, ndets_connected,  ndets_global_new, ndet_remote, ndet_min, ndet_max, ndet_tot, n_mc_sav
    real(rk) :: average_connections, nonzero_elems, ave_connections_min, ave_connections_max, nonzero_elems_min, nonzero_elems_max, nonzero_elems_tot
    character*32 fmt, filename
    real(rk), allocatable :: min_H_already_done(:) ! Min of H elements already done (i.e., all connected matrix elements exceeding this in magnitude have already been generated)
    real(rk),allocatable :: coeffs(:)
    real(rk),allocatable :: w(:),G0_np1(:,:,:),G0_nm1(:,:,:)
    real(rk),allocatable :: rdm_state_avg(:,:)
    logical :: use_var_active_space=.false.
    integer :: n_core_up,n_core_dn, n_virt_up,n_virt_dn
    logical :: active_only
    real(rk) :: delta_e_2pt_active,delta_e_2pt_virt
    logical :: is_con

    !MJO - det_map is used for the parallel run
    type(det_map) remote_det_map
    type(det_map_l) local_det_map

    write (6,'(/,''Performing HCI with eps_var_sched='',30es8.2e1,/,''and  eps_pt='',es12.4)') eps_var_sched, eps_pt ; call flush(6)
    if (present(hf_up))  write (6,*) "HCI starting det:",hf_up,hf_dn; call flush(6)
    write (6,'(/,''Computing'',i3,'' lowest-energy states'')') n_states ; call flush(6)
    allocate(energy(n_states))
    allocate(old_energy(n_states))

! Normally, when ncores=1, mpi is not used.
! use_mpi is used as a debugging tool to force mpi usage even when ncores=1.
    if (ncores>1) then
      use_mpi = .true.
    else
      use_mpi = .false.
!     use_mpi = .true. ! Use this to force program to use mpi even on 1 core (just for testing)
    endif

    ! If variational active space, assign active electrons and orbitals
    core_up = 0_ik
    core_dn = 0_ik
    virt_up = 0_ik
    virt_dn = 0_ik
    if (n_var_orbs>0) then
      write (6,*) "Using variational active space"; call flush(6)
      use_var_active_space=.true.
      if (var_orbs(1)==0) then ! automatically choose which orbitals are the active ones
        ! Set core orbital masks
        n_core_up = nup - n_var_e_up
        n_core_dn = ndn - n_var_e_dn
#ifdef NUM_ORBITALS_GT_127
        core_up = maskr_vec(n_core_up)
        core_dn = maskr_vec(n_core_dn)
#else
        core_up = maskr(n_core_up)
        core_dn = maskr(n_core_dn)
#endif
        ! Set virtual orbital masks
        n_virt_up = norb - n_core_up - n_var_orbs
        n_virt_dn = norb - n_core_dn - n_var_orbs
#ifdef NUM_ORBITALS_GT_127
        virt_up = maskr_vec(norb) - maskr_vec(n_core_up+n_var_orbs)
        virt_dn = maskr_vec(norb) - maskr_vec(n_core_dn+n_var_orbs)
#else
        virt_up = maskr(norb,ik) - maskr(n_core_up+n_var_orbs,ik)
        virt_dn = maskr(norb,ik) - maskr(n_core_dn+n_var_orbs,ik)
#endif
        write (6,*) "n_core_up,n_core_dn=",n_core_up,n_core_dn
        write (6,*) "core_up,core_dn=",core_up,core_dn
        write (6,*) "n_virt_up,n_virt_dn=",n_virt_up,n_virt_dn
        write (6,*) "virt_up,virt_dn=",virt_up,virt_dn
        call mpi_bsend(core_up)
        call mpi_bsend(core_dn)
        call mpi_bsend(virt_up)
        call mpi_bsend(virt_dn)
      else ! use the active orbitals that were input
        call mpi_stop("Hand-picked variational active orbs not implemented yet")
      endif
    else
      write (6,'(''Not using variational active space'')') ; call flush(6)
    endif


    ! If variational wavefn file for this eps_var exists, then read it
    eps_var=minval(eps_var_sched)
    write(fmt,'(es7.2e1)') eps_var
    write(filename,'(''wf_eps_var='',a)') trim(fmt)
    inquire(file=filename, exist=wf_file_exist)
    if(wf_file_exist) then !MJO TEMP PARALLEL
       write(6,'(/,''Reading variational wavefn from '',a)') trim(filename)
       call my_second(1, 'read variational wavefn')
       if (.not.use_mpi) then
          open(1, file=filename,form='unformatted',status='old')
          read(1) ndets_old
          ndets_global_old = ndets_old
          allocate(old_dets_up(ndets_old))
          allocate(old_dets_dn(ndets_old))
          allocate(old_wts(ndets_old,n_states))
          read(1) old_dets_up(1:ndets_old)
          read(1) old_dets_dn(1:ndets_old)
          read(1) old_wts(1:ndets_old,1:n_states)
          read(1) energy(1:n_states)
         !close(1)
          ndets_global_old=ndets_old
       else
         !MJO Only master core reads in
         if(master_core) then
           open(1, file=filename,form='unformatted',status='old')
           read(1) ndets_global_old
         endif
         call mpi_bsend(ndets_global_old)
         call shmem_allocate(remote_det_map%global_dets_up,int(ndets_global_old,i8b))
         call shmem_allocate(remote_det_map%global_dets_dn,int(ndets_global_old,i8b))
         call shmem_allocate(remote_det_map%remote_wts,int(ndets_global_old,i8b),int(n_states,i8b))
         if (master_core) then
           read(1) remote_det_map%global_dets_up(1:ndets_global_old)
           read(1) remote_det_map%global_dets_dn(1:ndets_global_old)
           read(1) remote_det_map%remote_wts(1:ndets_global_old,1:n_states)
           read(1) energy(1:n_states)
          !close(1)
         endif

         call mpi_bsend(energy)

         local_det_map%ndets_global_old = ndets_global_old
         !MJO Send everyone the old_wts, state by state
         !    FIXME Do all states at once! Write mpi_bsend_real_dparray2
         if (master_core_node) then
           do i=1,n_states
             call mpi_bsend_between_nodes(remote_det_map%remote_wts(:,i))
           enddo
         endif
         call mpi_barr_in_node()
         !MJO Send everyone the remote_det_map and assign the local_det_map pieces
         call mpi_distribute_remote_det_map(remote_det_map,local_det_map)
         ndets_old = local_det_map%ndets
         !Now, allocate or local dets list
         allocate(old_dets_up(ndets_old))
         allocate(old_dets_dn(ndets_old))
         allocate(old_wts(ndets_old,n_states))
         !Slice our our local determinants
         old_dets_up(:) = remote_det_map%global_dets_up(local_det_map%local_indices)
         old_dets_dn(:) = remote_det_map%global_dets_dn(local_det_map%local_indices)
         do i=1,n_states
           old_wts(:,i) = remote_det_map%remote_wts(local_det_map%local_indices,i)
         enddo
       endif
       call my_second(2, 'read variational wavefn')
       write (6,'(/,''Read final iteration''i9,'' ='',es9.2,'' dets_global, n_states, variational energy='',i3,10f16.6)') ndets_global_old, real(ndets_global_old), n_states, energy(1:n_states); call flush(6)

     else ! variational wavefn file for this eps_var does not exist

      write(6,'(/,''Did not find variational wavefn file '',a,'', so generate wavefn'')') trim(filename)

      allocate(old_dets_up(tdets))
      allocate(old_dets_dn(tdets))
      allocate(old_wts(tdets,n_states))

      ! Start with HF det
      ndets_old = 1
! For Hitesh's 3-band Hubbard model, the up and dn dets are not the same.    Take linear comb.
      if (present(hf_up)) then
        if (hf_up==hf_dn) then
          old_dets_up(1) = hf_up
          old_dets_dn(1) = hf_dn
          local_det_map%hf_up = hf_up
          local_det_map%hf_dn = hf_dn
        else
          !MJO For now, we only want to start with 1 determinant in chemistry,
          !    even if hf_up!=hf_dn. Note that this breaks the 3-band hubbard model
          ndets_old = 1
          old_dets_up(1) = hf_up
          old_dets_dn(1) = hf_dn

          ! ndets_old = 2
          ! old_dets_up(1) = hf_up
          ! old_dets_dn(1) = hf_dn
          ! old_dets_up(2) = hf_dn
          ! old_dets_dn(2) = hf_up
        endif
      else
#ifdef NUM_ORBITALS_GT_127
        old_dets_up(1) = maskr_vec(nup)
        old_dets_dn(1) = maskr_vec(ndn)
        local_det_map%hf_up = maskr_vec(nup)
        local_det_map%hf_dn = maskr_vec(ndn)
#else
        old_dets_up(1) = maskr(nup,ik)
        old_dets_dn(1) = maskr(ndn,ik)
        local_det_map%hf_up = maskr(nup,ik)
        local_det_map%hf_dn = maskr(ndn,ik)
#endif
      endif

      ndets_global_old = ndets_old
      core_id = get_det_owner(old_dets_up(1),old_dets_dn(1))
      if (core_id==whoami) then
        do i=1,n_states
          old_wts(i,i) = 1._rk
          energy(i) = 0._rk
          if (i<=ndets_old)  call hamiltonian(old_dets_up(i),old_dets_dn(i),old_dets_up(i),old_dets_dn(i),energy(i),is_con)
        enddo
        allocate(local_det_map%local_indices(1))
        local_det_map%local_indices(1) = 1
      else
        energy = 0
        ndets_old = 0
       !old_dets_up(1) = 0
       !old_dets_dn(1) = 0
      endif
      call mpi_allred(energy)

      write (6,'(''Iteration   0'','' eps1='',es7.1e1,'' ndets='',i9,'' ='',es9.2,'' energy='',10f16.6)') eps_var_sched(1), ndets_global_old, real(ndets_global_old), energy; call flush(6)
      old_energy = energy

      allocate(min_H_already_done(ndets_old))
      min_H_already_done(:) = 9.d99

      !MJO set the initial remote_det_map information
      !    we set everything to 0, because we want to treat HF as
      !    'generated' on the first iteration
      local_det_map%ndets_global = 0
      local_det_map%n_it         = 0
      local_det_map%hf_up = old_dets_up(1)
      local_det_map%hf_dn = old_dets_dn(1)
      local_det_map%ndets = 0
      local_det_map%average_connections = n_connected_dets_hf**0.5
      average_connections = n_connected_dets_hf**0.5
      allocate(local_det_map%n_det_remote(max_iters,0:ncores-1))
      local_det_map%n_det_remote = 0

      !MJO Allocate some initial small space for the global_dets_up routine
      call shmem_allocate(remote_det_map%global_dets_up,int(1,i8b))
      call shmem_allocate(remote_det_map%global_dets_dn,int(1,i8b))
      call shmem_allocate(remote_det_map%remote_wts,1_i8b,int(n_states,i8b))

      call my_second(1,'variational part')

!     n_eps_var_sched_steps=1
!     do i=2,size(eps_var_sched)
!       if(eps_var_sched(i).ne.eps_var_sched(i-1) n_eps_var_sched_steps=n_eps_var_sched_steps+1
!     enddo
!     write(6,'(/,''n_eps_var_sched_steps='',i2)') n_eps_var_sched_steps

!     iter_done=0
!     do i_eps_var_sched_steps=1,n_eps_var_sched_steps

!     do iter=iter_done+1,max_iters
      do iter=1,max_iters

         if(iter.le.size(eps_var_sched)) eps_var=eps_var_sched(iter)
         write(6,'(/,''HCI iter='',i3,'' eps_var='',es12.4)') iter, eps_var

         if (use_mpi) then !MJO TEMP PARALLEL
            ndets_global_old = local_det_map%ndets_global
         else
            ndets_global_old = ndets_old
         endif
        allocate(coeffs(max(ndets_old,1)))
        if (iter>1) then
          ! do i=1,ndets_global_old
          !   coeffs(i) = maxval(abs(old_wts(i,1:n_states)))
          ! enddo
          do i=1,ndets_old
           !coeffs(i) = sqrt(sum(old_wts(i,1:n_states)**2))
            coeffs(i) = maxval(abs(old_wts(i,1:n_states)))
          enddo
        else
          do i=1,ndets_old
            coeffs(i) = old_wts(i,1)
          enddo
        endif

        if (use_var_active_space) then
          if (use_mpi) then !MJO TEMP PARALLEL
            call get_next_det_list(n_states, eps_var, ndets_old, old_dets_up, old_dets_dn, coeffs, ndets_new, new_dets_up, new_dets_dn, min_H_already_done, remote_det_map, local_det_map,core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn)
            ndets_global_new = local_det_map%ndets_global
          else
            call get_next_det_list(n_states, eps_var, ndets_old, old_dets_up, old_dets_dn, coeffs, ndets_new, new_dets_up, new_dets_dn, min_H_already_done,core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn)
            ndets_global_new = ndets_new
          endif
        else
          if (use_mpi) then !MJO TEMP PARALLEL
            call get_next_det_list(n_states, eps_var, ndets_old, old_dets_up, old_dets_dn, coeffs, ndets_new, new_dets_up, new_dets_dn, min_H_already_done, remote_det_map, local_det_map)
            ndets_global_new = local_det_map%ndets_global
          else
            call get_next_det_list(n_states, eps_var, ndets_old, old_dets_up, old_dets_dn, coeffs, ndets_new, new_dets_up, new_dets_dn, min_H_already_done)
            ndets_global_new = ndets_new
          endif
        endif

        ! new dets should now be sorted by label
!        write(6,*) ndets_new,local_det_map%ndets
        ! do i=1,local_det_map%ndets
        !   write(6,*) local_det_map%local_indices(i)
        ! enddo
        ! tmp_dets_up(1:local_det_map%ndets) = remote_det_map%global_dets_up(local_det_map%local_indices)
        ! tmp_dets_dn(1:local_det_map%ndets) = remote_det_map%global_dets_dn(local_det_map%local_indices)
        ! do i=1,local_det_map%ndets
        !   write(6,*) new_dets_up(i),tmp_dets_up(i),new_dets_dn(i),tmp_dets_dn(i)
        ! enddo

        if(ndets_global_new==ndets_global_old) then
          write (6,'(/,''Cycling hci iteration because ndets_global_new==ndets_global_old'')')
          deallocate(coeffs)
          cycle
        endif

        ! 1st criterion for exiting HCI iterations: the number of new dets is less than 0.00001 of the number of old ones.
        if (ndets_global_new<=int(1.00001*ndets_global_old) .and. eps_var==eps_var_sched(size(eps_var_sched))) then
           write (6,'(/,''Exiting variational part because of number of new dets criterion'')')
           ndets_global_new = ndets_global_old
           write(6,*) 'ndets_new: ',ndets_global_new,'ndets_old',ndets_global_old
           if (use_mpi) then !MJO TEMP PARALLEL
              local_det_map%ndets_global = ndets_global_new
           endif
           exit
        endif

        ! Previously the determinants were getting reordered at each HCI iteration, so we had to find the existing dets and assign old_wts as the starting wts.
        ! Now the new dets. are all at the end, so it is simpler.
        !allocate(starting_wts(ndets_global_new,n_states))
        ! starting_wts(ndets_global_old+1:ndets_global_new,1:n_states) = 0._rk
        ! starting_wts(1:ndets_global_old,1:n_states) = old_wts(1:ndets_global_old,1:n_states)

        ! if (size(old_wts, 1) < ndets_global_new) then
        !   deallocate(old_wts)
        !   allocate(old_wts(ndets_global_new, n_states))
        ! endif

        ! if (iter==1) then
        !   starting_wts(1:ndets_global_new,1:n_states) = 0._rk
        !   old_wts(1:ndets_global_old,1:n_states) = 0
        !   core_id = get_det_owner(old_dets_up(1),old_dets_dn(1))
        !   !MJO Make sure HF, which might be on a different core, is the initial state
        !   old_wts(1+core_id,1) = 1._rk
        !   starting_wts(1+core_id,1) = 1._rk
        !   do i=2,n_states
        !     starting_wts(i+core_id,i) = 1._rk
        !   enddo
        ! endif

        allocate(starting_wts(ndets_new,n_states))
        starting_wts(ndets_old+1:ndets_new,1:n_states) = 0._rk
        starting_wts(1:ndets_old,1:n_states) = old_wts(1:ndets_old,1:n_states)

        if (size(old_wts, 1) < ndets_new) then
          deallocate(old_wts)
          allocate(old_wts(ndets_new, n_states))
        endif

        if (iter==1) then
          starting_wts(1:ndets_new,1:n_states) = 0._rk
          old_wts(1:ndets_old,1:n_states) = 0
          core_id = get_det_owner(old_dets_up(1),old_dets_dn(1))
          !MJO Make sure HF, which might be on a different core, is the initial state
          if (core_id.eq.whoami) then
            old_wts(1,1) = 1._rk
            starting_wts(1,1) = 1._rk
            do i=2,n_states
              starting_wts(i,i) = 1._rk
            enddo
          endif
        endif

        if (use_mpi) then !MJO TEMP PARALLEL
          if (use_hash_generation) then
            call iterative_diagonalize(ndets_new, n_states, new_dets_up, new_dets_dn, final_wts=old_wts, eigenvalues=energy, starting_wts=starting_wts, remote_det_map=remote_det_map,local_det_map=local_det_map,use_hash_generation=use_hash_generation)
          else
            call iterative_diagonalize(ndets_new, n_states, new_dets_up, new_dets_dn, final_wts=old_wts, eigenvalues=energy, starting_wts=starting_wts, remote_det_map=remote_det_map,local_det_map=local_det_map)
          endif
           average_connections = local_det_map%average_connections
        else
           call iterative_diagonalize(ndets_new, n_states, new_dets_up, new_dets_dn, final_wts=old_wts, eigenvalues=energy, starting_wts=starting_wts,average_connections=average_connections)
        endif
        deallocate(starting_wts,coeffs)
        write (6,'(''Iteration'',i4,'' eps1='',es7.1e1,'' ndets='',i9,'' ='',es9.2,'' energy='',10f16.6)') iter, eps_var, ndets_global_new, real(ndets_global_new), energy; call flush(6)
        iter_done=iter

        ndets_old = ndets_new

        if (size(old_dets_up) < ndets_new) then
          deallocate(old_dets_up)
          deallocate(old_dets_dn)
          allocate(old_dets_up(ndets_new))
          allocate(old_dets_dn(ndets_new))
        endif
        old_dets_up(1:ndets_old) = new_dets_up(1:ndets_old)
        old_dets_dn(1:ndets_old) = new_dets_dn(1:ndets_old)

        ! 2nd criterion for exiting HCI iterations
        if (maxval(abs(energy(:)-old_energy(:)))<1.e-5_rk .and. eps_var==eps_var_sched(size(eps_var_sched))) then
       !if (maxval(abs(energy(:)-old_energy(:)))<1.e-6_rk .and. eps_var==eps_var_sched(size(eps_var_sched))) then
           if (use_mpi) then
              ndets_global_old = local_det_map%ndets_global
              local_det_map%ndets_global_old = ndets_global_old
           else
              ndets_global_old = ndets_old
           endif
           write (6,'(/,''Exiting variational part because of energy criterion'')')
           old_energy = energy
           exit
        endif

        old_energy = energy

      enddo ! HCI iter

      call my_second(2,'variational part')
      if (iter>1.and.ndets_global_old>1) then
        write (6,*) "Deallocating Hamiltonian"
        deallocate(sparse_ham%indices,sparse_ham%nonzero_elements,sparse_ham%values)
        sparse_ham%ndet = 0
      endif

      !MJO Calculate the approximate load imbalance from determinants
      if (use_mpi) then !MJO TEMP PARALLEL
         write(6,*) 'Determinant distribution: '
         allocate(ave_connections_remote(0:ncores-1))
         ave_connections_remote         = 0
         ave_connections_remote(whoami) = local_det_map%average_connections
         call mpi_allred(ave_connections_remote)
         ndet_min=2147483647 ; ave_connections_min=9.d99; nonzero_elems_min=9.d99
         ndet_max=0 ; ave_connections_max=0; nonzero_elems_max=0
         ndet_tot=0 ; nonzero_elems_tot=0
         do i=0,ncores-1
            ndet_remote = sum(local_det_map%n_det_remote(:,i))
            nonzero_elems = ndet_remote*ave_connections_remote(i)
            ndet_tot=ndet_tot+ndet_remote
            nonzero_elems_tot=nonzero_elems_tot+nonzero_elems
            write(6,'(''Core:'',i6,'' ndet:'',i9,'' ave_connec:'',es11.4,'' nonzero_elems:'',es11.4)') i, ndet_remote, ave_connections_remote(i), nonzero_elems
            ndet_min=min(ndet_min,ndet_remote)
            ave_connections_min=min(ave_connections_min,ave_connections_remote(i))
            nonzero_elems_min=min(nonzero_elems_min,nonzero_elems)
            ndet_max=max(ndet_max,ndet_remote)
            ave_connections_max=max(ave_connections_max,ave_connections_remote(i))
            nonzero_elems_max=max(nonzero_elems_max,nonzero_elems)
         enddo
         deallocate(ave_connections_remote)
         write(6,'(''ndet_min,            ndet_max,            max/min, max/ave='',t60,2i10,9f9.2)') ndet_min, ndet_max, real(ndet_max,rk)/ndet_min, (real(ndet_max)*ncores)/ndet_tot
         write(6,'(''ave_connections_min, ave_connections_max, max/min='',t60,2es10.3,f9.2)') ave_connections_min, ave_connections_max, ave_connections_max/ave_connections_min
         write(6,'(''nonzero_elems_min,   nonzero_elems_max,   max/min, max/ave='',t60,2es10.3,9f9.2)') nonzero_elems_min, nonzero_elems_max, nonzero_elems_max/nonzero_elems_min, nonzero_elems_max*ncores/nonzero_elems_tot
      endif

      write (6,'(/,''Final Iteration'',i4,i9,'' ='',es9.2,'' dets, energy='',10f16.6)') iter_done, ndets_global_new, real(ndets_global_new), energy; call flush(6)
      ! Now, old_dets_up/dn, old_wts, old_energy are the variational wavefunction

     !write (6,'(''Final variational wavefunctions (20 sorted by label):'')')
     !if (use_mpi) then !MJO TEMP PARALLEL
     !   !MJO This print is not consisten between different numbers of cores because
     !   !    The global det list is not sorted
     !   do i=1,min(20,ndets_global_old)
     !      write (6,*) i,remote_det_map%global_dets_up(i),remote_det_map%global_dets_dn(i),remote_det_map%remote_wts(i,:)
     !   enddo
     !else
     !   do i=1,min(20,ndets_global_old)
     !      write (6,*) i,old_dets_up(i),old_dets_dn(i),old_wts(i,:)
     !   enddo
     !endif

! TODO: Add in write of n_states, iter_done, time_sym, z, eps_var_sched,

      ! Sort variational dets so that they can be binary searched

      write(6,'(/,''Sorting global dets list so that they can be binary searched'')') ; call flush(6)

      if (.not.use_mpi)  ndets_global_old=ndets_old
      if (.not.allocated(iorder))  allocate(iorder(ndets_global_old))
      if (.not.allocated(temp_i16_up))  allocate(temp_i16_up((ndets_global_old+1)/2))
      if (.not.allocated(temp_i16_dn))  allocate(temp_i16_dn((ndets_global_old+1)/2))
      if (.not.allocated(temp_i_2))  allocate(temp_i_2((ndets_global_old+1)/2))
      do j=1,ndets_global_old
        iorder(j)=j
      enddo

      if (use_mpi) then
        if (master_core_node) then
          call merge_sort2_up_dn(remote_det_map%global_dets_up(1:ndets_global_old),remote_det_map%global_dets_dn(1:ndets_global_old), iorder, ndets_global_old, temp_i16_up, temp_i16_dn, temp_i_2)
          remote_det_map%remote_wts(1:ndets_global_old,1:n_states) = remote_det_map%remote_wts(iorder(1:ndets_global_old),1:n_states)
        endif
        call mpi_barr_in_node()
      endif
      do j=1,ndets_old
        iorder(j)=j
      enddo
      call merge_sort2_up_dn(old_dets_up(1:ndets_old),old_dets_dn(1:ndets_old), iorder(1:ndets_old), ndets_old, temp_i16_up, temp_i16_dn, temp_i_2)
      old_wts(1:ndets_old,1:n_states) = old_wts(iorder(1:ndets_old),1:n_states)

      deallocate(iorder, temp_i16_up, temp_i16_dn, temp_i_2)

      ! Dump variational wavefn to file if it does not exist
      if(dump_wf_var .and. .not. wf_file_exist) then
        write(6,'(/,''Writing variational wavefn to '',a)') trim(filename)
        call my_second(1, 'dump variational wavefn')
        if(.not.use_mpi) then
           open(1, file=filename,form='unformatted',status='new')
           write(1) ndets_old
           write(1) old_dets_up(1:ndets_old)
           write(1) old_dets_dn(1:ndets_old)
           write(1) old_wts(1:ndets_old,1:n_states)
           write(1) energy(1:n_states)
           call flush(1)
          !close(1)
        elseif(master_core) then
           open(1, file=filename,form='unformatted',status='new')
           write(1) ndets_global_old
           write(1) remote_det_map%global_dets_up(1:ndets_global_old)
           write(1) remote_det_map%global_dets_dn(1:ndets_global_old)
           write(1) remote_det_map%remote_wts(1:ndets_global_old,1:n_states)
           write(1) energy(1:n_states)
           call flush(1)
          !close(1)
        endif ! ncores
        call my_second(2, 'dump variational wavefn')
      endif ! dump_wf_var

    endif ! read or compute variational wavefn

    call mem_avail(avail,avail2,avail3) ; write(6,'(/,''Before det. PT, avail, avail2, avail3='',t114,3f11.2)') avail, avail2, avail3
    if (n_max_connections <= 0) then
      n_max_connections=round_r(1e6*avail2/200,2) ! The 1e6 is for converting MB to B, and the 200 is 8 times an ill-defined factor
    endif
    !MJO Take a global average of all n_max_connections
    call mpi_allred(n_max_connections)
    n_max_connections = n_max_connections/ncores
    write(6,'(/,''Global average n_max_connections per node:'',i13,'' ='',es9.2)') int(n_max_connections,i8b), n_max_connections

    if (n_energy_batch.gt.0) then
      if(.not.use_mpi) then
        write(6,*) "ERROR: use_mpi must be true to use this routine"
      endif
      call energies_for_extrapolation(n_states,remote_det_map,local_det_map,n_energy_batch,eps_pt, eps_pt_big,target_error,n_mc,eps_var)
      stop
    endif

    ! Go back to determinant basis because lin combos of dets can potentially
    ! cause intruder state problems
    if (time_sym) then
      write (6,'(/,''Converting to determinant basis'')') ; call flush(6)
      if (use_mpi) then
        call convert_time_symmetrized_to_dets(ndets_old,old_dets_up,old_dets_dn,old_wts,n_states,remote_det_map,local_det_map)
        ndets_global_old = local_det_map%ndets_global_old
      else
        call convert_time_symmetrized_to_dets(ndets_old,old_dets_up,old_dets_dn,old_wts,n_states)
        ndets_global_old = ndets_old
      endif
      write (6,'(''From now on, no more time-reversal symmetry'',/)') ; call flush(6)
      time_sym = .false.
    endif
    ndets_global_new = ndets_global_old

    if (use_mpi) then
      !MJO Quick hack - compute all global diagonal elements
      call shmem_allocate(diag_elems_global,int(ndets_global_old,i8b))
      if (master_core_node) then
        do k=1,ndets_global_old
          call hamiltonian(remote_det_map%global_dets_up(k),remote_det_map%global_dets_dn(k),remote_det_map%global_dets_up(k),remote_det_map%global_dets_dn(k),diag_elems_global(k),is_con)
        enddo
      endif
      call mpi_barr_in_node()
    endif
    ! Quick hack: compute all local diag_elems now
    allocate(diag_elems(ndets_old))
    do i=1,ndets_old
      call hamiltonian(old_dets_up(i),old_dets_dn(i),old_dets_up(i),old_dets_dn(i),diag_elems(i),is_con)
    enddo

! Check available memory


! Calculate Natural Orbitals

    if (get_natorbs) then
      if (hamiltonian_type.ne.'chem')  stop "Natural orbitals only implemented for chem so far"
      if (use_pt) then
        write (6,'(/,''Getting natural orbitals from the lowest-order perturbative correction to the variational 1-RDM'')') ; call flush(6)
        write (6,*) "eps_pt=",eps_pt; call flush(6)
        if(eps_pt_big.le.0_rk) then
          eps_pt_big = eps_pt
          !MJO - Each core does a local estimate and then we reduce all of the values
          !      This gives us a much a larger sample and, thus, a more accurate estimate
          n_conn_estimate = estimate_n_connections(ndets_old,old_dets_up,old_dets_dn,old_wts(:,1),eps_pt)
          call mpi_allred(n_conn_estimate)
          write (6,*) "Estimated number of connections, maximum number of connections=",n_conn_estimate,n_max_connections; call flush(6)
          do while (n_conn_estimate>=nint(n_max_connections*nnodes,i8b))
            eps_pt_big = round_r(1.1*eps_pt_big*(real(n_conn_estimate)/(n_max_connections*nnodes))**.75,2)
            if(eps_pt_big >= eps_var) exit
            !MJO - Each core does a local estimate and then we reduce all of the values
            !      This gives us a much a larger sample and, thus, a more accurate estimate
            n_conn_estimate = estimate_n_connections(ndets_old, old_dets_up, old_dets_dn, old_wts(:,1), eps_pt_big)
            call mpi_allred(n_conn_estimate)
            write (6,'(''If eps_pt_big='',es11.4,'' estimated number of connections to variational wavefn='',i13,'' ='',es9.2)') eps_pt_big, n_conn_estimate, real(n_conn_estimate)
          enddo
        endif
        !Lowest order correction to 1RDM: <0|rho|0> + <1|rho|0> + <0|rho|1>   (leaves out <1|rho|1>)
        !(where |0> is variational wavefn and |1> is 1st order PT correction
        if (use_mpi) then
          stop "1rdm with pt not currently implemented in parallel! Try without pt."
        else
          call get_1rdm_with_pt(ndets_old, old_dets_up, old_dets_dn, old_wts(:,1), diag_elems, energy(1), eps_pt_big, rdm)
        endif
      else
        write (6,*) "Getting natural orbitals from the variational 1-RDM"; call flush(6)
        if (n_states>1) then
          write (6,*) "Using state-averaged 1-RDM from the lowest",n_states,"states"; call flush(6)
          allocate(rdm_state_avg(norb,norb))
          rdm_state_avg(:,:) = 0._rk
          do i=1,n_states
            write (6,*) "Computing 1-RDM of state",i; call flush(6)
            if (use_mpi) then
              call get_1rdm(ndets_old, remote_det_map%global_dets_up, remote_det_map%global_dets_dn, remote_det_map%remote_wts(:,i), rdm, remote_det_map, local_det_map, i)
            else
              call get_1rdm(ndets_old, old_dets_up, old_dets_dn, old_wts(:,i), rdm)
            endif
            rdm_state_avg = rdm_state_avg + rdm
          enddo ! i
          rdm = rdm_state_avg/n_states
          deallocate(rdm_state_avg)
        else
          if (use_mpi) then
            call get_1rdm(ndets_old, old_dets_up, old_dets_dn, old_wts(:,1), rdm, remote_det_map, local_det_map, 1)
          else
            call get_1rdm(ndets_old, old_dets_up, old_dets_dn, old_wts(:,1), rdm)
          endif
        endif
      endif ! use_pt
      write (6,'(''Done generating rdm'')') ; call flush(6)
      if (master_core) then
        !MJO Only master_core generates the natorb integrals,
        !    for now
        call generate_natorb_integrals(rdm)
      endif
      call mpi_barr()
      write (6,'(/,''Done generating natorbs; integrals in FCIDUMP.natorb'')') ; call flush(6)
      stop "Done generating natorbs; integrals in FCIDUMP.natorb"
    endif

    if (get_greens_function) then
      allocate(w(n_w))
      do i=1,n_w
        w(i) = w_min + (w_max-w_min)*(i-1._rk)/(n_w-1._rk)
      enddo
      write (6,*) "About to compute Green's function for",n_w,"frequencies:"
      write (6,*) w
      write (6,*) "E_0 =",energy(1)
      allocate(G0_np1(n_w,norb,norb))
      allocate(G0_nm1(n_w,norb,norb))
      call get_zeroth_order_variational_greens_function(ndets_old,old_dets_up,old_dets_dn,old_wts(:,1),energy(1),n_w,w,G0_np1,G0_nm1)
    endif

! Calculate PT correction

    ! Estimate of number of connections done separately for each state
    eps_pt_big_sav=eps_pt_big
    n_mc_sav=n_mc

    ! Save the read-in deterministic eps/pt here

    close(1)

    do i=1,n_states

      eps_pt_big = eps_pt_big_sav
      n_mc = n_mc_sav

      if (n_var_orbs>0) then ! Active space
        eps_pt_big_active = eps_pt_big_sav

        pt_big = 0._rk
        pt_diff = 0._rk
        pt_diff_std_dev = 0._rk

        write (6,*) "Computing PT correction within active space. This should be used for extrapolation of the total energy since it will be zero in the CAS+PT limit"; call flush(6)

        if (use_mpi) then
          call do_pt(ndets_old,old_dets_up(1:ndets_old),old_dets_dn(1:ndets_old),old_wts(1:ndets_old,i),diag_elems(1:ndets_old),ndets_global_old,remote_det_map%global_dets_up(1:ndets_global_old), remote_det_map%global_dets_dn(1:ndets_global_old), remote_det_map%remote_wts(1:ndets_global_old,i), diag_elems_global(1:ndets_global_old), energy(i),eps_var,eps_pt,eps_pt_big_active,target_error,n_mc,active_pt_big,active_pt_diff,active_pt_diff_std_dev,ndets_connected,remote_det_map,local_det_map,core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn,active_only_in=.true.,active_tot_energy=energy(i))
        else
          call do_pt(ndets_old,old_dets_up(1:ndets_old),old_dets_dn(1:ndets_old),old_wts(1:ndets_old,i),diag_elems(1:ndets_old),var_energy=energy(i),eps_var=eps_var,eps_pt=eps_pt,eps_pt_big=eps_pt_big_active,target_error=target_error,n_mc=n_mc,pt_big=active_pt_big,pt_diff=active_pt_diff,pt_diff_std_dev=active_pt_diff_std_dev,ndets_connected=ndets_connected,core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn,active_only_in=.true.,active_tot_energy=energy(i))
        endif

        write (6,*) "Computing PT correction outside of active space"; call flush(6)

        n_mc = n_mc_sav

        if (use_mpi) then
          call do_pt(ndets_old,old_dets_up(1:ndets_old),old_dets_dn(1:ndets_old),old_wts(1:ndets_old,i),diag_elems(1:ndets_old),ndets_global_old,remote_det_map%global_dets_up(1:ndets_global_old), remote_det_map%global_dets_dn(1:ndets_global_old), remote_det_map%remote_wts(1:ndets_global_old,i), diag_elems_global(1:ndets_global_old), energy(i),eps_var,eps_pt,eps_pt_big,target_error,n_mc,pt_big,pt_diff,pt_diff_std_dev,ndets_connected,remote_det_map,local_det_map,core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn,active_only_in=.false.,active_tot_energy=energy(i)+pt_big+pt_diff)
        else
          call do_pt(ndets_old,old_dets_up(1:ndets_old),old_dets_dn(1:ndets_old),old_wts(1:ndets_old,i),diag_elems(1:ndets_old),var_energy=energy(i),eps_var=eps_var,eps_pt=eps_pt,eps_pt_big=eps_pt_big,target_error=target_error,n_mc=n_mc,pt_big=pt_big,pt_diff=pt_diff,pt_diff_std_dev=pt_diff_std_dev,ndets_connected=ndets_connected,core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn,active_only_in=.false.,active_tot_energy=energy(i)+pt_big+pt_diff)
        endif

        ! Combine PT contributions
        pt_big = pt_big + active_pt_big
        pt_diff = pt_diff + active_pt_diff
        pt_diff_std_dev = sqrt(pt_diff_std_dev**2 + active_pt_diff_std_dev**2)

        write (6,'(/,''State'',i4,'':'')') i
        write (6,'(''Variational energy('',i1,'')='',t44,f15.9)') i, energy(i)
        if (active_pt_diff_std_dev==0._rk) then
          write (6,'(''Actv-sp. 2nd-order PT energy lowering('',i1,'')='',t44,f15.9)') i, active_pt_big
        else
          write (6,'(''Actv-sp. 2nd-order PT energy lowering('',i1,'')='',t44,f15.9,'' +-'',f12.9,'' ('',2f13.9,'')'')') i, active_pt_big+active_pt_diff, active_pt_diff_std_dev, active_pt_big, active_pt_diff
        endif
        if (pt_diff_std_dev==0._rk) then
          write (6,'(''Full-sp. 2nd-order PT energy lowering('',i1,'')='',t44,f15.9)') i, pt_big
          write (6,'(''Total energy('',i1,'')='',t44,f15.9)') i, energy(i)+pt_big
          write (6,'(''eps_var, eps_pt, ndets, ndets_connected(total), Variational, PT_actv, PT_full, Total Energies('',i1,'')='',2es8.1,i9,i11,f16.9,2f15.9,f16.9)') i, eps_var, eps_pt, ndets_global_old, ndets_connected, energy(i), active_pt_big, pt_big, energy(i)+pt_big
        else
          write (6,'(''Full-sp. 2nd-order PT energy lowering('',i1,'')='',t44,f15.9,'' +-'',f12.9,'' ('',2f13.9,'')'')') i, pt_big+pt_diff, pt_diff_std_dev, pt_big, pt_diff
          write (6,'(''Total energy('',i1,'')='',t44,f15.9,'' +-'',f12.9)') i, energy(i)+pt_big+pt_diff, pt_diff_std_dev
!         write (6,'(''ndets, ndets_connected(total), Variational, PT, Total Energies('',i1,'')='',i9,i11,f16.9,f15.9,f16.9,'' +-'',f12.9)') i, ndets_global_old, ndets_connected, energy(i), pt_big+pt_diff, energy(i)+pt_big+pt_diff, pt_diff_std_dev
          write (6,'(''eps_var, eps_pt, ndets, ndets_connected(total), Variational, PT_actv, PT_full, Total Energies('',i1,'')='',2es8.1,i9,i11,f16.9,2f15.9,f16.9,'' +-'',f12.9)') i, eps_var, eps_pt, ndets_global_old, ndets_connected, energy(i), active_pt_big+active_pt_diff, pt_big+pt_diff, energy(i)+pt_big+pt_diff, pt_diff_std_dev
        endif

      else ! No separation between active and inactive contributions

        if (use_mpi) then
          call do_pt(ndets_old,old_dets_up(1:ndets_old),old_dets_dn(1:ndets_old),old_wts(1:ndets_old,i),diag_elems(1:ndets_old),ndets_global_old,remote_det_map%global_dets_up(1:ndets_global_old), remote_det_map%global_dets_dn(1:ndets_global_old), remote_det_map%remote_wts(1:ndets_global_old,i), diag_elems_global(1:ndets_global_old), energy(i),eps_var,eps_pt,eps_pt_big,target_error,n_mc,pt_big,pt_diff,pt_diff_std_dev,ndets_connected,remote_det_map,local_det_map)
        else
          call do_pt(ndets_old,old_dets_up(1:ndets_old),old_dets_dn(1:ndets_old),old_wts(1:ndets_old,i),diag_elems(1:ndets_old),var_energy=energy(i),eps_var=eps_var,eps_pt=eps_pt,eps_pt_big=eps_pt_big,target_error=target_error,n_mc=n_mc,pt_big=pt_big,pt_diff=pt_diff,pt_diff_std_dev=pt_diff_std_dev,ndets_connected=ndets_connected)
        endif

        write (6,'(/,''State'',i4,'':'')') i
        write (6,'(''Variational energy('',i1,'')='',t36,f15.9)') i, energy(i)
        if (pt_diff_std_dev==0._rk) then
          write (6,'(''2nd-order PT energy lowering('',i1,'')='',t36,f15.9)') i, pt_big
          write (6,'(''Total energy('',i1,'')='',t36,f15.9)') i, energy(i)+pt_big
          write (6,'(''eps_var, eps_pt, ndets, ndets_connected(total), Variational, PT_actv, PT_full, Total Energies('',i1,'')='',2es8.1,i9,i11,f16.9,2f15.9,f16.9)') i, eps_var, eps_pt, ndets_global_old, ndets_connected, energy(i), pt_big, pt_big, energy(i)+pt_big
        else
          write (6,'(''2nd-order PT energy lowering('',i1,'')='',t36,f15.9,'' +-'',f12.9,'' ('',2f13.9,'')'')') i, pt_big+pt_diff, pt_diff_std_dev, pt_big, pt_diff
          write (6,'(''Total energy('',i1,'')='',t36,f15.9,'' +-'',f12.9)') i, energy(i)+pt_big+pt_diff, pt_diff_std_dev
          write (6,'(''eps_var, eps_pt, ndets, ndets_connected(total), Variational, PT_actv, PT_full, Total Energies('',i1,'')='',2es8.1,i9,i11,f16.9,2f15.9,f16.9,'' +-'',f12.9)') i, eps_var, eps_pt, ndets_global_old, ndets_connected, energy(i), pt_big+pt_diff, pt_big+pt_diff, energy(i)+pt_big+pt_diff, pt_diff_std_dev
        endif

        if(hamiltonian_type.eq.'heg') then
          if (pt_diff_std_dev==0._rk) then
            write(6,'(''Total energy (includ. Madelung)='',t36,f15.9)') energy(i)+pt_big+energy_madelung
            write(6,'(''Correlation energy ='',t36,f15.9)') energy(i)+pt_big-energy_hf
            write(6,'(''Correlation energy per electron ='',t36,f15.9)') (energy(i)+pt_big-energy_hf)/nelec
          else
            write(6,'(''Total energy (includ. Madelung)='',t36,f15.9,'' +-'',f12.9)') energy(i)+pt_big+pt_diff+energy_madelung, pt_diff_std_dev
            write(6,'(''Correlation energy ='',t36,f15.9,'' +-'',f12.9)') energy(i)+pt_big+pt_diff-energy_hf, pt_diff_std_dev
            write(6,'(''Correlation energy per electron ='',t36,f15.9)') (energy(i)+pt_big+pt_diff-energy_hf)/nelec
          endif
        endif

        call flush(6)

      endif ! Different active spaces
 
    enddo !n_states

  end subroutine perform_hci
!-------------------------------------------------------------------------------------

  subroutine get_next_det_list(n_states, eps_var, ndets_old, old_dets_up, old_dets_dn, old_wts, ndets_new, new_dets_up, new_dets_dn, min_H_already_done, remote_det_map,local_det_map,core_up,core_dn,virt_up,virt_dn)
  ! Gets next determinant list for the next iteration of Heat-bath CI
  ! A Holmes, Feb 2016
  ! Modified: A Holmes, 17 May 2017. Added optional arguments core_up/dn,
  !           virt_up/dn, for active space variational calculation

    use semistoch, only : find_doubly_excited
    use more_tools, only : binary_search
    use generic_sort, only : sort
    use mpi_routines, only : mpi_allgatherv_new_dets, ncores, det_map, det_map_l, master_core_node, mpi_barr_in_node

    integer,intent(in) :: n_states
    real(rk),intent(in) :: eps_var
    integer,intent(in) :: ndets_old
    integer,intent(out) :: ndets_new
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: old_dets_up(:),old_dets_dn(:)
    type(ik_vec),optional,intent(in) :: core_up,core_dn,virt_up,virt_dn
    type(ik_vec), allocatable, intent(out) :: new_dets_up(:),new_dets_dn(:)
    type(ik_vec),allocatable :: tmp_dets_up(:),tmp_dets_dn(:)
    type(ik_vec),allocatable :: tmp_dets_up2(:),tmp_dets_dn2(:)
    type(ik_vec),allocatable :: tmp_old_up(:),tmp_old_dn(:)
#else
    integer(ik),intent(in) :: old_dets_up(:),old_dets_dn(:)
    integer(ik),optional,intent(in) :: core_up,core_dn,virt_up,virt_dn
    integer(ik), allocatable, intent(out) :: new_dets_up(:),new_dets_dn(:)
    integer(ik),allocatable :: tmp_dets_up(:),tmp_dets_dn(:)
    integer(ik),allocatable :: tmp_dets_up2(:),tmp_dets_dn2(:)
    integer(ik),allocatable :: tmp_old_up(:),tmp_old_dn(:)
#endif
    real(rk),intent(in) :: old_wts(:)
    real(rk),allocatable,intent(inout) :: min_H_already_done(:)
    integer :: n_ref,n_det_remote_new(0:ncores-1),n_det_tmp
    integer(i8b) :: n_det_tmp2,ndets_global_old,ndets_global_new
    integer :: i,j,k,offset,ndets_tmp

    !MJO Det_map is used in parallel runs
    type(det_map), intent(inout),optional :: remote_det_map
    type(det_map_l), intent(inout),optional :: local_det_map

    call my_second(1,"get_next_det_list")

    ! Find dets i for which |H_ij*c_j| > eps_var for at least one j
    if (eps_var>=0._rk) then
      n_ref = ndets_old
      !MJO below is a little hackish way of breaking up the first iteration and the others
      !    They need to be treated differently because the first iteration has HF on a random
      !    core, and the other cores have no determinants. Later on, we need to slice out
      !    just our portion of the coeffs list
      if (present(core_up)) then
        if (n_ref == 1) then
          call find_doubly_excited(n_det=ndets_new, dets_up=tmp_dets_up, dets_dn=tmp_dets_dn, ref_up=old_dets_up(1:n_ref), ref_dn=old_dets_dn(1:n_ref), norb=norb, n_core_orb=n_core_orb, ref_coeffs=old_wts(1:n_ref), ninitiator=n_ref, eps_var_pt=eps_var, min_H_already_done=min_H_already_done(1:n_ref),core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn,active_only=.true.)
        elseif (n_ref == 0) then
           call find_doubly_excited(n_det=ndets_new, dets_up=tmp_dets_up, dets_dn=tmp_dets_dn, ref_up=old_dets_up(1:n_ref), ref_dn=old_dets_dn(1:n_ref), norb=norb, n_core_orb=n_core_orb, ref_coeffs=old_wts(1:n_ref), ninitiator=n_ref, eps_var_pt=eps_var, min_H_already_done=min_H_already_done(1:n_ref),core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn,active_only=.true.)
        else
           call find_doubly_excited(n_det=ndets_new, dets_up=tmp_dets_up, dets_dn=tmp_dets_dn, ref_up=old_dets_up(1:n_ref), ref_dn=old_dets_dn(1:n_ref), norb=norb, n_core_orb=n_core_orb, ref_coeffs=old_wts(1:n_ref), ninitiator=n_ref, eps_var_pt=eps_var, min_H_already_done=min_H_already_done(1:n_ref),core_up=core_up,core_dn=core_dn,virt_up=virt_up,virt_dn=virt_dn,active_only=.true.)
        endif
      else
        if (n_ref == 1) then
          call find_doubly_excited(n_det=ndets_new, dets_up=tmp_dets_up, dets_dn=tmp_dets_dn, ref_up=old_dets_up(1:n_ref), ref_dn=old_dets_dn(1:n_ref), norb=norb, n_core_orb=n_core_orb, ref_coeffs=old_wts(1:n_ref), ninitiator=n_ref, eps_var_pt=eps_var, min_H_already_done=min_H_already_done(1:n_ref))
        elseif (n_ref == 0) then
           call find_doubly_excited(n_det=ndets_new, dets_up=tmp_dets_up, dets_dn=tmp_dets_dn, ref_up=old_dets_up(1:n_ref), ref_dn=old_dets_dn(1:n_ref), norb=norb, n_core_orb=n_core_orb, ref_coeffs=old_wts(1:n_ref), ninitiator=n_ref, eps_var_pt=eps_var, min_H_already_done=min_H_already_done(1:n_ref))
        else
           call find_doubly_excited(n_det=ndets_new, dets_up=tmp_dets_up, dets_dn=tmp_dets_dn, ref_up=old_dets_up(1:n_ref), ref_dn=old_dets_dn(1:n_ref), norb=norb, n_core_orb=n_core_orb, ref_coeffs=old_wts(1:n_ref), ninitiator=n_ref, eps_var_pt=eps_var, min_H_already_done=min_H_already_done(1:n_ref))
        endif
      endif
    else
      stop "Must have eps_var >= 0!"
    endif

    ! Sort new dets, so that the first few are the old dets in the original order
    ! and the rest are the new dets in order by label

    if (allocated(new_dets_up)) then
      deallocate(new_dets_up)
      deallocate(new_dets_dn)
    endif
    allocate(new_dets_up(ndets_new))
    allocate(new_dets_dn(ndets_new))

    new_dets_up(1:ndets_old) = old_dets_up(1:ndets_old)
    new_dets_dn(1:ndets_old) = old_dets_dn(1:ndets_old)

    ! Sort old dets so that they can be binary searched
    allocate(tmp_old_up(ndets_old))
    allocate(tmp_old_dn(ndets_old))
    tmp_old_up(1:ndets_old) = old_dets_up(1:ndets_old)
    tmp_old_dn(1:ndets_old) = old_dets_dn(1:ndets_old)

    call sort(tmp_old_up,tmp_old_dn)

    j = ndets_old
    if (present(local_det_map)) then
       if (ndets_old==1.and.local_det_map%n_it.eq.0) then
          !MJO - if this is our first time, we want to include HF in the new determinants
          !      list to send around FIXME: Make this check more safe
          !      We do this by setting offset to 1, meaning we won't overwrite the
          !      first element of the array later on
          offset = 1
       else
          offset = 0
       endif
    else
       offset = 0
    endif

    allocate(tmp_dets_up2(size(tmp_dets_up)))
    allocate(tmp_dets_dn2(size(tmp_dets_dn)))

    if (offset==1) then
      tmp_dets_up2(:) = old_dets_up(1)
      tmp_dets_dn2(:) = old_dets_dn(1)
    endif

    do i=1,ndets_new
      call binary_search(tmp_dets_up(i),tmp_dets_dn(i),tmp_old_up(1:ndets_old),tmp_old_dn(1:ndets_old),k)
      if (k==0) then ! det not in old list
        j = j+1
        new_dets_up(j) = tmp_dets_up(i)
        new_dets_dn(j) = tmp_dets_dn(i)
        !MJO - reuse the first part of the tmp_dets_up array as our send buffer
        !AAH - Don't reuse because this assumed that the first det is old, which
        !is not always true
        tmp_dets_up2(j-ndets_old+offset) = tmp_dets_up(i)
        tmp_dets_dn2(j-ndets_old+offset) = tmp_dets_dn(i)
      endif
    enddo
    n_det_tmp = j-ndets_old+offset
    n_det_tmp2 = ndets_new

    !MJO - Get the list of new, unique dets from every core
    if(present(remote_det_map)) then
      call mpi_allgatherv_new_dets(tmp_dets_up2,tmp_dets_dn2,n_det_tmp,n_det_tmp2,remote_det_map,local_det_map)
    endif

!    do i=1,local_det_map%ndets_global
!      write(6,'(2B64)') remote_det_map%global_dets_up(i),remote_det_map%global_dets_dn(i)
!    enddo
!    deallocate(tmp_dets_up,tmp_dets_dn)
!    deallocate(tmp_old_up,tmp_old_dn)

    if (present(remote_det_map)) then
       ndets_global_old = local_det_map%ndets_global_old
       ndets_global_new = local_det_map%ndets_global
    else
       ndets_global_old = ndets_old
       ndets_global_new = ndets_new
    endif

    call alloc('min_H_already_done',min_H_already_done,int(ndets_new,i8b))
    min_H_already_done(1:ndets_old)=min(min_H_already_done(1:ndets_old),(eps_var/abs(old_wts(1:ndets_old))-1.d-14))
    min_H_already_done(ndets_old+1:ndets_new) = 9.d99
!   if(ndets_global_old.gt.size(old_wts)) then
!     write(6,'(''ndets_global_old,size(old_wts)='',9i10)') ndets_global_old,size(old_wts)
!   endif
!   write(6,'(''old_wts='',1000es9.1)') old_wts(1:ndets_global_old)

    write(6,'(''ndets_global='',i10)') ndets_global_new

     ! MJO Below prints the full sorted list, which is useful for debugging in parallel
    ! allocate(tmp_old_up(local_det_map%ndets_global))
    ! allocate(tmp_old_dn(local_det_map%ndets_global))
    ! tmp_old_up(1:local_det_map%ndets_global) = remote_det_map%global_dets_up(1:local_det_map%ndets_global)
    ! tmp_old_dn(1:local_det_map%npdets_global) = remote_det_map%global_dets_dn(1:local_det_map%ndets_global)

    ! call sort(tmp_old_up,tmp_old_dn)

    ! do i=1,local_det_map%ndets_global
    !   write(6,*) tmp_old_up(i),tmp_old_dn(i)
    ! enddo

    ! deallocate(tmp_old_up,tmp_old_dn)
    call my_second(2,"get_next_det_list")

  end subroutine get_next_det_list
!-------------------------------------------------------------------------------------

  subroutine iterative_diagonalize(ndets,n_states,dets_up,dets_dn,final_wts,eigenvalues,starting_wts,remote_det_map,local_det_map,average_connections,use_hash_generation)
  ! A Holmes, 10 Apr 2016

    use chemistry, only : diagonalize_sparse_hamiltonian_chem_excited,time_sym,find_important_connected_dets_chem,get_new_diag_elem
    use heg, only: diagonalize_sparse_hamiltonian_heg_excited
    use common_run, only : connected_dets_up,connected_dets_dn,connected_matrix_elements,max_connected_dets,diag_elem_info,connected_diag_elems_info
    use generic_sort, only : sort,sort_by_first_argument
    use common_selected_ci, only : sparse_ham
    use mpi_routines, only : det_map,det_map_l

    integer,intent(in) :: ndets,n_states
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(inout) :: dets_up(:),dets_dn(:)
#else
    integer(ik),intent(inout) :: dets_up(:),dets_dn(:)
#endif
    real(rk),intent(out) :: final_wts(:,:),eigenvalues(:)
    real(rk),optional,intent(in) :: starting_wts(:,:)
    real(rk), optional,intent(inout) :: average_connections
    type(det_map), optional,intent(inout) :: remote_det_map
    type(det_map_l), optional,intent(inout) :: local_det_map
    logical, optional, intent(in)           :: use_hash_generation
    real(rk) :: tau_out_for_diag_ham
    integer,allocatable :: tmp_ind_singles(:),tmp_info_singles(:)
    integer,allocatable :: tmp_ind_doubles(:),tmp_info_doubles(:)
    integer,allocatable :: ref_connection_indices(:)
    integer :: n_allocate_singles, n_allocate_doubles
    integer :: n_hash_singles, n_hash_doubles
    integer :: n_connected_dets
    integer :: i_ref,i,j

    if (hamiltonian_type.ne.'chem' .and. hamiltonian_type .ne. 'heg') then
      stop "Only chem and heg implemented so far!"
    endif

    if(present(remote_det_map)) then
      if (hamiltonian_type .eq. 'chem') then
        if (present(use_hash_generation)) then
          call diagonalize_sparse_hamiltonian_chem_excited(dets_up(1:ndets),dets_dn(1:ndets),ndets,n_states,final_wts(1:ndets,1:n_states),eigenvalues,time_sym,initial_vectors=starting_wts(1:ndets,1:n_states),sparse_ham=sparse_ham,remote_det_map=remote_det_map,local_det_map=local_det_map,use_hash=use_hash_generation)
        else
          call diagonalize_sparse_hamiltonian_chem_excited(dets_up(1:ndets),dets_dn(1:ndets),ndets,n_states,final_wts(1:ndets,1:n_states),eigenvalues,time_sym,initial_vectors=starting_wts(1:ndets,1:n_states),sparse_ham=sparse_ham,remote_det_map=remote_det_map,local_det_map=local_det_map)
        endif
      elseif (hamiltonian_type .eq. 'heg') then
        write(6,*) 'Parallel heg not yet implemented'
      endif
    else
      if (hamiltonian_type .eq. 'chem') then
        call diagonalize_sparse_hamiltonian_chem_excited(dets_up(1:ndets),dets_dn(1:ndets),ndets,n_states,final_wts(1:ndets,1:n_states),eigenvalues,time_sym,initial_vectors=starting_wts(1:ndets,1:n_states),sparse_ham=sparse_ham,average_connections=average_connections)
      elseif (hamiltonian_type .eq. 'heg') then
        call diagonalize_sparse_hamiltonian_heg_excited(dets_up(1:ndets),dets_dn(1:ndets),ndets,n_states,final_wts(1:ndets,1:n_states),eigenvalues,time_sym,initial_vectors=starting_wts(1:ndets,1:n_states),sparse_ham=sparse_ham)
      endif
    endif

   !call diagonalize_sparse_hamiltonian_chem(dets_up(1:ndets),dets_dn(1:ndets),ndets,final_wts(1:ndets,1:n_states),lowest_energy,tau_out_for_diag_ham,time_sym,sort_or_not=.false.,initial_vector=starting_wts(1:ndets,1:n_states),sparse_ham=sparse_ham)

  end subroutine iterative_diagonalize
!-------------------------------------------------------------------------------------

  subroutine second_order_pt(ndets,dets_up,dets_dn,wts,diag_elems,var_energy,eps_pt,delta_e_2pt,ndets_connected,remote_det_map,local_det_map,core_up,core_dn,virt_up,virt_dn,active_only)
    ! Compute and print the second order energy lowering
    ! A Holmes, Feb 2016
    ! Modified : A Holmes, 23 May 2017. Optional inputs core_up/dn, virt_up/dn,
    ! active_only for performing the PT either only in the active space or only
    ! outside the active space (active_only=true or false, respectively)

    use semistoch, only : find_doubly_excited,hamiltonian
    use mpi_routines, only : det_map,mpi_allred,det_map_l
    use chemistry, only : time_sym,hamiltonian_chem_time_sym
    use more_tools, only : binary_search

    type(det_map),optional,intent(in) :: remote_det_map
    type(det_map_l),optional,intent(in) :: local_det_map
    integer,intent(in) :: ndets
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
    type(ik_vec),optional,intent(in) :: core_up,core_dn,virt_up,virt_dn
    type(ik_vec),allocatable :: connected_dets_up(:),connected_dets_dn(:)
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
    integer(ik),optional,intent(in) :: core_up,core_dn,virt_up,virt_dn
    integer(ik),allocatable :: connected_dets_up(:),connected_dets_dn(:)
#endif
    real(rk),intent(in) :: wts(:),diag_elems(:),var_energy,eps_pt
    real(rk),intent(out) :: delta_e_2pt
    integer,intent(out) :: ndets_connected
    logical,optional,intent(in) :: active_only

    real(rk),allocatable :: connected_wts(:),new_diag_elems(:),e_mix_den(:)
    integer :: i, ndets_connected_global,m

    call my_second(1,'Deterministic 2nd-order PT correction')
#ifdef DEBUG
    do i=1,ndets-1
      if (dets_up(i)>dets_up(i+1).or.(dets_up(i)==dets_up(i+1).and.dets_dn(i)>dets_dn(i+1))) then
        write(6,*) 'i',ndets,i,dets_up(i),dets_up(i+1),dets_dn(i),dets_dn(i+1)
        write (6,*) "Bug! second_order_pt called with unsorted reference!"
        stop "DEBUG"
      endif
    enddo
#endif

    if (present(core_up)) then
      if (eps_pt>0._rk) then
         call find_doubly_excited(n_det=ndets_connected, dets_up=connected_dets_up, dets_dn=connected_dets_dn, ref_up=dets_up(1:ndets), ref_dn=dets_dn(1:ndets), norb=norb, n_core_orb=n_core_orb, ref_coeffs=wts(1:ndets), ninitiator=ndets, e_mix_num=connected_wts, e_mix_den=e_mix_den, eps_var_pt=eps_pt, ref_diag_elems=diag_elems(1:ndets), new_diag_elems=new_diag_elems, core_up=core_up, core_dn=core_dn, virt_up=virt_up, virt_dn=virt_dn, active_only=active_only)
      else ! use all connections
         call find_doubly_excited(n_det=ndets_connected, dets_up=connected_dets_up, dets_dn=connected_dets_dn, ref_up=dets_up(1:ndets), ref_dn=dets_dn(1:ndets), norb=norb, n_core_orb=n_core_orb, ref_coeffs=wts(1:ndets), ninitiator=ndets, e_mix_num=connected_wts, e_mix_den=e_mix_den, ref_diag_elems=diag_elems(1:ndets), new_diag_elems=new_diag_elems, core_up=core_up, core_dn=core_dn, virt_up=virt_up, virt_dn=virt_dn, active_only=active_only)
      endif
    else
      if (eps_pt>0._rk) then
         call find_doubly_excited(n_det=ndets_connected, dets_up=connected_dets_up, dets_dn=connected_dets_dn, ref_up=dets_up(1:ndets), ref_dn=dets_dn(1:ndets), norb=norb, n_core_orb=n_core_orb, ref_coeffs=wts(1:ndets), ninitiator=ndets, e_mix_num=connected_wts, e_mix_den=e_mix_den, eps_var_pt=eps_pt, ref_diag_elems=diag_elems(1:ndets), new_diag_elems=new_diag_elems)
      else ! use all connections
         call find_doubly_excited(n_det=ndets_connected, dets_up=connected_dets_up, dets_dn=connected_dets_dn, ref_up=dets_up(1:ndets), ref_dn=dets_dn(1:ndets), norb=norb, n_core_orb=n_core_orb, ref_coeffs=wts(1:ndets), ninitiator=ndets, e_mix_num=connected_wts, e_mix_den=e_mix_den, ref_diag_elems=diag_elems(1:ndets), new_diag_elems=new_diag_elems)
      endif
    endif

    write (6,'(''Number of connected dets on this core='',i13,'' ='',es9.2)') ndets_connected, real(ndets_connected) ; call flush(6)
!   write (6,'(''PT Contribution from the first 20 dets (sorted by label)='')')

    delta_e_2pt = 0._rk
    do i=1,ndets_connected
      call binary_search(connected_dets_up(i),connected_dets_dn(i),dets_up(1:ndets),dets_dn(1:ndets),m)
      if (m==0) then ! Not in variational space, so contributes to E_PT
     !if (e_mix_den(i)==0._rk) then
        if (time_sym) then ! Must recompute new diag elems, because they aren't a simple function of the reference elems and excitation orbitals anymore
          call hamiltonian_chem_time_sym(connected_dets_up(i),connected_dets_dn(i),connected_dets_up(i),connected_dets_dn(i),new_diag_elems(i))
        endif
        delta_e_2pt = delta_e_2pt + connected_wts(i)**2/(var_energy-new_diag_elems(i))
       !if (i<=20)  write (6,*) i,connected_dets_up(i),connected_dets_dn(i),connected_wts(i),new_diag_elems(i),connected_wts(i)**2/(var_energy-new_diag_elems(i)),delta_e_2pt
      endif
    enddo
    if (present(remote_det_map)) then !MJO TEMP PARALLEL
       call mpi_allred(delta_e_2pt)
    endif

    ndets_connected_global=ndets_connected
    call mpi_allred(ndets_connected_global)
    write (6,'(''Number of connected dets on all cores='',i13,'' ='',es9.2)') ndets_connected_global, real(ndets_connected_global) ; call flush(6)

    call my_second(2,'Deterministic 2nd-order PT correction')

  end subroutine second_order_pt
!-------------------------------------------------------------------------------------

!  subroutine second_order_pt_dtm_batches(ndets,dets_up,dets_dn,wts,diag_elems,var_energy,eps_pt,n_batches)
!    ! Compute and print the second order energy lowering using Cyrus's suggestion of dividing up the
!    ! reference into batches and exciting from all pairs of batches.
!    ! Not currently used
!    ! A Holmes, 20 June 2016
!    use types, only : i8b
!    use semistoch, only : find_doubly_excited,hamiltonian
!    use tools, only : merge_sort2_up_dn
!    use more_tools, only : binary_search
!    use heg, only : energy_hf, energy_madelung
!
!    integer,intent(in) :: ndets
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
!    type(ik_vec),allocatable :: connected_dets_up_i(:),connected_dets_dn_i(:)
!    type(ik_vec),allocatable :: connected_dets_up_j(:),connected_dets_dn_j(:)
!    type(ik_vec),allocatable :: sorted_dets_up(:),sorted_dets_dn(:)
!    type(ik_vec),allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
!#else
!    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
!    integer(ik),allocatable :: connected_dets_up_i(:),connected_dets_dn_i(:)
!    integer(ik),allocatable :: connected_dets_up_j(:),connected_dets_dn_j(:)
!    integer(ik),allocatable :: sorted_dets_up(:),sorted_dets_dn(:)
!    integer(ik),allocatable ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
!#endif
!    real(rk),intent(in) :: wts(:),diag_elems(:),var_energy,eps_pt
!    integer,intent(in) :: n_batches
!
!    real(rk) :: delta_e_2pt
!    real(rk),allocatable :: connected_wts_i(:),new_diag_elems_i(:),e_mix_den_i(:)
!    real(rk),allocatable :: connected_wts_j(:),new_diag_elems_j(:),e_mix_den_j(:)
!    integer :: i,j,ndets_connected,ndets_connected_i,ndets_connected_j
!    integer :: ibatch,jbatch
!    integer(i8b) :: k
!    integer,allocatable :: iorder(:),temp_i_2(:)
!
!    call my_second(1,'Deterministic batches 2nd-order PT correction')
!
!    allocate(sorted_dets_up(ndets))
!    allocate(sorted_dets_dn(ndets))
!
!    sorted_dets_up(:) = dets_up(1:ndets)
!    sorted_dets_dn(:) = dets_dn(1:ndets)
!
!    ! Sort by label for binary search
!    allocate(iorder(ndets))
!    allocate(temp_i16_up((ndets+1)/2))
!    allocate(temp_i16_dn((ndets+1)/2))
!    allocate(temp_i_2((ndets+1)/2))
!    do j=1,ndets
!      iorder(j)=j
!    enddo
!
!    call merge_sort2_up_dn(sorted_dets_up(1:ndets),sorted_dets_dn(1:ndets), iorder, ndets, temp_i16_up, temp_i16_dn, temp_i_2)
!
!    deallocate(iorder)
!    deallocate(temp_i16_up)
!    deallocate(temp_i16_dn)
!    deallocate(temp_i_2)
!
!    delta_e_2pt = 0._rk
!    ndets_connected = 0
!
!    ! Break up determinants into batches by dealing them out ("block cyclic distribution"),
!    ! since we don't want all of the big-wt variational dets in one batch because they have more connections exceeding eps2
!
!    do ibatch=1,n_batches
!
!      call find_doubly_excited(n_det=ndets_connected_i, dets_up=connected_dets_up_i, dets_dn=connected_dets_dn_i, ref_up=dets_up(ibatch:ndets:n_batches), ref_dn=dets_dn(ibatch:ndets:n_batches), norb=norb, n_core_orb=n_core_orb, ref_coeffs=wts(ibatch:ndets:n_batches), e_mix_num=connected_wts_i, e_mix_den=e_mix_den_i, eps_var_pt=eps_pt, ref_diag_elems=diag_elems(ibatch:ndets:n_batches), new_diag_elems=new_diag_elems_i)
!
!      do jbatch=ibatch,n_batches
!
!        call find_doubly_excited(n_det=ndets_connected_j, dets_up=connected_dets_up_j, dets_dn=connected_dets_dn_j, ref_up=dets_up(jbatch:ndets:n_batches), ref_dn=dets_dn(jbatch:ndets:n_batches), norb=norb, n_core_orb=n_core_orb, ref_coeffs=wts(jbatch:ndets:n_batches), e_mix_num=connected_wts_j, e_mix_den=e_mix_den_j, eps_var_pt=eps_pt, ref_diag_elems=diag_elems(jbatch:ndets:n_batches), new_diag_elems=new_diag_elems_j)
!
!
!        ! Find which connections k are in common in the i and j batches, since only terms of the form c_i H_ik H_kj c_j contribute
!        j=1
!        do i=1,ndets_connected_i
!          ! Increment j until it is no longer smaller than i
!          do while (connected_dets_up_i(i)>connected_dets_up_j(j).or.(connected_dets_up_i(i)==connected_dets_up_j(j).and.connected_dets_dn_i(i)>connected_dets_dn_j(j)))
!            if (j==ndets_connected_j)  exit
!            j = j+1
!          enddo
!          if (connected_dets_up_i(i)==connected_dets_up_j(j).and.connected_dets_dn_i(i)==connected_dets_dn_j(j)) then
!            if (n_batches>1) then
!              call binary_search(connected_dets_up_i(i),connected_dets_dn_i(i),sorted_dets_up(1:ndets),sorted_dets_dn(1:ndets),k)
!              if (k==0) then
!                if (ibatch==jbatch) then
!                  delta_e_2pt = delta_e_2pt + connected_wts_i(i)*connected_wts_j(j)/(var_energy-new_diag_elems_i(i))
!                else
!                  delta_e_2pt = delta_e_2pt + 2*connected_wts_i(i)*connected_wts_j(j)/(var_energy-new_diag_elems_i(i))
!                endif
!               !write (6,*) i,connected_dets_up_i(i),connected_dets_dn_i(i),connected_wts_i(i),new_diag_elems_i(i),connected_wts_i(i)**2/(var_energy-new_diag_elems_i(i)),delta_e_2pt
!              endif
!            else ! n_batches==1
!              if (e_mix_den_i(i)==0._rk) then
!                if (ibatch==jbatch) then
!                  delta_e_2pt = delta_e_2pt + connected_wts_i(i)*connected_wts_j(j)/(var_energy-new_diag_elems_i(i))
!                else
!                  delta_e_2pt = delta_e_2pt + 2*connected_wts_i(i)*connected_wts_j(j)/(var_energy-new_diag_elems_i(i))
!                endif
!               !write (6,*) i,connected_dets_up_i(i),connected_dets_dn_i(i),connected_wts_i(i),new_diag_elems_i(i),connected_wts_i(i)**2/(var_energy-new_diag_elems_i(i)),delta_e_2pt
!              endif
!            endif
!          endif
!        enddo
!
!      enddo
!
!      ndets_connected = ndets_connected + ndets_connected_i
!
!    enddo
!
!    write (6,'(''Variational energy='',t36,f15.9)') var_energy
!    write (6,'(''2nd-order PT energy lowering='',t36,f15.9)') delta_e_2pt
!    write (6,'(''Total energy='',t36,f15.9)') var_energy+delta_e_2pt
!    write (6,'(''ndets, ndets_connected(total), Variational, PT, Total Energies='',i9,i11,f16.9,f15.9,f16.9)') ndets, ndets_connected, var_energy, delta_e_2pt, var_energy+delta_e_2pt
!    if(hamiltonian_type.eq.'heg') then
!      write(6,'(''Total energy (includ. Madelung)='',t36,f15.9)') var_energy+delta_e_2pt+energy_madelung
!      write(6,'(''Correlation energy ='',t36,f15.9)') var_energy+delta_e_2pt-energy_hf
!      write(6,'(''Correlation energy per electron ='',t36,f15.9)') (var_energy+delta_e_2pt-energy_hf)/nelec
!    endif
!    call flush(6)
!
!    call my_second(2,'Deterministic batches 2nd-order PT correction')
!
!  end subroutine second_order_pt_dtm_batches
!-------------------------------------------------------------------------------------

  subroutine second_order_pt_alias(ndets,dets_up,dets_dn,wts,diag_elems,var_energy,eps_pt,n_mc,target_error,eps_pt_big,pt_energy,pt_energy_std_dev,ndets_connected,pt_big,core_up,core_dn,virt_up,virt_dn,active_only,active_tot_energy)
    ! Compute and print the second order energy lowering stochastically using Alias method
    ! Follows Cyrus' notes
    ! A Holmes, 15 June 2016
    use semistoch, only : find_doubly_excited,hamiltonian
    use tools, only : merge_sort2_up_dn, sort_and_merge_count_repeats, welford
    use more_tools, only : setup_alias,sample_alias,binary_search
    use mpi_routines, only : mpi_allred,ncores_on_node,whoami,mpi_barr,master_core_node,&
      &mpi_barr_in_node,shmem_allocate,shmem_deallocate,mpierr,nnodes,ncores,whoaminode
    use chemistry, only : time_sym,hamiltonian_chem_time_sym
    use generic_sort, only : sort,shell_sort_real_rank1_int_rank1
#ifdef MPI
    INCLUDE 'mpif.h'
#endif

    integer,intent(in) :: ndets
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
    type(ik_vec),optional,intent(in) :: core_up,core_dn,virt_up,virt_dn
    type(ik_vec),allocatable :: connected_dets_up(:),connected_dets_dn(:)
    type(ik_vec),allocatable :: sampled_dets_up(:),sampled_dets_dn(:)
    type(ik_vec),dimension(1) :: ref_up,ref_dn
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
    integer(ik),optional,intent(in) :: core_up,core_dn,virt_up,virt_dn
    integer(ik),allocatable :: connected_dets_up(:),connected_dets_dn(:)
    integer(ik),allocatable :: sampled_dets_up(:),sampled_dets_dn(:)
    integer(ik),dimension(1) :: ref_up,ref_dn
#endif
    real(rk),intent(in) :: wts(:),diag_elems(:),var_energy,eps_pt
    integer,intent(inout) :: n_mc
    real(rk),intent(in) :: target_error ! Target standard deviation in the energy (run until the sample error is less than this)
    real(rk),intent(in) :: eps_pt_big
    real(rk),intent(out) :: pt_energy,pt_energy_std_dev
    integer,intent(out) :: ndets_connected
    real(rk), intent(in) :: pt_big
    logical,optional,intent(in) :: active_only
    real(rk),optional,intent(in) :: active_tot_energy

    real(rk) :: delta_e_2pt,S_e2,var_e2 ! mean, Welford 'S', variance of 2nd order PT energy correction
    real(rk),allocatable :: new_diag_elems(:),e_mix_den(:)
    real(rk),allocatable :: term1(:),term2(:)
    real(rk),allocatable :: term1_big(:),term2_big(:)
    integer :: i
    real(rk),pointer,dimension(:) :: prob,q
    integer,pointer,dimension(:) :: smaller,larger,J
    real(rk) :: norm
    integer,allocatable :: samples(:),counts(:),my_samples(:),iorder(:)
    integer :: sample,max_samples,n_mc_diff
    real(rk),dimension(1) :: coeff,ref_diag_elem
    real(rk),allocatable :: coeffs(:),probs(:),w_over_p(:),sampled_diag_elems(:)
    real(rk) :: e_2pt_tmp
    integer :: batch_count,whoaminode_backwards
    real(rk) :: e_2pt_this_sample, avail, avail2, avail3, min_avail2
    integer :: m, k, n_mc_begin, first_sample, my_n_mc,my_start,my_end,remainder
    integer :: max_ndets_connected, j_node, k_core, my_node
    real(rk) :: target_variance

    call my_second(1,'Alias method 2nd-order PT correction')
#ifdef DEBUG
    do i=1,ndets-1
      if (dets_up(i)>dets_up(i+1).or.(dets_up(i)==dets_up(i+1).and.dets_dn(i)>dets_dn(i+1))) then
        write (6,*) "Bug! second_order_pt_alias called with unsorted reference!"
        stop "DEBUG"
      endif
    enddo
#endif
    delta_e_2pt = 0._rk
    batch_count = 1
    ! First, compute and store probability distribution
    call shmem_allocate(prob,int(ndets,i8b))
    call shmem_allocate(smaller,int(ndets,i8b))
    call shmem_allocate(larger,int(ndets,i8b))
    call shmem_allocate(J,int(ndets,i8b))
    call shmem_allocate(q,int(ndets,i8b))

    if(master_core_node) then
      norm = sum(abs(wts))
      prob = abs(wts)/norm
      call setup_alias(ndets,prob,smaller,larger,J,q)
    endif
    call mpi_barr_in_node()
    max_samples = 10**6

! If n_mc is not in input or n_mc<=0 in input, use an initial 10 samples to estimate how large n_mc can be.
    if(n_mc > 0) then
      n_mc_begin=n_mc
      if (n_mc.lt.ncores) then
        write(6,*) 'Warning! N_MC must be greater than or equal to 2*ncores! Setting N_MC=2*ncores'
        n_mc = 2*ncores
      endif
      first_sample=1
    else
      n_mc_begin=10*ncores
      n_mc=10*ncores
      first_sample=0
    endif

    allocate(samples(n_mc_begin))
    allocate(counts(n_mc_begin))
    allocate(iorder(n_mc_begin))
    allocate(my_samples(2*n_mc_begin/ncores))
    allocate(sampled_dets_up(min(n_mc_begin,ndets)))
    allocate(sampled_dets_dn(min(n_mc_begin,ndets)))
    allocate(sampled_diag_elems(min(n_mc_begin,ndets)))
    allocate(coeffs(min(n_mc_begin,ndets)))
    allocate(probs(min(n_mc_begin,ndets)))
    allocate(w_over_p(min(n_mc_begin,ndets))) ! counts / probs


    ! Sample max_samples samples or as many are needed for to achieve target_error.
    ! Use Alias method to pick each sample of size n_mc.

    write (6,'(''Performing stochastic PT until error bar drops below'',es10.2)') target_error ; call flush(6)
    write (6,'(''This is the difference between the PT correction for eps_pt and eps_pt_big'',/)') ; call flush(6)

    target_variance = target_error**2
    e_2pt_tmp = 0._rk
    do sample=first_sample,max_samples

     !call my_second(1, 'sample in second_order_pt_alias')

      e_2pt_this_sample = 0._rk

      do i=1,n_mc
        samples(i) = sample_alias(ndets,J,q)
        !MJO It is assumed that all cores will be drawing the same random sample here.
        !    If the are not (such as having different random number seeds), there will
        !    be a very subtle bug that is hard to find.
      enddo

      ! Now, combine all samples into distinct determinants, with wts counting number of copies of each !
      n_mc_diff = n_mc

      ! Sort and merge the list of integers samples(1:n_mc_diff) and their counts(1:n_mc_diff)
      call sort_and_merge_count_repeats(n_mc_diff,samples,counts)
      if (n_mc_diff<ncores) then
        write(6,*) 'Error: The sample size with duplicates removed needs to be'
        write(6,*) '       at least ncores! n_mc_diff: ',n_mc_diff,' ncores: ',ncores
        call flush(6)
        call mpi_barr()
        stop
      endif

     !call my_second(2, 'sort_and_merge_count_repeats')

      !MJO First, slice out the full sample
      sampled_dets_up(1:n_mc_diff) = dets_up(samples(1:n_mc_diff))
      sampled_dets_dn(1:n_mc_diff) = dets_dn(samples(1:n_mc_diff))
      coeffs(1:n_mc_diff) = wts(samples(1:n_mc_diff)) ! {c_i}
      probs(1:n_mc_diff) = prob(samples(1:n_mc_diff)) ! {p_i}
      sampled_diag_elems = diag_elems(samples(1:n_mc_diff)) ! {H_{ii}}
      w_over_p(1:n_mc_diff) = counts(1:n_mc_diff)/probs(1:n_mc_diff)

      do i=1,n_mc_diff
        iorder(i) = i
      enddo
      !MJO Sort the coeffs array and store sort in iorder array
      call shell_sort_real_rank1_int_rank1(coeffs(1:n_mc_diff),iorder(1:n_mc_diff))
      !Sort the rest of the arrays using iorder

      sampled_dets_up(1:n_mc_diff) = sampled_dets_up(iorder(1:n_mc_diff))
      sampled_dets_dn(1:n_mc_diff) = sampled_dets_dn(iorder(1:n_mc_diff))
      probs(1:n_mc_diff) = probs(iorder(1:n_mc_diff))
      sampled_diag_elems(1:n_mc_diff) = sampled_diag_elems(iorder(1:n_mc_diff)) ! {H_{ii}}
      w_over_p(1:n_mc_diff) = w_over_p(iorder(1:n_mc_diff))

      ! MJO Get our indices from the node interleaved distribution
      j_node = 0
      k_core = 0
      my_node = whoami/ncores_on_node !MJO Calculate node number using integer division
      ! MJO Calculate the 'backwards' whoami, i.e. with 16 ranks/node
      ! 0 => 15, 15 => 0, 1=>14, 14=>1, etc
      whoaminode_backwards = (ncores_on_node -1) - whoaminode

      do i=1,n_mc_diff
        if(mod(i,nnodes).eq.my_node) then
          !This sample belongs to my node
          j_node = j_node + 1
          !mod(j_node/ncores_on_node,2) is a way of checking whether we have
          !already assigned an even or odd number of dets per core. This will
          !allow for switching between the backwards and forwards sweep
          if(mod((j_node-1)/ncores_on_node,2).eq.0.and.mod(j_node,ncores_on_node).eq.whoaminode) then
            !MJO Forward sweep
            k_core = k_core + 1
            !MJO my_samples is an array of the indices of the array samples
            !    which belong to me
            my_samples(k_core) = i
          else if(mod((j_node-1)/ncores_on_node,2).eq.1.and.mod(j_node,ncores_on_node).eq.whoaminode_backwards) then
            !MJO Backwards sweep
            k_core = k_core + 1
            !MJO my_samples is an array of the indices of the array samples
            !    which belong to me
            my_samples(k_core) = i
          endif
        endif
      enddo

      my_n_mc = k_core

      sampled_dets_up(1:my_n_mc) = sampled_dets_up(my_samples(1:my_n_mc))
      sampled_dets_dn(1:my_n_mc) = sampled_dets_dn(my_samples(1:my_n_mc))
      coeffs(1:my_n_mc) = coeffs(my_samples(1:my_n_mc)) ! {c_i}
      probs(1:my_n_mc) = probs(my_samples(1:my_n_mc)) ! {p_i}
      sampled_diag_elems(1:my_n_mc) = sampled_diag_elems(my_samples(1:my_n_mc)) ! {H_{ii}}
      w_over_p(1:my_n_mc) = w_over_p(my_samples(1:my_n_mc))


      ! MJO Old contiguous distribution below
      ! my_n_mc = n_mc_diff/ncores

      ! ! MJO If we don't have perfect division we want to add some to get good load balancing
      ! remainder = mod(n_mc_diff,ncores)
      ! ! MJO Remainder is the number that needs to be added;
      ! !     we will add 1 to each of the final remainder cores
      ! if (whoami.gt.(ncores-remainder-1)) then
      !   my_n_mc  = my_n_mc + 1
      !   my_start = (my_n_mc-1)*(whoami) + 1 - (ncores-remainder-whoami)
      !   my_end   = my_start + my_n_mc - 1
      ! else
      !   my_start = my_n_mc*(whoami) + 1
      !   my_end   = my_start + my_n_mc - 1
      ! endif
      ! write(6,*) 'my_n_mc',my_n_mc
      ! sampled_dets_up(1:my_n_mc) = dets_up(samples(my_start:my_end))
      ! sampled_dets_dn(1:my_n_mc) = dets_dn(samples(my_start:my_end))
      ! coeffs(1:my_n_mc) = wts(samples(my_start:my_end)) ! {c_i}
      ! write(6,*) 'my_coeffs',(coeffs(1:my_n_mc))

      ! write(6,*) 'my_coeffs_sum',sum(coeffs(1:my_n_mc))
      ! write(6,*) 'my_det_indices',samples(my_start:my_end)
      ! probs(1:my_n_mc) = prob(samples(my_start:my_end)) ! {p_i}
      ! write(6,*) 'my_probs',(probs(1:my_n_mc))
      ! write(6,*) 'my_counts',counts(my_start:my_end)
      ! sampled_diag_elems = diag_elems(samples(my_start:my_end)) ! {H_{ii}}
      ! !MJO we don't need to slice probs again - we already sliced it once
      ! !    we have not sliced counts yet, though
      ! w_over_p(1:my_n_mc) = counts(my_start:my_end)/probs(1:my_n_mc)

      ! MJO - Old single core sampling code
      ! sampled_dets_up(1:n_mc_diff) = dets_up(samples(1:n_mc_diff))
      ! sampled_dets_dn(1:n_mc_diff) = dets_dn(samples(1:n_mc_diff))
      ! coeffs(1:n_mc_diff) = wts(samples(1:n_mc_diff)) ! {c_i}
      ! probs(1:n_mc_diff) = prob(samples(1:n_mc_diff)) ! {p_i}
      ! sampled_diag_elems = diag_elems(samples(1:n_mc_diff)) ! {H_{ii}}
      ! w_over_p(1:n_mc_diff) = counts(1:n_mc_diff)/probs(1:n_mc_diff)


      if (present(core_up)) then
        call find_doubly_excited(n_det=ndets_connected, dets_up=connected_dets_up, dets_dn=connected_dets_dn, ref_up=sampled_dets_up(1:my_n_mc), ref_dn=sampled_dets_dn(1:my_n_mc), norb=norb, n_core_orb=n_core_orb, ref_coeffs=coeffs(1:my_n_mc), e_mix_num=term1, e_mix_den=term2, term1_big=term1_big, term2_big=term2_big, eps_var_pt=eps_pt, eps_var_pt_big=eps_pt_big, ref_diag_elems=sampled_diag_elems(1:my_n_mc), new_diag_elems=new_diag_elems, n_mc=n_mc, w_over_p=w_over_p(1:my_n_mc), core_up=core_up, core_dn=core_dn, virt_up=virt_up, virt_dn=virt_dn, active_only=active_only)
      else
        call find_doubly_excited(n_det=ndets_connected, dets_up=connected_dets_up, dets_dn=connected_dets_dn, ref_up=sampled_dets_up(1:my_n_mc), ref_dn=sampled_dets_dn(1:my_n_mc), norb=norb, n_core_orb=n_core_orb, ref_coeffs=coeffs(1:my_n_mc), e_mix_num=term1, e_mix_den=term2, term1_big=term1_big, term2_big=term2_big, eps_var_pt=eps_pt, eps_var_pt_big=eps_pt_big, ref_diag_elems=sampled_diag_elems(1:my_n_mc), new_diag_elems=new_diag_elems, n_mc=n_mc, w_over_p=w_over_p(1:my_n_mc))
      endif

      call my_second(2, 'find_doubly_excited')

      ! If n_mc was <=0 in input, then pick n_mc based on available memory
      if(sample.eq.0) then
        call mem_avail(avail,avail2,avail3) ; write(6,'(''after find_doubly_excited in second_order_pt_alias, avail, avail2,avail3='',9f11.2)') avail, avail2, avail3
#ifdef MPI
        call MPI_ALLREDUCE(avail2,min_avail2,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpierr)
        call MPI_ALLREDUCE(ndets_connected,max_ndets_connected,1,MPI_INTEGER4,MPI_MAX,MPI_COMM_WORLD,mpierr)
#else
        min_avail2 = avail2
        max_ndets_connected = ndets_connected
#endif
        if (time_sym) then
          n_mc=round_i(int(0.5e6*n_mc_begin/ncores*min_avail2/real(real(max_ndets_connected)*200*ncores_on_node)),2)*ncores ! The 1e6 is for converting MB to B, and the 200 is 8 times an ill-defined factor
        else
          n_mc=round_i(int(0.5e6*n_mc_begin/ncores*min_avail2/real(real(max_ndets_connected)*200*ncores_on_node)),2)*ncores ! The 1e6 is for converting MB to B, and the 200 is 8 times an ill-defined factor
        endif
        !MJO ncores_on_node so that we don't overallocate in parallel. Multiply by ncores because we
        !   want the global population size; use min_avail2 and max_ndets_connected to be conservative
        !   with our guess
        !   Also, since we multiplied n_mc_begin by ncores, we need to divide here to get
        !   the number per core that gave us the min_avail2 memory and max_ndets_connected
        write(6,'(''n_mc_begin, avail2, ndets_connected, ncores_on_node'',i5,es9.2,i10,i4)') n_mc_begin, avail2, ndets_connected, ncores_on_node
        write(6,'(''min_avail2, max_ndets_connected'',es9.2,i10)') min_avail2, max_ndets_connected
        if (n_mc.lt.ncores) then
          write(6,*) 'Warning! N_MC must be greater than or equal to 2*ncores! Setting N_MC=2*ncores'
          n_mc = 2*ncores
        endif
        write (6,'(/,''Computing PT energy using single-list stochastic method with computed N_MC='',i6,'' determinants per sample'',/)') n_mc ; call flush(6)
        deallocate(samples, counts, sampled_dets_up, sampled_dets_dn, sampled_diag_elems, coeffs, probs, w_over_p, iorder,my_samples)
        allocate(samples(n_mc))
        allocate(counts(n_mc))
        allocate(iorder(min(n_mc,ndets)))
        allocate(my_samples(2*min(n_mc,ndets)/ncores))
        allocate(sampled_dets_up(min(n_mc,ndets)))
        allocate(sampled_dets_dn(min(n_mc,ndets)))
        allocate(sampled_diag_elems(min(n_mc,ndets)))
        allocate(coeffs(min(n_mc,ndets)))
        allocate(probs(min(n_mc,ndets)))
        allocate(w_over_p(min(n_mc,ndets))) ! counts / probs
        cycle
      endif

      ! term1(k) = sum_i^{N_MC_diff} H_{ki} c_i w_i / p_i
      ! term2(k) = sum_i^{N_MC_diff} (H_{ki} c_i)**2 * ( (n_mc-1) * w_i / p_i - w_i**2 / p_i**2 )

      do k=1,ndets_connected
#ifdef DEBUG
      if (time_sym.and.connected_dets_up(k)>connected_dets_dn(k)) then
       write (6,*) "Bug! Time reversal shouldn't allow up > dn! k=",k; call flush(6)
       stop "DEBUG"
      endif
#endif
        ! Throw out k values corresponding to dets in variational space!!
        call binary_search(connected_dets_up(k),connected_dets_dn(k),dets_up(1:ndets),dets_dn(1:ndets),m)
        if (m==0) then ! Not in variational space, so contributes to E_PT
          if (time_sym) then ! Must recompute new diag elems, because they aren't a simple formula anymore
            call hamiltonian_chem_time_sym(connected_dets_up(k),connected_dets_dn(k),connected_dets_up(k),connected_dets_dn(k),new_diag_elems(k))
          endif
          e_2pt_this_sample = e_2pt_this_sample + 1._rk/(var_energy-new_diag_elems(k)) * (term1(k)**2 + term2(k) - term1_big(k)**2 - term2_big(k))
          !write(6,'(''k, term1(k), term2(k), term1_big(k), term2_big(k), e_2pt_this_sample'',i7,9es12.4)') k, term1(k), term2(k), term1_big(k), term2_big(k), e_2pt_this_sample
         endif
      enddo ! k

     !call my_second(2, 'binary_search')

      if (.false.) then !Don't batch for now, maybe implement later
        !MJO Previously, we were using independent, local populations.
        !    Now, we use global populations. Batching was for the local
        !    populations. The code is left here for now; should be removed
        !    when global populations are thoroughly vetted
        if (mod(sample,batch_count)==0) then
          e_2pt_this_sample = e_2pt_this_sample / (n_mc*real(n_mc - 1))
          e_2pt_tmp = e_2pt_tmp + e_2pt_this_sample
          e_2pt_tmp = e_2pt_tmp/batch_count
          call mpi_allred(e_2pt_tmp)
          e_2pt_tmp = e_2pt_tmp / ncores
          !            e_2pt_this_sample = e_2pt_this_sample / (n_mc*real(n_mc-1))
          call welford(sample/batch_count,e_2pt_tmp,delta_e_2pt,S_e2,var_e2)
          e_2pt_tmp = 0
        else
          e_2pt_this_sample = e_2pt_this_sample / (n_mc*real(n_mc - 1))
          e_2pt_tmp = e_2pt_tmp + e_2pt_this_sample
        endif
      else
        call mpi_allred(e_2pt_this_sample)
        e_2pt_this_sample = e_2pt_this_sample / (n_mc*real(n_mc-1))
        call welford(sample,e_2pt_this_sample,delta_e_2pt,S_e2,var_e2)
      endif
!      e_2pt_this_sample = e_2pt_this_sample / (n_mc*real(n_mc-1))

!      call welford(sample,e_2pt_this_sample,delta_e_2pt,S_e2,var_e2)
      write (6,*)
      if (present(active_tot_energy)) then
        write (6,'(''Sample, E_2pt_now, E_2pt estimate, total energy='',i6,f15.9,f12.8,f15.8,'' +-'',f12.8)') sample, e_2pt_this_sample, delta_e_2pt, active_tot_energy + pt_big + delta_e_2pt, sqrt(var_e2) ; call flush(6)
      else
        write (6,'(''Sample, E_2pt_now, E_2pt estimate, total energy='',i6,f15.9,f12.8,f15.8,'' +-'',f12.8)') sample, e_2pt_this_sample, delta_e_2pt, var_energy + pt_big + delta_e_2pt, sqrt(var_e2) ; call flush(6)
      endif
      write (6,*)
!     write (6,'(''Total energy = '',f15.9,'' +- '',f12.9)') var_energy + pt_big + delta_e_2pt, sqrt(var_e2)

      if (sample.ge.10*batch_count .and. var_e2<target_variance)  exit

    enddo ! sample

    call my_second(2,'Alias method 2nd-order PT correction')

    pt_energy = delta_e_2pt
    pt_energy_std_dev = sqrt(var_e2)

    call shmem_deallocate(prob)
    call shmem_deallocate(smaller)
    call shmem_deallocate(larger)
    call shmem_deallocate(J)
    call shmem_deallocate(q)

  end subroutine second_order_pt_alias
!-------------------------------------------------------------------------------------

!  subroutine second_order_pt_multi_index_hashing(ndets,dets_up,dets_dn,wts,diag_elems,var_energy,eps_pt,iorder,temp_i16_up,temp_i16_dn,temp_i_2)
!    ! Compute and print the second order energy lowering
!    ! Not currently used
!    use types, only : i8b
!    use semistoch, only : find_doubly_excited,hamiltonian
!    use chemistry, only : find_important_connected_dets_chem,get_new_diag_elem
!    use heg, only : energy_hf, energy_madelung
!    use common_run, only : connected_dets_up,connected_dets_dn,connected_matrix_elements,max_connected_dets,diag_elem_info,connected_diag_elems_info
!
!    integer,intent(in) :: ndets
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
!    type(ik_vec),intent(inout) ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
!    type(ik_vec),allocatable :: potential_connected_dets_up(:),potential_connected_dets_dn(:)
!#else
!    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
!    integer(ik),intent(inout) ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
!    integer(ik),allocatable :: potential_connected_dets_up(:),potential_connected_dets_dn(:)
!#endif
!    real(rk),intent(in) :: wts(:),var_energy,eps_pt
!    integer,intent(inout) :: iorder(:),temp_i_2(:)
!
!    real(rk) :: delta_e_2pt
!    integer :: i,j,k
!    integer :: n_hash
!    real(rk),allocatable :: diag_elems(:)
!    type(neighbor_table) :: neighbors,ref_dets
!    real(rk) :: H_ii,H_ij
!    real(rk) :: old_sum
!    integer :: n_connected_dets,n_connections
!    real(rk),allocatable :: potential_connected_coeffs(:)
!    logical :: is_con
!    integer(i8b) :: n_potential, n_actual
!
!    n_potential = 0_i8b
!    n_actual = 0_i8b
!
!    call my_second(1,'multi_index_hashing 2nd-order PT correction')
!
!    n_hash = 5 ! 1 more than the number of orbitals that change occupancy in a double excitation
!
!    call init_hash_fn ! Sets up hash function
!
!    allocate(potential_connected_dets_up(n_hash*ndets))
!    allocate(potential_connected_dets_dn(n_hash*ndets))
!    allocate(potential_connected_coeffs(n_hash*ndets))
!
!    ! First, create neighbor tables, for putting reference dets into, so that connected dets can quickly check for potential connections
!    call create_neighbor_tables(ndets,dets_up,dets_dn,n_hash,n_bins,neighbors) !,orb_group_mat)
!    write (6,'(''Done creating neighbor tables'')')
!
!    ! Note: "neighbors" has been allocated but is currently empty!
!    ! This is because no reference dets have been visited yet
!
!    ! Also, create a table of just the ref_dets, to enable us to quickly check whether a connected det should be thrown out
!    ! (since the PT expression sums over only dets not in reference)
!    call create_neighbor_tables(ndets,dets_up,dets_dn,1,n_bins,ref_dets)
!
!    ! This table should be filled up with all references!
!    call fill_neighbor_tables(ndets,dets_up,dets_dn,wts,ref_dets)
!
!
!    ! Next, loop over all references:
!    ! For each, compute diagonal elements, tthen generate all connections exceeding eps_pt
!    ! For each connection, query neighbors to find potential connections
!    ! After all connections for a given reference det have been looped over,
!    ! add that reference det to neighbors
!
!    delta_e_2pt = 0._rk
!
!    do i=1,ndets
!
!      call find_important_connected_dets_chem(dets_up(i), dets_dn(i), eps_pt/abs(wts(i)), n_connected_dets, connected_dets_up, connected_dets_dn,connected_matrix_elements,diag_elems(i),connected_diag_elems_info) ! connected_diag_elems_info is needed for computing diagonal elements in O(N) time, but don't use more than once per connected det!
!
!      if (n_connected_dets==0)  cycle
!
!      do j=2,n_connected_dets
!
!        ! Check by a single hash table whether it is a reference det
!        if (is_in_reference(connected_dets_up(j),connected_dets_dn(j),ref_dets))  cycle
!
!        call get_potential_connections(connected_dets_up(j),connected_dets_dn(j),neighbors,n_connections,potential_connected_dets_up,potential_connected_dets_dn,potential_connected_coeffs,iorder,temp_i16_up,temp_i16_dn,temp_i_2) !,orb_group_mat)
!
!        n_potential = n_potential + int(n_connections,i8b)
!
!        old_sum = 0._rk
!        if (n_connections>0) then
!          if (is_connected(connected_dets_up(j),connected_dets_dn(j),potential_connected_dets_up(1),potential_connected_dets_dn(1))) then
!            n_actual = n_actual + 1_i8b
!            call hamiltonian(connected_dets_up(j),connected_dets_dn(j),potential_connected_dets_up(1),potential_connected_dets_dn(1),H_ij,is_con)
!            if (abs(H_ij*potential_connected_coeffs(1))>eps_pt)  old_sum = H_ij*potential_connected_coeffs(1)
!          endif
!          do k=2,n_connections
!            if (potential_connected_dets_up(k)==potential_connected_dets_up(k-1).and.potential_connected_dets_dn(k)==potential_connected_dets_dn(k-1))  cycle
!            if (.not.is_connected(connected_dets_up(j),connected_dets_dn(j),potential_connected_dets_up(k),potential_connected_dets_dn(k)))  cycle
!            n_actual = n_actual + 1_i8b
!            call hamiltonian(connected_dets_up(j),connected_dets_dn(j),potential_connected_dets_up(k),potential_connected_dets_dn(k),H_ij,is_con)
!            if (abs(H_ij*potential_connected_coeffs(k))>eps_pt)  old_sum = old_sum + H_ij*potential_connected_coeffs(k)
!          enddo ! potential connections (in reference det list)
!        endif
!
!        ! get H_ii
!        if (connected_diag_elems_info(j)%old_diag_elem>1.e50_rk) then ! single excitation; compute the easy way for now
!          call hamiltonian(connected_dets_up(j),connected_dets_dn(j),connected_dets_up(j),connected_dets_dn(j),H_ii,is_con)
!        else
!          call get_new_diag_elem(connected_diag_elems_info(j),connected_dets_up(j),connected_dets_dn(j),H_ii)
!        endif
!
!        ! Update contribution from determinant j
!        call hamiltonian(connected_dets_up(j),connected_dets_dn(j),dets_up(i),dets_dn(i),H_ij,is_con)
!        delta_e_2pt = delta_e_2pt + (-(old_sum**2) + (old_sum + H_ij*wts(i))**2) / (var_energy - H_ii)
!
!      enddo ! n_connected_dets
!
!      ! Add reference i to neighbors
!      call add_det_to_neighbor_tables(dets_up(i),dets_dn(i),wts(i),neighbors) !,orb_group_mat)
!
!    enddo ! ndets
!
!    write (6,*) "n_actual=",n_actual,"n_potential=",n_potential,"success rate=",real(n_actual)/real(n_potential)
!
!    write (6,'(''Variational energy='',t36,f15.9)') var_energy
!    write (6,'(''2nd-order PT energy lowering='',t36,f15.9)') delta_e_2pt
!    write (6,'(''Total energy='',t36,f15.9)') var_energy+delta_e_2pt
!    if(hamiltonian_type.eq.'heg') then
!      write(6,'(''Total energy (includ. Madelung)='',t36,f15.9)') var_energy+delta_e_2pt+energy_madelung
!      write(6,'(''Correlation energy ='',t36,f15.9)') var_energy+delta_e_2pt-energy_hf
!      write(6,'(''Correlation energy per electron ='',t36,f15.9)') (var_energy+delta_e_2pt-energy_hf)/nelec
!    endif
!    call flush(6)
!
!    call my_second(2,'multi_index_hashing 2nd-order PT correction')
!
!  end subroutine second_order_pt_multi_index_hashing



  subroutine energies_for_extrapolation(n_states,remote_det_map,local_det_map,n_energy_batch,eps_pt,eps_pt_big,target_error,n_mc,eps_var)
    use chemistry, only : time_sym,z
    use tools, only : merge_sort2_up_dn
    use mpi_routines, only : master_core_node, master_core, det_map, det_map_l, mpi_distribute_remote_det_map, shmem_allocate, shmem_deallocate, shmem_reallocate, mpi_barr_in_node, mpi_barr
    use semistoch, only : hamiltonian
    use generic_sort, only : sort,shell_sort_real_rank1_int_rank1
    integer,allocatable :: iorder(:),temp_i_2(:)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),allocatable :: new_dets_up(:), new_dets_dn(:)
    type(ik_vec),allocatable :: temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#else
    integer(ik),allocatable :: new_dets_up(:), new_dets_dn(:)
    integer(ik),allocatable :: temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
#endif
    integer, intent(in) :: n_states,n_energy_batch
    integer, intent(inout) :: n_mc
    real(rk), intent(in)   :: eps_pt, target_error, eps_var
    real(rk), intent(inout) :: eps_pt_big
    type(det_map), intent(inout) :: remote_det_map
    type(det_map_l), intent(in) :: local_det_map
    integer :: ndets_connected
    integer ndets_local, batch_size, ndets_symmetrized
    logical :: is_con, time_sym_old
    integer i,j,k,i_state,ndets_global_old,ndets_full
    real(rk),allocatable :: old_wts(:,:),starting_wts(:,:), energy(:)
    real(rk),allocatable :: diag_elems(:)
    real(rk), pointer, dimension(:) :: diag_elems_global(:)
    real(rk) :: pt_diff, pt_big, pt_diff_std_dev

    !MJO - det_map is used for the parallel run
    type(det_map) fake_remote_det_map
    type(det_map_l) fake_local_det_map


    ndets_full = local_det_map%ndets_global_old
    time_sym_old = time_sym
    ! Sort global list
    !MJO Use 2*ndets_full to allow for the expanded symmetrized wavefunction, used before PT
    if (.not.allocated(iorder))  allocate(iorder(2*ndets_full))
    if (.not.allocated(temp_i16_up))  allocate(temp_i16_up((ndets_full+1)))
    if (.not.allocated(temp_i16_dn))  allocate(temp_i16_dn((ndets_full+1)))
    if (.not.allocated(temp_i_2))  allocate(temp_i_2((ndets_full+1)))
    do j=1,ndets_full
      iorder(j)=j
    enddo

    write(6,'(/,''Sorting global dets list to create truncated wavefunctions'')') ; call flush(6)

    if (master_core_node) then
      call shell_sort_real_rank1_int_rank1(remote_det_map%remote_wts(1:ndets_full,1),iorder(1:ndets_full))
      remote_det_map%global_dets_up(1:ndets_full) = remote_det_map%global_dets_up(iorder(1:ndets_full))
      remote_det_map%global_dets_dn(1:ndets_full) = remote_det_map%global_dets_dn(iorder(1:ndets_full))
      if (n_states>1) then
        do i=2,n_states
          remote_det_map%remote_wts(1:ndets_full,i) = remote_det_map%remote_wts(iorder(1:ndets_full),i)
        enddo
      endif
    endif


    allocate(new_dets_up(local_det_map%ndets )) ! Allocate enough space for your whole sectionn
    allocate(new_dets_dn(local_det_map%ndets )) ! Allocate enough space for your whole section
    allocate(starting_wts(local_det_map%ndets,n_states )) ! Allocate enough space for your whole section
    allocate(old_wts(local_det_map%ndets,n_states )) ! Allocate enough space for your whole section
    allocate(energy(n_states))
    !MJO Use 3*ndets to allow for the expanded symmetrized wavefunction, used before PT
    !    Maybe, in parallel, you will end up with more than 2x determinants because of hashing
    allocate(diag_elems(3*local_det_map%ndets))
    batch_size = ndets_full/n_energy_batch
    call shmem_allocate(fake_remote_det_map%global_dets_up,int(batch_size,i8b))
    call shmem_allocate(fake_remote_det_map%global_dets_dn,int(batch_size,i8b))
    call shmem_allocate(fake_remote_det_map%remote_wts,int(batch_size,i8b),int(n_states,i8b))
    !MJO Use 2*ndets_full to allow for the expanded symmetrized wavefunction, used before PT
    call shmem_allocate(diag_elems_global,int(2*ndets_full,i8b))

    do i=1,n_energy_batch+1 !MJO +1 so that we always go up to ndets_full
      ndets_global_old = batch_size * i
      if (ndets_global_old.gt.ndets_full) then
        ndets_global_old = ndets_full
      endif
      call shmem_deallocate(fake_remote_det_map%remote_wts)
      call shmem_deallocate(fake_remote_det_map%global_dets_up)
      call shmem_deallocate(fake_remote_det_map%global_dets_dn)
      call shmem_allocate(fake_remote_det_map%remote_wts,int(ndets_global_old,i8b),int(n_states,i8b))
      call shmem_allocate(fake_remote_det_map%global_dets_up,int(ndets_global_old,i8b))
      call shmem_allocate(fake_remote_det_map%global_dets_dn,int(ndets_global_old,i8b))

      if (master_core_node) then
        fake_remote_det_map%global_dets_up(1:ndets_global_old) = remote_det_map%global_dets_up(1:ndets_global_old)
        fake_remote_det_map%global_dets_dn(1:ndets_global_old) = remote_det_map%global_dets_dn(1:ndets_global_old)
        fake_remote_det_map%remote_wts(1:ndets_global_old,1:n_states) = remote_det_map%remote_wts(1:ndets_global_old,1:n_states)
      endif

      fake_local_det_map%ndets_global_old = ndets_global_old

      call mpi_distribute_remote_det_map(fake_remote_det_map,fake_local_det_map)

      fake_local_det_map%ndets_global_old = batch_size * (i-1)
      fake_local_det_map%ndets_global = ndets_global_old
      deallocate(new_dets_up,new_dets_dn,starting_wts,old_wts)
      allocate(new_dets_up(fake_local_det_map%ndets )) ! Allocate enough space for your whole sectionn
      allocate(new_dets_dn(fake_local_det_map%ndets )) ! Allocate enough space for your whole section
      allocate(starting_wts(fake_local_det_map%ndets,n_states )) ! Allocate enough space for your whole section
      allocate(old_wts(fake_local_det_map%ndets,n_states )) ! Allocate enough space for your whole section

      new_dets_up(:) = fake_remote_det_map%global_dets_up(fake_local_det_map%local_indices)
      new_dets_dn(:) = fake_remote_det_map%global_dets_dn(fake_local_det_map%local_indices)

      do j=1,n_states
        old_wts(:,j)      = fake_remote_det_map%remote_wts(fake_local_det_map%local_indices,j)
        starting_wts(:,j) = fake_remote_det_map%remote_wts(fake_local_det_map%local_indices,j)
      enddo

      ndets_local = fake_local_det_map%ndets

      call iterative_diagonalize(ndets_local, n_states, new_dets_up, new_dets_dn, final_wts=old_wts, eigenvalues=energy, starting_wts=starting_wts, remote_det_map=fake_remote_det_map,local_det_map=fake_local_det_map)
      ndets_symmetrized = batch_size*i
      fake_local_det_map%ndets_global_old = batch_size * i
      if (fake_local_det_map%ndets_global_old.gt.ndets_full) then
        fake_local_det_map%ndets_global_old = ndets_full
      endif

      if (time_sym) then
        write (6,'(/,''Converting to determinant basis'')') ; call flush(6)
        if (use_mpi) then
          call convert_time_symmetrized_to_dets(ndets_local,new_dets_up,new_dets_dn,old_wts,n_states,fake_remote_det_map,fake_local_det_map)
          ndets_global_old = fake_local_det_map%ndets_global_old
        else
          call convert_time_symmetrized_to_dets(ndets_local,new_dets_up,new_dets_dn,old_wts,n_states)
          ndets_global_old = ndets_local
        endif
        write (6,'(''From now on, no more time-reversal symmetry'',/)') ; call flush(6)
        time_sym = .false.
      endif

      do j=1,ndets_global_old
        iorder(j)=j
      enddo
      ! Sort variational dets so that they can be binary searched
      if (master_core_node) then
        call merge_sort2_up_dn(fake_remote_det_map%global_dets_up(1:ndets_global_old),fake_remote_det_map%global_dets_dn(1:ndets_global_old), iorder, ndets_global_old, temp_i16_up, temp_i16_dn, temp_i_2)
        fake_remote_det_map%remote_wts(1:ndets_global_old,1:n_states) = fake_remote_det_map%remote_wts(iorder(1:ndets_global_old),1:n_states)
      endif
      call mpi_barr_in_node()
      do j=1,ndets_local
        iorder(j)=j
      enddo
      call merge_sort2_up_dn(new_dets_up(1:ndets_local),new_dets_dn(1:ndets_local), iorder(1:ndets_local), ndets_local, temp_i16_up, temp_i16_dn, temp_i_2)
      old_wts(1:ndets_local,1:n_states) = old_wts(iorder(1:ndets_local),1:n_states)

      if (use_mpi) then
        !MJO Quick hack - compute all global diagonal elements
        if (master_core_node) then
          do k=1,ndets_global_old
            call hamiltonian(fake_remote_det_map%global_dets_up(k),fake_remote_det_map%global_dets_dn(k),fake_remote_det_map%global_dets_up(k),fake_remote_det_map%global_dets_dn(k),diag_elems_global(k),is_con)
          enddo
        endif
        call mpi_barr_in_node()
      endif
      ! Quick hack: compute all local diag_elems now

      do i_state=1,ndets_local
        call hamiltonian(new_dets_up(i_state),new_dets_dn(i_state),new_dets_up(i_state),new_dets_dn(i_state),diag_elems(i_state),is_con)
      enddo

      do i_state=1,n_states
        pt_big = 0._rk
        pt_diff = 0._rk
        pt_diff_std_dev = 0._rk

        call do_pt(ndets_local,new_dets_up(1:ndets_local),new_dets_dn(1:ndets_local),old_wts(1:ndets_local,i_state),diag_elems(1:ndets_local),ndets_global_old,fake_remote_det_map%global_dets_up(1:ndets_global_old), fake_remote_det_map%global_dets_dn(1:ndets_global_old), fake_remote_det_map%remote_wts(1:ndets_global_old,i_state), diag_elems_global(1:ndets_global_old), energy(i_state),eps_var,eps_pt,eps_pt_big,target_error,n_mc,pt_big,pt_diff,pt_diff_std_dev,ndets_connected,fake_remote_det_map,fake_local_det_map)
        write (6,'(/,''State'',i4,'':'')') i_state
        write (6,'(''Variational energy('',i1,'')='',t36,f15.9)') i_state, energy(i_state)
        if (pt_diff_std_dev==0._rk) then
          write (6,'(''2nd-order PT energy lowering('',i1,'')='',t36,f15.9)') i_state, pt_big
          write (6,'(''Total energy('',i2,'')='',t36,f15.9)') i_state, energy(i_state)+pt_big
          write (6,'(''eps_var, eps_pt, ndets, ndets_connected(total), Variational, PT_actv, PT_full, Total Energies('',i1,'')='',2e9.2,i9,i11,f16.9,2f15.9,f16.9)') i_state, eps_var, eps_pt, ndets_symmetrized, ndets_connected, energy(i_state), pt_big, pt_big, energy(i_state)+pt_big
        else
          write (6,'(''2nd-order PT energy lowering('',i1,'')='',t36,f15.9,'' +-'',f12.9,'' ('',2f13.9,'')'')') i_state, pt_big+pt_diff, pt_diff_std_dev, pt_big, pt_diff
          write (6,'(''Total energy('',i1,'')='',t36,f15.9,'' +-'',f12.9)') i_state, energy(i_state)+pt_big+pt_diff, pt_diff_std_dev
          write (6,'(''eps_var, eps_pt, ndets, ndets_connected(total), Variational, PT_actv, PT_full, Total Energies('',i1,'')='',2e9.2,i9,i11,f16.9,2f15.9,f16.9,'' +-'',f12.9)') i_state, eps_var, eps_pt, ndets_symmetrized, ndets_connected, energy(i_state), pt_big+pt_diff, pt_big+pt_diff, energy(i_state)+pt_big+pt_diff, pt_diff_std_dev
        endif
      enddo
      write(6,*) 'ndets: ',ndets_global_old,'energy: ',energy(1); call flush(6)
      deallocate(fake_local_det_map%local_indices)
      !MJO Reset time_sym
      time_sym = time_sym_old
    enddo

    call shmem_deallocate(fake_remote_det_map%global_dets_up)
    call shmem_deallocate(fake_remote_det_map%global_dets_dn)
    call shmem_deallocate(fake_remote_det_map%remote_wts)
    stop
  end subroutine energies_for_extrapolation


!-------------------------------------------------------------------------------------

!  subroutine second_order_pt_ah(ndets,dets_up,dets_dn,wts,diag_elems,var_energy,eps_pt)
!    ! Compute and print the second order energy lowering
!    ! For each reference, generate all connections exceeding eps_pt.
!    ! Then, loop over all connections to see if any are connections to references already passed through; update E_2pt appropriately
!    ! When done generating connections to a reference dets, if the number of connections is >0, store all determinants
!    ! that can be constructed by removing the most-excited electron that was just excited.
!    ! When checking whether a connection is connected to any previous reference dets, loop over all of its electrons,
!    ! and look up each resulting determinant to find any potential connections.
!
!    ! A Holmes, 4 Apr 2016
!
!    use semistoch, only : find_doubly_excited,hamiltonian
!    use common_run, only : connected_dets_up,connected_dets_dn,connected_matrix_elements,max_connected_dets,diag_elem_info,connected_diag_elems_info
!    use chemistry, only : find_important_connected_dets_chem,find_important_connected_dets_1e_removed,find_important_singly_connected_dets_1e_removed,find_important_doubly_connected_dets_2e_removed,get_new_diag_elem
!    use types, only : i8b
!    use generic_sort, only : sort
!    use heg, only : energy_hf, energy_madelung
!
!    integer,intent(in) :: ndets
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
!#else
!    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
!#endif
!    real(rk),intent(in) :: wts(:),var_energy,eps_pt
!
!    real(rk) :: delta_e_2pt
!    integer :: i,j,k
!    real(rk),allocatable :: diag_elems(:)
!    type(neighbor_table) :: neighbors,ref_dets
!    real(rk) :: H_ii,H_ij
!    real(rk) :: old_sum
!    integer :: n_connected_dets,n_connections
!    real(rk),allocatable :: potential_connected_coeffs(:)
!    logical :: is_con
!    integer(i8b) :: n_potential, n_actual
!    integer :: i_ref
!    integer,allocatable :: tmp_ind_singles(:),tmp_info_singles(:)
!    integer,allocatable :: tmp_ind_doubles(:),tmp_info_doubles(:)
!   !integer,allocatable :: tmp_ind(:),tmp_info(:)
!    integer,allocatable :: ref_connection_indices(:)
!    integer :: n_allocate_singles, n_allocate_doubles
!    integer :: n_hash_singles, n_hash_doubles
!   !integer :: n_allocate
!   !integer :: n_hash ! number of dets placed in hash table
!    integer :: n_sum_terms
!
!    call my_second(1,'AH 2nd-order PT correction')
!
!    call init_hash_fn ! Sets up hash function
!    n_potential = 0_i8b
!    n_actual = 0_i8b
!
!    ! First, create hash tables
!    n_allocate_singles = ndets*nelec
!    n_allocate_doubles = ndets*nelec ! will probably need to be reallocated
!    allocate(hash_ind_singles(n_allocate_singles))
!    allocate(hash_info_singles(n_allocate_singles))
!    allocate(hash_ind_doubles(n_allocate_doubles))
!    allocate(hash_info_doubles(n_allocate_doubles))
!    n_hash_singles = 0
!    n_hash_doubles = 0
!
!    ! Single excitations
!    do i_ref=1,ndets
!      call find_important_singly_connected_dets_1e_removed(dets_up(i_ref), dets_dn(i_ref), n_connected_dets, connected_dets_up, connected_dets_dn)
!      ! Now, connected_dets represents all those 1-electron-removed singly excited dets connected to reference
!      ! These are the places in the hash table where the index i_ref needs to be stored.
!      do j=1,n_connected_dets
!        n_hash_singles = n_hash_singles + 1
!        hash_ind_singles(n_hash_singles) = hash_det(connected_dets_up(j),connected_dets_dn(j)) ! index of hash table
!       !hash_info(n_hash) = i_ref ! what is stored in hash table
!      enddo ! j
!    enddo ! i_ref
!
!    ! Double excitations
!    do i_ref=1,ndets
!      call find_important_doubly_connected_dets_2e_removed(dets_up(i_ref), dets_dn(i_ref), eps_pt/abs(wts(i_ref)), n_connected_dets, connected_dets_up, connected_dets_dn)
!      ! Now, connected_dets represents all those 2-electron-removed doubly excited dets connected to reference by more than eps_pt
!      ! These are the places in the hash table where the index i_ref needs to be stored.
!      if (n_hash_doubles+n_connected_dets>n_allocate_doubles) then ! reallocate
!        write (6,*) "Reallocating hash tables from size",n_allocate_doubles,"to size", 2*n_allocate_doubles; call flush(6)
!        allocate(tmp_ind_doubles(n_allocate_doubles))
!        allocate(tmp_info_doubles(n_allocate_doubles))
!        tmp_ind_doubles(:) = hash_ind_doubles(:)
!        tmp_info_doubles(:) = hash_info_doubles(:)
!        deallocate(hash_ind_doubles,hash_info_doubles)
!        n_allocate_doubles = 2*n_allocate_doubles
!        allocate(hash_ind_doubles(n_allocate_doubles))
!        allocate(hash_info_doubles(n_allocate_doubles))
!        hash_ind_doubles(1:size(tmp_ind_doubles)) = tmp_ind_doubles(:)
!        hash_info_doubles(1:size(tmp_info_doubles)) = tmp_info_doubles(:)
!        deallocate(tmp_ind_doubles,tmp_info_doubles)
!      endif
!      do j=1,n_connected_dets
!        n_hash_doubles = n_hash_doubles + 1
!        hash_ind_doubles(n_hash_doubles) = hash_det(connected_dets_up(j),connected_dets_dn(j)) ! index of hash table
!       !hash_info(n_hash) = i_ref ! what is stored in hash table
!      enddo ! j
!    enddo ! i_ref
!
!    ! Now, sort the hash tables by hash_ind
!    call sort(n_hash_singles,hash_ind_singles)
!    call sort(n_hash_doubles,hash_ind_doubles)
!
!    ! Now, store the start_index and end_index for each hash_ind
!    allocate(start_index_singles(n_bins))
!    start_index_singles(:) = 0
!    start_index_singles(hash_ind_singles(1)) = 1
!    do i=2,n_hash_singles
!      if (hash_ind_singles(i).ne.hash_ind_singles(i-1)) then
!        start_index_singles(hash_ind_singles(i)) = i
!      endif
!    enddo
!    deallocate(hash_ind_singles)
!
!    allocate(start_index_doubles(n_bins))
!    start_index_doubles(:) = 0
!    start_index_doubles(hash_ind_doubles(1)) = 1
!    do i=2,n_hash_doubles
!      if (hash_ind_doubles(i).ne.hash_ind_doubles(i-1)) then
!        start_index_doubles(hash_ind_doubles(i)) = i
!      endif
!    enddo
!    deallocate(hash_ind_doubles)
!
!    ! At this point, hash_info(start_index(i):end_index(i)), where i=hash_det(det_up,det_dn), are the
!    ! reference dets connected to det_up/dn
!
!    allocate(n_per_bin_singles(n_bins))
!    n_per_bin_singles(:) = 0
!    hash_info_singles(:) = 0 ! Because we only want to check whether a connection is connected to a reference that has already been looped over!
!
!    allocate(n_per_bin_doubles(n_bins))
!    n_per_bin_doubles(:) = 0
!    hash_info_doubles(:) = 0 ! Because we only want to check whether a connection is connected to a reference that has already been looped over!
!
!    ! Also, create a table of just the ref_dets, to enable us to quickly check whether a connected det should be thrown out
!    ! (since the PT expression sums over only dets not in reference)
!    call create_neighbor_tables(ndets,dets_up,dets_dn,1,n_bins,ref_dets)
!
!    ! This table should be filled up with all references!
!    call fill_neighbor_tables(ndets,dets_up,dets_dn,wts,ref_dets)
!
!    ! Next, loop over all references:
!    ! For each, compute diagonal elements, then generate all connections exceeding eps_pt
!    ! For each connection, loop over all of the electrons; hash each determinant obtained by removing one electron
!    ! This yields the set of reference determinants connected to the query
!    ! After all connections for a given reference det have been looped over,
!    ! add that reference det to neighbors
!
!    allocate(ref_connection_indices(ndets))
!
!    delta_e_2pt = 0._rk
!
!    max_sum_terms = 0
!
!    do i=1,ndets
!
!      call find_important_connected_dets_chem(dets_up(i), dets_dn(i), eps_pt/abs(wts(i)), n_connected_dets, connected_dets_up, connected_dets_dn,connected_matrix_elements,diag_elems(i),connected_diag_elems_info) ! connected_diag_elems_info is needed for computing diagonal elements in O(N) time, but don't use more than once per connected det!
!
!      if (n_connected_dets==0)  cycle
!
!      do j=2,n_connected_dets ! skip 1 because it is the diagonal element
!
!        ! Check by a single hash table whether it is a reference det
!        if (is_in_reference(connected_dets_up(j),connected_dets_dn(j),ref_dets))  cycle
!
!        call get_reference_connections_less_memory(ndets,connected_dets_up(j),connected_dets_dn(j),n_connections,ref_connection_indices)
!       !call get_reference_connections(connected_dets_up(j),connected_dets_dn(j),n_connections,ref_connection_indices)
!
!        n_potential = n_potential + n_connections
!
!        old_sum = 0._rk
!        n_sum_terms = 0
!
!        if (n_connections>0) then
!          do k=1,n_connections
!            ! Following check has to be in because of hash collisions
!            if (.not.is_connected(connected_dets_up(j),connected_dets_dn(j),dets_up(ref_connection_indices(k)),dets_dn(ref_connection_indices(k))))  cycle
!            n_actual = n_actual + 1_i8b
!            call hamiltonian(connected_dets_up(j),connected_dets_dn(j),dets_up(ref_connection_indices(k)),dets_dn(ref_connection_indices(k)),H_ij,is_con)
!            if (abs(H_ij*wts(ref_connection_indices(k)))>eps_pt) then
!              old_sum = old_sum + H_ij*wts(ref_connection_indices(k))
!              n_sum_terms = n_sum_terms + 1
!            endif
!          enddo ! reference connections (in reference det list)
!        endif
!
!        if (n_sum_terms>max_sum_terms)  max_sum_terms = n_sum_terms
!
!        ! get H_ii
!        if (connected_diag_elems_info(j)%old_diag_elem>1.e50_rk) then ! single excitation; compute the easy way for now
!          call hamiltonian(connected_dets_up(j),connected_dets_dn(j),connected_dets_up(j),connected_dets_dn(j),H_ii,is_con)
!        else
!          call get_new_diag_elem(connected_diag_elems_info(j),connected_dets_up(j),connected_dets_dn(j),H_ii)
!        endif
!
!        ! Update contribution from determinant j
!        call hamiltonian(connected_dets_up(j),connected_dets_dn(j),dets_up(i),dets_dn(i),H_ij,is_con)
!        delta_e_2pt = delta_e_2pt + (-(old_sum**2) + (old_sum + H_ij*wts(i))**2) / (var_energy - H_ii)
!
!      enddo ! j=2,n_connected_dets
!
!      ! Add reference i to hash table
!      call add_reference_to_hash_table(i,dets_up(i),dets_dn(i),wts(i),eps_pt)
!
!    enddo ! i=1,ndets
!
!    write (6,*) "Max terms in one sum=",max_sum_terms
!
!    write (6,*) "n_actual=",n_actual,"n_potential=",n_potential,"success rate=",real(n_actual)/real(n_potential)
!
!    write (6,'(''Variational energy='',t36,f15.9)') var_energy
!    write (6,'(''2nd-order PT energy lowering='',t36,f15.9)') delta_e_2pt
!    write (6,'(''Total energy='',t36,f15.9)') var_energy+delta_e_2pt
!    if(hamiltonian_type.eq.'heg') then
!      write(6,'(''Total energy (includ. Madelung)='',t36,f15.9)') var_energy+delta_e_2pt+energy_madelung
!      write(6,'(''Correlation energy ='',t36,f15.9)') var_energy+delta_e_2pt-energy_hf
!      write(6,'(''Correlation energy per electron ='',t36,f15.9)') (var_energy+delta_e_2pt-energy_hf)/nelec
!    endif
!    call flush(6)
!
!    call my_second(2,'AH 2nd-order PT correction')
!
!  end subroutine second_order_pt_ah
!-------------------------------------------------------------------------------------

!  subroutine create_neighbor_tables(n_det,dets_up,dets_dn,n_hash,n_bins,neighbors,orb_group_mat)
!  ! For multi-index hashing (not currently being used)
!  ! Create a structure to be used for approximate string matching.
!  ! It should be empty for computing 2PT energy, so
!  ! when used for constructing H quickly, call a separate
!  ! subroutine to fill it.
!  ! orb_group_mat(i,j) is the j'th orbital in the i'th group; if not
!  ! present, then just divide up the orbitals into equal-sized groups
!  ! by putting every 5th one in a different group
!  ! A Holmes, 22 Mar 2016
!
!    integer,intent(in) :: n_det
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
!#else
!    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
!#endif
!    integer,intent(in) :: n_hash,n_bins ! want n_hash = 5, n_bins ~ n_det
!    type(neighbor_table),intent(out) :: neighbors
!    integer,optional,intent(in) :: orb_group_mat(:,:)
!
!    integer,allocatable :: bin_counts(:,:)
!    integer :: i,j
!    integer,allocatable :: hashes(:)
!
!    neighbors%n_hash_functions = n_hash
!    neighbors%n_bins_per_hash = n_bins
!
!    allocate(hashes(n_hash))
!
!    allocate(neighbors%is_empty(n_hash,n_bins))
!    allocate(neighbors%bin_counts(n_hash,n_bins))
!    allocate(neighbors%start_index(n_hash,n_bins))
!    allocate(neighbors%end_index(n_hash,n_bins))
!    allocate(neighbors%dets_up(n_hash,n_det))
!    allocate(neighbors%dets_dn(n_hash,n_det))
!    allocate(neighbors%coeffs(n_hash,n_det))
!
!    ! First, count number of dets that hash into each bin
!    allocate(bin_counts(n_hash,n_bins))
!    bin_counts(:,:) = 0
!    if (present(orb_group_mat)) then
!      do i=1,n_det
!        call get_hashes(dets_up(i),dets_dn(i),n_hash,n_bins,hashes,orb_group_mat)
!        do j=1,n_hash
!          bin_counts(j,hashes(j)) = bin_counts(j,hashes(j)) + 1
!        enddo
!      enddo
!    else
!      do i=1,n_det
!        call get_hashes(dets_up(i),dets_dn(i),n_hash,n_bins,hashes)
!        do j=1,n_hash
!          bin_counts(j,hashes(j)) = bin_counts(j,hashes(j)) + 1
!        enddo
!      enddo
!    endif
!
!    ! Then, "allocate" appropriately, i.e., determine which indices correspond to which bins of hash table
!    do j=1,n_hash
!      neighbors%is_empty(j,1) = (bin_counts(j,1)==0)
!      neighbors%start_index(j,1) = 1
!      neighbors%end_index(j,1) = bin_counts(j,1)
!      do i=2,n_bins
!        neighbors%is_empty(j,i) = (bin_counts(j,i)==0)
!        neighbors%start_index(j,i) = neighbors%end_index(j,i-1) + 1
!        neighbors%end_index(j,i) = neighbors%end_index(j,i-1) + bin_counts(j,i)
!      enddo
!    enddo
!
!    neighbors%is_empty(:,:) = .true. ! Because the neighbors table is currently empty! (call fill_neighbor_tables afterwards if needed)
!    neighbors%bin_counts(:,:) = 0
!
!  end subroutine create_neighbor_tables
!!-------------------------------------------------------------------------------------
!
!  subroutine fill_neighbor_tables(n_det,dets_up,dets_dn,coeffs,neighbors,orb_group_mat)
!  ! For multi-index hashing (not currently being used)
!  ! Fill neighbor tables with dets_up/dn
!  ! Assumes neighbor tables have already been created with create_neighbor_tables
!  ! A Holmes, 22 Mar 2016
!
!    integer,intent(in) :: n_det
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
!#else
!    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
!#endif
!    real(rk),intent(in) :: coeffs(:)
!    type(neighbor_table),intent(inout) :: neighbors
!    integer,optional,intent(in) :: orb_group_mat(:,:)
!
!    integer,allocatable :: bin_counts(:,:)
!    integer :: i,j,ind
!    integer :: n_hash
!    integer :: n_bins
!    integer,allocatable :: hashes(:)
!
!    n_hash = neighbors%n_hash_functions
!    n_bins = neighbors%n_bins_per_hash
!
!    allocate(hashes(n_hash))
!
!    allocate(bin_counts(n_hash,n_bins))
!    bin_counts(:,:) = 0
!    if (present(orb_group_mat)) then
!      do i=1,n_det
!        call get_hashes(dets_up(i),dets_dn(i),n_hash,n_bins,hashes,orb_group_mat)
!        do j=1,n_hash
!          ind = neighbors%start_index(j,hashes(j)) + bin_counts(j,hashes(j))
!          neighbors%dets_up(j,ind) = dets_up(i)
!          neighbors%dets_dn(j,ind) = dets_dn(i)
!          neighbors%coeffs(j,ind) = coeffs(i)
!          neighbors%is_empty(j,hashes(j)) = .false.
!          bin_counts(j,hashes(j)) = bin_counts(j,hashes(j)) + 1
!        enddo
!      enddo
!    else
!      do i=1,n_det
!        call get_hashes(dets_up(i),dets_dn(i),n_hash,n_bins,hashes)
!        do j=1,n_hash
!          ind = neighbors%start_index(j,hashes(j)) + bin_counts(j,hashes(j))
!          neighbors%dets_up(j,ind) = dets_up(i)
!          neighbors%dets_dn(j,ind) = dets_dn(i)
!          neighbors%coeffs(j,ind) = coeffs(i)
!          neighbors%is_empty(j,hashes(j)) = .false.
!          bin_counts(j,hashes(j)) = bin_counts(j,hashes(j)) + 1
!        enddo
!      enddo
!    endif
!
!    neighbors%bin_counts = bin_counts
!
!  end subroutine fill_neighbor_tables
!!-------------------------------------------------------------------------------------
!
!  subroutine add_det_to_neighbor_tables(det_up,det_dn,coeff,neighbors,orb_group_mat)
!  ! For multi-index hashing (not currently being used)
!  ! Add a single det to neighbor tables
!  ! Assumes this is one of the dets that was known about when
!  ! create_neighbor_tables was called, and that this det has
!  ! not yet been placed in the table!
!  ! A Holmes, 22 Mar 2016
!
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),intent(in) :: det_up,det_dn
!#else
!    integer(ik),intent(in) :: det_up,det_dn
!#endif
!    real(rk),intent(in) :: coeff
!    type(neighbor_table),intent(inout) :: neighbors
!    integer,optional,intent(in) :: orb_group_mat(:,:)
!
!    integer :: n_hash
!    integer :: n_bins
!    integer :: ind,j
!    integer,allocatable :: hashes(:)
!
!    n_hash = neighbors%n_hash_functions
!    n_bins = neighbors%n_bins_per_hash
!
!    allocate(hashes(n_hash))
!
!    if (present(orb_group_mat)) then
!      call get_hashes(det_up,det_dn,n_hash,n_bins,hashes,orb_group_mat)
!      do j=1,n_hash
!        ind = neighbors%start_index(j,hashes(j)) + neighbors%bin_counts(j,hashes(j))
!        neighbors%dets_up(j,ind) = det_up
!        neighbors%dets_dn(j,ind) = det_dn
!        neighbors%coeffs(j,ind) = coeff
!        neighbors%is_empty(j,hashes(j)) = .false.
!        neighbors%bin_counts(j,hashes(j)) = neighbors%bin_counts(j,hashes(j)) + 1
!      enddo
!    else
!      call get_hashes(det_up,det_dn,n_hash,n_bins,hashes)
!      do j=1,n_hash
!        ind = neighbors%start_index(j,hashes(j)) + neighbors%bin_counts(j,hashes(j))
!        neighbors%dets_up(j,ind) = det_up
!        neighbors%dets_dn(j,ind) = det_dn
!        neighbors%coeffs(j,ind) = coeff
!        neighbors%is_empty(j,hashes(j)) = .false.
!        neighbors%bin_counts(j,hashes(j)) = neighbors%bin_counts(j,hashes(j)) + 1
!      enddo
!    endif
!
!  end subroutine add_det_to_neighbor_tables
!-------------------------------------------------------------------------------------

!  subroutine get_hashes(det_up,det_dn,n_hash,n_bins,hashes,orb_group_mat)
!  ! Get hashes of all the subsets of the orbitals of det_up/dn.
!  ! If orb_group_mat not present, then assume the orbitals alternate their groups
!  ! A Holmes, 22 Mar 2016
!
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),intent(in) :: det_up,det_dn
!    type(ik_vec),allocatable :: pieces(:)
!#else
!    integer(ik),intent(in) :: det_up,det_dn
!    integer(ik),allocatable :: pieces(:)
!#endif
!    integer,intent(in) :: n_hash,n_bins
!    integer,intent(out) :: hashes(:)
!    integer,optional,intent(in) :: orb_group_mat(:,:)
!
!    integer :: i,k,bit
!
!    ! First, divide bits up into different groups
!    allocate(pieces(n_hash))
!    pieces(:) = 0_ik
!    bit = 0
!    if (present(orb_group_mat)) then
!      stop "orb_group_mat not implemented yet!"
!    else
!      do i=1,2*norb
!        k = mod(i-1,n_hash)+1
!        if (k==n_hash)  bit = bit + 1
!        if (i<=norb) then
!          if (btest(det_up,i-1))  pieces(k) = ibset(pieces(k),bit)
!        else
!          if (btest(det_dn,i-norb-1))  pieces(k) = ibset(pieces(k),bit)
!        endif
!      enddo
!    endif
!
!    ! Then, apply hash function to each group
!    do i=1,n_hash
!      call hash(pieces(i),n_bins,hashes(i),RandomHash1)
!    enddo
!
!  end subroutine get_hashes
!-------------------------------------------------------------------------------------

!  integer function hash_det(det_up,det_dn)
!  ! Hash the determinant det_up/dn, into a hash table with n_bins bins
!  ! A Holmes, 5 Apr 2016
!
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),intent(in) :: det_up,det_dn
!#else
!    integer(ik),intent(in) :: det_up,det_dn
!#endif
!
!    call hash_2_dets(det_up,det_dn,n_bins,hash_det,RandomHash1)
!
!  end function hash_det
!-------------------------------------------------------------------------------------

!  subroutine get_potential_connections(det_up,det_dn,neighbors,n_connections,connected_dets_up,connected_dets_dn,connected_coeffs,iorder,temp_i16_up,temp_i16_dn,temp_i_2,orb_group_mat)
!  ! For multi-index hashing (not currently being used)
!  ! Query the neighbor tables (which may not be completely filled in!) with a query
!  ! determinant, det_up/dn. Return a set of (potential) connected determinants and
!  ! their corresponding coefficients
!  ! Returns them sorted by det label, but there may be repeats
!  ! A Holmes, 22 Mar 2016
!
!    use tools, only : merge_sort2_up_dn
!
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),intent(in) :: det_up,det_dn
!    type(ik_vec),intent(out) :: connected_dets_up(:),connected_dets_dn(:)
!    type(ik_vec),intent(inout) ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
!#else
!    integer(ik),intent(in) :: det_up,det_dn
!    integer(ik),intent(out) :: connected_dets_up(:),connected_dets_dn(:)
!    integer(ik),intent(inout) ::temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
!#endif
!    type(neighbor_table),intent(in) :: neighbors
!    integer,intent(out) :: n_connections
!    real(rk),intent(out) :: connected_coeffs(:)
!    integer,intent(inout) :: iorder(:),temp_i_2(:)
!    integer,optional,intent(in) :: orb_group_mat(:,:)
!
!    integer :: i,ind
!    integer,allocatable :: hashes(:)
!    integer :: n_hash
!    integer :: n_bins
!
!    n_hash = neighbors%n_hash_functions
!    n_bins = neighbors%n_bins_per_hash
!
!    allocate(hashes(n_hash))
!
!    n_connections = 0
!
!    if (present(orb_group_mat)) then
!      call get_hashes(det_up,det_dn,n_hash,n_bins,hashes,orb_group_mat)
!    else
!      call get_hashes(det_up,det_dn,n_hash,n_bins,hashes)
!    endif
!
!    do i=1,n_hash
!      ind = hashes(i)
!      if (neighbors%is_empty(i,ind))  cycle
!      connected_dets_up(n_connections+1:n_connections+neighbors%bin_counts(i,ind)) = neighbors%dets_up(i,neighbors%start_index(i,ind):neighbors%start_index(i,ind)+neighbors%bin_counts(i,ind)-1)
!      connected_dets_dn(n_connections+1:n_connections+neighbors%bin_counts(i,ind)) = neighbors%dets_dn(i,neighbors%start_index(i,ind):neighbors%start_index(i,ind)+neighbors%bin_counts(i,ind)-1)
!      connected_coeffs(n_connections+1:n_connections+neighbors%bin_counts(i,ind)) = neighbors%coeffs(i,neighbors%start_index(i,ind):neighbors%start_index(i,ind)+neighbors%bin_counts(i,ind)-1)
!      n_connections = n_connections + neighbors%bin_counts(i,ind)
!    enddo
!
!    if (n_connections==0)  return
!
!    ! Finally, sort and merge to remove repeats!
!    do i=1,n_connections
!      iorder(i) = i
!    enddo
!    call merge_sort2_up_dn(connected_dets_up(1:n_connections),connected_dets_dn(1:n_connections), iorder, n_connections, temp_i16_up, temp_i16_dn, temp_i_2)
!    connected_coeffs(1:n_connections) = connected_coeffs(iorder(1:n_connections))
! ! MERGE! (TODO)
!
!  end subroutine get_potential_connections
!!-------------------------------------------------------------------------------------

!  subroutine get_reference_connections(det_up,det_dn,n_connections,ref_connection_indices)
!  ! Return the set of reference determinants connected to det_up/dn
!  ! All determinants will be either single or double excitations (unless there are hash collisions)
!  ! Time complexity is O(N_elec + N_connections*(1+N_hash_collisions))
!  ! A Holmes, 5 Apr 2016
!
!    use chemistry, only : occ_up,occ_dn
!    use tools, only : sort_and_merge
!
!    integer,intent(out) :: n_connections
!    integer,intent(out) :: ref_connection_indices(:)
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),intent(in) :: det_up,det_dn
!    type(ik_vec) :: tmp
!#else
!    integer(ik),intent(in) :: det_up,det_dn
!    integer(ik) :: tmp
!#endif
!
!    integer :: i_elec,ind,j
!
!    ! Get occ_up, occ_dn in O(N_elec) time
!    tmp = det_up
!    do i_elec=1,nup
!      occ_up(i_elec) = trailz(tmp) + 1
!      tmp = ibclr(tmp,occ_up(i_elec)-1)
!    enddo
!    tmp = det_dn
!    do i_elec=1,ndn
!      occ_dn(i_elec) = trailz(tmp) + 1
!      tmp = ibclr(tmp,occ_dn(i_elec)-1)
!    enddo
!
!    ! Loop over occupied electron lists; for each, remove that electron and hash to look up connections
!    n_connections = 0
!    do i_elec=1,nup
!      ind = hash_det(ibclr(det_up,occ_up(i_elec)-1),det_dn)
!      if (n_per_bin(ind)==0)  cycle
!      do j=1,n_per_bin(ind)
!        n_connections = n_connections + 1
!        ref_connection_indices(n_connections) = hash_info(start_index(ind)+j-1)
!      enddo ! j
!    enddo ! i_elec
!    do i_elec=1,ndn
!      ind = hash_det(det_up,ibclr(det_dn,occ_dn(i_elec)-1))
!      if (n_per_bin(ind)==0)  cycle
!      do j=1,n_per_bin(ind)
!        n_connections = n_connections + 1
!        ref_connection_indices(n_connections) = hash_info(start_index(ind)+j-1)
!      enddo ! j
!    enddo ! i_elec
!
!    ! Finally, we must sort and merge because there can be up to
!    ! N_elec repeats of each single excitation and up to
!    ! 2 repeats of each double excitation!
!
!    call sort_and_merge(n_connections,ref_connection_indices)
!
!  end subroutine get_reference_connections
!-------------------------------------------------------------------------------------

!  subroutine get_reference_connections_less_memory(n_ref,det_up,det_dn,n_connections,ref_connection_indices)
!  ! For multi-index hashing (not currently being used)
!  ! Works like get_reference_connections, but uses a factor of N_orb less memory
!  ! In worst case, takes a factor of N_elec more time, but since it avoids the problem of having up to
!  ! N_elec duplicates that get_reference_connections has, it should not be too much slower in practice
!  ! It also avoids the need to sort and merge because there are no duplicates
!  ! A Holmes, 6 Apr 2016
!
!    use chemistry, only : occ_up,occ_dn
!
!    integer,intent(in) :: n_ref
!    integer,intent(out) :: n_connections
!    integer,intent(out) :: ref_connection_indices(:)
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),intent(in) :: det_up,det_dn
!    type(ik_vec) :: tmp
!#else
!    integer(ik),intent(in) :: det_up,det_dn
!    integer(ik) :: tmp
!#endif
!
!    integer :: i_elec
!    integer :: n_single, n_double
!
!    ! Get occ_up, occ_dn in O(N_elec) time
!    tmp = det_up
!    do i_elec=1,nup
!      occ_up(i_elec) = trailz(tmp) + 1
!      tmp = ibclr(tmp,occ_up(i_elec)-1)
!    enddo
!    tmp = det_dn
!    do i_elec=1,ndn
!      occ_dn(i_elec) = trailz(tmp) + 1
!      tmp = ibclr(tmp,occ_dn(i_elec)-1)
!    enddo
!
!    call get_singly_connected_ref_dets(det_up,det_dn,n_single,ref_connection_indices,occ_up,occ_dn,hash_info_singles,n_per_bin_singles,start_index_singles)
!    call get_doubly_connected_ref_dets(det_up,det_dn,n_double,ref_connection_indices(n_single+1:n_ref),occ_up,occ_dn,hash_info_doubles,n_per_bin_doubles,start_index_doubles)
!   !call get_singly_connected_ref_dets(det_up,det_dn,n_single,ref_connection_indices,occ_up,occ_dn)
!   !call get_doubly_connected_ref_dets(det_up,det_dn,n_double,ref_connection_indices(n_single+1:n_ref),occ_up,occ_dn)
!
!    n_connections = n_single + n_double
!
!  end subroutine get_reference_connections_less_memory
!!-------------------------------------------------------------------------------------

!  subroutine get_singly_connected_ref_dets(det_up,det_dn,n_connections,ref_connection_indices,occ_up,occ_dn,hash_info_singles,n_per_bin_singles,start_index_singles)
!  ! For partial connections (look up single excitations in a list)
!  ! Return the set of singly excited reference determinants connected to det_up/dn
!  ! All determinants will be single excitations (unless there are hash collisions)
!  ! There should be no repeats, so no need to sort and merge
!  ! Time complexity is O(N_elec + N_connections*(1+N_hash_collisions))
!  ! A Holmes, 6 Apr 2016
!
!    integer,intent(out) :: n_connections
!    integer,intent(out) :: ref_connection_indices(:)
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),intent(in) :: det_up,det_dn
!    type(ik_vec) :: tmp
!#else
!    integer(ik),intent(in) :: det_up,det_dn
!    integer(ik) :: tmp
!#endif
!    integer,intent(in) :: occ_up(:),occ_dn(:)
!    integer,intent(in) :: hash_info_singles(:),n_per_bin_singles(:),start_index_singles(:)
!
!    integer :: i_elec,ind,j
!
!    ! Loop over occupied electron lists; for each, remove that electron and hash to look up connections
!    n_connections = 0
!    do i_elec=1,nup
!      ind = hash_det(ibclr(det_up,occ_up(i_elec)-1),det_dn)
!      if (n_per_bin_singles(ind)==0)  cycle
!      do j=1,n_per_bin_singles(ind)
!        n_connections = n_connections + 1
!        ref_connection_indices(n_connections) = hash_info_singles(start_index_singles(ind)+j-1)
!      enddo ! j
!    enddo ! i_elec
!    do i_elec=1,ndn
!      ind = hash_det(det_up,ibclr(det_dn,occ_dn(i_elec)-1))
!      if (n_per_bin_singles(ind)==0)  cycle
!      do j=1,n_per_bin_singles(ind)
!        n_connections = n_connections + 1
!        ref_connection_indices(n_connections) = hash_info_singles(start_index_singles(ind)+j-1)
!      enddo ! j
!    enddo ! i_elec
!
!  end subroutine get_singly_connected_ref_dets
!!-------------------------------------------------------------------------------------
!
!  subroutine get_doubly_connected_ref_dets(det_up,det_dn,n_connections,ref_connection_indices,occ_up,occ_dn,hash_info_doubles,n_per_bin_doubles,start_index_doubles)
!  ! For partial connections (look up double excitations in a list)
!  ! Return the set of doubly excited reference determinants connected to det_up/dn
!  ! All determinants will be double excitations (unless there are hash collisions)
!  ! There should be no repeats, so no need to sort and merge
!  ! Time complexity is O(N_elec^2 + N_connections*(1+N_hash_collisions))
!  ! Warning: also returns N_elec copies of each singly connected ref det,
!  !          which must be filtered out after calling this routine!
!  ! A Holmes, 6 Apr 2016
!
!    use chemistry, only : pairs_e1,pairs_e2
!
!    integer,intent(out) :: n_connections
!    integer,intent(out) :: ref_connection_indices(:)
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),intent(in) :: det_up,det_dn
!    type(ik_vec) :: tmp, new_up, new_dn
!#else
!    integer(ik),intent(in) :: det_up,det_dn
!    integer(ik) :: tmp,new_up,new_dn
!#endif
!    integer,intent(in) :: occ_up(:),occ_dn(:)
!    integer,intent(in) :: hash_info_doubles(:),n_per_bin_doubles(:),start_index_doubles(:)
!
!    integer :: i_elec,ind,j
!    integer :: ipair,npairs
!    integer :: p,q
!
!    npairs = 0
!    ! Same spin, up
!    do p=1,nup-1
!      do q=p+1,nup
!        npairs = npairs + 1
!        pairs_e1(npairs) = occ_up(p)
!        pairs_e2(npairs) = occ_up(q)
!      enddo
!    enddo
!    ! Same spin, dn
!    do p=1,ndn-1
!      do q=p+1,ndn
!        npairs = npairs + 1
!        pairs_e1(npairs) = occ_dn(p)+norb
!        pairs_e2(npairs) = occ_dn(q)+norb
!      enddo
!    enddo
!    ! Opposite spin
!    do p=1,nup
!      do q=1,ndn
!        npairs = npairs + 1
!        pairs_e1(npairs) = occ_up(p)
!        pairs_e2(npairs) = occ_dn(q)+norb
!      enddo
!    enddo
!
!    ! Loop over occupied electron pairs; for each, remove that electron and hash to look up connections
!
!    n_connections = 0
!    do ipair=1,npairs
!
!      p = pairs_e1(ipair)
!      q = pairs_e2(ipair)
!
!      ! Commenting out for now, since I think we want to retrieve all partial connections that have been stored:
!     !epair = int(combine_2_indices(p,q))
!     !if (dtm_hb(epair,1)%absH<eps_var)  exit ! If no double excitations that start from p,q are larger than eps_var, then skip this p,q pair
!
!      ! Generate new determinants by removing p and q
!      new_up = det_up
!      new_dn = det_dn
!      if (p<=norb) then
!        new_up = ibclr(new_up,p-1)
!      else
!        new_dn = ibclr(new_dn,p-norb-1)
!      endif
!      if (q<=norb) then
!        new_up = ibclr(new_up,q-1)
!      else
!        new_dn = ibclr(new_dn,q-norb-1)
!      endif
!
!      ind = hash_det(new_up,new_dn)
!
!      if (n_per_bin_doubles(ind)==0)  cycle
!
!      do j=1,n_per_bin_doubles(ind)
!        n_connections = n_connections + 1
!        ref_connection_indices(n_connections) = hash_info_doubles(start_index_doubles(ind)+j-1)
!      enddo ! j
!
!    enddo ! ipair (p,q)
!
!  end subroutine get_doubly_connected_ref_dets
!!-------------------------------------------------------------------------------------

  logical function is_connected(up1,dn1,up2,dn2)
  ! Returns true if the excitation level between dets 1 and 2 is <=2 and >0
  ! false otherwise
  ! A Holmes, 22 Mar 2016

    use chemistry, only : excitation_level

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: up1,dn1,up2,dn2
#else
    integer(ik),intent(in) :: up1,dn1,up2,dn2
#endif

    integer :: excite_level

    call excitation_level(up1,dn1,up2,dn2,excite_level)
    is_connected = (excite_level>0)

  end function is_connected
!-------------------------------------------------------------------------------------

  subroutine init_hash_fn
    ! These are random numbers:
    RandomHash1 = (/4735613, 313959, 650067,2918183,2123508,1824383, 368641,4665125, 889838,1653649, 567626, 700470,2072070,3265229,1218533,1765295,1278102,2587196,3192811,4079255,3370036,5073942,1910990,2852420,2908100,1267516,1016462,3757537,4184599, 304348,2123593,4425848,3581496,1791107,3497983,1391199,2198178,2947481,3044171,4776221,2298252,3308124,4255082, 756103,3464666,3571942,3888499,1721390,4148643, 221259,1768488,1669008,1336151, 899703,5072037,2064427,2291053,2385335,3980263,3427668,4590425,3130633,1795017,1298250, 437459,1125424,1810093,5031931,1982167,4491159,1664669,3584140,1229843,1109613,2897092,2456647,2870772,5057078,2580274, 246617,3489605, 743439,2495387,3268699, 945352,2387627,4946247,1554925,1258363,1796338,3509947,3914421,2572562, 907556, 455356,2301574,4622378,2526352,4954820, 850178, 535417, 286038,3498667,3241142,2524137,2937664,2249144,4947862,3411414,4485867,4610744,  17755,4813498,4686397,4337008, 178525,1955778,2093500,3448090,4787430,3377083,  43503,4225437, 761386,5010808,3790117,1625775,5031659,2308473,4029366,1280477,4861902,3954799,2985978,2901621,3061381,1543058, 653181,4504563,2731597,2634156,2414947,1399327, 425706,2280240,2920945,1939056,3003820,2989859,2873210,5001290,3667487,1452012,2119338,4194025,1593590,3594528,3385764,4118903, 308089, 170759,1138565, 903090,4696165,1358916,2896407,1462530,3696410,2300452,3165629,2134285,3437968,1946733,2601552,1185301,3922477,2544259, 260349,2558668,3904972,4254275,3685299,1264351,4273640,2627883,2919797,4986566,1057643,2066436,2962802, 610070,3901649,1381809,2342696,3827048, 138016,1993221,3701099, 658366, 871615,1086936,2323667,1609455,3283476,3220993,4161146,1102928,4279104,1346917,5038082,4026472,1083272,5080232,1091946,4863898,1349339,3686856,2421638, 616475, 147825,  36701,4707242,2266338,3872717, 651830, 404626, 498554,3907313, 855440,2487375,2716318,4383594,1890043,4229224, 506837,1815062,2767149,1773485, 614411,2654510,4654238,5059003,  41437,3902708,1384039,1358399,4635825,3555781,3241332,4683214,1461543,3618113,4767352,5028089,2733360,2020587/)
  end subroutine init_hash_fn
!-------------------------------------------------------------------------------------

  subroutine hash(det,range,hashx,RandomHash)
  ! hash function borrowed from mpi_routines, which I believe is
  ! a Pearson hash
  ! returns hashx, which is between 1 and range

    implicit none
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det
#else
    integer(ik), intent(in) :: det
#endif
    integer, intent(in) :: range,RandomHash(0:255)
    integer, intent(out) :: hashx
    integer(ik) :: acc
    integer(ik) :: test_nk
    integer :: test_int,i,j
    integer :: val
   !integer, parameter :: n_bits_nk = bit_size(test_nk)
   !integer, parameter :: n_bits = bit_size(test_int)

    acc = 0
#ifdef NUM_ORBITALS_GT_127
    do i=1,num_words    !run over integers defining the bit-string
        do j = 0, norb -1, 8
       !do j = 0, n_bits_nk -1, 8
            val = int(iand(ishft(det%v(i)*368296850_ik,-j), int(255,ik)),kind(test_int)) ! This change from Matt Otten has better load balancing
            !val is now a number between 0 -> 255
            !Ensure that RandomHash has a mapping for the 0th element too
            !1099511628211_ik = ibset(0_ik,27)
            acc = (1099511628211_ik * acc) + (RandomHash(val) * (i + j))
        enddo
    enddo
#else
    do j = 0, norb -1, 8
   !do j = 0, n_bits_nk -1, 8
        val = int(iand(ishft(det*368296850_ik,-j), int(255,ik)),kind(test_int)) ! This change from Matt Otten has better load balancing
        !val is now a number between 0 -> 255
        !Ensure that RandomHash has a mapping for the 0th element too
        !1099511628211_ik = ibset(0_ik,27)
        acc = (1099511628211_ik * acc) + (RandomHash(val) * (i + j))
    enddo
#endif
    hashx = int(abs(mod(acc,int(range,ik))),kind(test_int))+1

  end subroutine hash
!-------------------------------------------------------------------------------------

  subroutine hash_2_dets(det1,det2,range,hashx,RandomHash)
  ! Borrowed from mpi_routines.f90

    implicit none
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: det1,det2
#else
    integer(ik), intent(in) :: det1,det2
#endif
    integer, intent(in) :: range,RandomHash(0:255)
    integer, intent(out) :: hashx
    integer(ik) :: acc
    integer(ik) :: test_nk
    integer :: test_int,i,j
    integer :: val
   !integer, parameter :: n_bits_nk = bit_size(test_nk)

    acc = 0
    i = 1
    do j = 0, norb -1, 8
   !do j = 0, n_bits_nk -1, 8
        val = int(iand(ishft(368296850_ik*det1,-j), int(255,ik)),kind(test_int)) ! This change from Matt Otten has better load balancing
        !val is now a number between 0 -> 255
        !Ensure that RandomHash has a mapping for the 0th element too
        !1099511628211_ik = ibset(0_ik,27)
        acc = (1099511628211_ik * acc) + (RandomHash(val) * (i + j))
    enddo
    i = 2
    do j = 0, norb -1, 8
   !do j = 0, n_bits_nk -1, 8
        val = int(iand(ishft(368296850_ik*det2,-j), int(255,ik)),kind(test_int)) ! This change from Matt Otten has better load balancing
        !val is now a number between 0 -> 255
        !Ensure that RandomHash has a mapping for the 0th element too
        !1099511628211_ik = ibset(0_ik,27)
        acc = (1099511628211_ik * acc) + (RandomHash(val) * (i + j))
    enddo
    hashx = int(abs(mod(acc,int(range,ik))),kind(test_int)) + 1

  end subroutine hash_2_dets
!-------------------------------------------------------------------------------------

!  logical function is_in_reference(det_up,det_dn,ref_dets)
!  ! Check whether det is in ref_dets in O(1+N_hash_collisions) time
!  ! A Holmes, 23 Mar 2016
!
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),intent(in) :: det_up,det_dn
!#else
!    integer(ik), intent(in) :: det_up,det_dn
!#endif
!    type(neighbor_table),intent(in) :: ref_dets
!
!    integer :: n_bins
!    integer :: hashes(1)
!    integer :: i,ind
!
!    n_bins = ref_dets%n_bins_per_hash
!
!    call get_hashes(det_up,det_dn,1,n_bins,hashes)
!
!    ind = hashes(1)
!    if (ref_dets%is_empty(1,ind)) then
!      is_in_reference = .false.
!      return
!    endif
!
!    do i=1,ref_dets%bin_counts(1,ind) ! Only iterates over more than 1 element if there are hash collisions!
!      if (det_up==ref_dets%dets_up(1,ref_dets%start_index(1,ind)+i-1).and.det_dn==ref_dets%dets_dn(1,ref_dets%start_index(1,ind)+i-1)) then
!        is_in_reference = .true.
!        return
!      endif
!    enddo
!
!    is_in_reference = .false.
!
!  end function is_in_reference
!-------------------------------------------------------------------------------------

!  subroutine add_reference_to_hash_table(i_ref,det_up,det_dn,wt,eps_pt)
!  ! Add the reference det_up/dn to the hash table
!  ! A Holmes, 5 Apr 2016
!
!    use common_run, only : connected_dets_up, connected_dets_dn
!    use chemistry, only : find_important_connected_dets_1e_removed,find_important_singly_connected_dets_1e_removed,find_important_doubly_connected_dets_2e_removed
!
!    integer,intent(in) :: i_ref ! because we only want to store the index of the reference, not the actual configuration
!    real(rk),intent(in) :: wt,eps_pt
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),intent(in) :: det_up,det_dn
!    type(ik_vec) :: tmp
!#else
!    integer(ik), intent(in) :: det_up,det_dn
!    integer(ik) :: tmp
!#endif
!
!    integer :: ind,i_elec,j
!    integer :: n_connected_dets
!
!      call find_important_singly_connected_dets_1e_removed(det_up, det_dn, n_connected_dets, connected_dets_up, connected_dets_dn)
!
!      do j=1,n_connected_dets
!        ind = hash_det(connected_dets_up(j),connected_dets_dn(j))
!        hash_info_singles(start_index_singles(ind)+n_per_bin_singles(ind)) = i_ref
!        n_per_bin_singles(ind) = n_per_bin_singles(ind) + 1
!      enddo ! j
!
!      call find_important_doubly_connected_dets_2e_removed(det_up, det_dn, eps_pt/abs(wt), n_connected_dets, connected_dets_up, connected_dets_dn)
!
!      do j=1,n_connected_dets
!        ind = hash_det(connected_dets_up(j),connected_dets_dn(j))
!        hash_info_doubles(start_index_doubles(ind)+n_per_bin_doubles(ind)) = i_ref
!        n_per_bin_doubles(ind) = n_per_bin_doubles(ind) + 1
!      enddo ! j
!
!     !call find_important_connected_dets_1e_removed(det_up, det_dn, eps_pt/abs(wt), n_connected_dets, connected_dets_up, connected_dets_dn)
!
!     !do j=1,n_connected_dets
!     !  ind = hash_det(connected_dets_up(j),connected_dets_dn(j))
!     !  hash_info(start_index(ind)+n_per_bin(ind)) = i_ref
!     !  n_per_bin(ind) = n_per_bin(ind) + 1
!     !enddo ! j
!
!  end subroutine add_reference_to_hash_table

!=====================================================================================================================
! subroutine matrix_lanczos_partial_connections(n,dets_up,dets_dn,lowest_eigenvector,lowest_eigenvalue,hash_info_singles,n_per_bin_singles,start_index_singles,hash_info_doubles,n_per_bin_doubles,start_index_doubles,highest_eigenvalue,second_lowest_eigenvalue,initial_vector)
! ! A Holmes, 10 Apr 2016. Similar to Hitesh's Lanczos routine, but for the case where instead of storing the
! !                        Hamiltonian, we only store the partial connections (i.e., all ways to remove 1 or 2
! !                        electrons) of all the reference determinants
!!=====================================================================================================================
!
! use chemistry, only : occ_up,occ_dn
! use common_ham, only : hamiltonian_type,nelec
! use common_run, only : max_connected_dets
! use semistoch, only : hamiltonian
!
! implicit none
!
! ! Dummy
! integer,intent(in)          :: n
!#ifdef NUM_ORBITALS_GT_127
! type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
! type(ik_vec) :: tmp
!#else
! integer(ik),intent(in) :: dets_up(:),dets_dn(:)
! integer(ik) :: tmp
!#endif
! real(rk),intent(out)        :: lowest_eigenvector(:)
! real(rk),intent(out)        :: lowest_eigenvalue
! integer,intent(in)          :: hash_info_singles(:),n_per_bin_singles(:),start_index_singles(:)
! integer,intent(in)          :: hash_info_doubles(:),n_per_bin_doubles(:),start_index_doubles(:)
! real(rk),intent(out),optional :: highest_eigenvalue,second_lowest_eigenvalue
! real(rk),intent(in),optional :: initial_vector(:)
!
! ! Local
! integer                    :: i,it,i_elec,j
!!real(rk)                   :: rannyu
! real(rk)                   :: energy_shift
! real(rk)                   :: norm,norm_inv
! real(rk),allocatable       :: w(:),v(:,:)
! real(rk),allocatable       :: alphas(:),betas(:)
! integer                    :: iterations
! real(rk)                   :: lowest_eigenvalue_prev
! logical                    :: converged=.false.
! integer                    :: len_work,info
! real(rk),allocatable       :: work(:),eigenvalues(:),tridiag(:,:)
! real(rk)                   :: lanczos_epsilon=1.e-10_rk
! integer,allocatable        :: ref_con_indices(:)
! integer                    :: n_con
! real(rk)                   :: H_ij
! logical                    :: is_con
!
!  allocate(ref_con_indices(max_connected_dets+(nelec*(nelec+1))/2)) ! Because there can be up to nelec(nelec-1)/2 copies of the diagonal element coming from doubles and up to nelec copies coming from singles
!  iterations=100          ! User option
!  iterations=min(n,iterations)
!  allocate (v(n,iterations+1))
!  allocate (w(n))
!  allocate(alphas(iterations+1))
!  allocate(betas(iterations+1))
!  w(:)=0._rk
!
!  if (present(initial_vector)) then
!    norm = 1._rk/sqrt(dot_product(initial_vector,initial_vector))
!    v(:,1) = norm*initial_vector(:)
!  else
!    v(:,1)=0
!    v(1,1)=1
!  endif
!
!  energy_shift=0._rk
!  betas(1)=0._rk
!  allocate (tridiag(iterations,iterations))
!
!  allocate(eigenvalues(iterations))
!  len_work = 3*iterations-1
!  allocate(work(len_work))
!
!  converged=.false.
!
!  write(6,'(/,''Executing matrix_lanczos_partial_connections in more_tools.f90'')')
!  if (n>1) then
!      do it=1,iterations
!         ! w(:) = H*v(:,it)
!         do i=1,n
!
!           ! Get occ_up, occ_dn in O(N_elec) time
!           tmp = dets_up(i)
!           do i_elec=1,nup
!             occ_up(i_elec) = trailz(tmp) + 1
!             tmp = ibclr(tmp,occ_up(i_elec)-1)
!           enddo
!           tmp = dets_dn(i)
!           do i_elec=1,ndn
!             occ_dn(i_elec) = trailz(tmp) + 1
!             tmp = ibclr(tmp,occ_dn(i_elec)-1)
!           enddo
!
!           ! Diagonal element
!           call hamiltonian(dets_up(i),dets_dn(i),dets_up(i),dets_dn(i),H_ij,is_con)
!           w(i) = w(i) + H_ij*v(i,it)
!
!           ! Single excitations
!           call get_singly_connected_ref_dets(dets_up(i),dets_dn(i),n_con,ref_con_indices,occ_up,occ_dn,hash_info_singles,n_per_bin_singles,start_index_singles)
!           do j=1,n_con
!             if (i==ref_con_indices(j))  cycle ! already did diagonal element
!             call hamiltonian(dets_up(i),dets_dn(i),dets_up(ref_con_indices(j)),dets_dn(ref_con_indices(j)),H_ij,is_con)
!             if (is_con)  w(ref_con_indices(j)) = w(ref_con_indices(j)) + H_ij*v(i,it)
!           enddo
!
!           ! Double excitations
!           call get_doubly_connected_ref_dets(dets_up(i),dets_dn(i),n_con,ref_con_indices,occ_up,occ_dn,hash_info_doubles,n_per_bin_doubles,start_index_doubles)
!           do j=1,n_con
!             if (i==ref_con_indices(j))  cycle ! already did diagonal element
!             ! Now filter out single excitations:
!             if (dets_up(i)==dets_up(ref_con_indices(j))) then
!               if (popcnt(iand(dets_dn(i),not(dets_dn(ref_con_indices(j)))))==1)  cycle
!             elseif (dets_dn(i)==dets_dn(ref_con_indices(j))) then
!               if (popcnt(iand(dets_up(i),not(dets_up(ref_con_indices(j)))))==1)  cycle
!             endif
!             call hamiltonian(dets_up(i),dets_dn(i),dets_up(ref_con_indices(j)),dets_dn(ref_con_indices(j)),H_ij,is_con)
!             if (is_con)  w(ref_con_indices(j)) = w(ref_con_indices(j)) + H_ij*v(i,it)
!           enddo
!
!         enddo ! i
!
!         if (it .gt. 1) w(:)=w(:)-betas(it)*v(:,it-1)
!         alphas(it)=dot_product(w,v(:,it))
!         w(:)=w(:)-alphas(it)*v(:,it)
!         norm=dot_product(w,w)
!         if (norm<(1.e-12_rk))  converged=.true.
!         betas(it+1)=norm**(0.5_rk)
!         norm_inv=1._rk/betas(it+1)
!         v(:,it+1)=w(:)*norm_inv
!         w(:)=v(:,it+1)
!         do i=1,it        ! Reorthogonalization
!             norm=dot_product(v(:,it+1),v(:,i))
!             call flush(6)
!             w(:)=w(:)-norm*v(:,i)
!         enddo
!         v(:,it+1)=w(:)
!         w(:)=0._rk
!         norm=dot_product(v(:,it+1),v(:,it+1))
!         norm_inv=1._rk/(norm**(0.5_rk))
!         v(:,it+1)=v(:,it+1)*norm_inv
!         tridiag(:,:)=0._rk
!
!         eigenvalues(:)=0._rk
!
!         do i=1,it
!             tridiag(i,i)=alphas(i)
!             if (i<it) then
!                 tridiag(i,i+1)=betas(i+1)
!                 tridiag(i+1,i)=betas(i+1)
!             endif
!         enddo
!
!         !diagonalize with lapack routine
!         len_work = 3*it-1
!
!         call dsyev('V', 'U', it, tridiag(1:it,1:it), it, eigenvalues, work, len_work, info)
!
!         !deallocate(tridiag) !deallocate(work)
!         lowest_eigenvalue=eigenvalues(1)
!         if (present(highest_eigenvalue))  highest_eigenvalue=eigenvalues(it)
!         if (present(second_lowest_eigenvalue).and.it>1)  second_lowest_eigenvalue=eigenvalues(2)
!         !call print_real_matrix(size(eigenvalues),1,eigenvalues)
!         if (it.gt.1 .and. abs(lowest_eigenvalue-lowest_eigenvalue_prev)<lanczos_epsilon) then
!             converged=.true.
!             exit
!         else
!             lowest_eigenvalue_prev=lowest_eigenvalue
!             write(6,'(''Iteration, Eigenvalue='',i3,f15.9)') it, lowest_eigenvalue
!             call flush(6)
!         endif
!         if (converged)  exit
!
!      enddo
!
!      it=min(it,iterations)
!      write(6,'(''matrix_lanczos_partial_connections: n, Lowest eigenvalue ='',i10, f16.10)') n, lowest_eigenvalue
!      call flush(6)
!
!      v(:,1)=matmul(v(:,1:it),tridiag(1:it,1))
!
!      if (allocated(eigenvalues)) deallocate(eigenvalues)
!      if (allocated(tridiag))     deallocate(tridiag)
!      if (allocated(work))        deallocate(work)
!
!  else
!      write (6,'(''Diagonalization attempted with n=1'')')
!      call hamiltonian(dets_up(1),dets_dn(1),dets_up(1),dets_dn(1),lowest_eigenvalue,is_con)
!  endif
!
!  lowest_eigenvector(1:n)=v(1:n,1)
!
!  end subroutine matrix_lanczos_partial_connections
!-------------------------------------------------------------------------------------

  subroutine get_1rdm(ndet, dets_up_in, dets_dn_in, coeffs_in, rdm, remote_det_map, local_det_map, state)
  ! Gets the 1-RDM of the variational wavefunction
  ! A Holmes, 20 Jun 2016

    use tools, only : merge_sort2_up_dn,permutation_factor
    use more_tools, only : get_occ_orbs,binary_search,print_real_matrix
    use chemistry, only : orbital_symmetries, time_sym, z
    use mpi_routines, only : det_map, det_map_l, mpi_allred, master_core_node,&
      & mpi_barr_in_node,shmem_allocate,shmem_deallocate

    type(det_map),optional,intent(in) :: remote_det_map
    type(det_map_l),optional,intent(in) :: local_det_map
    integer,intent(in) :: ndet
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up_in(:),dets_dn_in(:)
    type(ik_vec),pointer,dimension(:) :: dets_up,dets_dn
    type(ik_vec),allocatable :: temp_i16_up(:),temp_i16_dn(:)
#else
    integer(ik),intent(in) :: dets_up_in(:),dets_dn_in(:)
    integer(ik),pointer,dimension(:) :: dets_up,dets_dn
    integer(ik),allocatable :: temp_i16_up(:),temp_i16_dn(:)
#endif
    real(rk),intent(in) :: coeffs_in(:)
    real(rk),intent(out) :: rdm(norb,norb)
    integer,optional,intent(in) :: state ! only need for mpi

    real(rk),pointer,dimension(:) :: coeffs
    integer,allocatable :: iorder(:),temp_i_2(:)
    integer :: idet,i_elec,p,r,ndets_global
    integer :: j,k
    integer,allocatable :: occ_up(:),occ_dn(:)
    integer(i8b) :: n_allocate
    integer :: i

    call my_second(1,'1rdm')

    ! If time-reversal symmetry is used, re-express wavefunction is basis of determinants 
    ! (rather than linear combinations of determinants) first, before computing 1-RDM

    if (present(remote_det_map)) then
      ndets_global = local_det_map%ndets_global_old
      if (time_sym) then
        n_allocate = ndets_global
        do i=1,ndets_global
          if (remote_det_map%global_dets_up(i).ne.remote_det_map%global_dets_dn(i))  n_allocate = n_allocate + 1
        enddo
      else
        n_allocate = ndets_global
      endif
      write (6,*) "mpi: n_allocate=",n_allocate; call flush(6)
      call shmem_allocate(dets_up,int(n_allocate,i8b))
      call shmem_allocate(dets_dn,int(n_allocate,i8b))
      call shmem_allocate(coeffs,int(n_allocate,i8b))

      if (master_core_node) then
        if (time_sym) then
          j = 0
          do i=1,ndets_global
            j = j + 1
            dets_up(j) = remote_det_map%global_dets_up(i)
            dets_dn(j) = remote_det_map%global_dets_dn(i)
            if (dets_up(j)==dets_dn(j)) then
              coeffs(j) = remote_det_map%remote_wts(i,state)
            else
              coeffs(j) = inv_sqrt2*(remote_det_map%remote_wts(i,state))
              j = j + 1
              dets_up(j) = remote_det_map%global_dets_dn(i)
              dets_dn(j) = remote_det_map%global_dets_up(i)
              coeffs(j) = z*inv_sqrt2*(remote_det_map%remote_wts(i,state))
            endif
          enddo
          ndets_global = n_allocate
        else
          dets_up(1:ndets_global) = remote_det_map%global_dets_up(1:ndets_global)
          dets_dn(1:ndets_global) = remote_det_map%global_dets_dn(1:ndets_global)
          coeffs(1:ndets_global) = remote_det_map%remote_wts(1:ndets_global,state)
        endif
      endif
      call mpi_barr_in_node()

    else ! serial code

      if (time_sym) then
        n_allocate = ndet
        do i=1,ndet
          if (dets_up_in(i).ne.dets_dn_in(i))  n_allocate = n_allocate + 1
        enddo
      else
        n_allocate = ndet
      endif
      write (6,*) "n_allocate=",n_allocate; call flush(6)
      allocate(dets_up(n_allocate))
      allocate(dets_dn(n_allocate))
      allocate(coeffs(n_allocate))

      if (time_sym) then
        j = 0
        do i=1,ndet
          j = j + 1
          dets_up(j) = dets_up_in(i)
          dets_dn(j) = dets_dn_in(i)
          if (dets_up(j)==dets_dn(j)) then
            coeffs(j) = coeffs_in(i)
          else
            coeffs(j) = inv_sqrt2*coeffs_in(i)
            j = j + 1
            dets_up(j) = dets_dn_in(i)
            dets_dn(j) = dets_up_in(i)
            coeffs(j) = z*inv_sqrt2*coeffs_in(i)
          endif
        enddo
        ndets_global = n_allocate
      else
        dets_up(1:ndet) = dets_up_in(1:ndet)
        dets_dn(1:ndet) = dets_dn_in(1:ndet)
        coeffs(1:ndet) = coeffs_in(1:ndet)
        ndets_global = ndet
      endif
    endif

    if (master_core_node) then
      ! Sort by label for binary search

      allocate(iorder(ndets_global))
      allocate(temp_i16_up((ndets_global+1)/2))
      allocate(temp_i16_dn((ndets_global+1)/2))
      allocate(temp_i_2((ndets_global+1)/2))
      do j=1,ndets_global
        iorder(j)=j
      enddo
      call merge_sort2_up_dn(dets_up(1:ndets_global),dets_dn(1:ndets_global), iorder, ndets_global, temp_i16_up, temp_i16_dn, temp_i_2)
      coeffs(1:ndets_global) = coeffs(iorder(1:ndets_global))

      deallocate(iorder)
      deallocate(temp_i16_up)
      deallocate(temp_i16_dn)
      deallocate(temp_i_2)
    endif
    call mpi_barr_in_node()

    allocate(occ_up(nup))
    allocate(occ_dn(ndn))

    write (6,'(/,''Computing 1-RDM'')')
    rdm(:,:) = 0._rk

    !MJO Loop through local determinants
    do idet=1,ndets_global
    !do idet=1,ndet
      call get_occ_orbs(dets_up(idet),dets_dn(idet),occ_up,occ_dn)
      do i_elec=1,nup
        p = occ_up(i_elec)
        do r=1,norb
          if (orbital_symmetries(p).ne.orbital_symmetries(r))  cycle
          if (p.ne.r.and.btest(dets_up(idet),r-1))  cycle
          ! MJO But search in global list
          call binary_search(ibset(ibclr(dets_up(idet),p-1),r-1),dets_dn(idet),dets_up,dets_dn,k)
          if (k>0) then
           !rdm(p,r) = rdm(p,r) + coeffs(idet)*coeffs(k)
            rdm(p,r) = rdm(p,r) + permutation_factor(dets_up(idet),dets_up(k))*coeffs(idet)*coeffs(k)
          endif
        enddo ! r
      enddo ! i_elec
      do i_elec=1,ndn
        p = occ_dn(i_elec)
        do r=1,norb
          if (orbital_symmetries(p).ne.orbital_symmetries(r))  cycle
          if (p.ne.r.and.btest(dets_dn(idet),r-1))  cycle
          call binary_search(dets_up(idet),ibset(ibclr(dets_dn(idet),p-1),r-1),dets_up,dets_dn,k)
          if (k>0) then
           !rdm(p,r) = rdm(p,r) + coeffs(idet)*coeffs(k)
            rdm(p,r) = rdm(p,r) + permutation_factor(dets_dn(idet),dets_dn(k))*coeffs(idet)*coeffs(k)
          endif
        enddo ! r
      enddo ! i_elec
    enddo ! idet


   !call mpi_allred(rdm)
    write (6,'(/,''1-RDM:'')')
    call print_real_matrix(norb,norb,rdm)

   !write (6,*) "Writing 1-RDM to file rdm_file"
   !open(8, file='rdm_file',status='new')

   !do p=1,norb
   !  write(8,'(1000es13.5)') rdm(p,:)
   !enddo ! p

   !close(8)

    !write (6,*) "Done computing and printing 1-RDM"
    if (present(remote_det_map)) then
      call shmem_deallocate(dets_up)
      call shmem_deallocate(dets_dn)
      call shmem_deallocate(coeffs)
    endif

    call my_second(2,'1rdm')

  end subroutine get_1rdm

    subroutine get_1rdm_with_pt(ndet, dets_up_in, dets_dn_in, coeffs_in, diag_elems, var_energy, eps_pt_big, rdm)
  ! Gets the 1-RDM to lowest nonzero order in PT, i.e., if psi = psi0+psi1+..., then
  ! <psi|rho|psi> ~ <psi0|rho|psi0> + 2 <psi0|rho|psi1>
  ! Just like in the HCI PT energy correction, only uses terms in the numerator for which |H_{ij} c_j| > eps_pt_big
  ! A Holmes, 27 Jun 2016

    use tools, only : merge_sort2_up_dn,permutation_factor
    use more_tools, only : get_occ_orbs,binary_search,print_real_matrix
    use semistoch, only : find_doubly_excited,hamiltonian
   !use common_run, only : connected_dets_up,connected_dets_dn

    integer,intent(in) :: ndet
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up_in(:),dets_dn_in(:)
    type(ik_vec),allocatable :: dets_up(:),dets_dn(:)
    type(ik_vec),allocatable :: temp_i16_up(:),temp_i16_dn(:)
    type(ik_vec),allocatable :: connected_dets_up(:),connected_dets_dn(:)
#else
    integer(ik),intent(in) :: dets_up_in(:),dets_dn_in(:)
    integer(ik),allocatable :: dets_up(:),dets_dn(:)
    integer(ik),allocatable :: temp_i16_up(:),temp_i16_dn(:)
    integer(ik),allocatable :: connected_dets_up(:),connected_dets_dn(:)
#endif
    real(rk),intent(in) :: coeffs_in(:),diag_elems(:),var_energy,eps_pt_big
   !integer,intent(in) :: n_batches
    real(rk),intent(out) :: rdm(:,:)

    real(rk),allocatable :: coeffs(:)
    integer,allocatable :: iorder(:),temp_i_2(:)
    integer :: idet,i_elec,p,r
    integer :: j,k,k2
    integer,allocatable :: occ_up(:),occ_dn(:)
    real(rk) :: norm ! because 1st-order corrected wavefunction is not normalized
    integer :: i !batch,i
    real(rk),allocatable :: e_mix_den(:)
    real(rk),allocatable :: connected_wts(:),new_diag_elems(:)
    integer :: ndets_connected

    allocate(dets_up(ndet))
    allocate(dets_dn(ndet))
    allocate(coeffs(ndet))
    dets_up(1:ndet) = dets_up_in(1:ndet)
    dets_dn(1:ndet) = dets_dn_in(1:ndet)
    coeffs(1:ndet) = coeffs_in(1:ndet)

    ! Sort by label for binary search
    allocate(iorder(ndet))
    allocate(temp_i16_up((ndet+1)/2))
    allocate(temp_i16_dn((ndet+1)/2))
    allocate(temp_i_2((ndet+1)/2))
    do j=1,ndet
      iorder(j)=j
    enddo

    call merge_sort2_up_dn(dets_up(1:ndet),dets_dn(1:ndet), iorder, ndet, temp_i16_up, temp_i16_dn, temp_i_2)
    coeffs(1:ndet) = coeffs(iorder(1:ndet))

    deallocate(iorder)
    deallocate(temp_i16_up)
    deallocate(temp_i16_dn)
    deallocate(temp_i_2)

    allocate(occ_up(nup))
    allocate(occ_dn(ndn))

    write (6,'(/,''Computing 1-RDM to lowest nonzero order in PT'',/)')


    rdm(:,:) = 0._rk

    call my_second(1,'call to find_doubly_excited in 1rdm')
    call find_doubly_excited(n_det=ndets_connected, dets_up=connected_dets_up, dets_dn=connected_dets_dn, ref_up=dets_up(1:ndet), ref_dn=dets_dn(1:ndet), norb=norb, n_core_orb=n_core_orb, ref_coeffs=coeffs(1:ndet), e_mix_num=connected_wts, e_mix_den=e_mix_den, eps_var_pt=eps_pt_big, ref_diag_elems=diag_elems(1:ndet), new_diag_elems=new_diag_elems)

    call my_second(2,'call to find_doubly_excited in 1rdm')
    write(6,'(/,''Actual ndets_connected in get_1rdm_with_pt'',i13,'' ='',es9.2)') ndets_connected, real(ndets_connected)

    ! Now loop over single excitations from psi0, binary searching each one

    do idet=1,ndet
      call get_occ_orbs(dets_up(idet),dets_dn(idet),occ_up,occ_dn)
      do i_elec=1,nup
        p = occ_up(i_elec)
        do r=1,norb
          if (p.ne.r.and.btest(dets_up(idet),r-1))  cycle
          if (p==r) then
            rdm(p,r) = rdm(p,r) + coeffs(idet)**2
            cycle
          endif
          call binary_search(ibset(ibclr(dets_up(idet),p-1),r-1),dets_dn(idet),connected_dets_up(1:ndets_connected),connected_dets_dn(1:ndets_connected),k)
          if (k>0) then
            call binary_search(ibset(ibclr(dets_up(idet),p-1),r-1),dets_dn(idet),dets_up(1:ndet),dets_dn(1:ndet),k2)
            if (k2>0) then ! connected to variational det
             !rdm(p,r) = rdm(p,r) + coeffs(idet)*coeffs(k2)
              rdm(p,r) = rdm(p,r) + permutation_factor(dets_up(idet),dets_up(k2))*coeffs(idet)*coeffs(k2)
            else ! connected to PT det
              ! Have to do both because <psi0|rho_{pr}|psi1> != <psi1|rho_{pr}|psi0>
             !rdm(p,r) = rdm(p,r) + coeffs(idet)*connected_wts(k)/(var_energy-new_diag_elems(k))
             !rdm(r,p) = rdm(r,p) + coeffs(idet)*connected_wts(k)/(var_energy-new_diag_elems(k))
              rdm(p,r) = rdm(p,r) + permutation_factor(dets_up(idet),connected_dets_up(k))*coeffs(idet)*connected_wts(k)/(var_energy-new_diag_elems(k))
              rdm(r,p) = rdm(r,p) + permutation_factor(dets_up(idet),connected_dets_up(k))*coeffs(idet)*connected_wts(k)/(var_energy-new_diag_elems(k))
            endif
          endif
        enddo ! r
      enddo ! i_elec
      do i_elec=1,ndn
        p = occ_dn(i_elec)
        do r=1,norb
          if (p.ne.r.and.btest(dets_dn(idet),r-1))  cycle
          if (p==r) then
            rdm(p,r) = rdm(p,r) + coeffs(idet)**2
            cycle
          endif
          call binary_search(dets_up(idet),ibset(ibclr(dets_dn(idet),p-1),r-1),connected_dets_up(1:ndets_connected),connected_dets_dn(1:ndets_connected),k)
          if (k>0) then
            call binary_search(dets_up(idet),ibset(ibclr(dets_dn(idet),p-1),r-1),dets_up(1:ndet),dets_dn(1:ndet),k2)
            if (k2>0) then ! connected to variational det
             !rdm(p,r) = rdm(p,r) + coeffs(idet)*coeffs(k2)
              rdm(p,r) = rdm(p,r) + permutation_factor(dets_dn(idet),dets_dn(k2))*coeffs(idet)*coeffs(k2)
            else ! connected to connected det
              ! Have to do both because <psi0|rho_{pr}|psi1> != <psi1|rho_{pr}|psi0>
             !rdm(p,r) = rdm(p,r) + coeffs(idet)*connected_wts(k)/(var_energy-new_diag_elems(k))
             !rdm(r,p) = rdm(r,p) + coeffs(idet)*connected_wts(k)/(var_energy-new_diag_elems(k))
              rdm(p,r) = rdm(p,r) + permutation_factor(dets_dn(idet),connected_dets_dn(k))*coeffs(idet)*connected_wts(k)/(var_energy-new_diag_elems(k))
              rdm(r,p) = rdm(r,p) + permutation_factor(dets_dn(idet),connected_dets_dn(k))*coeffs(idet)*connected_wts(k)/(var_energy-new_diag_elems(k))
            endif
          endif
        enddo ! r
      enddo ! i_elec

    enddo ! idet

   !enddo

    write (6,'(/,''PT-corrected 1-RDM:'')')
    call print_real_matrix(norb,norb,rdm)
   !do p=1,norb
   !  write(6,'(1000es16.8)') rdm(p,:)
   !enddo
   !write (6,*)

   !write (6,*) "Writing 1-RDM to file rdm_file"
   !open(8, file='rdm_file',status='new')

   !do p=1,norb
   !  write(8,'(1000es13.5)') rdm(p,:)
   !enddo ! p

   !close(8)

   !write (6,*) "Done computing and printing 1-RDM"

  end subroutine get_1rdm_with_pt
!-------------------------------------------------------------------------------------

  subroutine generate_natorb_integrals(rdm)
  ! Generates a new FCIDUMP file, FCIDUMP.natorb, with integrals transformed by rotating to the natural orbital basis
  ! A Holmes, 17 Jan 2017
    use more_tools, only : print_real_matrix
    use chemistry, only : integral_value,orbital_symmetries,n_group_elements,nuclear_nuclear_energy,n_group_elements,orb_order,orb_order_inv,point_group
    use common_selected_ci, only : eps_var_sched
    use generic_sort, only : sort_by_first_argument
    use mpi_routines, only : master_core

    real(rk),intent(inout) :: rdm(norb,norb)

    integer :: len_work,info
    real(rk),allocatable :: work(:),eigenvalues(:)
    real(rk),allocatable :: new_integrals(:,:,:,:),tmp_integrals(:,:,:,:)
    integer :: i,j,k,n,p,q,r,s
    integer,allocatable :: new_orbital_symmetries(:)
    integer,allocatable :: ints(:)
    integer,allocatable :: inds(:,:),n_in_group(:)
    integer :: irrep
    real(rk),allocatable :: tmp_rdm(:,:),tmp_eigenvalues(:)
    character*32 fmt, filename

    call my_second(1,'generate_natorb_integrals')

    ! First, find eigenvectors of rdm

    ! The following simple approach does not work because redundant rotations can mess up symmetry labels:
   !allocate(eigenvalues(norb))
   !len_work = 3*norb-1
   !allocate(work(len_work))
   !call dsyev('V', 'U', norb, rdm(1:norb,1:norb), norb, eigenvalues, work, len_work, info)
   !write (6,*) "Eigenvalues,eigenvectors of 1RDM:"
   !rdm(1:norb,1:norb) = rdm(1:norb,norb:1:-1)
   !eigenvalues(1:norb) = eigenvalues(norb:1:-1)
   !call print_real_matrix(norb,norb,rdm)
   !write (6,*) eigenvalues

    ! Instead, must perform diagonalization in each irrep separately:
    allocate(inds(n_group_elements,norb))
    allocate(n_in_group(n_group_elements))

    do irrep=1,n_group_elements
      n_in_group(irrep) = 0
      do i=1,norb
        if (orbital_symmetries(i)==irrep) then
          n_in_group(irrep) = n_in_group(irrep) + 1
          inds(irrep,n_in_group(irrep)) = i
        endif
      enddo
    enddo ! irrep

    allocate(tmp_rdm(norb,norb))
    allocate(eigenvalues(norb))
    allocate(tmp_eigenvalues(norb))
    allocate(work(3*norb-1))
    do irrep=1,n_group_elements
      if (n_in_group(irrep)==0)  cycle
      n = n_in_group(irrep)
      len_work = 3*n_in_group(irrep)-1
      write (6,*) "About to call dsyev on orbitals",inds(irrep,1:n); call flush(6)
      do i=1,n
        do j=1,n
          tmp_rdm(i,j) = rdm(inds(irrep,i),inds(irrep,j))
        enddo
      enddo
      call dsyev('V', 'U', n, tmp_rdm(1:n,1:n), n, tmp_eigenvalues(1:n), work(1:n), len_work, info)
      ! Reverse because want in decreasing order of eigenvalue (occupation
      ! number)
      do i=1,n
        eigenvalues(inds(irrep,i)) = tmp_eigenvalues(n-i+1)
        do j=1,n
          rdm(inds(irrep,i),inds(irrep,j)) = tmp_rdm(i,n-j+1)
     !do i=1,n
     !  eigenvalues(inds(irrep,i)) = tmp_eigenvalues(i)
     !  do j=1,n
     !    rdm(inds(irrep,i),inds(irrep,j)) = tmp_rdm(i,j)
        enddo
      enddo
    enddo ! irrep

    ! Reorder based on natural orbital occupation number
   !tmp_eigenvalues(1:norb) = eigenvalues(1:norb)
   !allocate(ints(norb))
   !do i=1,norb
   !  ints(i) = i
   !enddo
   !call sort_by_first_argument(norb,tmp_eigenvalues,ints)
   !eigenvalues(1:norb) = eigenvalues(ints(norb:1:-1))
   !rdm(:,1:norb) = rdm(:,ints(norb:1:-1))
   !allocate(new_orbital_symmetries(norb))
   !new_orbital_symmetries(1:norb) = orbital_symmetries(ints(norb:1:-1))
   !deallocate(ints)

    allocate(new_orbital_symmetries(norb))
    new_orbital_symmetries(1:norb) = orbital_symmetries(1:norb)
    write (6,'(''New orbital symmetries:'',1000i2)') new_orbital_symmetries(orb_order_inv(1:norb)) ; call flush(6)
    if (hamiltonian_type.eq.'chem') then
      if (point_group.eq.'dih') then
        ! Convert to Sandeep's irrep labels
        do i=1,norb
          n = new_orbital_symmetries(i)
          if (n<=2)  cycle
          if (mod(n+1,4)<=1) then ! + Lz
            if (mod(n,2)==1) then
              new_orbital_symmetries(i) = (n-1)/2+4
            else
              new_orbital_symmetries(i) = (n-1)/2+5
            endif
          else ! - Lz
            if (mod(n,2)==1) then
              new_orbital_symmetries(i) = -(n-3)/2-4
            else
              new_orbital_symmetries(i) = -(n-3)/2-5
            endif
          endif
        enddo ! i
        write (6,'(''With Lz labels, new orbital symmetries:'',1000i2)') new_orbital_symmetries(orb_order_inv(1:norb)) ; call flush(6)
      endif
    endif

    write (6,'(/,''Eigenvectors of 1RDM:'')')
    call print_real_matrix(norb,norb,rdm)
    write (6,'(''Eigenvalues of 1RDM:'')')
    write (6,'(10es11.4)') eigenvalues

    write (6,*) "New orbital symmetries:",new_orbital_symmetries; call flush(6)

    write (6,'(/,''Rotating integrals to new natural orbitals'')') ; call flush(6)

    ! Now, rotate integrals
    allocate(new_integrals(norb,norb,norb+1,norb+1))
    allocate(tmp_integrals(norb,norb,norb+1,norb+1))

    ! Two-body integrals
    do p=1,norb
      do q=1,norb
        do r=1,norb
          do s=1,norb
            tmp_integrals(p,q,r,s) = integral_value(p,q,r,s)
          enddo
        enddo
      enddo
    enddo

    new_integrals(:,:,:,:) = 0._rk

    do p=1,norb
      do q=1,norb
        do r=1,norb
          do s=1,norb
            new_integrals(p,q,r,s) = dot_product(rdm(:,p),tmp_integrals(:,q,r,s))
           !do k=1,norb
           ! !if (orbital_symmetries(k).ne.orbital_symmetries(p))  cycle
           !  new_integrals(p,q,r,s) = new_integrals(p,q,r,s) + rdm(k,p)*tmp_integrals(k,q,r,s)
           !enddo
          enddo
        enddo
      enddo
    enddo

    tmp_integrals = new_integrals
    new_integrals(:,:,:,:) = 0._rk

    do p=1,norb
      do q=1,norb
        do r=1,norb
          do s=1,norb
            new_integrals(p,q,r,s) = dot_product(rdm(:,q),tmp_integrals(p,:,r,s))
           !do k=1,norb
           ! !if (orbital_symmetries(k).ne.orbital_symmetries(q))  cycle
           !  new_integrals(p,q,r,s) = new_integrals(p,q,r,s) + rdm(k,q)*tmp_integrals(p,k,r,s)
           !enddo
          enddo
        enddo
      enddo
    enddo

    tmp_integrals = new_integrals
    new_integrals(:,:,:,:) = 0._rk

    do p=1,norb
      do q=1,norb
        do r=1,norb
          do s=1,norb
            new_integrals(p,q,r,s) = dot_product(rdm(:,r),tmp_integrals(p,q,1:norb,s))
           !do k=1,norb
           ! !if (orbital_symmetries(k).ne.orbital_symmetries(r))  cycle
           !  new_integrals(p,q,r,s) = new_integrals(p,q,r,s) + rdm(k,r)*tmp_integrals(p,q,k,s)
           !enddo
          enddo
        enddo
      enddo
    enddo

    tmp_integrals = new_integrals
    new_integrals(:,:,:,:) = 0._rk

    do p=1,norb
      do q=1,norb
        do r=1,norb
          do s=1,norb
            new_integrals(p,q,r,s) = dot_product(rdm(:,s),tmp_integrals(p,q,r,1:norb))
           !do k=1,norb
           ! !if (orbital_symmetries(k).ne.orbital_symmetries(s))  cycle
           !  new_integrals(p,q,r,s) = new_integrals(p,q,r,s) + rdm(k,s)*tmp_integrals(p,q,r,k)
           !enddo
          enddo
        enddo
      enddo
    enddo

    ! One-body integrals
    do p=1,norb
      do q=1,norb
        tmp_integrals(p,q,norb+1,norb+1) = integral_value(p,q,norb+1,norb+1)
      enddo
    enddo

    new_integrals(:,:,norb+1,norb+1) = 0._rk

    do p=1,norb
      do q=1,norb
        do k=1,norb
          new_integrals(p,q,norb+1,norb+1) = new_integrals(p,q,norb+1,norb+1) + rdm(k,p)*tmp_integrals(k,q,norb+1,norb+1)
        enddo
      enddo
    enddo

    tmp_integrals(:,:,norb+1,norb+1) = new_integrals(:,:,norb+1,norb+1)
    new_integrals(:,:,norb+1,norb+1) = 0._rk

    do p=1,norb
      do q=1,norb
        do k=1,norb
          new_integrals(p,q,norb+1,norb+1) = new_integrals(p,q,norb+1,norb+1) + rdm(k,q)*tmp_integrals(p,k,norb+1,norb+1)
        enddo
      enddo
    enddo

    ! Finally, dump into file FCIDUMP_eps_var...
    if(master_core) then
      write(fmt,'(es7.2e1)') minval(eps_var_sched)
      write(filename,'(''FCIDUMP_eps_var='',a)') trim(fmt)
      open(88, file=filename, status='new')
      write (6,'(''Done computing new integrals. Dumping new integrals into file '',a)') filename ; call flush(6)

      ! Header
      write (88,'(''&FCI NORB='',i5,'', NELEC='',i5,'', MS2=0'')') norb,nelec
      write (88,'(''ORBSYM='',100000i3)') new_orbital_symmetries(orb_order_inv(1:norb))
      write (88,'(''ISYM=1'')')
      allocate(ints(max(nup,ndn)))
      do i=1,max(nup,ndn)
        ints(i) = i
      enddo
     !write(fmt,'(i5,''i3'')') nup
     !write (88,'(' // fmt // ', '' hf_up'')') ints(1:nup)
     !write(fmt,'(i5,''i3'')') ndn
     !write (88,'(' // fmt // ', '' hf_dn'')') ints(1:ndn)
      write (88,'(''&END'')')
      deallocate(ints)

      ! Two-body integrals
      do p=1,norb
        do q=1,p
          do r=1,p
            do s=1,r
              if (p==r.and.q<s)  cycle
              if (abs(new_integrals(p,q,r,s))>1e-8_rk) then
                write(88,'(es19.12e2,4i4)')  new_integrals(p,q,r,s),orb_order(p),orb_order(q),orb_order(r),orb_order(s)
              endif
            enddo
          enddo
        enddo
      enddo

      ! One-body integrals
      do p=1,norb
        do q=1,p
          if (abs(new_integrals(p,q,norb+1,norb+1))>1e-8_rk) then
            write(88,'(es19.12e2,4i4)')  new_integrals(p,q,norb+1,norb+1),orb_order(p),orb_order(q),0,0
          endif
        enddo
      enddo

      ! Nuclear-nuclear energy
      write(88,'(es19.12e2,4i4)')  nuclear_nuclear_energy,0,0,0,0

      close(88)
    endif ! master_core

    call my_second(2,'generate_natorb_integrals')

  end subroutine generate_natorb_integrals


  subroutine get_zeroth_order_variational_greens_function(ndet,dets_up_in,dets_dn_in,coeffs_in,e0,n_w,w,G0_np1,G0_nm1)

  ! Get the zeroth-order Green's function of the variational wavefunction
  ! By zeroth-order, I mean that the zeroth-order Hamiltonian in the (N_el+/-1)-space is diagonal
  ! This is an efficient choice when there are many energies w to evaluate the Green's function for,
  ! but more accurate choices of H0 can be chosen if we can afford to iteratively compute the inverse of a matrix
  ! many times

  ! The Green's function G0 is divided into two pieces, for N+1 and N-1 electrons (G0_np1 and G0_nm1, respectively):
  !   G0(w,p,q) = G0_np1(w,p,q) + G0_nm1(w,p,q)
  !   G0_np1(w,p,q) = < var | a_p 1/(w-(H_0-E_0)) a_q^\dagger | var > (where H_0 is the N+1 electron H0)
  !   G0_nm1(w,p,q) = < var | a_p^\dagger 1/(w-(E_0-H_0)) a_q | var > (where H_0 is the N-1 electron H0)

  ! Time scaling of this implementation is N_el*N_orb*N_var*log(N_var)

  ! A Holmes, 27 Jan 2017

    use tools, only : merge_sort2_up_dn
    use more_tools, only : get_occ_orbs,binary_search,print_real_matrix
    use chemistry, only : hamiltonian_chem, time_sym

    integer,intent(in) :: ndet
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up_in(:),dets_dn_in(:)
    type(ik_vec),allocatable :: dets_up(:),dets_dn(:)
    type(ik_vec),allocatable :: temp_i16_up(:),temp_i16_dn(:)
    type(ik_vec) :: tmp_up,tmp_dn
#else
    integer(ik),intent(in) :: dets_up_in(:),dets_dn_in(:)
    integer(ik),allocatable :: dets_up(:),dets_dn(:)
    integer(ik),allocatable :: temp_i16_up(:),temp_i16_dn(:)
    integer(ik) :: tmp_up,tmp_dn
#endif
    real(rk),intent(in) :: coeffs_in(:)
    real(rk),intent(in) :: e0
    integer,intent(in) :: n_w
    real(rk),intent(in) :: w(:)
    real(rk),intent(out) :: G0_np1(:,:,:),G0_nm1(:,:,:)

    real(rk),allocatable :: coeffs(:),rho(:,:)
    integer,allocatable :: iorder(:),temp_i_2(:)
    integer :: idet,i_elec,p,q
    integer :: j,k,i_w
    integer,allocatable :: occ_up(:),occ_dn(:)
    real(rk) :: h_ii

    if (time_sym)  stop "Time sym not implemented for Green's functions yet"

    allocate(dets_up(ndet))
    allocate(dets_dn(ndet))
    allocate(coeffs(ndet))
    dets_up(1:ndet) = dets_up_in(1:ndet)
    dets_dn(1:ndet) = dets_dn_in(1:ndet)
    coeffs(1:ndet) = coeffs_in(1:ndet)

    ! Sort by label for binary search
    allocate(iorder(ndet))
    allocate(temp_i16_up((ndet+1)/2))
    allocate(temp_i16_dn((ndet+1)/2))
    allocate(temp_i_2((ndet+1)/2))
    do j=1,ndet
      iorder(j)=j
    enddo

    call merge_sort2_up_dn(dets_up(1:ndet),dets_dn(1:ndet), iorder, ndet, temp_i16_up, temp_i16_dn, temp_i_2)
    coeffs(1:ndet) = coeffs(iorder(1:ndet))

    deallocate(iorder)
    deallocate(temp_i16_up)
    deallocate(temp_i16_dn)
    deallocate(temp_i_2)

    allocate(occ_up(nup+1))
    allocate(occ_dn(ndn+1))

    G0_np1(:,:,:) = 0._rk
    G0_nm1(:,:,:) = 0._rk

    do idet=1,ndet

      ! G0_N+1

      ! up electrons
      tmp_dn = dets_dn(idet)
      do q=1,norb
        if (btest(dets_up(idet),q-1))  cycle
        tmp_up = ibset(dets_up(idet),q-1)
        call hamiltonian_chem(tmp_up,tmp_dn,tmp_up,tmp_dn,0,h_ii) ! N+1 electron Hamiltonian
        call get_occ_orbs(tmp_up,tmp_dn,occ_up(1:nup+1),occ_dn(1:ndn))
        do i_elec=1,nup+1
          p = occ_up(i_elec)
          ! binary search on dets_up,dn
          call binary_search(ibclr(tmp_up,p-1),tmp_dn,dets_up(1:ndet),dets_dn(1:ndet),k)
          if (k>0) then ! found in variational wf
            do i_w=1,n_w
              G0_np1(i_w,p,q) = G0_np1(i_w,p,q) + coeffs(idet)*coeffs(k) * 1._rk/(w(i_w)-(h_ii-e0))
            enddo ! w
          endif ! found in variational wf
        enddo ! p (annihilation)
      enddo ! q (creation)

      ! dn electrons
      tmp_up = dets_up(idet)
      do q=1,norb
        if (btest(dets_dn(idet),q-1))  cycle
        tmp_dn = ibset(dets_dn(idet),q-1)
        call hamiltonian_chem(tmp_up,tmp_dn,tmp_up,tmp_dn,0,h_ii) ! N+1 electron Hamiltonian
        call get_occ_orbs(tmp_up,tmp_dn,occ_up(1:nup),occ_dn(1:ndn+1))
        do i_elec=1,ndn+1
          p = occ_dn(i_elec)
          ! binary search on dets_up,dn
          call binary_search(tmp_up,ibclr(tmp_dn,p-1),dets_up(1:ndet),dets_dn(1:ndet),k)
          if (k>0) then ! found in variational wf
            do i_w=1,n_w
              G0_np1(i_w,p,q) = G0_np1(i_w,p,q) + coeffs(idet)*coeffs(k) * 1._rk/(w(i_w)-(h_ii-e0))
            enddo ! w
          endif ! found in variational wf
        enddo ! p (annihilation)
      enddo ! q (creation)

      ! G0_N-1
      call get_occ_orbs(dets_up(idet),dets_dn(idet),occ_up(1:nup),occ_dn(1:ndn))

      ! up electrons
      tmp_dn = dets_dn(idet)
      do i_elec=1,nup
        q = occ_up(i_elec)
        tmp_up = ibclr(dets_up(idet),q-1)
        call hamiltonian_chem(tmp_up,tmp_dn,tmp_up,tmp_dn,0,h_ii) ! N-1 electron Hamiltonian
        do p=1,norb
          if (p.ne.q.and.btest(dets_up(idet),p-1))  cycle
          ! binary search on dets_up,dn
          call binary_search(ibset(tmp_up,p-1),tmp_dn,dets_up(1:ndet),dets_dn(1:ndet),k)
          if (k>0) then ! found in variational wf
            do i_w=1,n_w
              G0_nm1(i_w,p,q) = G0_nm1(i_w,p,q) + coeffs(idet)*coeffs(k) * 1._rk/(w(i_w)-(e0-h_ii))
            enddo ! w
          endif ! found in variational wf
        enddo ! p (creation)
      enddo ! q (annihilation)

      ! dn electrons
      tmp_up = dets_up(idet)
      do i_elec=1,ndn
        q = occ_dn(i_elec)
        tmp_dn = ibclr(dets_dn(idet),q-1)
        call hamiltonian_chem(tmp_up,tmp_dn,tmp_up,tmp_dn,0,h_ii) ! N-1 electron Hamiltonian
        do p=1,norb
          if (p.ne.q.and.btest(dets_dn(idet),p-1))  cycle
          ! binary search on dets_up,dn
          call binary_search(tmp_up,ibset(tmp_dn,p-1),dets_up(1:ndet),dets_dn(1:ndet),k)
          if (k>0) then ! found in variational wf
            do i_w=1,n_w
              G0_nm1(i_w,p,q) = G0_nm1(i_w,p,q) + coeffs(idet)*coeffs(k) * 1._rk/(w(i_w)-(e0-h_ii))
            enddo ! w
          endif ! found in variational wf
        enddo ! p (creation)
      enddo ! q (annihilation)

    enddo ! idet

    allocate(rho(n_w,2))
    rho(:,:) = 0._rk

    ! Green's functons

    open(60, file='G0_N+1', status='new')
    write (60,*) "w,p,q,G0_{N+1}(w,p,q)"
    do i_w=1,n_w
      do q=1,norb
        do p=1,norb
          if (abs(G0_np1(i_w,p,q))>1.e-8_rk)  write (60,*) w(i_w),p,q,G0_np1(i_w,p,q)
          if (p==q)  rho(i_w,1) = rho(i_w,1) + G0_np1(i_w,p,q)
        enddo ! p
      enddo ! q
    enddo ! w
    call flush(60)
    close(60)

    open(61, file='G0_N-1', status='new')
    write (61,*) "w,p,q,G0_{N-1}(w,p,q)"
    do i_w=1,n_w
      do q=1,norb
        do p=1,norb
          if (abs(G0_nm1(i_w,p,q))>1.e-8_rk)  write (61,*) w(i_w),p,q,G0_nm1(i_w,p,q)
          if (p==q)  rho(i_w,2) = rho(i_w,2) + G0_nm1(i_w,p,q)
        enddo ! p
      enddo ! q
    enddo ! w
    call flush(61)
    close(61)

    open(62, file='rho', status='new')
    call flush(62)
    write (62,*) "Density of states: w, rho_{N+1}(w),rho_{N-1}(w),rho_tot(w)"
    do i_w=1,n_w
      write (62,*) w(i_w),rho(i_w,1),rho(i_w,2),rho(i_w,1)+rho(i_w,2)
    enddo
    close(62)


  end subroutine get_zeroth_order_variational_greens_function
!--------------------------------------------------------------------------------

!  subroutine write_wf(eps_pt, eps_pt_big, target_error, n_states, hf_up, hf_dn)
!
!    use chemistry, only : orbital_energies,time_sym,z
!    use types, only : i8b
!    use semistoch, only : hamiltonian, estimate_n_connections
!    use tools, only : merge_sort2_up_dn
!    use heg, only : energy_hf, energy_madelung
!    use common_psi_t, only : trial_wf_iters
!    use common_run, only : n_connected_dets_hf
!    use common_selected_ci, only : eps_var_sched, dump_wf_var, sparse_ham, n_mc, get_natorbs, get_greens_function, n_w, w_min, w_max, use_pt, eps_pt_big_energy, n_var_orbs, var_orbs, n_var_e_up, n_var_e_dn
!    use mpi_routines, only : master_core_node, master_core, init_hash_owners, get_det_owner, whoami, &
!         &det_map, mpi_barr, det_map_l, mpi_barr_in_node, shmem_allocate, shmem_deallocate, sharedComm, mpierr, &
!         &nnodes, ncores, ncores_on_node, mpi_allred, mpi_bsend, mpi_distribute_remote_det_map, mpi_bsend_between_nodes, &
!         get_det_owner, mpi_stop
!#ifdef MPI
!    INCLUDE 'mpif.h'
!#endif
!
!    real(rk),intent(in)    :: eps_pt, target_error
!    real(rk)               :: eps_var
!    real(rk),intent(inout) :: eps_pt_big
!
!    integer :: iter, iter_done=0, max_iters=50
!#ifdef NUM_ORBITALS_GT_127
!    type(ik_vec),optional, intent(in) :: hf_up, hf_dn
!    type(ik_vec),allocatable :: old_dets_up(:), old_dets_dn(:), new_dets_up(:), new_dets_dn(:)
!    type(ik_vec),allocatable :: temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
!    type(ik_vec) :: core_up,core_dn,virt_up,virt_dn
!
!#else
!    integer(ik),optional, intent(in) :: hf_up, hf_dn
!    integer(ik),allocatable :: old_dets_up(:), old_dets_dn(:), new_dets_up(:), new_dets_dn(:)
!    integer(ik),allocatable :: temp_i16_up(:), temp_i16_dn(:) ! For walk_dets_up, walk_dets_dn. Needed for merge_sort2_up_dn
!    integer(ik) :: core_up,core_dn,virt_up,virt_dn
!#endif
!    integer,intent(in) :: n_states ! Number of low lying states to compute (1 for just ground state, 2 for ground and 1st excited, etc)
!    real(rk),allocatable :: old_wts(:,:),starting_wts(:,:)
!    integer :: ndets_old,ndets_new
!    real(rk),allocatable :: energy(:),old_energy(:)
!    logical :: is_con
!    integer :: i,j,k
!    integer,allocatable :: iorder(:),temp_i_2(:)
!!   integer :: tdets=9999999 ! max number of dets.  Warning: this should be redone in a more robust way.
!    real(rk),allocatable :: diag_elems(:),ave_connections_remote(:)
!    integer :: tdets=9999 ! Initial number of dets.
!    integer(i8b) :: n_conn_estimate
!    integer :: n_batches, ierr=0
!    real(rk) :: rdm(norb,norb)
!    real(rk) :: delta_e_2pt
!    real(rk) :: pt_big, pt_diff, pt_diff_std_dev, eps_pt_big_sav, eps_pt_big_read=9.d99, pt_big_read, avail, avail2, avail3
!    integer :: core_id, ndets_connected, ndets_global_old, ndets_global_new, ndet_remote, ndet_min, ndet_max, ndet_tot, n_mc_sav
!    real(rk) :: average_connections, nonzero_elems, ave_connections_min, ave_connections_max, nonzero_elems_min, nonzero_elems_max, nonzero_elems_tot
!    character*32 fmt, filename
!    logical wf_file_exist
!    real(rk), allocatable :: min_H_already_done(:) ! Min of H elements already done (i.e., all connected matrix elements exceeding this in magnitude have already been generated)
!    real(rk), pointer, dimension(:),contiguous :: diag_elems_global
!    real(rk),allocatable :: coeffs(:)
!    real(rk),allocatable :: w(:),G0_np1(:,:,:),G0_nm1(:,:,:)
!    real(rk),allocatable :: rdm_state_avg(:,:)
!    logical :: use_var_active_space=.false.
!    integer :: n_core_up,n_core_dn, n_virt_up,n_virt_dn
!
!    !MJO - det_map is used for the parallel run
!    type(det_map) remote_det_map
!    type(det_map_l) local_det_map
!
!      ! Dump variational wavefn to file if it does not exist
!      if(dump_wf_var .and. .not. wf_file_exist) then
!        write(6,'(/,''Writing variational wavefn to '',a)') trim(filename)
!        call my_second(1, 'dump variational wavefn')
!        if(.not.use_mpi) then
!           open(1, file=filename,form='unformatted',status='new')
!           write(1) ndets_old
!           write(1) old_dets_up(1:ndets_old)
!           write(1) old_dets_dn(1:ndets_old)
!           write(1) old_wts(1:ndets_old,1:n_states)
!           write(1) energy(1:n_states)
!           call flush(1)
!          !close(1)
!        elseif(master_core) then
!           open(1, file=filename,form='unformatted',status='new')
!           write(1) ndets_global_old
!           write(1) remote_det_map%global_dets_up(1:ndets_global_old)
!           write(1) remote_det_map%global_dets_dn(1:ndets_global_old)
!           write(1) remote_det_map%remote_wts(1:ndets_global_old,1:n_states)
!           write(1) energy(1:n_states)
!           call flush(1)
!          !close(1)
!        endif ! ncores
!        call my_second(2, 'dump variational wavefn')
!      endif ! dump_wf_var
!
!
!  end subroutine write_wf

  subroutine do_pt(ndets,dets_up,dets_dn,wts,diag_elems,ndets_global,dets_up_global,dets_dn_global,wts_global,diag_elems_global,var_energy,eps_var,eps_pt,eps_pt_big,target_error,n_mc,pt_big,pt_diff,pt_diff_std_dev,ndets_connected,remote_det_map,local_det_map,core_up,core_dn,virt_up,virt_dn,active_only_in,active_tot_energy)
  ! Performs uncontracted second-order Epstein-Nesbet perturbation theory using
  ! the HCI screened sum (using eps_pt)
  ! Sums the smallest contributions stochastically if necessary
  ! Skips reading of deterministic PT from wf file for now (FIXME)
  ! Moved to its own subroutine by A Holmes on 25 May 2017

    use semistoch, only : find_doubly_excited,hamiltonian
    use mpi_routines, only : det_map,mpi_allred,det_map_l,master_core,master_core_node,nnodes,mpi_barr_in_node
    use chemistry, only : time_sym,hamiltonian_chem
    use semistoch, only : hamiltonian, estimate_n_connections
    use common_selected_ci, only : dump_wf_var,eps_pt_big_energy
    use tools, only : merge_sort2_up_dn

    integer,intent(in) :: ndets
    integer,optional,intent(in) :: ndets_global
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:) ! inout because they can be sorted in this routine if MPI not used
    type(ik_vec),optional,intent(in) :: dets_up_global(:),dets_dn_global(:) ! inout because they can be sorted in this routine if MPI not used
    type(ik_vec),optional,intent(in) :: core_up,core_dn,virt_up,virt_dn
    type(ik_vec),allocatable :: connected_dets_up(:),connected_dets_dn(:)
    type(ik_vec),allocatable :: temp_i16_up(:),temp_i16_dn(:)
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
    integer(ik),optional,intent(in) :: dets_up_global(:),dets_dn_global(:) ! inout because they can be sorted in this routine if MPI not used
    integer(ik),optional,intent(in) :: core_up,core_dn,virt_up,virt_dn
    integer(ik),allocatable :: connected_dets_up(:),connected_dets_dn(:)
    integer(ik),allocatable :: temp_i16_up(:),temp_i16_dn(:)
#endif
    real(rk),intent(in) :: wts(:),diag_elems(:)
    real(rk),optional,intent(in) :: wts_global(:),diag_elems_global(:)
    real(rk),intent(in) :: var_energy,eps_var,eps_pt
    real(rk),intent(inout) :: eps_pt_big
    real(rk),intent(in) :: target_error
    integer,intent(inout) :: n_mc
    real(rk),intent(out) :: pt_big,pt_diff,pt_diff_std_dev
    integer,intent(out) :: ndets_connected
   !character*32,intent(in) :: filename
    type(det_map),optional,intent(in) :: remote_det_map
    type(det_map_l),optional,intent(in) :: local_det_map
    logical,optional,intent(in) :: active_only_in
    real(rk),optional,intent(in) :: active_tot_energy

    real(rk),allocatable :: connected_wts(:),new_diag_elems(:),e_mix_den(:)
    integer :: i,j,k, ndets_connected_global
    integer,allocatable :: iorder(:),temp_i_2(:)
    real(rk) :: eps_pt_big_read=9.d99, pt_big_read, eps_pt_big_bsend
    integer(i8b) :: n_conn_estimate
    logical :: is_con
    integer :: ierr=0
    logical :: active_only = .true.
    logical :: wf_file_exit

    ! First, see if we have any deterministic PT info from wf file
   !inquire(file=filename, exist=wf_file_exist)

    if (present(active_only_in))  active_only = active_only_in

    !MJO - Each core does a local estimate and then we reduce all of the values
    write(6,'(/,''Estimating # of connections of all var. dets. by sampling connections of local list and summing (duplicates are not merged)'')')
    if (present(core_up)) then
      n_conn_estimate = estimate_n_connections(ndets,dets_up,dets_dn,wts(:),eps_pt,core_up,core_dn,virt_up,virt_dn,active_only)
    else
      n_conn_estimate = estimate_n_connections(ndets,dets_up,dets_dn,wts(:),eps_pt)
    endif
    call mpi_allred(n_conn_estimate)

    write (6,'(/,''If eps_pt='',es15.4,'' estimated number of connections to variational wavefn='',i13,'' ='',es9.2)') eps_pt, n_conn_estimate, real(n_conn_estimate) ; call flush(6)

    if (n_conn_estimate<nint(n_max_connections*nnodes,i8b) .and. eps_pt_big <= 0) then ! If you uncomment this line, then it will use eps_pt, if it can, only if a +ve eps_pt_big is not inputted.  This is useful for testing.
   !if (n_conn_estimate<nint(n_max_connections*nnodes,i8b)                      ) then ! If you uncomment this line, then it will use eps_pt, if it can, regardless of the value of eps_pt_big.  This is more efficient.
      ! Can be done deterministically with just 1 batch
      write (6,'(/,''Computing PT energy using straightforward deterministic method'')') ; call flush(6)

      if (present(core_up)) then
        if (use_mpi) then!MJO TEMP PARALLEL
          call second_order_pt(ndets, dets_up(1:ndets), dets_dn(1:ndets), wts(1:ndets), diag_elems(1:ndets), var_energy, eps_pt, pt_big, ndets_connected, remote_det_map, local_det_map, core_up=core_up, core_dn=core_dn, virt_up=virt_up, virt_dn=virt_dn, active_only=active_only)
        else
          call second_order_pt(ndets, dets_up(1:ndets), dets_dn(1:ndets), wts(1:ndets), diag_elems(1:ndets), var_energy, eps_pt, pt_big, ndets_connected, core_up=core_up, core_dn=core_dn, virt_up=virt_up, virt_dn=virt_dn, active_only=active_only)
        endif
      else
        if (use_mpi) then!MJO TEMP PARALLEL
          call second_order_pt(ndets, dets_up(1:ndets), dets_dn(1:ndets), wts(1:ndets), diag_elems(1:ndets), var_energy, eps_pt, pt_big, ndets_connected, remote_det_map, local_det_map)
        else
          call second_order_pt(ndets, dets_up(1:ndets), dets_dn(1:ndets), wts(1:ndets), diag_elems(1:ndets), var_energy, eps_pt, pt_big, ndets_connected)
        endif
      endif

      call mpi_allred(ndets_connected) ! Energy gets allreduced in second_order_pt, but ndets_connected does not (actually it does, but is not returned, so it has to be done again here - fix?)

      write (6,'(/,''PT_correction, eps_pt, ndets_connected for fully deterministic run='',f15.9,es12.4,i13,'' ='',es9.2)') pt_big, eps_pt, ndets_connected, real(ndets_connected) ; call flush(6)

      pt_diff = 0._rk
      pt_diff_std_dev = 0._rk

      return

    else ! Semistochastic: deterministic 1 batch with eps_pt_big and stochastic with both eps_pt_big and eps_pt


      ! If eps_pt_big is positive, use it.  Otherwise find an eps_pt_big value for which single-batch deterministic PT is possible
      if(eps_pt_big.le.0_rk) then
        eps_pt_big = eps_pt
        do while (n_conn_estimate>=nint(n_max_connections*nnodes,i8b))
          eps_pt_big = round_r(1.1*eps_pt_big*(real(n_conn_estimate)/(n_max_connections*nnodes))**.75,2)

          if(eps_pt_big >= eps_var) then
            write(6,'(''eps_pt_big>= eps_var, eps_pt_big, eps_var='',2es10.2)') eps_pt_big, eps_var
            exit
          endif
          !MJO - Each core does a local estimate and then we reduce all of the values
          if (present(core_up)) then
            n_conn_estimate = estimate_n_connections(ndets,dets_up,dets_dn,wts(:),eps_pt_big,core_up,core_dn,virt_up,virt_dn,active_only)
          else
            n_conn_estimate = estimate_n_connections(ndets,dets_up,dets_dn,wts(:),eps_pt_big)
          endif
          call mpi_allred(n_conn_estimate)
          write (6,'(''If eps_pt_big='',es11.4,'' estimated number of connections to variational wavefn='',i13,'' ='',es9.2)') eps_pt_big, n_conn_estimate, real(n_conn_estimate) ; call flush(6)
        enddo
      else
        !MJO - Each core does a local estimate and then we reduce all of the values
        if (present(core_up)) then
          n_conn_estimate = estimate_n_connections(ndets,dets_up,dets_dn,wts(:),eps_pt_big,core_up,core_dn,virt_up,virt_dn,active_only)
        else
          n_conn_estimate = estimate_n_connections(ndets,dets_up,dets_dn,wts(:),eps_pt_big)
        endif
        call mpi_allred(n_conn_estimate)
        write (6,'(''If eps_pt_big='',es11.4,'' estimated number of connections to variational wavefn='',i13,'' ='',es9.2)') eps_pt_big, n_conn_estimate, real(n_conn_estimate) ; call flush(6)
      endif ! eps_pt_big.le.0_rk

      ! If active spaces are different, then stochastic PT (not semistochastic) requires eps_pt_big=inf, not just eps_pt_big>=eps_var!
      if (.not.active_only) then
!       write(6,'(/,''Different variational and PT spaces, so resetting eps_pt_big from'',es9.2,'' to'',es9.2)') eps_pt_big, 3*eps_pt_big
!       eps_pt_big=3*eps_pt_big
        if (eps_pt_big>=eps_var) then
          eps_pt_big = 1.e99_rk
          write (6,*) "Active spaces for variational and PT stages are different. Setting eps_pt_big=",eps_pt_big ; call flush(6)
        endif
      endif

     !if(wf_file_exist .and. master_core) then
     !  read(1,iostat=ierr) eps_pt_big_read, pt_big_read
     !endif

     !call mpi_bsend(eps_pt_big_read)
     !call mpi_bsend(pt_big_read)
     !call mpi_bsend(ierr)
      if(wf_file_exist .and. ierr.eq.0 .and. eps_pt_big_read.le.eps_pt_big) then ! If eps_pt_big read in is smaller than the one for this run, just use the read in value.
        eps_pt_big=eps_pt_big_read
        pt_big=pt_big_read
        write(6,'(/,''PT_correction, eps_pt_big for short deterministic run read in='',f15.9,es12.4)') pt_big, eps_pt_big
      else
        ! Perform deterministic PT with eps2_big first
        if(eps_pt_big < eps_var) then
          if (eps_pt_big_energy.gt.0) then ! eps_pt_big_energy is initialized to 1 and will always be negative if it is read in
            write (6,'(/,''Computing deterministic part of PT first with eps_pt_big, estimated number of connections='',es12.4,i13,'' ='',es9.2)') eps_pt_big, n_conn_estimate, real(n_conn_estimate)
            call flush(6)
            if (present(core_up)) then
              if(use_mpi) then !MJO TEMP PARALLLEL
                call second_order_pt(ndets, dets_up(1:ndets), dets_dn(1:ndets), wts(1:ndets), diag_elems(1:ndets), var_energy, eps_pt_big, pt_big, ndets_connected,remote_det_map,local_det_map, core_up=core_up, core_dn=core_dn, virt_up=virt_up, virt_dn=virt_dn, active_only=active_only)
              else
                call second_order_pt(ndets, dets_up(1:ndets), dets_dn(1:ndets), wts(1:ndets), diag_elems(1:ndets), var_energy, eps_pt_big, pt_big, ndets_connected, core_up=core_up, core_dn=core_dn, virt_up=virt_up, virt_dn=virt_dn, active_only=active_only)
              endif
            else
              if(use_mpi) then !MJO TEMP PARALLLEL
                call second_order_pt(ndets, dets_up(1:ndets), dets_dn(1:ndets), wts(1:ndets), diag_elems(1:ndets), var_energy, eps_pt_big, pt_big, ndets_connected,remote_det_map,local_det_map)
              else
                call second_order_pt(ndets, dets_up(1:ndets), dets_dn(1:ndets), wts(1:ndets), diag_elems(1:ndets), var_energy, eps_pt_big, pt_big, ndets_connected)
              endif
            endif
            call mpi_allred(ndets_connected)
            write (6,'(/,''PT_correction, eps_pt_big, ndets_connected for short deterministic run='',f15.9,es12.4,i13,'' ='',es9.2)') pt_big, eps_pt_big, ndets_connected, real(ndets_connected)
          else
            pt_big = eps_pt_big_energy
            write (6,'(/,''Used input PT_correction with eps_pt_big for short deterministic run='',f15.9,es12.4)') pt_big, eps_pt_big
          endif
        else
          pt_big=0
          write (6,'(''Short deterministic run not performed because eps_pt_big > eps_var either in the input or because of insufficient memory'')')
        endif ! eps_pt_big < eps_var
       !if(dump_wf_var .and. master_core) then
       !  if(wf_file_exist .and. eps_pt_big_read.gt.eps_pt_big) backspace(1)
       !  write(1) eps_pt_big, pt_big
       !  call flush(1)
       !endif
      endif ! eps_pt_big_read<=eps_pt_big
     !close(1)

      ! Now compute the difference in PT energy, for eps_pt and eps_pt_big, using single-list stochastic PT, with variational dets sampled using the Alias method.

      if(n_mc > 0) then
         write (6,'(/,''Computing PT energy using single-list stochastic method with N_MC='',i8,'' determinants per sample'')') n_mc
      else
         write (6,'(/,''Number of dets per sample, N_MC, will be set as big as available memory allows'')')
      endif
      call flush(6)

      if (present(core_up)) then
        if (use_mpi) then
          call second_order_pt_alias(ndets_global, dets_up_global(1:ndets_global), dets_dn_global(1:ndets_global), wts_global(1:ndets_global), diag_elems_global(1:ndets_global), var_energy, eps_pt, n_mc, target_error, eps_pt_big, pt_diff, pt_diff_std_dev, ndets_connected, pt_big, core_up=core_up, core_dn=core_dn, virt_up=virt_up, virt_dn=virt_dn, active_only=active_only, active_tot_energy=active_tot_energy)
        else
          call second_order_pt_alias(ndets, dets_up(1:ndets), dets_dn(1:ndets), wts(1:ndets), diag_elems(1:ndets), var_energy, eps_pt, n_mc, target_error, eps_pt_big, pt_diff, pt_diff_std_dev, ndets_connected, pt_big, core_up=core_up, core_dn=core_dn, virt_up=virt_up, virt_dn=virt_dn, active_only=active_only, active_tot_energy=active_tot_energy)
        endif
      else
        if (use_mpi) then
          call second_order_pt_alias(ndets_global, dets_up_global(1:ndets_global), dets_dn_global(1:ndets_global), wts_global(1:ndets_global), diag_elems_global(1:ndets_global), var_energy, eps_pt, n_mc, target_error, eps_pt_big, pt_diff, pt_diff_std_dev, ndets_connected, pt_big)
        else
          call second_order_pt_alias(ndets, dets_up(1:ndets), dets_dn(1:ndets), wts(1:ndets), diag_elems(1:ndets), var_energy, eps_pt, n_mc, target_error, eps_pt_big, pt_diff, pt_diff_std_dev, ndets_connected, pt_big)
        endif
      endif
      call mpi_allred(ndets_connected)

    endif ! Semistochastic

  end subroutine do_pt


  subroutine convert_time_symmetrized_to_dets(ndet,dets_up,dets_dn,coeffs,n_states,remote_det_map,local_det_map)
  ! If time_symmetry is used, convert to determinant basis (rather than basis of
  ! symmetrized linear combinations of dets)
  ! A Holmes, 4 Jun 2017

  use chemistry, only : time_sym,z
  use mpi_routines, only : ncores, det_map, det_map_l, master_core_node, mpi_barr_in_node, mpi_distribute_remote_det_map,shmem_allocate,shmem_deallocate, mpi_bsend_between_nodes
  use tools, only : merge_sort2_up_dn

    type(det_map),optional,intent(inout) :: remote_det_map
    type(det_map_l),optional,intent(inout) :: local_det_map
    integer,intent(inout) :: ndet
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),allocatable,dimension(:) :: dets_up,dets_dn
    type(ik_vec),allocatable :: temp_i16_up(:),temp_i16_dn(:)
    type(ik_vec),allocatable :: dets_up_out(:),dets_dn_out(:)
#else
    integer(ik),allocatable,dimension(:) :: dets_up,dets_dn
    integer(ik),allocatable :: temp_i16_up(:),temp_i16_dn(:)
    integer(ik),allocatable :: dets_up_out(:),dets_dn_out(:)
#endif
    real(rk),allocatable,intent(inout) :: coeffs(:,:)
    integer,intent(in) :: n_states ! only need for mpi

    real(rk),allocatable :: coeffs_out(:,:)
    integer,allocatable :: iorder(:),temp_i_2(:)
    integer :: ndets_global, n_allocate
    integer :: i,j

    if (.not.time_sym)  return

    if (present(remote_det_map)) then

      ndets_global = local_det_map%ndets_global_old
      n_allocate = ndets_global
      do i=1,ndets_global
        if (remote_det_map%global_dets_up(i).ne.remote_det_map%global_dets_dn(i))  n_allocate = n_allocate + 1
      enddo
      write (6,'(''Converting from'',i10,'' time-reversal symmetrized linear combinations to'',i10,'' individual determinants'')') ndets_global ,n_allocate ; call flush(6)

      deallocate(dets_up,dets_dn,coeffs)
      if (master_core_node) then
        allocate(dets_up(int(n_allocate,i8b)))
        allocate(dets_dn(int(n_allocate,i8b)))
        allocate(coeffs(int(n_allocate,i8b),n_states))
      endif

      if (master_core_node) then
        j = 0
        do i=1,ndets_global
          j = j + 1
          dets_up(j) = remote_det_map%global_dets_up(i)
          dets_dn(j) = remote_det_map%global_dets_dn(i)
          if (dets_up(j)==dets_dn(j)) then
            coeffs(j,:) = remote_det_map%remote_wts(i,:)
          else
            coeffs(j,:) = inv_sqrt2*(remote_det_map%remote_wts(i,:))
            j = j + 1
            dets_up(j) = remote_det_map%global_dets_dn(i)
            dets_dn(j) = remote_det_map%global_dets_up(i)
            coeffs(j,:) = z*inv_sqrt2*(remote_det_map%remote_wts(i,:))
          endif
        enddo

        ndets_global = n_allocate
        local_det_map%ndets_global = n_allocate

        ! Sort by label for binary search

        allocate(iorder(ndets_global))
        allocate(temp_i16_up((ndets_global+1)/2))
        allocate(temp_i16_dn((ndets_global+1)/2))
        allocate(temp_i_2((ndets_global+1)/2))
        do j=1,ndets_global
          iorder(j)=j
        enddo
        call merge_sort2_up_dn(dets_up(1:ndets_global),dets_dn(1:ndets_global), iorder, ndets_global, temp_i16_up, temp_i16_dn, temp_i_2)
        do i=1,n_states
          coeffs(1:ndets_global,i) = coeffs(iorder(1:ndets_global),i)
        enddo

        deallocate(iorder)
        deallocate(temp_i16_up)
        deallocate(temp_i16_dn)
        deallocate(temp_i_2)

      endif

      call mpi_barr_in_node()

      call shmem_deallocate(remote_det_map%global_dets_up)
      call shmem_deallocate(remote_det_map%global_dets_dn)
      call shmem_deallocate(remote_det_map%remote_wts)

      call shmem_allocate(remote_det_map%global_dets_up,int(n_allocate,i8b))
      if (master_core_node)  remote_det_map%global_dets_up(1:n_allocate) = dets_up(1:n_allocate)
      if (master_core_node)  deallocate(dets_up)

      call shmem_allocate(remote_det_map%global_dets_dn,int(n_allocate,i8b))
      if (master_core_node)  remote_det_map%global_dets_dn(1:n_allocate) = dets_dn(1:n_allocate)
      if (master_core_node)  deallocate(dets_dn)

      call shmem_allocate(remote_det_map%remote_wts,int(n_allocate,i8b),int(n_states,i8b))
      if (master_core_node)  remote_det_map%remote_wts(1:n_allocate,1:n_states) = coeffs(1:n_allocate,1:n_states)
      if (master_core_node)  deallocate(coeffs)

      ! Now distribute

      !MJO Send everyone the old_wts, state by state
      !    FIXME Do all states at once! Write mpi_bsend_real_dparray2
      if (master_core_node) then
        do i=1,n_states
          call mpi_bsend_between_nodes(remote_det_map%remote_wts(:,i))
        enddo
      endif
      call mpi_barr_in_node()
      !MJO Send everyone the remote_det_map and assign the local_det_map pieces
      local_det_map%ndets_global_old = n_allocate
      deallocate(local_det_map%local_indices)
      call mpi_distribute_remote_det_map(remote_det_map,local_det_map)

      ndet = local_det_map%ndets

      !Now, allocate or local dets list
      allocate(dets_up(ndet))
      allocate(dets_dn(ndet))
      allocate(coeffs(ndet,n_states))

      !Slice our our local determinants
      dets_up(:) = remote_det_map%global_dets_up(local_det_map%local_indices)
      dets_dn(:) = remote_det_map%global_dets_dn(local_det_map%local_indices)
      do i=1,n_states
        coeffs(:,i) = remote_det_map%remote_wts(local_det_map%local_indices,i)
      enddo

    else ! serial code

      n_allocate = ndet
      do i=1,ndet
        if (dets_up(i).ne.dets_dn(i))  n_allocate = n_allocate + 1
      enddo
      write (6,'(''Converting from'',i10,'' time-reversal symmetrized linear combinations to'',i10,'' individual determinants'')') ndet, n_allocate ; call flush(6)

      allocate(dets_up_out(n_allocate))
      allocate(dets_dn_out(n_allocate))
      allocate(coeffs_out(n_allocate,n_states))

      j = 0
      do i=1,ndet
        j = j + 1
        dets_up_out(j) = dets_up(i)
        dets_dn_out(j) = dets_dn(i)
        if (dets_up_out(j)==dets_dn_out(j)) then
          coeffs_out(j,:) = coeffs(i,:)
        else
          coeffs_out(j,:) = inv_sqrt2*coeffs(i,:)
          j = j + 1
          dets_up_out(j) = dets_dn(i)
          dets_dn_out(j) = dets_up(i)
          coeffs_out(j,:) = z*inv_sqrt2*coeffs(i,:)
        endif
      enddo

      ndet = n_allocate
      deallocate(dets_up,dets_dn,coeffs)

      allocate(dets_up(ndet))
      dets_up(1:ndet) = dets_up_out(1:ndet)
      deallocate(dets_up_out)

      allocate(dets_dn(ndet))
      dets_dn(1:ndet) = dets_dn_out(1:ndet)
      deallocate(dets_dn_out)

      allocate(coeffs(ndet,n_states))
      coeffs(1:ndet,1:n_states) = coeffs_out(1:ndet,1:n_states)
      deallocate(coeffs_out)

        ! Sort by label for binary search

        allocate(iorder(ndet))
        allocate(temp_i16_up((ndet+1)/2))
        allocate(temp_i16_dn((ndet+1)/2))
        allocate(temp_i_2((ndet+1)/2))
        do j=1,ndet
          iorder(j)=j
        enddo
        call merge_sort2_up_dn(dets_up(1:ndet),dets_dn(1:ndet), iorder, ndet, temp_i16_up, temp_i16_dn, temp_i_2)
        do i=1,n_states
          coeffs(1:ndet,i) = coeffs(iorder(1:ndet),i)
        enddo

        deallocate(iorder)
        deallocate(temp_i16_up)
        deallocate(temp_i16_dn)
        deallocate(temp_i_2)

    endif

  end subroutine convert_time_symmetrized_to_dets


end module hci
