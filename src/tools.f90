module tools
  use types, only: ik, ik_vec, i8b, i16b, rk
  use common_run, only : diag_elem_info ! type
#ifdef NUM_ORBITALS_GT_127
  use overload
#endif
  implicit none
  save
  private
  public ::  n_choose_k, random_int, count_bits_naive, count_bits_imp1, round_r, round_i, alloc, &!subroutines
           & merge_sort2_up_dn, merge_sort_real, merge2_up_dn, merge_original_with_spawned3, &
           & get_free_memory, verbose_allocate, permutation_factor, &
           & permutation_factor2, count_excitations, print_excitation_levels_and_wts, &
           & matrix_inversion, do_once, do_n_times, choose_det_with_prob_prop_to_abs_wt, &
           & sort_and_merge, sort_and_merge_count_repeats, welford

  interface merge_original_with_spawned3
     module procedure merge_original_with_spawned3_nowts, merge_original_with_spawned3_nowts_replace, merge_original_with_spawned3_num_denom, merge_original_with_spawned3_num_denom_diag_elems_info, merge_original_with_spawned3_num, merge_original_with_spawned3_5items, merge_original_with_spawned3_num_denom_replace, merge_original_with_spawned3_num_denom_diag_elems_info_replace, merge_original_with_spawned3_5items_replace, merge_original_with_spawned3_add_wts
  end interface merge_original_with_spawned3

  interface sort_and_merge
    module procedure sort_and_merge_int, sort_and_merge_det, sort_and_merge_int_rk
  end interface sort_and_merge

! merge_original_with_spawned3 can now be called with any of the following lists of inputs:
!   (old_up,old_dn,nwalk,new_up,new_dn)  -  Merge dets, allocate output to size nwalk
!   (nwalk,dets_up,dets_dn)  -  Merge dets, replace input (don't reallocate)
!   (old_up,old_dn,nwalk,new_up,new_dn,old_e_mix_num,old_e_mix_den,e_mix_num,e_mix_den)  -  Merge dets and local energies, allocate output to size nwalk
!   (old_up,old_dn,nwalk,new_up,new_dn,old_e_mix_num,old_e_mix_den,old_diag_elems_info,e_mix_num,e_mix_den,diag_elems_info)  -  Merge dets and local energies and diag_elems_info, allocate output to size nwalk
!   (old_up,old_dn,nwalk,new_up,new_dn,old_e_mix_num,e_mix_num)  -  Merge dets and local energies, allocate output to size nwalk
!   (nwalk,dets_up,dets_dn,e_mix_num,e_mix_den)  -  Merge dets and local energies, replace input (don't reallocate)
!   (nwalk,dets_up,dets_dn,e_mix_num,e_mix_den,diag_elems_info)  -  Merge dets and local energies and diag_elems_info, replace input (don't reallocate)
!   (dets_up,dets_dn,walk_wt,nwalk)  -  Merge dets and weights, replace input by reallocating it to size nwalk (only used for chemistry with time_sym on)

  interface n_choose_k
     module procedure n_choose_k_int4, n_choose_k_int16
  end interface

  interface verbose_allocate
     module procedure verbose_allocate_int,verbose_allocate_ik,verbose_allocate_rk,verbose_allocate_rk2,verbose_allocate_int_test
  end interface

  interface alloc
     module procedure alloc_i, alloc_ik, alloc_rk_1, alloc_rk_2, alloc_diag_elem_info
  end interface

contains
  !===========================================================================
  function n_choose_k_int4(n,k)
    !---------------------------------------------------------------------------
    ! Description : Generate binomial coefficient
    ! Created     : F. Petruzielo, 5 Nov 2010
    !---------------------------------------------------------------------------
    implicit none

    !dummy arguments
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer(i16b) :: n_choose_k_int4

    !local variables
    integer :: i
    real(rk) :: log_n_factorial, log_k_factorial, log_n_minus_k_factorial

    if (n < k) then
      write(6,'(''n < k in binomial coefficient!, n,k='',9i5)') n,k
      stop 'n < k in binomial coefficient n_choose_k!'
    endif
    log_n_factorial = 0._rk
    log_k_factorial = 0._rk
    log_n_minus_k_factorial = 0._rk
    do i = 2, n
       log_n_factorial = log_n_factorial + log(i*1._rk)
    enddo

    do i = 2, k
       log_k_factorial = log_k_factorial + log(i*1._rk)
    enddo

    do i = 2, n-k
       log_n_minus_k_factorial = log_n_minus_k_factorial + log(i*1._rk)
    enddo

    n_choose_k_int4 = int(exp(log_n_factorial - log_k_factorial - log_n_minus_k_factorial)+0.5,i16b)

  end function n_choose_k_int4
  !===========================================================================
  !===========================================================================
  function n_choose_k_int16(n,k)
    !---------------------------------------------------------------------------
    ! Description : Generate binomial coefficient
    ! Created     : F. Petruzielo, 5 Nov 2010
    !---------------------------------------------------------------------------
    implicit none

    !dummy arguments
    integer(ik), intent(in) :: n
    integer(ik), intent(in) :: k
    integer(i16b) :: n_choose_k_int16

    !local variables
    integer(ik) :: i
    real(rk) :: log_n_factorial, log_k_factorial, log_n_minus_k_factorial

    if (n < k) then
      write(6,'(''n < k in binomial coefficient!, n,k='',9i5)') n,k
      stop 'n < k in binomial coefficient n_choose_k!'
    endif
    log_n_factorial = 0._rk
    log_k_factorial = 0._rk
    log_n_minus_k_factorial = 0._rk
    do i = 2, n
       log_n_factorial = log_n_factorial + log(i*1._rk)
    enddo

    do i = 2, k
       log_k_factorial = log_k_factorial + log(i*1._rk)
    enddo

    do i = 2, n-k
       log_n_minus_k_factorial = log_n_minus_k_factorial + log(i*1._rk)
    enddo

    n_choose_k_int16 = int(exp(log_n_factorial - log_k_factorial - log_n_minus_k_factorial)+0.5,i16b)

  end function n_choose_k_int16
  !===========================================================================

  !===========================================================================
  function random_int(n)
    !---------------------------------------------------------------------------
    ! Description : Generate random integer in range 1,...,n
    !
    ! Created     : F. Petruzielo, 2 Nov 2010
    !---------------------------------------------------------------------------
    implicit none

    !dummy arguments
    integer, intent(in) :: n
    integer :: random_int

    !local variables
     real(rk) :: rannyu
    !real :: r

   !call random_number(r) !replace this with call to Dr. Umrigar random number generator
    random_int = int(n * rannyu()) + 1

  end function random_int
  !===========================================================================
    subroutine count_bits_naive(bit_pattern, n_bits, n_set)
    ! Description:  Naive bit counting
    ! Created    :  Frank Petruzielo , added in sqmc by Hitesh on April 18 2012
    implicit none
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: bit_pattern
#else
    integer(ik), intent(in) :: bit_pattern
#endif
    integer, intent(in) :: n_bits
    integer, intent(out) :: n_set
    integer :: i

    n_set = 0
    do i = 1, n_bits
       if (btest(bit_pattern,i-1)) then
          n_set = n_set + 1
       endif
    enddo

  end subroutine count_bits_naive

  !===========================================================================
  subroutine count_bits_imp1(bit_pattern, n_set)
    ! Description:  Intelligent bit counting
    ! Created    :  Frank Petruzielo , added in sqmc by Hitesh on April 18 2012

    implicit none
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(inout) :: bit_pattern
#else
    integer(ik), intent(inout) :: bit_pattern
#endif
    integer, intent(out) :: n_set
    !From wikipedia: The above implementations have the best worst-case behavior of any known algorithm.
    !However, when a value is expected to have few nonzero bits, it may instead be more efficient to use algorithms that count these bits one at a time.
    !As Wegner (1960) described,[4] the bitwise and of x with x − 1 differs from x only in zeroing out the least significant nonzero bit: subtracting 1 changes the rightmost string of 0s to 1s, and changes the rightmost 1 to a 0.
    !If x originally had n bits that were 1, then after only n iterations of this operation, x will be reduced to zero.
    !The following implementation is based on this principle.
    !This is better when most bits in x are 0
    !It uses 3 arithmetic operations and one comparison/branch per "1" bit in x.
    n_set = 0
    do while (bit_pattern .ne. 0)
       n_set = n_set + 1
       bit_pattern = iand(bit_pattern, bit_pattern - 1)
    enddo

  end subroutine count_bits_imp1

!------------------------------------------------------------------------------------------------------------
recursive subroutine merge_sort2_up_dn(key_up, key_dn, iorder, nwalk, temp_i16_up, temp_i16_dn, temp_i_2)
! ===========================================================================================================
! Description   : Merge sort in ascending order.
!               : It sorts nwalk items in key.  Index is then used outside this routine to sort auxilliary items.
!               : temp_i16 is an auxilliary array of length (nwalk+1)/2 of same type as key.
!               : temp_i_2 is an auxilliary array of length (nwalk+1)/2 of same type as iorder.
!               : In the present version key is type(ik_vec) and iorder is integer.
!               : temp_i16 and temp_i_2 are passed in to avoid overhead of automatic arrays or allocations.
! Adapted by    : Cyrus Umrigar
! Modified by   : A Holmes, 29 Jan 2015. Change sorting order so that
!                 HF comes first even if using norb = number of bits
!                 allowed in signed integer type.
! -----------------------------------------------------------------------------------------------------------
  use types, only: ik
  implicit none
  integer, intent(inout) :: iorder(nwalk)
  integer, intent(in) :: nwalk
  integer,     dimension((nwalk+1)/2), intent (out) :: temp_i_2                 ! temporary array actually neither in nor out
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), intent(inout) :: key_up(nwalk), key_dn(nwalk)
  type(ik_vec), dimension((nwalk+1)/2), intent (out) :: temp_i16_up, temp_i16_dn ! temporary array actually neither in nor out
! Local variables
  type(ik_vec)  t0
#else
  integer(ik), intent(inout) :: key_up(nwalk), key_dn(nwalk)
  integer(ik), dimension((nwalk+1)/2), intent (out) :: temp_i16_up, temp_i16_dn ! temporary array actually neither in nor out
! Local variables
  integer(ik)  t0
#endif
  integer :: na,nb,i

  if (nwalk < 2) return
  if (nwalk == 2) then
    if (abs(key_up(1)) > abs(key_up(2))) then ! If key_up's in wrong order, sort on key_up
      t0 = key_up(1)
      key_up(1) = key_up(2)
      key_up(2) = t0
      t0 = key_dn(1)
      key_dn(1) = key_dn(2)
      key_dn(2) = t0
      i = iorder(1)
      iorder(1) = iorder(2)
      iorder(2) = i
    endif
    if (key_up(1) == key_up(2) .and. abs(key_dn(1)) > abs(key_dn(2))) then ! If key_up's are equal, sort on key_dn
      t0 = key_dn(1)
      key_dn(1) = key_dn(2)
      key_dn(2) = t0
      i = iorder(1)
      iorder(1) = iorder(2)
      iorder(2) = i
    endif
    return
  endif

  na=(nwalk+1)/2
  nb=nwalk-na
  call merge_sort2_up_dn(key_up, key_dn, iorder, na, temp_i16_up, temp_i16_dn, temp_i_2)
  call merge_sort2_up_dn(key_up(na+1), key_dn(na+1), iorder(na+1), nb, temp_i16_up, temp_i16_dn, temp_i_2)

  if (abs(key_up(na)) > abs(key_up(na+1)) .or. (key_up(na) == key_up(na+1) .and. abs(key_dn(na)) > abs(key_dn(na+1)))) then
    temp_i16_up(1:na)=key_up(1:na)
    temp_i16_dn(1:na)=key_dn(1:na)
    temp_i_2(1:na)=iorder(1:na)
    call merge2_up_dn(temp_i16_up, temp_i16_dn, temp_i_2, na, key_up(na+1), key_dn(na+1), iorder(na+1), nb, key_up, key_dn, iorder, nwalk)
  endif

  return
end subroutine merge_sort2_up_dn
!-----------------------------------------------------------------------------------------------

recursive subroutine merge_sort_real(key_real, iorder, nelts, temp_real, temp_i_2)
! ===========================================================================================================
! Description   : Merge sort in ascending order.
!               : It sorts nwalk items in key which is real.  Index is then used outside this routine to sort auxilliary items.
!               : temp_real is an auxilliary array of length (nwalk+1)/2 of same type as key.
!               : temp_i_2 is an auxilliary array of length (nwalk+1)/2 of same type as iorder.
!               : In the present version key is real(rk) and iorder is integer.
!               : temp_i16 and temp_i_2 are passed in to avoid overhead of automatic arrays or allocations.
! Adapted by    : Cyrus Umrigar
! -----------------------------------------------------------------------------------------------------------
  use types, only: ik
  implicit none
  real(rk), intent(inout) :: key_real(nelts)
  integer, intent(inout) :: iorder(nelts)
  integer, intent(in) :: nelts
  real(rk), dimension((nelts+1)/2), intent (out) :: temp_real
  integer,  dimension((nelts+1)/2), intent (out) :: temp_i_2     ! temporary array actually neither in nor out
! Local variables
  real(rk)  t0
  integer :: na,nb,i

  if (nelts < 2) return
  if (nelts == 2) then
    if (key_real(1) > key_real(2)) then ! If key_real's in wrong order, sort on key_real
      t0 = key_real(1)
      key_real(1) = key_real(2)
      key_real(2) = t0
      i = iorder(1)
      iorder(1) = iorder(2)
      iorder(2) = i
    endif
    return
  endif

  na=(nelts+1)/2
  nb=nelts-na
  call merge_sort_real(key_real, iorder, na, temp_real, temp_i_2)
  call merge_sort_real(key_real(na+1), iorder(na+1), nb, temp_real, temp_i_2)

  if (key_real(na) > key_real(na+1)) then
    temp_real(1:na)=key_real(1:na)
    temp_i_2(1:na)=iorder(1:na)
    call merge_real(temp_real, temp_i_2, na, key_real(na+1), iorder(na+1), nb, key_real, iorder, nelts)
  endif

  return
end subroutine merge_sort_real
!-----------------------------------------------------------------------------------------------


subroutine merge2_up_dn(a_up,a_dn,a2,na, b_up,b_dn,b2,nb, c_up,c_dn,c2,nc)
! ===============================================================================================
! Description   : Called by merge_sort2_up_dn
! Adapted by    : Cyrus Umrigar
! Comments      : Pulled out from semistoch_mod.f90 and chemistry.f90 and put in tools.f90 for general purpose use
! -----------------------------------------------------------------------------------------------
  use types, only: ik
  implicit none
  integer, intent(in) :: na,nb,nc                                      ! na+nb = nc
  integer    , intent(inout) :: a2(na), c2(nc)                         ! b2   overlays c2(na+1:nc)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), intent(inout) :: a_up(na), a_dn(na), c_up(nc), c_dn(nc) ! b_up overlays c_up(na+1:nc), and similarly for b_dn, c_dn
  type(ik_vec), intent(in)    :: b_up(nb), b_dn(nb)
#else
  integer(ik), intent(inout) :: a_up(na), a_dn(na), c_up(nc), c_dn(nc) ! b_up overlays c_up(na+1:nc), and similarly for b_dn, c_dn
  integer(ik), intent(in)    :: b_up(nb), b_dn(nb)
#endif
  integer    , intent(in)    :: b2(nb)
  integer :: i,j,k

  i = 1; j = 1; k = 1;
  do while(i <= abs(na) .and. j <= abs(nb))
    if (abs(a_up(i)) < abs(b_up(j)) .or. (a_up(i) == b_up(j) .and. abs(a_dn(i)) <= abs(b_dn(j)))) then
      c_up(k) = a_up(i)
      c_dn(k) = a_dn(i)
      c2(k) = a2(i)
      i = i+1
    else
      c_up(k) = b_up(j)
      c_dn(k) = b_dn(j)
      c2(k) = b2(j)
      j = j+1
    endif
    k = k + 1
  enddo

  do while (i <= abs(na))
    c_up(k) = a_up(i)
    c_dn(k) = a_dn(i)
    c2(k) = a2(i)
    i = i + 1
    k = k + 1
  enddo

  return
end subroutine merge2_up_dn
!----------------------------------------------------------------------------------------------------------------------

subroutine merge_real(a,a2,na,b,b2,nb,c,c2,nc)
! ===============================================================================================
! Description   : Called by merge_sort2_up_dn
! Adapted by    : Cyrus Umrigar
! Comments      : Pulled out from semistoch_mod.f90 and chemistry.f90 and put in tools.f90 for general purpose use
! -----------------------------------------------------------------------------------------------
  use types, only: ik
  implicit none
  integer, intent(in) :: na,nb,nc                                      ! na+nb = nc
  real(rk), intent(inout) :: a(na), c(nc) ! b_up overlays c_up(na+1:nc), and similarly for b_dn, c_dn
  integer    , intent(inout) :: a2(na), c2(nc)                         ! b2   overlays c2(na+1:nc)
  real(rk), intent(in)    :: b(nb)
  integer    , intent(in)    :: b2(nb)
  integer :: i,j,k

  i = 1; j = 1; k = 1;
  do while(i <= na .and. j <= nb)
    if (a(i) < b(j)) then
      c(k) = a(i)
      c2(k) = a2(i)
      i = i+1
    else
      c(k) = b(j)
      c2(k) = b2(j)
      j = j+1
    endif
    k = k + 1
  enddo

  do while (i <= na)
    c(k) = a(i)
    c2(k) = a2(i)
    i = i + 1
    k = k + 1
  enddo

  return
end subroutine merge_real
!----------------------------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
 subroutine merge_original_with_spawned3_nowts(old_up, old_dn, nwalk, new_up, new_dn)
! ==============================================================================
! Description   : Merge walkers
!                 Creates new list of dets, allocated to size nwalk
! Author        : Cyrus Umrigar
! Modified      : Adam Holmes, Pulled out of chemistry.f90 by HJC
! ------------------------------------------------------------------------------
  use types, only: ik
  implicit none
  integer, intent(inout) :: nwalk
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: old_up(nwalk),old_dn(nwalk)
  type(ik_vec),intent(out) :: new_up(:),new_dn(:)
  type(ik_vec),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#else
  integer(ik),intent(in) :: old_up(nwalk),old_dn(nwalk)
  integer(ik),intent(out) :: new_up(:),new_dn(:)
  integer(ik),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#endif
  integer iwalk,nshift

  allocate(walk_dets_up(nwalk))
  allocate(walk_dets_dn(nwalk))

  walk_dets_up(1:nwalk) = old_up(1:nwalk)
  walk_dets_dn(1:nwalk) = old_dn(1:nwalk)

  nshift=0
  do iwalk=2,nwalk
    if(walk_dets_up(iwalk).eq.walk_dets_up(iwalk-nshift-1) .and. walk_dets_dn(iwalk).eq.walk_dets_dn(iwalk-nshift-1)) then ! This walker is on same det as previous walker
      nshift=nshift+1
    else                                                                ! This walker is on different det from previous walker. Complete working on previous det before working on this one.
! Done working on previous det.  Now move position of current one.
      walk_dets_up(iwalk-nshift)=walk_dets_up(iwalk)
      walk_dets_dn(iwalk-nshift)=walk_dets_dn(iwalk)
    endif
  enddo
! Fix up last det, just as we did above after "else".
  nwalk=nwalk-nshift

  new_up(1:nwalk) = walk_dets_up(1:nwalk)
  new_dn(1:nwalk) = walk_dets_dn(1:nwalk)

  deallocate(walk_dets_up,walk_dets_dn)

 end subroutine merge_original_with_spawned3_nowts
!-------------------------------------------------------------------------------
 subroutine merge_original_with_spawned3_nowts_replace(nwalk, dets_up, dets_dn)
! ==============================================================================
! Description   : Merge walkers
!                 Replaces current list of dets (keeps same size)
! Author        : Cyrus Umrigar
! Modified      : Adam Holmes, Pulled out of chemistry.f90 by HJC
! ------------------------------------------------------------------------------
  use types, only: ik
  implicit none
  integer, intent(inout) :: nwalk
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(inout) :: dets_up(:),dets_dn(:)
  type(ik_vec),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#else
  integer(ik),intent(inout) :: dets_up(:),dets_dn(:)
  integer(ik),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#endif
  integer iwalk,nshift

  allocate(walk_dets_up(nwalk))
  allocate(walk_dets_dn(nwalk))

  walk_dets_up(1:nwalk) = dets_up(1:nwalk)
  walk_dets_dn(1:nwalk) = dets_dn(1:nwalk)

  nshift=0
  do iwalk=2,nwalk
    if(walk_dets_up(iwalk).eq.walk_dets_up(iwalk-nshift-1) .and. walk_dets_dn(iwalk).eq.walk_dets_dn(iwalk-nshift-1)) then ! This walker is on same det as previous walker
      nshift=nshift+1
    else                                                                ! This walker is on different det from previous walker. Complete working on previous det before working on this one.
! Done working on previous det.  Now move position of current one.
      walk_dets_up(iwalk-nshift)=walk_dets_up(iwalk)
      walk_dets_dn(iwalk-nshift)=walk_dets_dn(iwalk)
    endif
  enddo
! Fix up last det, just as we did above after "else".
  nwalk=nwalk-nshift

  dets_up(1:nwalk) = walk_dets_up(1:nwalk)
  dets_dn(1:nwalk) = walk_dets_dn(1:nwalk)

  deallocate(walk_dets_up,walk_dets_dn)

 end subroutine merge_original_with_spawned3_nowts_replace
!-------------------------------------------------------------------------------
subroutine merge_original_with_spawned3_num_denom(old_up, old_dn, nwalk, new_up, new_dn, old_e_mix_num, old_e_mix_den, e_mix_num, e_mix_den)
! =====================================================================================================================
! Description   : Merge walkers
! Author        : Cyrus Umrigar
! Modified      : A Holmes
!                 Creates new lists allocated to size nwalk for e_mix_num,den as well as dets
! Comments      : Pulled out from semistoch_mod.f90 by HJC and used in a common area
! Last edited   : A Holmes, 8 Aug 2012. Made e_mix_den,old_e_mix_den optional.
! ---------------------------------------------------------------------------------------------------------------------

  implicit none
  integer, intent(inout) :: nwalk
  real(rk),intent(in) :: old_e_mix_num(nwalk)
  real(rk),intent(in) :: old_e_mix_den(nwalk)
  real(rk),allocatable,intent(out) :: e_mix_num(:)
  real(rk),allocatable,intent(out) :: e_mix_den(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: old_up(nwalk),old_dn(nwalk)
  type(ik_vec),allocatable,intent(out) :: new_up(:),new_dn(:)
  type(ik_vec),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#else
  integer(ik),intent(in) :: old_up(nwalk),old_dn(nwalk)
  integer(ik),allocatable,intent(out) :: new_up(:),new_dn(:)
  integer(ik),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#endif
  integer iwalk,nshift
  real(rk),allocatable :: tmp_num(:),tmp_den(:)

  allocate(walk_dets_up(nwalk))
  allocate(walk_dets_dn(nwalk))

  walk_dets_up = old_up
  walk_dets_dn = old_dn

  allocate(tmp_num(nwalk))
  tmp_num = old_e_mix_num
  allocate(tmp_den(nwalk))
  tmp_den = old_e_mix_den

  nshift=0
  do iwalk=2,nwalk
    if(walk_dets_up(iwalk).eq.walk_dets_up(iwalk-nshift-1) .and. walk_dets_dn(iwalk).eq.walk_dets_dn(iwalk-nshift-1)) then ! This walker is on same det as previous walker
      nshift=nshift+1
      tmp_num(iwalk-nshift)=tmp_num(iwalk-nshift)+tmp_num(iwalk)
      tmp_den(iwalk-nshift)=tmp_den(iwalk-nshift)+tmp_den(iwalk)
    else                                                                ! This walker is on different det from previous walker. Complete working on previous det before working on this one.
! Done working on previous det.  Now move position of current one.
      walk_dets_up(iwalk-nshift)=walk_dets_up(iwalk)
      walk_dets_dn(iwalk-nshift)=walk_dets_dn(iwalk)
      tmp_num(iwalk-nshift)=tmp_num(iwalk)
      tmp_den(iwalk-nshift)=tmp_den(iwalk)
    endif
  enddo
! Fix up last det, just as we did above after "else".
  nwalk=nwalk-nshift

  allocate(new_up(nwalk))
  allocate(new_dn(nwalk))

  new_up(1:nwalk) = walk_dets_up(1:nwalk)
  new_dn(1:nwalk) = walk_dets_dn(1:nwalk)

  deallocate(walk_dets_up,walk_dets_dn)

  allocate(e_mix_num(nwalk))
  e_mix_num(1:nwalk) = tmp_num(1:nwalk)
  deallocate(tmp_num)
  allocate(e_mix_den(nwalk))
  e_mix_den(1:nwalk) = tmp_den(1:nwalk)
  deallocate(tmp_den)

end subroutine merge_original_with_spawned3_num_denom
!---------------------------------------------------------------------------------------------------------------------

subroutine merge_original_with_spawned3_num_denom_diag_elems_info(old_up, old_dn, nwalk, new_up, new_dn, old_e_mix_num, old_e_mix_den, old_diag_elems_info, e_mix_num, e_mix_den, diag_elems_info)
! =====================================================================================================================
! Description   : Merge walkers
! Author        : Cyrus Umrigar
! Modified      : A Holmes
!                 Creates new lists allocated to size nwalk for e_mix_num,den as well as dets
! Comments      : Pulled out from semistoch_mod.f90 by HJC and used in a common area
! Last edited   : A Holmes, 10 Mar 2016. Added diag_elems_info
! ---------------------------------------------------------------------------------------------------------------------

  use common_run, only : diag_elem_info

  implicit none
  integer, intent(inout) :: nwalk
  real(rk),intent(in) :: old_e_mix_num(nwalk)
  real(rk),intent(in) :: old_e_mix_den(nwalk)
  type(diag_elem_info),intent(in) :: old_diag_elems_info(:)
  real(rk),allocatable,intent(out) :: e_mix_num(:)
  real(rk),allocatable,intent(out) :: e_mix_den(:)
  type(diag_elem_info),allocatable,intent(out) :: diag_elems_info(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: old_up(nwalk),old_dn(nwalk)
  type(ik_vec),allocatable,intent(out) :: new_up(:),new_dn(:)
  type(ik_vec),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#else
  integer(ik),intent(in) :: old_up(nwalk),old_dn(nwalk)
  integer(ik),allocatable,intent(out) :: new_up(:),new_dn(:)
  integer(ik),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#endif
  integer iwalk,nshift
  real(rk),allocatable :: tmp_num(:),tmp_den(:)
  type(diag_elem_info),allocatable :: tmp_diag_elems_info(:)

  allocate(walk_dets_up(nwalk))
  allocate(walk_dets_dn(nwalk))

  walk_dets_up = old_up
  walk_dets_dn = old_dn

  allocate(tmp_num(nwalk))
  tmp_num = old_e_mix_num
  allocate(tmp_den(nwalk))
  tmp_den = old_e_mix_den
  allocate(tmp_diag_elems_info(nwalk))
  tmp_diag_elems_info = old_diag_elems_info

  nshift=0
  do iwalk=2,nwalk
    if(walk_dets_up(iwalk).eq.walk_dets_up(iwalk-nshift-1) .and. walk_dets_dn(iwalk).eq.walk_dets_dn(iwalk-nshift-1)) then ! This walker is on same det as previous walker
      nshift=nshift+1
      tmp_num(iwalk-nshift)=tmp_num(iwalk-nshift)+tmp_num(iwalk)
      tmp_den(iwalk-nshift)=tmp_den(iwalk-nshift)+tmp_den(iwalk)
      tmp_diag_elems_info(iwalk-nshift)=tmp_diag_elems_info(iwalk)
    else                                                                ! This walker is on different det from previous walker. Complete working on previous det before working on this one.
! Done working on previous det.  Now move position of current one.
      walk_dets_up(iwalk-nshift)=walk_dets_up(iwalk)
      walk_dets_dn(iwalk-nshift)=walk_dets_dn(iwalk)
      tmp_num(iwalk-nshift)=tmp_num(iwalk)
      tmp_den(iwalk-nshift)=tmp_den(iwalk)
      tmp_diag_elems_info(iwalk-nshift)=tmp_diag_elems_info(iwalk)
    endif
  enddo
! Fix up last det, just as we did above after "else".
  nwalk=nwalk-nshift

  allocate(new_up(nwalk))
  allocate(new_dn(nwalk))

  new_up(1:nwalk) = walk_dets_up(1:nwalk)
  new_dn(1:nwalk) = walk_dets_dn(1:nwalk)

  deallocate(walk_dets_up,walk_dets_dn)

  allocate(e_mix_num(nwalk))
  e_mix_num(1:nwalk) = tmp_num(1:nwalk)
  deallocate(tmp_num)
  allocate(e_mix_den(nwalk))
  e_mix_den(1:nwalk) = tmp_den(1:nwalk)
  deallocate(tmp_den)
  allocate(diag_elems_info(nwalk))
  diag_elems_info(1:nwalk) = tmp_diag_elems_info(1:nwalk)
  deallocate(tmp_diag_elems_info)

end subroutine merge_original_with_spawned3_num_denom_diag_elems_info
!---------------------------------------------------------------------------------------------------------------------

subroutine merge_original_with_spawned3_5items(old_up, old_dn, nwalk, new_up, new_dn, old_e_mix_num, old_e_mix_den, old_diag_elems_info, old_term1_big, old_term2_big, e_mix_num, e_mix_den, diag_elems_info, term1_big, term2_big)
! =====================================================================================================================
! Description   : Merge walkers
! Author        : Cyrus Umrigar
! Modified      : A Holmes
!                 Creates new lists allocated to size nwalk for e_mix_num,den as well as dets
! Comments      : Pulled out from semistoch_mod.f90 by HJC and used in a common area
! Last edited   : A Holmes, 10 Mar 2016. Added diag_elems_info
! ---------------------------------------------------------------------------------------------------------------------

  use common_run, only : diag_elem_info

  implicit none
  integer, intent(inout) :: nwalk
  real(rk),intent(in) :: old_e_mix_num(nwalk)
  real(rk),intent(in) :: old_e_mix_den(nwalk)
  type(diag_elem_info),intent(in) :: old_diag_elems_info(:)
  real(rk),allocatable,intent(out) :: e_mix_num(:)
  real(rk),allocatable,intent(out) :: e_mix_den(:)
  type(diag_elem_info),allocatable,intent(out) :: diag_elems_info(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: old_up(nwalk),old_dn(nwalk)
  type(ik_vec),allocatable,intent(out) :: new_up(:),new_dn(:)
  type(ik_vec),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#else
  integer(ik),intent(in) :: old_up(nwalk),old_dn(nwalk)
  integer(ik),allocatable,intent(out) :: new_up(:),new_dn(:)
  integer(ik),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#endif
  real(rk),intent(in) :: old_term1_big(:)
  real(rk),intent(in) :: old_term2_big(:)
  real(rk),allocatable,intent(out) :: term1_big(:)
  real(rk),allocatable,intent(out) :: term2_big(:)
  integer iwalk,nshift
  real(rk),allocatable :: tmp_num(:),tmp_den(:)
  real(rk),allocatable :: tmp_1big(:),tmp_2big(:)
  type(diag_elem_info),allocatable :: tmp_diag_elems_info(:)

  allocate(walk_dets_up(nwalk))
  allocate(walk_dets_dn(nwalk))

  walk_dets_up = old_up
  walk_dets_dn = old_dn

  allocate(tmp_num(nwalk))
  tmp_num = old_e_mix_num
  allocate(tmp_den(nwalk))
  tmp_den = old_e_mix_den
  allocate(tmp_diag_elems_info(nwalk))
  tmp_diag_elems_info = old_diag_elems_info
  allocate(tmp_1big(nwalk))
  tmp_1big = old_term1_big
  allocate(tmp_2big(nwalk))
  tmp_2big = old_term2_big

  nshift=0
  do iwalk=2,nwalk
    if(walk_dets_up(iwalk).eq.walk_dets_up(iwalk-nshift-1) .and. walk_dets_dn(iwalk).eq.walk_dets_dn(iwalk-nshift-1)) then ! This walker is on same det as previous walker
      nshift=nshift+1
      tmp_num(iwalk-nshift)=tmp_num(iwalk-nshift)+tmp_num(iwalk)
      tmp_den(iwalk-nshift)=tmp_den(iwalk-nshift)+tmp_den(iwalk)
      tmp_diag_elems_info(iwalk-nshift)=tmp_diag_elems_info(iwalk)
      tmp_1big(iwalk-nshift)=tmp_1big(iwalk-nshift)+tmp_1big(iwalk)
      tmp_2big(iwalk-nshift)=tmp_2big(iwalk-nshift)+tmp_2big(iwalk)
    else                                                                ! This walker is on different det from previous walker. Complete working on previous det before working on this one.
! Done working on previous det.  Now move position of current one.
      walk_dets_up(iwalk-nshift)=walk_dets_up(iwalk)
      walk_dets_dn(iwalk-nshift)=walk_dets_dn(iwalk)
      tmp_num(iwalk-nshift)=tmp_num(iwalk)
      tmp_den(iwalk-nshift)=tmp_den(iwalk)
      tmp_diag_elems_info(iwalk-nshift)=tmp_diag_elems_info(iwalk)
      tmp_1big(iwalk-nshift)=tmp_1big(iwalk)
      tmp_2big(iwalk-nshift)=tmp_2big(iwalk)
    endif
  enddo
! Fix up last det, just as we did above after "else".
  nwalk=nwalk-nshift

  allocate(new_up(nwalk))
  allocate(new_dn(nwalk))

  new_up(1:nwalk) = walk_dets_up(1:nwalk)
  new_dn(1:nwalk) = walk_dets_dn(1:nwalk)

  deallocate(walk_dets_up,walk_dets_dn)

  allocate(e_mix_num(nwalk))
  e_mix_num(1:nwalk) = tmp_num(1:nwalk)
  deallocate(tmp_num)
  allocate(e_mix_den(nwalk))
  e_mix_den(1:nwalk) = tmp_den(1:nwalk)
  deallocate(tmp_den)
  allocate(diag_elems_info(nwalk))
  diag_elems_info(1:nwalk) = tmp_diag_elems_info(1:nwalk)
  deallocate(tmp_diag_elems_info)
  allocate(term1_big(nwalk))
  term1_big(1:nwalk) = tmp_1big(1:nwalk)
  deallocate(tmp_1big)
  allocate(term2_big(nwalk))
  term2_big(1:nwalk) = tmp_2big(1:nwalk)
  deallocate(tmp_2big)

end subroutine merge_original_with_spawned3_5items
!---------------------------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine merge_original_with_spawned3_num(old_up, old_dn, nwalk, new_up, new_dn, old_e_mix_num, e_mix_num)
! =====================================================================================================================
! Description   : Merge walkers
! Author        : Cyrus Umrigar
! Modified      : A Holmes
!                 Creates new lists allocated to size nwalk for e_mix_num,den as well as dets
! Comments      : Pulled out from semistoch_mod.f90 by HJC and used in a common area
! Last edited   : A Holmes, 8 Aug 2012. Made e_mix_den,old_e_mix_den optional.
! ---------------------------------------------------------------------------------------------------------------------

  implicit none
  integer, intent(inout) :: nwalk
  real(rk),intent(in) :: old_e_mix_num(nwalk)
  real(rk),allocatable,intent(out) :: e_mix_num(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: old_up(nwalk),old_dn(nwalk)
  type(ik_vec),allocatable,intent(out) :: new_up(:),new_dn(:)
  type(ik_vec),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#else
  integer(ik),intent(in) :: old_up(nwalk),old_dn(nwalk)
  integer(ik),allocatable,intent(out) :: new_up(:),new_dn(:)
  integer(ik),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#endif
  integer iwalk,nshift
  real(rk),allocatable :: tmp_num(:)

  allocate(walk_dets_up(nwalk))
  allocate(walk_dets_dn(nwalk))

  walk_dets_up = old_up
  walk_dets_dn = old_dn

  allocate(tmp_num(nwalk))
  tmp_num = old_e_mix_num

  nshift=0
  do iwalk=2,nwalk
    if(walk_dets_up(iwalk).eq.walk_dets_up(iwalk-nshift-1) .and. walk_dets_dn(iwalk).eq.walk_dets_dn(iwalk-nshift-1)) then ! This walker is on same det as previous walker
      nshift=nshift+1
      tmp_num(iwalk-nshift)=tmp_num(iwalk-nshift)+tmp_num(iwalk)
    else                                                                ! This walker is on different det from previous walker. Complete working on previous det before working on this one.
! Done working on previous det.  Now move position of current one.
      walk_dets_up(iwalk-nshift)=walk_dets_up(iwalk)
      walk_dets_dn(iwalk-nshift)=walk_dets_dn(iwalk)
      tmp_num(iwalk-nshift)=tmp_num(iwalk)
    endif
  enddo
! Fix up last det, just as we did above after "else".
  nwalk=nwalk-nshift

  allocate(new_up(nwalk))
  allocate(new_dn(nwalk))

  new_up(1:nwalk) = walk_dets_up(1:nwalk)
  new_dn(1:nwalk) = walk_dets_dn(1:nwalk)

  deallocate(walk_dets_up,walk_dets_dn)

  allocate(e_mix_num(nwalk))
  e_mix_num(1:nwalk) = tmp_num(1:nwalk)
  deallocate(tmp_num)

end subroutine merge_original_with_spawned3_num
!---------------------------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine merge_original_with_spawned3_num_denom_replace(nwalk, dets_up, dets_dn, e_mix_num, e_mix_den)
! =====================================================================================================================
! Description   : Merge walkers
! Author        : Cyrus Umrigar
! Modified      : A Holmes
!                 Replaces lists of dets and local energies (does not reallocate)
! Comments      : Pulled out from semistoch_mod.f90 by HJC and used in a common area
! Last edited   : A Holmes, 8 Aug 2012. Made e_mix_den optional.
! ---------------------------------------------------------------------------------------------------------------------

  implicit none
  integer, intent(inout) :: nwalk
  real(rk),intent(inout) :: e_mix_num(:)
  real(rk),optional,intent(inout) :: e_mix_den(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(inout) :: dets_up(:),dets_dn(:)
  type(ik_vec),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#else
  integer(ik),intent(inout) :: dets_up(:),dets_dn(:)
  integer(ik),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#endif
  integer iwalk,nshift
  real(rk),allocatable :: tmp_num(:),tmp_den(:)

  allocate(walk_dets_up(nwalk))
  allocate(walk_dets_dn(nwalk))

  walk_dets_up = dets_up(1:nwalk)
  walk_dets_dn = dets_dn(1:nwalk)

  allocate(tmp_num(nwalk))
  tmp_num = e_mix_num(1:nwalk)
  if (present(e_mix_den)) then
    allocate(tmp_den(nwalk))
    tmp_den = e_mix_den(1:nwalk)
  endif

  nshift=0
  do iwalk=2,nwalk
    if(walk_dets_up(iwalk).eq.walk_dets_up(iwalk-nshift-1) .and. walk_dets_dn(iwalk).eq.walk_dets_dn(iwalk-nshift-1)) then ! This walker is on same det as previous walker
      nshift=nshift+1
      tmp_num(iwalk-nshift)=tmp_num(iwalk-nshift)+tmp_num(iwalk)
      if (present(e_mix_den))  tmp_den(iwalk-nshift)=tmp_den(iwalk-nshift)+tmp_den(iwalk)
    else                                                                ! This walker is on different det from previous walker. Complete working on previous det before working on this one.
! Done working on previous det.  Now move position of current one.
      walk_dets_up(iwalk-nshift)=walk_dets_up(iwalk)
      walk_dets_dn(iwalk-nshift)=walk_dets_dn(iwalk)
      tmp_num(iwalk-nshift)=tmp_num(iwalk)
      if (present(e_mix_den))  tmp_den(iwalk-nshift)=tmp_den(iwalk)
    endif
  enddo
! Fix up last det, just as we did above after "else".
  nwalk=nwalk-nshift

  dets_up(1:nwalk) = walk_dets_up(1:nwalk)
  dets_dn(1:nwalk) = walk_dets_dn(1:nwalk)

  deallocate(walk_dets_up,walk_dets_dn)

  e_mix_num(1:nwalk) = tmp_num(1:nwalk)
  deallocate(tmp_num)
  if (present(e_mix_den)) then
    e_mix_den(1:nwalk) = tmp_den(1:nwalk)
    deallocate(tmp_den)
  endif

end subroutine merge_original_with_spawned3_num_denom_replace
!---------------------------------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine merge_original_with_spawned3_num_denom_diag_elems_info_replace(nwalk, dets_up, dets_dn, e_mix_num, e_mix_den, diag_elems_info)
! =====================================================================================================================
! Description   : Merge walkers
! Author        : Cyrus Umrigar
! Modified      : A Holmes
!                 Replaces lists of dets and local energies (does not reallocate)
! Comments      : Pulled out from semistoch_mod.f90 by HJC and used in a common area
! Last edited   : A Holmes, 8 Aug 2012. Made e_mix_den optional.
! ---------------------------------------------------------------------------------------------------------------------

  use common_run, only : diag_elem_info

  implicit none
  integer, intent(inout) :: nwalk
  real(rk),intent(inout) :: e_mix_num(:)
  real(rk),intent(inout) :: e_mix_den(:)
  type(diag_elem_info),intent(inout) :: diag_elems_info(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(inout) :: dets_up(:),dets_dn(:)
  type(ik_vec),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#else
  integer(ik),intent(inout) :: dets_up(:),dets_dn(:)
  integer(ik),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#endif
  integer iwalk,nshift
  real(rk),allocatable :: tmp_num(:),tmp_den(:)
  type(diag_elem_info),allocatable :: tmp_diag_elems_info(:)

  allocate(walk_dets_up(nwalk))
  allocate(walk_dets_dn(nwalk))

  walk_dets_up = dets_up(1:nwalk)
  walk_dets_dn = dets_dn(1:nwalk)

  allocate(tmp_num(nwalk))
  allocate(tmp_den(nwalk))
  allocate(tmp_diag_elems_info(nwalk))
  tmp_num = e_mix_num(1:nwalk)
  tmp_den = e_mix_den(1:nwalk)
  tmp_diag_elems_info = diag_elems_info(1:nwalk)

  nshift=0
  do iwalk=2,nwalk
    if(walk_dets_up(iwalk).eq.walk_dets_up(iwalk-nshift-1) .and. walk_dets_dn(iwalk).eq.walk_dets_dn(iwalk-nshift-1)) then ! This walker is on same det as previous walker
      nshift=nshift+1
      tmp_num(iwalk-nshift)=tmp_num(iwalk-nshift)+tmp_num(iwalk)
      tmp_den(iwalk-nshift)=tmp_den(iwalk-nshift)+tmp_den(iwalk)
      tmp_diag_elems_info(iwalk-nshift)=tmp_diag_elems_info(iwalk)
    else                                                                ! This walker is on different det from previous walker. Complete working on previous det before working on this one.
! Done working on previous det.  Now move position of current one.
      walk_dets_up(iwalk-nshift)=walk_dets_up(iwalk)
      walk_dets_dn(iwalk-nshift)=walk_dets_dn(iwalk)
      tmp_num(iwalk-nshift)=tmp_num(iwalk)
      tmp_den(iwalk-nshift)=tmp_den(iwalk)
      tmp_diag_elems_info(iwalk-nshift)=tmp_diag_elems_info(iwalk)
    endif
  enddo
! Fix up last det, just as we did above after "else".
  nwalk=nwalk-nshift

  dets_up(1:nwalk) = walk_dets_up(1:nwalk)
  dets_dn(1:nwalk) = walk_dets_dn(1:nwalk)

  deallocate(walk_dets_up,walk_dets_dn)

  e_mix_num(1:nwalk) = tmp_num(1:nwalk)
  e_mix_den(1:nwalk) = tmp_den(1:nwalk)
  diag_elems_info(1:nwalk) = tmp_diag_elems_info(1:nwalk)

  deallocate(tmp_num)
  deallocate(tmp_den)
  deallocate(tmp_diag_elems_info)

end subroutine merge_original_with_spawned3_num_denom_diag_elems_info_replace
!---------------------------------------------------------------------------------------------------------------------

subroutine merge_original_with_spawned3_5items_replace(nwalk, dets_up, dets_dn, e_mix_num, e_mix_den, diag_elems_info, term1_big, term2_big)
! =====================================================================================================================
! Description   : Merge walkers
! Author        : Cyrus Umrigar
! Modified      : A Holmes
!                 Replaces lists of dets and local energies (does not reallocate)
! Comments      : Pulled out from semistoch_mod.f90 by HJC and used in a common area
! Last edited   : A Holmes, 8 Aug 2012. Made e_mix_den optional.
! ---------------------------------------------------------------------------------------------------------------------

  use common_run, only : diag_elem_info

  implicit none
  integer, intent(inout) :: nwalk
  real(rk),intent(inout) :: e_mix_num(:)
  real(rk),intent(inout) :: e_mix_den(:)
  type(diag_elem_info),intent(inout) :: diag_elems_info(:)
  real(rk),intent(inout) :: term1_big(:)
  real(rk),intent(inout) :: term2_big(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(inout) :: dets_up(:),dets_dn(:)
  type(ik_vec),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#else
  integer(ik),intent(inout) :: dets_up(:),dets_dn(:)
  integer(ik),allocatable :: walk_dets_up(:),walk_dets_dn(:)
#endif
  integer iwalk,nshift
  real(rk),allocatable :: tmp_num(:),tmp_den(:)
  real(rk),allocatable :: tmp_1big(:),tmp_2big(:)
  type(diag_elem_info),allocatable :: tmp_diag_elems_info(:)

  allocate(walk_dets_up(nwalk))
  allocate(walk_dets_dn(nwalk))

  walk_dets_up = dets_up(1:nwalk)
  walk_dets_dn = dets_dn(1:nwalk)

  allocate(tmp_num(nwalk))
  allocate(tmp_den(nwalk))
  allocate(tmp_diag_elems_info(nwalk))
  allocate(tmp_1big(nwalk))
  allocate(tmp_2big(nwalk))

  tmp_num = e_mix_num(1:nwalk)
  tmp_den = e_mix_den(1:nwalk)
  tmp_diag_elems_info = diag_elems_info(1:nwalk)
  tmp_1big = term1_big(1:nwalk)
  tmp_2big = term2_big(1:nwalk)

  nshift=0
  do iwalk=2,nwalk
    if(walk_dets_up(iwalk).eq.walk_dets_up(iwalk-nshift-1) .and. walk_dets_dn(iwalk).eq.walk_dets_dn(iwalk-nshift-1)) then ! This walker is on same det as previous walker
      nshift=nshift+1
      tmp_num(iwalk-nshift)=tmp_num(iwalk-nshift)+tmp_num(iwalk)
      tmp_den(iwalk-nshift)=tmp_den(iwalk-nshift)+tmp_den(iwalk)
      tmp_diag_elems_info(iwalk-nshift)=tmp_diag_elems_info(iwalk)
      tmp_1big(iwalk-nshift)=tmp_1big(iwalk-nshift)+tmp_1big(iwalk)
      tmp_2big(iwalk-nshift)=tmp_2big(iwalk-nshift)+tmp_2big(iwalk)
    else                                                                ! This walker is on different det from previous walker. Complete working on previous det before working on this one.
! Done working on previous det.  Now move position of current one.
      walk_dets_up(iwalk-nshift)=walk_dets_up(iwalk)
      walk_dets_dn(iwalk-nshift)=walk_dets_dn(iwalk)
      tmp_num(iwalk-nshift)=tmp_num(iwalk)
      tmp_den(iwalk-nshift)=tmp_den(iwalk)
      tmp_diag_elems_info(iwalk-nshift)=tmp_diag_elems_info(iwalk)
      tmp_1big(iwalk-nshift)=tmp_1big(iwalk)
      tmp_2big(iwalk-nshift)=tmp_2big(iwalk)
    endif
  enddo
! Fix up last det, just as we did above after "else".
  nwalk=nwalk-nshift

  dets_up(1:nwalk) = walk_dets_up(1:nwalk)
  dets_dn(1:nwalk) = walk_dets_dn(1:nwalk)

  deallocate(walk_dets_up,walk_dets_dn)

  e_mix_num(1:nwalk) = tmp_num(1:nwalk)
  e_mix_den(1:nwalk) = tmp_den(1:nwalk)
  diag_elems_info(1:nwalk) = tmp_diag_elems_info(1:nwalk)
  term1_big(1:nwalk) = tmp_1big(1:nwalk)
  term2_big(1:nwalk) = tmp_2big(1:nwalk)

  deallocate(tmp_num)
  deallocate(tmp_den)
  deallocate(tmp_diag_elems_info)
  deallocate(tmp_1big)
  deallocate(tmp_2big)

end subroutine merge_original_with_spawned3_5items_replace
!---------------------------------------------------------------------------------------------------------------------

subroutine merge_original_with_spawned3_add_wts(dets_up, dets_dn, walk_wt, nwalk)
!=====================================================================================================================
! Description   : Merge walkers
! Author        : Cyrus Umrigari/HJC
! Comments      : Modified by HJC for replacing walker list with a shortened list with added weights
!               : Primarily used when one wants to add up matrix elements of determinants related by symmetry
!                 Reallocates dets and weights!
! Last edited   : A Holmes, 2 Aug 2012. Make arguments unallocatable.
! ---------------------------------------------------------------------------------------------------------------------
  use types, only: ik
  implicit none
  integer, intent(inout)                :: nwalk
  real(rk),intent(inout)    :: walk_wt(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(inout) :: dets_up(:),dets_dn(:)
  type(ik_vec),allocatable               :: walk_dets_up(:), walk_dets_dn(:)
#else
  integer(ik),intent(inout) :: dets_up(:),dets_dn(:)
  integer(ik),allocatable               :: walk_dets_up(:), walk_dets_dn(:)
#endif
  real(rk),allocatable                  :: tmp_wt(:)
  integer                               :: iwalk,nshift

  allocate(walk_dets_up(size(dets_up)))
  allocate(walk_dets_dn(size(dets_dn)))
  allocate(tmp_wt(size(walk_wt)))

  walk_dets_up = dets_up
  walk_dets_dn = dets_dn
  tmp_wt = walk_wt

  nshift=0
  do iwalk=2,nwalk
    if(walk_dets_up(iwalk).eq.walk_dets_up(iwalk-nshift-1) .and. walk_dets_dn(iwalk).eq.walk_dets_dn(iwalk-nshift-1)) then ! This walker is on same det as previous walker
      nshift=nshift+1
      tmp_wt(iwalk-nshift)=tmp_wt(iwalk-nshift)+tmp_wt(iwalk)
    else                                                                ! This walker is on different det from previous walker. Complete working on previous det before working on this one.
! Done working on previous det.  Now move position of current one.
      walk_dets_up(iwalk-nshift)=walk_dets_up(iwalk)
      walk_dets_dn(iwalk-nshift)=walk_dets_dn(iwalk)
      tmp_wt(iwalk-nshift)=tmp_wt(iwalk)
    endif
  enddo
! Fix up last det, just as we did above after "else".
  nwalk=nwalk-nshift

 !deallocate (dets_up,dets_dn)
 !allocate(dets_up(nwalk))
 !allocate(dets_dn(nwalk))

  dets_up(1:nwalk) = walk_dets_up(1:nwalk)
  dets_dn(1:nwalk) = walk_dets_dn(1:nwalk)

 !deallocate(walk_wt)
 !allocate(walk_wt(nwalk))
  walk_wt(1:nwalk) = tmp_wt(1:nwalk)

  end subroutine merge_original_with_spawned3_add_wts
!---------------------------------------------------------------------------------------------------------------------

function get_free_memory()
!=====================================================================================================================
! Description   : Returns free physical memory in bytes
!               : Does not work right.  Use mem_avail.f90 instead.
! Author        : A Holmes, 18 May 2012
! ---------------------------------------------------------------------------------------------------------------------
  implicit none
  real :: get_free_memory

  call system("free -b | grep Mem | awk '{print $4}' > tmp.sD2fiTHTLO")
  open(7, file = 'tmp.sD2fiTHTLO',status='old')
  read(7,*) get_free_memory
  close(7)
  call system("rm -f tmp.sD2fiTHTLO")

end function get_free_memory
!---------------------------------------------------------------------------------------------------------------------

subroutine verbose_allocate_int(variable,length)
!=====================================================================================================================
! Description   : Allocates, but prints the size you're allocating, the memory used in allocation, and remaining
!                 free memory before and after allocation.
! Author        : A Holmes, 22 May 2012
! ---------------------------------------------------------------------------------------------------------------------
  implicit none
  integer,intent(inout),allocatable :: variable(:)
  integer,intent(in) :: length
  real :: memory,memory_used

  memory = get_free_memory()
  write(6,'(/,''Attempting to allocate integer array of size '',i12,''; current physical memory ='',es10.3,'' bytes.'')') length,memory
  call flush(6)
  allocate(variable(length))
  memory = get_free_memory()
  memory_used = length*real(storage_size(variable))/8
  if (nint(memory_used) < 10000) then
    write(6,'(''Allocate successful! Used '',i4,'' bytes; remaining physical memory ='',es10.3,'' bytes.'')') nint(memory_used),memory
  else
    write(6,'(''Allocate successful! Used'',es10.3,'' bytes; remaining physical memory ='',es10.3,'' bytes.'')') memory_used,memory
  endif
  call flush(6)

end subroutine verbose_allocate_int
!---------------------------------------------------------------------------------------------------------------------

subroutine verbose_allocate_ik(variable,length)
!=====================================================================================================================
! Description   : Allocates, but prints the size you're allocating, the memory used in allocation, and remaining
!                 free memory before and after allocation.
! Author        : A Holmes, 22 May 2012
! ---------------------------------------------------------------------------------------------------------------------
  implicit none
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(inout),allocatable :: variable(:)
#else
  integer(ik),intent(inout),allocatable :: variable(:)
#endif
  integer,intent(in) :: length
  real :: memory,memory_used

  memory = get_free_memory()
#ifdef NUM_ORBITALS_GT_127
  write(6,'(/,''Attempting to allocate type(ik_vec) array of size '',i12,''; current physical memory ='',es10.3,'' bytes.'')') length,memory
#else
  write(6,'(/,''Attempting to allocate integer(ik) array of size '',i12,''; current physical memory ='',es10.3,'' bytes.'')') length,memory
#endif
  call flush(6)
  allocate(variable(length))
  memory = get_free_memory()
  memory_used = length*real(storage_size(variable))/8
  if (nint(memory_used) < 10000) then
    write(6,'(''Allocate successful! Used '',i4,'' bytes; remaining physical memory ='',es10.3,'' bytes.'')') nint(memory_used),memory
  else
    write(6,'(''Allocate successful! Used'',es10.3,'' bytes; remaining physical memory ='',es10.3,'' bytes.'')') memory_used,memory
  endif
  call flush(6)

end subroutine verbose_allocate_ik
!---------------------------------------------------------------------------------------------------------------------

subroutine verbose_allocate_rk(variable,length)
!=====================================================================================================================
! Description   : Allocates, but prints the size you're allocating, the memory used in allocation, and remaining
!                 free memory before and after allocation.
! Author        : A Holmes, 22 May 2012
! ---------------------------------------------------------------------------------------------------------------------
  implicit none
  real(rk),intent(inout),allocatable :: variable(:)
  integer,intent(in) :: length
  real :: memory,memory_used

  memory = get_free_memory()
  write(6,'(/,''Attempting to allocate real(rk) array of size '',i12,''; current physical memory ='',es10.3,'' bytes.'')') length,memory
  call flush(6)
  allocate(variable(length))
  memory = get_free_memory()
  memory_used = length*real(storage_size(variable))/8
  if (nint(memory_used) < 10000) then
    write(6,'(''Allocate successful! Used '',i4,'' bytes; remaining physical memory ='',es10.3,'' bytes.'')') nint(memory_used),memory
  else
    write(6,'(''Allocate successful! Used'',es10.3,'' bytes; remaining physical memory ='',es10.3,'' bytes.'')') memory_used,memory
  endif
  call flush(6)

end subroutine verbose_allocate_rk
!---------------------------------------------------------------------------------------------------------------------

subroutine verbose_allocate_rk2(variable,length1,length2)
!=====================================================================================================================
! Description   : Allocates, but prints the size you're allocating, the memory used in allocation, and remaining
!                 free memory before and after allocation.
! Author        : A Holmes, 22 May 2012
! ---------------------------------------------------------------------------------------------------------------------
  implicit none
  real(rk),intent(inout),allocatable :: variable(:,:)
  integer,intent(in) :: length1,length2
  real :: memory,memory_used

  memory = get_free_memory()
  write(6,'(/,''Attempting to allocate 2D real(rk) array of size '',i12,''; current physical memory ='',es10.3,'' bytes.'')') length1*length2,memory
  call flush(6)
  allocate(variable(length1,length2))
  memory = get_free_memory()
  memory_used = length1*length2*real(storage_size(variable))/8
  if (nint(memory_used) < 10000) then
    write(6,'(''Allocate successful! Used '',i4,'' bytes; remaining physical memory ='',es10.3,'' bytes.'')') nint(memory_used),memory
  else
    write(6,'(''Allocate successful! Used'',es10.3,'' bytes; remaining physical memory ='',es10.3,'' bytes.'')') memory_used,memory
  endif
  call flush(6)

end subroutine verbose_allocate_rk2
!---------------------------------------------------------------------------------------------------------------------

subroutine verbose_allocate_int_test(variable,length,ratio)
!=====================================================================================================================
! Description   : Allocates, but prints the size you're allocating, the memory used in allocation, and remaining
!                 free memory before and after allocation.
!                 Outputs the ratio of the size of the allocated variable to the difference in free memory before and after allocation
!                 Only used for testing
! Author        : A Holmes, 22 May 2012
! ---------------------------------------------------------------------------------------------------------------------
  implicit none
  integer,intent(inout),allocatable :: variable(:)
  integer,intent(in) :: length
  real(rk),intent(out) :: ratio
  real :: memory,memory_used

  memory = get_free_memory()
  allocate(variable(length))
  ratio = memory
  memory = get_free_memory()
  memory_used = length*real(storage_size(variable))/8
  ratio = memory_used/(ratio - memory)
  write (6,*) "size",length,"ratio",ratio
  call flush(6)

end subroutine verbose_allocate_int_test
!---------------------------------------------------------------------------------------------------------------------

function permutation_factor(det1,det2)
!=====================================================================================================================
! Description   : For two bit-packed dets that are one excitation apart, calculates permutation factor as follows:
!                 Count number of electrons between them and return 1 if that number is even, -1 if odd.
! Author        : A Holmes, 22 May 2012
! ---------------------------------------------------------------------------------------------------------------------
  use types, only : i8b
  implicit none
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: det1,det2
  type(ik_vec) :: diff
#else
  integer(ik),intent(in) :: det1,det2
  integer(ik) :: diff
#endif
  integer :: permutation_factor

  if (ik==i8b) then
    if (det1<0.and.det2>0) then
      diff = iand(det1,det1-det2)
    elseif (det1>0.and.det2<0) then
      diff = iand(det2,det2-det1)
    else
      if (det1>det2) then
        diff = iand(det1,det1-det2)
      else
        diff = iand(det2,det2-det1)
      endif
    endif
  else
    if (det1>det2) then
      diff = iand(det1,det1-det2)
    else
      diff = iand(det2,det2-det1)
    endif
  endif

  if (poppar(diff)==0) then
    permutation_factor = 1
  else
    permutation_factor = -1
  endif

end function permutation_factor
!---------------------------------------------------------------------------------------------------------------------
! ===============================================================================================


subroutine permutation_factor2(det_i,det_j,gamma,first_i_bit,second_i_bit,first_j_bit,second_j_bit,nosign)
!=====================================================================================================================
! Description   : For two bit-packed dets that are two excitations apart, calculates permutation factor.
!                 Also returns "from" and "to" sites for the excitations.
! Author        : A Holmes, 25 May 2012
! ---------------------------------------------------------------------------------------------------------------------
  use types, only : ik

  implicit none
  integer,intent(out) :: gamma,first_i_bit,second_i_bit,first_j_bit,second_j_bit
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: det_i,det_j
  type(ik_vec) :: diff
#else
  integer(ik),intent(in) :: det_i,det_j
  integer(ik) :: diff
#endif
  logical,optional,intent(in) :: nosign

  ! First find the ordered electrons that do the hopping, i.e., the electrons
  ! present in one det but not the other. These are called {first_i_bit, second_i_bit}
  ! for det_i, and are similarly named for det_j.

  diff = iand(det_i,not(det_j))
  first_i_bit = trailz(diff)
  second_i_bit = trailz(ibclr(diff,first_i_bit))

  diff = iand(det_j,not(det_i))
  first_j_bit = trailz(diff)
  second_j_bit = trailz(ibclr(diff,first_j_bit))

  if (present(nosign)) return

  ! Count the number of electrons that are hopped over exactly once by the electrons that are exciting.
  ! This is done by first identifying which orbitals are hopped over by the first excitation,
  ! ieor(maskr(first_i_bit,ik),maskr(first_j_bit,ik)), and which ones are hopped over by the second
  ! excitation, ieor(maskr(second_i_bit,ik),maskr(second_j_bit,ik)), and taking an "exclusive or"
  ! between those two bitmasks (in case they overlap). Next, take an "and" between the resulting
  ! bitmask and dets i and j, to convert the "orbitals hopped over" into "electrons hopped over."
  ! Both dets i and j are used here so that the electrons that do the hopping are not accidentally
  ! included. Finally, if the number of electrons hopped over once is even, gamma=1; else, gamma=-1.

#ifdef NUM_ORBITALS_GT_127
  diff = iand(det_i,iand(det_j,ieor(ieor(maskr_vec(first_i_bit),maskr_vec(first_j_bit)),ieor(maskr_vec(second_i_bit),maskr_vec(second_j_bit)))))
#else
  diff = iand(det_i,iand(det_j,ieor(ieor(maskr(first_i_bit,ik),maskr(first_j_bit,ik)),ieor(maskr(second_i_bit,ik),maskr(second_j_bit,ik)))))
#endif

  if (poppar(diff)==0) then
    gamma = 1
  else
    gamma = -1
  endif

end subroutine permutation_factor2
!---------------------------------------------------------------------------------------------------------------------

function count_excitations(det1,det2)
!=====================================================================================================================
! Description   : For two bit-packed dets, counts number of excitations using efficient bit counting.
!               : It does this by counting the number of bits that have changed from 1 to 0
! Author        : A Holmes, 23 May 2012
! ---------------------------------------------------------------------------------------------------------------------
  implicit none
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: det1,det2
#else
  integer(ik),intent(in) :: det1,det2
#endif
  integer :: count_excitations

  count_excitations = popcnt(iand(det1,not(det2)))
  !count_excitations = popcnt(ieor(det1,det2))/2 ! This should give the same

end function count_excitations
!---------------------------------------------------------------------------------------------------------------------

subroutine print_excitation_levels_and_wts(n_det,dets_up,dets_dn,lowest_eigenvector,norb,orbital_symmetries,hf_up_remote,hf_dn_remote)
!=====================================================================================================================
! Description   : Calculate number of dets with various excitation levels and the sum of their squared wts.
! Author        : Cyrus Umrigar, 3 Jul 2012. (moved to this subroutine from chemistry.f90 by A Holmes)
! ---------------------------------------------------------------------------------------------------------------------

    implicit none

    integer,intent(in) :: n_det
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
    type(ik_vec), optional,intent(in) :: hf_up_remote,hf_dn_remote
    type(ik_vec) :: hf_up,hf_dn
#else
    integer(ik),intent(in) :: dets_up(:),dets_dn(:)
    integer(ik), optional,intent(in) :: hf_up_remote,hf_dn_remote
    integer(ik) :: hf_up,hf_dn
#endif
    real(rk),intent(in) :: lowest_eigenvector(:)
    integer,intent(in) :: norb
    integer,optional,intent(in) :: orbital_symmetries(:)

    integer :: i,iorb
    real(rk), allocatable :: sum_w2_orb(:)
    integer :: excite_level
    integer, parameter :: max_excite_level=20 ! Used for printing statistics only
    integer :: popul_excit_level(0:max_excite_level)
    real(rk) ::sum_w2_excit_level(0:max_excite_level)

    if (present(hf_up_remote).and.present(hf_dn_remote)) then
      !MJO In parallel, hf is not guaranteed to be the first determinant, so we set
      !    it explicitly
      hf_up = hf_up_remote 
      hf_dn = hf_dn_remote 
    else
      hf_up = dets_up(1)
      hf_dn = dets_dn(1)
    endif

    popul_excit_level(1:max_excite_level)=0
    sum_w2_excit_level(1:max_excite_level)=0
    do i=1,n_det
      excite_level=min(popcnt(iand(dets_up(i),not(hf_up)))+popcnt(iand(dets_dn(i),not(hf_dn))),max_excite_level)
      popul_excit_level(excite_level)=popul_excit_level(excite_level)+1
      sum_w2_excit_level(excite_level)=sum_w2_excit_level(excite_level)+lowest_eigenvector(i)**2
    enddo
    write(6,'(/,''Excit_level  # dets  sum_w2_excit_level'')')
    do excite_level=0,max_excite_level
      if(sum_w2_excit_level(excite_level).ne.0_rk) write(6,'(i6,i10,f14.8)') excite_level, popul_excit_level(excite_level), sum_w2_excit_level(excite_level)
    enddo
! Calculate sum of squared wts on each orbital
    allocate(sum_w2_orb(norb))
    sum_w2_orb(1:norb)=0
    do i=1,n_det
        ! Total occupation being measured - instead of up and down separately
       do iorb=1,norb
        if(btest(dets_up(i),iorb-1)) sum_w2_orb(iorb)=sum_w2_orb(iorb)+lowest_eigenvector(i)**2
        if(btest(dets_dn(i),iorb-1)) sum_w2_orb(iorb)=sum_w2_orb(iorb)+lowest_eigenvector(i)**2
       enddo
    enddo
    if (present(orbital_symmetries)) then
      write(6,'(/,''Orbital orb_sym sum_w2_orb'')')
      do iorb=1,norb
        if(sum_w2_orb(iorb).ne.0_rk) write(6,'(2i6,f12.6)') iorb, orbital_symmetries(iorb), sum_w2_orb(iorb)
      enddo
    else
      write(6,'(/,''Orbital sum_w2_orb'')')
      do iorb=1,norb
        if(sum_w2_orb(iorb).ne.0_rk) write(6,'(i6,f12.6)') iorb, sum_w2_orb(iorb)
      enddo
    endif
    deallocate(sum_w2_orb)
    call flush(6)

end subroutine print_excitation_levels_and_wts
!---------------------------------------------------------------------------------------------------------------------

subroutine matrix_inversion(n,A,x,b,B0,ipiv,work)
  ! Solves Ax=b for x by using lapack, by finding A^-1 first, which is not the most efficient way.
  ! B0 must be same size as A; ipiv must be same size as x and b.
  ! A Holmes, 22 Aug 2012.

  integer,intent(in) :: n
  real(rk),intent(in) :: A(:,:)
  real(rk),intent(out) :: x(:)
  real(rk),intent(in) :: b(:)
  real(rk),intent(out) :: B0(:,:)
  integer,intent(out) :: ipiv(:)
  real(rk),intent(out) :: work(:)

  integer :: info,lwork

  B0 = A
  lwork = size(work)
  call dgetrf(n,n,B0,n,ipiv,info)
  call dgetri(n,B0,n,ipiv,work,lwork,info)
  x = matmul(B0,b)

end subroutine matrix_inversion
!---------------------------------------------------------------------------------------------------------------------

logical function do_once()
  ! Function that returns true the first time it is called and false otherwise.
  ! Only used for debugging.
  ! A Holmes, 9 Jan 2013.

  use common_run, only : first_time

  do_once = first_time
  if (first_time)  first_time = .false.

end function do_once

logical function do_n_times(n)
  ! Function that returns true the first n times it is called and false otherwise.
  ! Only used for debugging.
  ! A Holmes, 25 Jan 2015.

  use common_run, only : debug_counter
  integer,intent(in) :: n

  debug_counter = debug_counter + 1
  do_n_times = (debug_counter<=n)

end function do_n_times


subroutine choose_det_with_prob_prop_to_abs_wt(iwalk,tot_wt)
  ! Given nwalk walkers of weight walk_wt(1:nwalk), select one
  ! of them with probability proportional to its absolute wt.
  ! Time taken scales as nwalk.
  ! A Holmes, 15 Aug 2013.
  use common_walk, only : walk_wt,nwalk
  integer,intent(out) :: iwalk
  real(rk),intent(out) :: tot_wt
  real(rk) :: rand_wt,tmp_wt,rannyu
  integer :: i

  tot_wt = sum(abs(walk_wt(1:nwalk)))
  rand_wt = tot_wt*rannyu()
  tmp_wt = 0._rk

  do i=1,nwalk
    tmp_wt = tmp_wt + abs(walk_wt(i))
    if (tmp_wt > rand_wt) then
      iwalk = i
      return
    endif
  enddo
  iwalk = nwalk

end subroutine choose_det_with_prob_prop_to_abs_wt
!---------------------------------------------------------------------------------------------------------------------


  subroutine sort_and_merge_count_repeats(n,arr,counts)
    ! Sorts and merges dets list to remove duplicates
    ! A Holmes, 29 Sep 2016

    use generic_sort, only : sort

    integer,intent(inout) :: n
    integer,intent(inout) :: arr(:)
    integer,intent(out) :: counts(:)

    integer :: i,nshift

    counts(:) = 1

    call sort(n,arr(1:n))

    nshift = 0
    do i=2,n
      if (arr(i)==arr(i-1)) then
        nshift = nshift + 1
        counts(i) = counts(i) + counts(i-1)
      endif
      arr(i-nshift) = arr(i)
      counts(i-nshift) = counts(i)
    enddo

    n = n - nshift

  end subroutine sort_and_merge_count_repeats
!---------------------------------------------------------------------------------------------------------------------

  subroutine sort_and_merge_int(n,arr)
    ! Sorts and merges dets list to remove duplicates
    ! A Holmes, 5 Apr 2016

    use generic_sort, only : sort

    integer,intent(inout) :: n
    integer,intent(inout) :: arr(:)

    integer :: i,nshift

    call sort(n,arr(1:n))

    nshift = 0
    do i=2,n
      if (arr(i)==arr(i-1)) then
        nshift = nshift + 1
      endif
      arr(i-nshift) = arr(i)
    enddo

    n = n - nshift

  end subroutine sort_and_merge_int
!---------------------------------------------------------------------------------------------------------------------

  subroutine sort_and_merge_int_rk(n,arr,arr_rk)
    ! Sorts and merges dets list to remove duplicates
    ! A Holmes, 5 Apr 2016

    use generic_sort, only : sort_by_first_argument

    integer,intent(inout) :: n
    integer,intent(inout) :: arr(:)
    real(rk),intent(inout) :: arr_rk(:)

    integer :: i,nshift

    call sort_by_first_argument(n,arr(1:n),arr_rk(1:n))

    nshift = 0
    do i=2,n
      if (arr(i)==arr(i-1)) then
        nshift = nshift + 1
      endif
      arr(i-nshift) = arr(i)
      arr_rk(i-nshift) = arr_rk(i)
    enddo

    n = n - nshift

  end subroutine sort_and_merge_int_rk
!---------------------------------------------------------------------------------------------------------------------

  subroutine sort_and_merge_det(n_dets,dets_up,dets_dn)
    ! Sorts and merges dets list to remove duplicates
    ! A Holmes, 4 Apr 2016

    integer,intent(inout) :: n_dets
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(inout) :: dets_up(:),dets_dn(:)
#else
    integer(ik),intent(inout) :: dets_up(:),dets_dn(:)
#endif

    integer :: i,nshift

    call quick_sort(1, n_dets)

    nshift = 0
    do i=2,n_dets
      if (dets_up(i)==dets_up(i-1).and.dets_dn(i)==dets_dn(i-1)) then
        nshift = nshift + 1
      endif
      dets_up(i-nshift) = dets_up(i)
      dets_dn(i-nshift) = dets_dn(i)
    enddo

    n_dets = n_dets - nshift

    contains

      recursive subroutine quick_sort(left_end, right_end)

        integer, intent(in) :: left_end, right_end

        !     local variables
        integer             :: i, j
#ifdef NUM_ORBITALS_GT_127
        type(ik_vec)        :: reference_up, temp_up
        type(ik_vec)        :: reference_dn, temp_dn
#else
        integer(ik)         :: reference_up, temp_up
        integer(ik)         :: reference_dn, temp_dn
#endif
        integer, parameter  :: max_simple_sort_size = 6

        if (right_end<left_end+max_simple_sort_size) then
          ! use interchange sort for small arrs
          call interchange_sort(left_end, right_end)
        else
          ! use partition ("quick") sort
          reference_up = dets_up((left_end+ right_end)/2)
          reference_dn = dets_dn((left_end+ right_end)/2)
          i = left_end- 1; j = right_end + 1
          do
            ! scan dets from left end until element >= reference is found
            do
              i = i + 1
              if (dets_up(i) > reference_up.or.(dets_up(i)==reference_up.and.dets_dn(i) >= reference_dn)) exit
            enddo
            ! scan arr from right end until element <= reference is found
            do
              j = j - 1
              if (dets_up(j) < reference_up.or.(dets_up(j)==reference_up.and.dets_dn(j) <= reference_dn)) exit
            enddo
            if (i < j) then
              ! swap two out-of-order elements
              temp_up = dets_up(i); dets_up(i) = dets_up(j); dets_up(j) = temp_up
              temp_dn = dets_dn(i); dets_dn(i) = dets_dn(j); dets_dn(j) = temp_dn
            elseif (i == j) then
              i = i + 1
              exit
            else
              exit
            endif
          enddo
          if (left_end<j)  call quick_sort(left_end, j)
          if (i<right_end)  call quick_sort(i, right_end)
        endif

      end subroutine quick_sort
!---------------------------------------------------------------------------------------------------------------------

      subroutine interchange_sort(left_end, right_end)
        integer, intent(in) :: left_end, right_end
        integer             :: i, j
#ifdef NUM_ORBITALS_GT_127
        type(ik_vec)           :: temp
#else
        integer(ik)            :: temp
#endif

        do i = left_end, right_end - 1
          do j = i+1, right_end
            if (dets_up(i)>dets_up(j).or.(dets_up(i)==dets_up(j).and.dets_dn(i)>dets_dn(j))) then
              temp = dets_up(i); dets_up(i) = dets_up(j); dets_up(j) = temp
              temp = dets_dn(i); dets_dn(i) = dets_dn(j); dets_dn(j) = temp
            endif
          enddo
        enddo
      end subroutine interchange_sort

  end subroutine sort_and_merge_det
!---------------------------------------------------------------------------------------------------------------------

  subroutine welford(n,x,M,S,var)
  ! Welford's method for single-pass mean and variance
  ! A Holmes, 30 Sep 2016

    integer,intent(in) :: n ! Data point number (1,2,3,...)
    real(rk),intent(in) :: x ! Data point
    real(rk),intent(inout) :: M ! Mean (initialize to 0)
    real(rk),intent(inout) :: S ! Variance times n-1 (initialize to 0)
    real(rk),intent(out) :: var ! Variance

    real(rk) :: oldM

      oldM = M
      M = M + (x-M)/n
      S = S + (x-M)*(x-oldM)
      var = S/real(n-1)/n

  end subroutine welford
!---------------------------------------------------------------------------------------------------------------------

function round_r(r,n)
!==============================================================
! Rounds a positive real(rk) number to n significant digits
! Author: Cyrus Umrigar
!==============================================================
integer, intent(in)  :: n
real(rk), intent(in) :: r
real(rk)             :: round_r, ten_to_log10i
if(r.lt.0) stop 'r < 0'
ten_to_log10i = 10.d0**(floor(log10(r))-n+1)
round_r=(int(r/ten_to_log10i))*ten_to_log10i
return
end function round_r
!---------------------------------------------------------------------------------------------------------------------

function round_i(i,n)
!==============================================================
! Rounds a positive 4-byte integer down to n significant digits
! Author: Cyrus Umrigar
!==============================================================
integer, intent(in) :: i,n
integer round_i, ten_to_log10i
if(i.lt.0) stop 'i < 0'
ten_to_log10i = max(1,10**(floor(log10(real(i)))-n+1))
round_i=(i/ten_to_log10i)*ten_to_log10i
return
end function round_i
!---------------------------------------------------------------------------------------------------------------------

  subroutine alloc_i(array_name, array, dim1)
!---------------------------------------------------------------------------
! Allocate or reallocate an ik or ik_vec 1-index array
! Author : Cyrus Umrigar, Jan 2017
!---------------------------------------------------------------------------
  implicit none

  character(*), intent(in)          :: array_name
  integer(i8b), intent(in)          :: dim1
  integer,intent(inout),allocatable :: array(:)

  integer ierr
  integer(i8b) array_dim
  integer, allocatable              :: array_tmp(:)
  real(rk) avail,avail2,avail3

! allocate array if not already allocated
  if(.not. allocated(array)) then

    allocate(array(dim1), stat=ierr)
    if(ierr /= 0) then
      write(6,'(''alloc_i: allocation for array'',a,'' dim'',i10,'' failed'')') trim(array_name), dim1
      stop 'allocation failed in alloc_i'
    endif

!   array(:) = 0.d0

  else

! resize array if already allocated with different dimension
    array_dim = size(array)

    if(dim1 < array_dim) then
      array=array(1:dim1)

    else

      allocate(array_tmp(dim1), stat=ierr)
      if(ierr /= 0) then
        write(6,'(''alloc_i: allocation for array'',a,'' dim'',i10,'' failed'')') trim(array_name), dim1
        stop 'allocation failed in alloc_i'
      endif

      array_tmp(1:array_dim)=array
      call move_alloc(array_tmp,array)
    endif

  endif

  call mem_avail(avail,avail2,avail3)
  write(6,'(/,''After allocating '',a,'' of size'',i10,'', Memory avail (MB)='',t117,9f9.2)') array_name, dim1, avail, avail2, avail3

  end subroutine alloc_i

!===========================================================================

  subroutine alloc_ik(array_name, array, dim1)
!---------------------------------------------------------------------------
! Allocate or reallocate an ik or ik_vec 1-index array
! Author : Cyrus Umrigar, Jan 2017
!---------------------------------------------------------------------------
  implicit none

  character(*), intent(in)       :: array_name
  integer(i8b), intent(in)       :: dim1
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(inout),allocatable :: array(:)
  type(ik_vec),              allocatable :: array_tmp(:)
#else
  integer(ik),intent(inout),allocatable  :: array(:)
  integer(ik),              allocatable  :: array_tmp(:)
#endif

  integer ierr
  integer(i8b) array_dim
  real(rk) avail,avail2,avail3

! allocate array if not already allocated
  if(.not. allocated(array)) then

    allocate(array(dim1), stat=ierr)
    if(ierr /= 0) then
      write(6,'(''alloc_ik: allocation for array'',a,'' dim'',i10,'' failed'')') trim(array_name), dim1
      stop 'allocation failed in alloc_ik'
    endif

!   array(:) = 0.d0

  else

! resize array if already allocated with different dimension
    array_dim = size(array)

    if(dim1 < array_dim) then
      array=array(1:dim1)

    else

      allocate(array_tmp(dim1), stat=ierr)
      if(ierr /= 0) then
        write(6,'(''alloc_ik: allocation for array'',a,'' dim'',i10,'' failed'')') trim(array_name), dim1
        stop 'allocation failed in alloc_ik'
      endif

      array_tmp(1:array_dim)=array
      call move_alloc(array_tmp,array)
    endif

  endif

  call mem_avail(avail,avail2,avail3)
  write(6,'(/,''After allocating '',a,'' of size'',i10,'', Memory avail (MB)='',t117,9f9.2)') array_name, dim1, avail, avail2, avail3

  end subroutine alloc_ik

!===========================================================================

  subroutine alloc_rk_1(array_name, array, dim1)
!---------------------------------------------------------------------------
! Allocate or reallocate a kind rk (double precision), 1-index array
! Author : Cyrus Umrigar, Jan 2017
!---------------------------------------------------------------------------
  implicit none

  character(*), intent(in)        :: array_name
  integer(i8b), intent(in)        :: dim1
  real(rk), allocatable           :: array(:)

  integer ierr
  integer(i8b) array_dim
  real(rk), allocatable           :: array_tmp(:)
  real(rk) avail,avail2,avail3

! allocate array if not already allocated
  if(.not. allocated(array)) then

    allocate(array(dim1), stat=ierr)
    if(ierr /= 0) then
      write(6,'(''alloc_rk_1: allocation for array'',a,'' dim'',i10,'' failed'')') trim(array_name), dim1
      stop 'allocation failed in alloc_rk_1'
    endif

!   array(:) = 0.d0

  else

! resize array if already allocated with different dimension
    array_dim = size(array)

    if(dim1 < array_dim) then
      array=array(1:dim1)

    else

      allocate(array_tmp(dim1), stat=ierr)
      if(ierr /= 0) then
        write(6,'(''alloc_rk_1: allocation for array'',a,'' dim'',i10,'' failed'')') trim(array_name), dim1
        stop 'allocation failed in alloc_rk_1'
      endif

      array_tmp(1:array_dim)=array
      call move_alloc(array_tmp,array)
    endif

  endif

  call mem_avail(avail,avail2,avail3)
  write(6,'(/,''After allocating '',a,'' of size'',i10,'', Memory avail (MB)='',t117,9f9.2)') array_name, dim1, avail, avail2, avail3

  end subroutine alloc_rk_1

!===========================================================================
  subroutine alloc_rk_2(array_name, array, dim1, dim2)
!---------------------------------------------------------------------------
! Allocate or reallocate a double precision, 2-index array
! Author : Cyrus Umrigar, Jan 2017
!---------------------------------------------------------------------------
  implicit none

  character(len=*), intent(in)    :: array_name
  integer(i8b), intent(in)        :: dim1, dim2
  real(rk), allocatable           :: array(:,:)

  integer ierr
  integer(i8b) array_dim1, array_dim2

  real(rk), allocatable           :: array_tmp(:,:)
  real(rk) avail,avail2,avail3

! allocate array if not already allocated
  if(.not. allocated(array)) then

    allocate(array(dim1,dim2), stat=ierr)
    if(ierr /= 0) then
      write(6,'(''alloc_rk_2: allocation for array'',a,'' dim1,dim2'',2i10,'' failed'')') trim(array_name), dim1,dim2
      stop 'allocation failed in alloc_rk_2'
    endif

!   array(:,:) = 0.d0

  else

! resize array if already allocated with different dimension
    array_dim1 = size(array,1)
    array_dim2 = size(array,2)

    if(dim1 <= array_dim1 .and. dim2 <= array_dim2) then
      array=array(1:dim1,1:dim2)

    else

      allocate(array_tmp(dim1,dim2), stat=ierr)
      if(ierr /= 0) then
        write(6,'(''alloc_rk_2: allocation for array'',a,'' dim1,dim2'',2i10,'' failed'')') trim(array_name), dim1,dim2
        stop 'allocation failed in alloc_rk_2'
      endif

      array_tmp(1:array_dim1,1:array_dim2)=array
      call move_alloc(array_tmp,array)
    endif

  endif

  call mem_avail(avail,avail2,avail3)
  write(6,'(/,''After allocating '',a,'', of size'',2i10,'', Memory avail (MB)='',t117,9f9.2)') array_name, dim1, dim2, avail, avail2, avail3
  call flush(6)

  end subroutine alloc_rk_2
!===========================================================================

  subroutine alloc_diag_elem_info(array_name, array, dim1)
!---------------------------------------------------------------------------
! Allocate or reallocate an alloc_diag_elem_info type 1-index array
! Author : Cyrus Umrigar, Jan 2017
!---------------------------------------------------------------------------
  implicit none

  character(*), intent(in)       :: array_name
  integer(i8b), intent(in)       :: dim1
  type(diag_elem_info),intent(inout),allocatable :: array(:)
  type(diag_elem_info),              allocatable :: array_tmp(:)

  integer ierr
  integer(i8b) array_dim
  real(rk) avail,avail2,avail3

! allocate array if not already allocated
  if(.not. allocated(array)) then

    allocate(array(dim1), stat=ierr)
    if(ierr /= 0) then
      write(6,'(''alloc_ik: allocation for array'',a,'' dim'',i10,'' failed'')') trim(array_name), dim1
      stop 'allocation failed in alloc_ik'
    endif

!   array(:) = 0.d0

  else

! resize array if already allocated with different dimension
    array_dim = size(array)

    if(dim1 < array_dim) then
      array=array(1:dim1)

    else

      allocate(array_tmp(dim1), stat=ierr)
      if(ierr /= 0) then
        write(6,'(''alloc_ik: allocation for array'',a,'' dim'',i10,'' failed'')') trim(array_name), dim1
        stop 'allocation failed in alloc_ik'
      endif

      array_tmp(1:array_dim)=array
      call move_alloc(array_tmp,array)
    endif

  endif

  call mem_avail(avail,avail2,avail3)
  write(6,'(/,''After allocating '',a,'' of size'',i10,'', Memory avail (MB)='',t117,9f9.2)') array_name, dim1, avail, avail2, avail3

  end subroutine alloc_diag_elem_info

end module tools
