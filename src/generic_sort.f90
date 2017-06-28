module generic_sort
  use overload
  use types, only: rk, ik, ik_vec
  implicit none
  save
  private
  public :: sort,quick_sort,sort_by_first_argument
  public :: shell_sort_real_rank1_int_rank1
  !setup interface for ease of use and to allow sort to be
  !changed without editing the rest of the code
  interface sort
#ifdef NUM_ORBITALS_GT_127
    module procedure shell_sort_int, shell_sort_simple_int, shell_sort_int_int, shell_sort_ik_vec, shell_sort_real_rank2, shell_sort_real_rank1_int_rank1_int_rank1,shell_sort_real_rank1_ik_vec_rank1_ik_vec_rank1,shell_sort_real_rank1_int_rank1
#else
     module procedure shell_sort_int, shell_sort_simple_int, shell_sort_int_int, shell_sort_real_rank2, shell_sort_real_rank1_int_rank1_int_rank1,shell_sort_real_rank1_ik_vec_rank1_ik_vec_rank1
#endif
  end interface sort

  interface quick_sort
    module procedure quick_sort_int, quick_sort_i8b
  end interface quick_sort

  interface sort_by_first_argument
    module procedure sort_by_first_argument_int_int, sort_by_first_argument_ik, sort_by_first_argument_int_rk, sort_by_first_argument_rk_int
  end interface sort_by_first_argument

contains

  !=============================================================
  subroutine sort_by_first_argument_int_int(n,list_int_sort_by,list_int)
    !---------------------------------------------------------------------------
    ! Description : sort lists such that list_int_sort_by is in increasing order using shell sort O(nlog^2 n)
    !               Note shell sort is not stable sort.
    !
    ! Created     : H.J.Changlani- borrowed from F. Petruzielo, 28 Oct 2010 (modified from rosetta code)
    ! Modified    : A Holmes, Apr 2016. Removed unneeded arguments
    !---------------------------------------------------------------------------
   
    implicit none

    !dummy arguments
    integer,intent(in)     :: n
    integer, intent(inout) :: list_int_sort_by(n)
    integer, intent(inout) :: list_int(n)

    !local variables
    integer :: i, j, increment
    integer :: temp_int
    integer :: temp_int2

    increment = size(list_int_sort_by) / 2
    do while (increment > 0)
       do i = increment+1, size(list_int_sort_by)
          j = i
          temp_int = list_int_sort_by(i)
          temp_int2 = list_int(i)
          do while (j >= increment+1)
             if (list_int_sort_by(j-increment) <= temp_int) exit
             list_int_sort_by(j) = list_int_sort_by(j-increment)
             list_int(j) = list_int(j-increment)
             j = j - increment
          enddo
          list_int_sort_by(j) = temp_int
          list_int(j) = temp_int2
       enddo
       if (increment == 2) then
          increment = 1
       else
          increment = increment * 5 / 11
       endif
    enddo

   end subroutine sort_by_first_argument_int_int

  subroutine sort_by_first_argument_int_rk(n,list_rk_sort_by,list_rk)
    !---------------------------------------------------------------------------
    ! Description : sort lists such that list_rk_sort_by is in increasing order using shell sort O(nlog^2 n)
    !               Note shell sort is not stable sort.
    !
    ! Created     : H.J.Changlani- borrowed from F. Petruzielo, 28 Oct 2010 (modified from rosetta code)
    ! Modified    : A Holmes, Apr 2016. Removed unneeded arguments
    !---------------------------------------------------------------------------
   
    implicit none

    !dummy arguments
    integer,intent(in)     :: n
    integer, intent(inout) :: list_rk_sort_by(n)
    real(rk), intent(inout) :: list_rk(n)

    !local variables
    integer :: i, j, increment
    integer :: temp_int
    real(rk) :: temp_int2

    increment = size(list_rk_sort_by) / 2
    do while (increment > 0)
       do i = increment+1, size(list_rk_sort_by)
          j = i
          temp_int = list_rk_sort_by(i)
          temp_int2 = list_rk(i)
          do while (j >= increment+1)
             if (list_rk_sort_by(j-increment) <= temp_int) exit
             list_rk_sort_by(j) = list_rk_sort_by(j-increment)
             list_rk(j) = list_rk(j-increment)
             j = j - increment
          enddo
          list_rk_sort_by(j) = temp_int
          list_rk(j) = temp_int2
       enddo
       if (increment == 2) then
          increment = 1
       else
          increment = increment * 5 / 11
       endif
    enddo

   end subroutine sort_by_first_argument_int_rk

  subroutine sort_by_first_argument_rk_int(n,list_rk_sort_by,list_rk)
    !---------------------------------------------------------------------------
    ! Description : sort lists such that list_rk_sort_by is in increasing order using shell sort O(nlog^2 n)
    !               Note shell sort is not stable sort.
    !
    ! Created     : H.J.Changlani- borrowed from F. Petruzielo, 28 Oct 2010 (modified from rosetta code)
    ! Modified    : A Holmes, Apr 2016. Removed unneeded arguments
    !---------------------------------------------------------------------------
   
    implicit none

    !dummy arguments
    integer,intent(in)     :: n
    real(rk), intent(inout) :: list_rk_sort_by(n)
    integer, intent(inout) :: list_rk(n)

    !local variables
    integer :: i, j, increment
    real(rk) :: temp_int
    integer :: temp_int2

    increment = size(list_rk_sort_by) / 2
    do while (increment > 0)
       do i = increment+1, size(list_rk_sort_by)
          j = i
          temp_int = list_rk_sort_by(i)
          temp_int2 = list_rk(i)
          do while (j >= increment+1)
             if (list_rk_sort_by(j-increment) <= temp_int) exit
             list_rk_sort_by(j) = list_rk_sort_by(j-increment)
             list_rk(j) = list_rk(j-increment)
             j = j - increment
          enddo
          list_rk_sort_by(j) = temp_int
          list_rk(j) = temp_int2
       enddo
       if (increment == 2) then
          increment = 1
       else
          increment = increment * 5 / 11
       endif
    enddo

   end subroutine sort_by_first_argument_rk_int

  !=============================================================
  subroutine sort_by_first_argument_ik(n,list_int_sort_by,list_int)
    !---------------------------------------------------------------------------
    ! Description : sort lists such that list_int_sort_by is in increasing order using shell sort O(nlog^2 n)
    !               Note shell sort is not stable sort.
    !
    ! Created     : H.J.Changlani- borrowed from F. Petruzielo, 28 Oct 2010 (modified from rosetta code)
    ! Modified    : A Holmes, Apr 2016. Removed unneeded arguments
    !---------------------------------------------------------------------------
   
    implicit none

    !dummy arguments
    integer,intent(in)     :: n
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(inout) :: list_int_sort_by(n)
#else
    integer(ik), intent(inout) :: list_int_sort_by(n)
#endif
    integer, intent(inout) :: list_int(n)

    !local variables
    integer :: i, j, increment
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec) :: temp_int
#else
    integer(ik) :: temp_int
#endif
    integer :: temp_int2

    increment = size(list_int_sort_by) / 2
    do while (increment > 0)
       do i = increment+1, size(list_int_sort_by)
          j = i
          temp_int = list_int_sort_by(i)
          temp_int2 = list_int(i)
          do while (j >= increment+1)
             if (list_int_sort_by(j-increment) <= temp_int) exit
             list_int_sort_by(j) = list_int_sort_by(j-increment)
             list_int(j) = list_int(j-increment)
             j = j - increment
          enddo
          list_int_sort_by(j) = temp_int
          list_int(j) = temp_int2
       enddo
       if (increment == 2) then
          increment = 1
       else
          increment = increment * 5 / 11
       endif
    enddo

   end subroutine sort_by_first_argument_ik

  !=============================================================
  subroutine shell_sort_simple_int(n,list_int)
    !---------------------------------------------------------------------------
    ! Description : sort list_real by magnitude using shell sort O(nlog^2 n)
    !               Note shell sort is not stable sort.
    !
    ! Created     : H.J.Changlani- borrowed from F. Petruzielo, 28 Oct 2010 (modified from rosetta code)
    !---------------------------------------------------------------------------
   
    implicit none

    !dummy arguments
    integer,intent(in)     :: n
    integer, intent(inout) :: list_int(n)

    !local variables
    integer :: i, j, increment
    integer :: temp_int

    increment = size(list_int) / 2
    do while (increment > 0)
       do i = increment+1, size(list_int)
          j = i
          temp_int = list_int(i)
          do while (j >= increment+1)
             if (list_int(j-increment) >= temp_int) exit
             list_int(j) = list_int(j-increment)
             j = j - increment
          enddo
          list_int(j) = temp_int
       enddo
       if (increment == 2) then
          increment = 1
       else
          increment = increment * 5 / 11
       endif
    enddo

   end subroutine shell_sort_simple_int
  
  !=============================================================
  subroutine shell_sort_int(n,list_int,num_ops)
    !---------------------------------------------------------------------------
    ! Description : sort list_real by magnitude using shell sort O(nlog^2 n)
    !               Note shell sort is not stable sort.
    !
    ! Created     : H.J.Changlani- borrowed from F. Petruzielo, 28 Oct 2010 (modified from rosetta code)
    !---------------------------------------------------------------------------
   
    implicit none

    !dummy arguments
    integer,intent(in)     :: n
    integer, intent(inout) :: list_int(n)
    integer, intent(out)   :: num_ops

    !local variables
    integer :: i, j, increment
    integer :: temp_int

    num_ops=0

    increment = size(list_int) / 2
    do while (increment > 0)
       do i = increment+1, size(list_int)
          j = i
          temp_int = list_int(i)
          do while (j >= increment+1)
             if (list_int(j-increment) >= temp_int) exit
             list_int(j) = list_int(j-increment)
             j = j - increment
             num_ops=num_ops+1
          enddo
          list_int(j) = temp_int
       enddo
       if (increment == 2) then
          increment = 1
       else
          increment = increment * 5 / 11
       endif
    enddo

   end subroutine shell_sort_int
  
  !=============================================================
  recursive subroutine quick_sort_int(left_in,right,list_int,list_real)
    !---------------------------------------------------------------------------
    ! Description : sort list_int by magnitude using quick sort O(nlog n)
    !               Real numbers come with ints, and order changes in corresponding way
    !
    ! Created     : A Holmes, 20 Mar 2014
    !---------------------------------------------------------------------------
   
    implicit none

    !dummy arguments
    integer,intent(in),optional :: left_in
    integer,intent(in)     :: right
    integer, intent(inout) :: list_int(:)
    real(rk),intent(inout) :: list_real(:)

    !local variables
    integer :: i,j,left
    integer :: tmp,pivot
    real(rk) :: tmp_real

    if (present(left_in)) then
      left = left_in
    else
      left = 1
    endif

    if (left==right)  return ! already sorted (asked to sort 1-element list)

    if (left==right-1) then ! 2-element list
      if (list_int(left)>list_int(right)) then
        tmp = list_int(left)
        list_int(left) = list_int(right)
        list_int(right) = tmp
        tmp_real = list_real(left)
        list_real(left) = list_real(right)
        list_real(right) = tmp_real
      endif
      return
    endif

    pivot = (left+right)/2

    ! Here, move elements less than pivot to left of pivot, more than pivot to right of pivot
    i = left-1
    do j=left,pivot-1
      i = i+1
      if (list_int(i)>list_int(pivot)) then
        tmp = list_int(i)
        list_int(i:pivot-1) = list_int(i+1:pivot)
        list_int(pivot) = tmp
        tmp_real = list_real(i)
        list_real(i:pivot-1) = list_real(i+1:pivot)
        list_real(pivot) = tmp_real
        pivot = pivot - 1
        i = i-1
      endif
    enddo
    i = right+1
    do j=right,pivot+1,-1
      i = i-1
      if (list_int(i)<list_int(pivot)) then
        tmp = list_int(i)
        list_int(pivot+1:i) = list_int(pivot:i-1)
        list_int(pivot) = tmp
        tmp_real = list_real(i)
        list_real(pivot+1:i) = list_real(pivot:i-1)
        list_real(pivot) = tmp_real
        pivot = pivot + 1
        i = i+1
      endif
    enddo

    if (left < pivot-1)  call quick_sort_int(left,pivot-1,list_int,list_real)
    if (right > pivot+1)  call quick_sort_int(pivot+1,right,list_int,list_real)

   end subroutine quick_sort_int
  
  recursive subroutine quick_sort_i8b(left_in,right,list_int,list_real)
    !---------------------------------------------------------------------------
    ! Description : sort list_int by magnitude using quick sort O(nlog n)
    !               Real numbers come with ints, and order changes in corresponding way
    !
    ! Created     : A Holmes, 20 Mar 2014
    !---------------------------------------------------------------------------

    use types, only : i8b
   
    implicit none

    !dummy arguments
    integer(i8b),intent(in),optional :: left_in
    integer(i8b),intent(in)     :: right
    integer(i8b), intent(inout) :: list_int(:)
    real(rk),intent(inout) :: list_real(:)

    !local variables
    integer(i8b) :: i,j,left
    integer(i8b) :: tmp,pivot
    real(rk) :: tmp_real

    if (present(left_in)) then
      left = left_in
    else
      left = 1
    endif

    if (left==right)  return ! already sorted (asked to sort 1-element list)

    if (left==right-1) then ! 2-element list
      if (list_int(left)>list_int(right)) then
        tmp = list_int(left)
        list_int(left) = list_int(right)
        list_int(right) = tmp
        tmp_real = list_real(left)
        list_real(left) = list_real(right)
        list_real(right) = tmp_real
      endif
      return
    endif

    pivot = (left+right)/2

    ! Here, move elements less than pivot to left of pivot, more than pivot to right of pivot
    i = left-1
    do j=left,pivot-1
      i = i+1
      if (list_int(i)>list_int(pivot)) then
        tmp = list_int(i)
        list_int(i:pivot-1) = list_int(i+1:pivot)
        list_int(pivot) = tmp
        tmp_real = list_real(i)
        list_real(i:pivot-1) = list_real(i+1:pivot)
        list_real(pivot) = tmp_real
        pivot = pivot - 1
        i = i-1
      endif
    enddo
    i = right+1
    do j=right,pivot+1,-1
      i = i-1
      if (list_int(i)<list_int(pivot)) then
        tmp = list_int(i)
        list_int(pivot+1:i) = list_int(pivot:i-1)
        list_int(pivot) = tmp
        tmp_real = list_real(i)
        list_real(pivot+1:i) = list_real(pivot:i-1)
        list_real(pivot) = tmp_real
        pivot = pivot + 1
        i = i+1
      endif
    enddo

    if (left < pivot-1)  call quick_sort_i8b(left,pivot-1,list_int,list_real)
    if (right > pivot+1)  call quick_sort_i8b(pivot+1,right,list_int,list_real)

   end subroutine quick_sort_i8b
  
  !=============================================================
  subroutine shell_sort_int_int(list_int1,list_int2)
    !---------------------------------------------------------------------------
    ! Description : sort list_real by magnitude using shell sort O(nlog^2 n)
    !               Note shell sort is not stable sort.
    !
    ! Created     : H.J.Changlani- borrowed from F. Petruzielo, 28 Oct 2010 (modified from rosetta code)
    !---------------------------------------------------------------------------
   
    implicit none

    !dummy arguments
    integer(ik), intent(inout) :: list_int1(:),list_int2(:)

    !local variables
    integer     :: i, j, increment
    integer(ik) :: temp_int1,temp_int2

    increment = size(list_int1) / 2
    do while (increment > 0)
       do i = increment+1, size(list_int1)
          j = i
          temp_int1 = list_int1(i)
          temp_int2 = list_int2(i)
          do while (j >= increment+1)
             if (list_int1(j-increment) < temp_int1) exit
             if ((list_int1(j-increment) .eq. temp_int1) .and. list_int2(j-increment) < temp_int2) exit
             list_int1(j) = list_int1(j-increment)
             list_int2(j) = list_int2(j-increment)
             j = j - increment
          enddo
          list_int1(j) = temp_int1
          list_int2(j) = temp_int2
       enddo
       if (increment == 2) then
          increment = 1
       else
          increment = increment * 5 / 11
       endif
    enddo

   end subroutine shell_sort_int_int
 
#ifdef NUM_ORBITALS_GT_127
  !=============================================================
  subroutine shell_sort_ik_vec(list_int1,list_int2)
    !---------------------------------------------------------------------------
    ! Description : sort list_real by magnitude using shell sort O(nlog^2 n)
    !               Note shell sort is not stable sort.
    !
    ! Created     : H.J.Changlani- borrowed from F. Petruzielo, 28 Oct 2010 (modified from rosetta code)
    !---------------------------------------------------------------------------
    use types, only : ik_vec
    use overload
    implicit none

    !dummy arguments
    type(ik_vec), intent(inout) :: list_int1(:),list_int2(:)

    !local variables
    integer     :: i, j, increment
    type(ik_vec) :: temp_int1,temp_int2

    increment = size(list_int1) / 2
    do while (increment > 0)
       do i = increment+1, size(list_int1)
          j = i
          temp_int1 = list_int1(i)
          temp_int2 = list_int2(i)
          do while (j >= increment+1)
             if (list_int1(j-increment) < temp_int1) exit
             if ((list_int1(j-increment) .eq. temp_int1) .and. list_int2(j-increment) < temp_int2) exit
             list_int1(j) = list_int1(j-increment)
             list_int2(j) = list_int2(j-increment)
             j = j - increment
          enddo
          list_int1(j) = temp_int1
          list_int2(j) = temp_int2
       enddo
       if (increment == 2) then
          increment = 1
       else
          increment = increment * 5 / 11
       endif
    enddo

   end subroutine shell_sort_ik_vec
 
  !=============================================================
#endif
 
  !=============================================================
  subroutine shell_sort_real_rank2(list_real)
    !---------------------------------------------------------------------------
    ! Description : sort list_real by magnitude using shell sort O(nlog^2 n)
    !               Note shell sort is not stable sort.
    !
    ! Created     : F. Petruzielo, 28 Oct 2010 (modified from rosetta code)
    !---------------------------------------------------------------------------
    implicit none

    !dummy arguments
    real(rk), intent(inout) :: list_real(:,:)

    !local variables
    integer :: i, j, increment
    real(rk), allocatable :: temp(:)

    allocate(temp(size(list_real, 1)))

    increment = size(list_real, 2) / 2
    do while (increment > 0)
       do i = increment+1, size(list_real, 2)
          j = i
          temp(:) = list_real(:,i)
          do while (j >= increment+1)
             if (sum(list_real(:,j-increment)**2) <= sum(temp(:)**2)) exit
             list_real(:,j) = list_real(:, j-increment)
             j = j - increment
          enddo
          list_real(:,j) = temp(:)
       enddo
       if (increment == 2) then
          increment = 1
       else
          increment = increment * 5 / 11
       endif
    enddo

  end subroutine shell_sort_real_rank2

  !==================================================================================================

  !===================================================================================================
  subroutine shell_sort_real_rank1_int_rank1_int_rank1(list_real, list_int1, list_int2,consider_sign)
    !------------------------------------------------------------------------------------------------
    ! Description : sort list_real by magnitude (greatest to least) using shell sort O(nlog^2 n)
    !               sort list_int1 and list_int2 into the same order
    !
    ! Created     : F. Petruzielo, 6 Oct 2010 (modified from rosetta code)
    ! Edited      : H.J. Changlani, Feb 9 2012 to include option to sort by considering sign as well
    !               instead of just magnitude
    !-----------------------------------------------------------------------------------------------
    implicit none

    !dummy arguments
    real(rk), intent(inout) :: list_real(:)
    integer(ik), intent(inout) :: list_int1(:)
    integer(ik), intent(inout) :: list_int2(:)
    logical, optional  :: consider_sign

    !local variables
    integer :: i, j, increment
    real(rk) :: temp_real
    integer(ik) :: temp_int1, temp_int2
    logical     :: sign_on

    sign_on=.false.                       ! Default is to sort by magnitude

    if (present(consider_sign)) sign_on=consider_sign

    if (sign_on .eqv. .false.) then       ! Frank's code untouched
        increment = size(list_real) / 2
        do while (increment > 0)
           do i = increment+1, size(list_real)
              j = i
              temp_real = list_real(i)
              temp_int1 = list_int1(i)
              temp_int2 = list_int2(i)
              do while (j >= increment+1)
                 if (abs(list_real(j-increment)) >= abs(temp_real)) exit
                 list_real(j) = list_real(j-increment)
                 list_int1(j) = list_int1(j-increment)
                 list_int2(j) = list_int2(j-increment)
                 j = j - increment
              enddo
              list_real(j) = temp_real
              list_int1(j) = temp_int1
              list_int2(j) = temp_int2
           enddo
           if (increment == 2) then
              increment = 1
           else
              increment = increment * 5 / 11
           endif
        enddo

        !Switch signs if necessary to make largest in magnitude coef be positive
        if(list_real(1).lt.0) list_real=-list_real
    
     else                                ! Sort considering sign as well - HJC
                                         ! Minor change - just removed abs function in sorting metric                                             
        increment = size(list_real) / 2
        do while (increment > 0)
           do i = increment+1, size(list_real)
              j = i
              temp_real = list_real(i)
              temp_int1 = list_int1(i)
              temp_int2 = list_int2(i)
              do while (j >= increment+1)
                 if ((list_real(j-increment)) >= (temp_real)) exit
                 list_real(j) = list_real(j-increment)
                 list_int1(j) = list_int1(j-increment)
                 list_int2(j) = list_int2(j-increment)
                 j = j - increment
              enddo
              list_real(j) = temp_real
              list_int1(j) = temp_int1
              list_int2(j) = temp_int2
           enddo
           if (increment == 2) then
              increment = 1
           else
              increment = increment * 5 / 11
           endif
        enddo
     endif


  end subroutine shell_sort_real_rank1_int_rank1_int_rank1
  !=============================================================

    !===================================================================================================
  subroutine shell_sort_real_rank1_int_rank1(list_real, list_int1, consider_sign)
    !------------------------------------------------------------------------------------------------
    ! Description : sort list_real by magnitude (greatest to least) using shell sort O(nlog^2 n)
    !               sort list_int1 into the same order
    !
    ! Created     : F. Petruzielo, 6 Oct 2010 (modified from rosetta code)
    ! Edited      : H.J. Changlani, Feb 9 2012 to include option to sort by considering sign as well
    !               instead of just magnitude
    !-----------------------------------------------------------------------------------------------
    implicit none

    !dummy arguments
    real(rk), intent(inout) :: list_real(:)
    integer, intent(inout) :: list_int1(:)
    logical, optional  :: consider_sign

    !local variables
    integer :: i, j, increment
    real(rk) :: temp_real
    integer :: temp_int1
    logical     :: sign_on

    sign_on=.false.                       ! Default is to sort by magnitude

    if (present(consider_sign)) sign_on=consider_sign

    if (sign_on .eqv. .false.) then       ! Frank's code untouched
        increment = size(list_real) / 2
        do while (increment > 0)
           do i = increment+1, size(list_real)
              j = i
              temp_real = list_real(i)
              temp_int1 = list_int1(i)
              do while (j >= increment+1)
                 if (abs(list_real(j-increment)) >= abs(temp_real)) exit
                 list_real(j) = list_real(j-increment)
                 list_int1(j) = list_int1(j-increment)
                 j = j - increment
              enddo
              list_real(j) = temp_real
              list_int1(j) = temp_int1
           enddo
           if (increment == 2) then
              increment = 1
           else
              increment = increment * 5 / 11
           endif
        enddo

        !Switch signs if necessary to make largest in magnitude coef be positive
        if(list_real(1).lt.0) list_real=-list_real

     else                                ! Sort considering sign as well - HJC
                                         ! Minor change - just removed abs function in sorting metri
        increment = size(list_real) / 2
        do while (increment > 0)
           do i = increment+1, size(list_real)
              j = i
              temp_real = list_real(i)
              temp_int1 = list_int1(i)
              do while (j >= increment+1)
                 if ((list_real(j-increment)) >= (temp_real)) exit
                 list_real(j) = list_real(j-increment)
                 list_int1(j) = list_int1(j-increment)
                 j = j - increment
              enddo
              list_real(j) = temp_real
              list_int1(j) = temp_int1
           enddo
           if (increment == 2) then
              increment = 1
           else
              increment = increment * 5 / 11
           endif
        enddo
     endif


  end subroutine shell_sort_real_rank1_int_rank1
  !=============================================================


  !===================================================================================================
  subroutine shell_sort_real_rank1_ik_vec_rank1_ik_vec_rank1(list_real, list_int1, list_int2,consider_sign)
    !------------------------------------------------------------------------------------------------
    ! Description : sort list_real by magnitude (greatest to least) using shell sort O(nlog^2 n)
    !               sort list_int1 and list_int2 into the same order
    !
    ! Created     : F. Petruzielo, 6 Oct 2010 (modified from rosetta code)
    ! Edited      : H.J. Changlani, Feb 9 2012 to include option to sort by considering sign as well
    !               instead of just magnitude
    !-----------------------------------------------------------------------------------------------
    use overload
    implicit none

    !dummy arguments
    real(rk), intent(inout) :: list_real(:)
    type(ik_vec), intent(inout) :: list_int1(:)
    type(ik_vec), intent(inout) :: list_int2(:)
    logical, optional  :: consider_sign

    !local variables
    integer :: i, j, increment
    real(rk) :: temp_real
    type(ik_vec) :: temp_int1, temp_int2
    logical     :: sign_on

    sign_on=.false.                       ! Default is to sort by magnitude

    if (present(consider_sign)) sign_on=consider_sign

    if (sign_on .eqv. .false.) then       ! Frank's code untouched
        increment = size(list_real) / 2
        do while (increment > 0)
           do i = increment+1, size(list_real)
              j = i
              temp_real = list_real(i)
              temp_int1 = list_int1(i)
              temp_int2 = list_int2(i)
              do while (j >= increment+1)
                 if (abs(list_real(j-increment)) >= abs(temp_real)) exit
                 list_real(j) = list_real(j-increment)
                 list_int1(j) = list_int1(j-increment)
                 list_int2(j) = list_int2(j-increment)
                 j = j - increment
              enddo
              list_real(j) = temp_real
              list_int1(j) = temp_int1
              list_int2(j) = temp_int2
           enddo
           if (increment == 2) then
              increment = 1
           else
              increment = increment * 5 / 11
           endif
        enddo

        !Switch signs if necessary to make largest in magnitude coef be positive
        if(list_real(1).lt.0) list_real=-list_real
    
     else                                ! Sort considering sign as well - HJC
                                         ! Minor change - just removed abs function in sorting metric                                             
        increment = size(list_real) / 2
        do while (increment > 0)
           do i = increment+1, size(list_real)
              j = i
              temp_real = list_real(i)
              temp_int1 = list_int1(i)
              temp_int2 = list_int2(i)
              do while (j >= increment+1)
                 if ((list_real(j-increment)) >= (temp_real)) exit
                 list_real(j) = list_real(j-increment)
                 list_int1(j) = list_int1(j-increment)
                 list_int2(j) = list_int2(j-increment)
                 j = j - increment
              enddo
              list_real(j) = temp_real
              list_int1(j) = temp_int1
              list_int2(j) = temp_int2
           enddo
           if (increment == 2) then
              increment = 1
           else
              increment = increment * 5 / 11
           endif
        enddo
     endif


  end subroutine shell_sort_real_rank1_ik_vec_rank1_ik_vec_rank1
  !=============================================================

end module generic_sort
