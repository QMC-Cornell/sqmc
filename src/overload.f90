module overload
! Module that overloads operators and intrinsic functions that previously worked on integer(ik), so that they
! will also work on type(ik_vec), a vector of integer(ik)'s.
! Throughout this module, assume that any integer(ik)'s are less than or equal to huge(ik).
! A Holmes, 22 Apr 2013.

#ifdef NUM_ORBITALS_GT_127

  use types, only : rk, i4b, ik, ik_vec, num_words, bits_per_word
  implicit none

  interface assignment (=)
    module procedure assign_from_vec, assign_from_int, assign_from_ik, assign_from_rk, &
      & assign_vec_from_int, assign_vec_from_ik, assign_array_from_ik, assign_to_ik
  end interface

  interface operator (==)
    module procedure eq_vec, eq_int, eq_ik
  end interface

  interface operator (/=)
    module procedure ne_vec, ne_int, ne_ik
  end interface

  interface operator (<)
    module procedure lt_vec, lt_int, lt_ik
  end interface

  interface operator (>)
    module procedure gt_vec, gt_int, gt_ik
  end interface

  interface operator (<=)
    module procedure le_vec, le_int, le_ik
  end interface

  interface operator (>=)
    module procedure ge_vec, ge_int, ge_ik
  end interface

  interface operator (+)
    module procedure det_plus_det, det_plus_int, int_plus_det, det_plus_ik, ik_plus_det
  end interface

  interface operator (-)
    module procedure negative_det, det_minus_det, det_minus_int, det_minus_ik
  end interface

  interface operator (*)
    module procedure det_times_rk, int_times_det, ik_times_det
  end interface

  interface operator (/)
    module procedure det_div_det
  end interface

  interface abs
   module procedure abs_vec
  end interface

  interface btest
    module procedure btest_vec
  end interface

  interface ibset
    module procedure ibset_vec
  end interface

  interface ibclr
    module procedure ibclr_vec
  end interface

  interface trailz
    module procedure trailz_vec
  end interface

 !interface maskr
 !  module procedure maskr_vec  ! This is the only one which was not overloaded.
 !end interface

  interface iand
    module procedure iand_vec
  end interface

  interface ior
    module procedure ior_vec
  end interface

  interface ieor
    module procedure ieor_vec
  end interface

  interface not
    module procedure not_vec
  end interface

  interface ishft
    module procedure ishft_vec
  end interface

  interface popcnt
    module procedure popcnt_vec
  end interface

  interface poppar
    module procedure poppar_vec
  end interface

  interface max
    module procedure max_vec
  end interface

  interface maxval
    module procedure maxval_vec
  end interface

! interface disp
!   module procedure disp_int, disp_vec
! end interface

! interface comp
!   module procedure compare_int_to_vec,compare_ik_to_vec,compare_int_to_int,compare_rk_to_rk,compare_log_to_log
! end interface

contains

  subroutine assign_from_vec(det1,det2)
  ! A Holmes, 22 Apr 2013.
    type(ik_vec),intent(out) :: det1
    type(ik_vec),intent(in) :: det2
    integer :: i

    do i=1,num_words
      det1%v(i) = det2%v(i)
    enddo

  end subroutine assign_from_vec

  subroutine assign_from_int(det,itgr)
  ! A Holmes, 22 Apr 2013.
    integer,intent(in) :: itgr
    type(ik_vec),intent(out) :: det

    det%v(:) = 0
    det%v(1) = itgr

  end subroutine assign_from_int

  subroutine assign_from_ik(det,itgr)
  ! A Holmes, 22 Apr 2013.
    integer(ik),intent(in) :: itgr
    type(ik_vec),intent(out) :: det

    det%v(:) = 0
    det%v(1) = itgr

  end subroutine assign_from_ik

  subroutine assign_from_rk(det,x)
  ! A Holmes, 22 Apr 2013.
    real(rk),intent(in) :: x
    type(ik_vec),intent(out) :: det
    real(rk) :: y
    integer :: i

    det = 0
    y = x
    do i=num_words*bits_per_word,1,-1
      if (y>2._rk**(i-1)) then
        det = ibset(det,i-1)
        y = y - 2._rk**(i-1)
      endif
    enddo

  end subroutine assign_from_rk

  subroutine assign_vec_from_int(det,itgr)
  ! A Holmes, 24 Apr 2013.
    type(ik_vec),intent(out) :: det(:)
    integer,intent(in) :: itgr
    integer :: i

    do i=1,size(det)
      det(i)%v(:) = 0
      det(i)%v(1) = itgr
    enddo

  end subroutine assign_vec_from_int

  subroutine assign_vec_from_ik(det,itgr)
  ! A Holmes, 24 Apr 2013.
    type(ik_vec),intent(out) :: det(:)
    integer(ik),intent(in) :: itgr
    integer :: i

    do i=1,size(det)
      det(i)%v(:) = 0
      det(i)%v(1) = itgr
    enddo
    
  end subroutine assign_vec_from_ik

  subroutine assign_array_from_ik(det,itgr)
  ! A Holmes, 24 Apr 2013.
    type(ik_vec),intent(out) :: det(:,:)
    integer(ik),intent(in) :: itgr
    integer :: i,j

    do i=1,size(det,1)
      do j=1,size(det,2)
        det(i,j)%v(:) = 0
        det(i,j)%v(1) = itgr
      enddo
    enddo
    
  end subroutine assign_array_from_ik

  subroutine assign_to_ik(itgr,det)
  ! A Holmes, 7 May 2013.
    integer(ik),intent(out) :: itgr
    type(ik_vec),intent(in) :: det

    itgr = det%v(1)

  end subroutine

  logical function eq_vec(det1,det2)
  ! A Holmes, 22 Apr 2013.
    type(ik_vec),intent(in) :: det1,det2
    integer :: i

    eq_vec = .false.
    do i=1,num_words
      if (det1%v(i) .ne. det2%v(i))  return
    enddo
    eq_vec = .true.

  end function eq_vec

  logical function eq_int(det,itgr)
  ! A Holmes, 22 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer,intent(in) :: itgr
    integer :: i

    eq_int = .false.
    if (det%v(1) .ne. itgr)  return
    do i=2,num_words
      if (det%v(i) .ne. 0_ik)  return
    enddo
    eq_int = .true.

  end function eq_int

  logical function eq_ik(det,itgr)
  ! A Holmes, 22 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer(ik),intent(in) :: itgr
    integer :: i

    eq_ik = .false.
    if (det%v(1) .ne. itgr)  return
    do i=2,num_words
      if (det%v(i) .ne. 0_ik)  return
    enddo
    eq_ik = .true.

  end function eq_ik

  logical function ne_vec(det1,det2)
  ! A Holmes, 22 Apr 2013.
    type(ik_vec),intent(in) :: det1,det2
    integer :: i

    ne_vec = .true.
    do i=1,num_words
      if (det1%v(i) .ne. det2%v(i))  return
    enddo
    ne_vec = .false.

  end function ne_vec

  logical function ne_int(det,itgr)
  ! A Holmes, 22 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer,intent(in) :: itgr
    integer :: i

    ne_int = .true.
    if (det%v(1) .ne. itgr)  return
    do i=2,num_words
      if (det%v(i) .ne. 0_ik)  return
    enddo
    ne_int = .false.

  end function ne_int

  logical function ne_ik(det,itgr)
  ! A Holmes, 22 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer(ik),intent(in) :: itgr
    integer :: i

    ne_ik = .true.
    if (det%v(1) .ne. itgr)  return
    do i=2,num_words
      if (det%v(i) .ne. 0_ik)  return
    enddo
    ne_ik = .false.

  end function ne_ik

  function lt_vec(det1,det2) result (answer)
  ! Returns true if det1 is less than det2, false otherwise
  ! A Holmes, 2 Apr 2013.
    type(ik_vec),intent(in) :: det1,det2
    logical :: answer
    integer :: i

    do i=num_words,1,-1
      if (det1%v(i).lt.det2%v(i)) then
        answer = .true.; return
      elseif (det1%v(i).gt.det2%v(i)) then
        answer = .false.; return
      endif
    enddo
    answer = .false. ! dets are equal

  end function lt_vec

  logical function lt_int(det,itgr)
  ! A Holmes, 22 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer,intent(in) :: itgr
    integer :: i

    lt_int = .false.
    if (det%v(1) .ge. itgr)  return
    do i=2,num_words
      if (det%v(i) .ne. 0_ik)  return
    enddo
    lt_int = .true.

  end function lt_int

  logical function lt_ik(det,itgr)
  ! A Holmes, 22 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer(ik),intent(in) :: itgr
    integer :: i

    lt_ik = .false.
    if (det%v(1) .ge. itgr)  return
    do i=2,num_words
      if (det%v(i) .ne. 0_ik)  return
    enddo
    lt_ik = .true.

  end function lt_ik

  logical function gt_vec(det1,det2)
  ! Returns true if det1 is greater than det2, false otherwise
  ! A Holmes, 2 Apr 2013.
    type(ik_vec),intent(in) :: det1,det2
    integer :: i
 
    do i=num_words,1,-1
      if (det1%v(i).gt.det2%v(i)) then
        gt_vec = .true.; return
      elseif (det1%v(i).lt.det2%v(i)) then
        gt_vec = .false.; return
      endif
    enddo
    gt_vec = .false. ! dets are equal

  end function gt_vec

  logical function gt_int(det,itgr)
  ! A Holmes, 22 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer,intent(in) :: itgr
    integer :: i

    gt_int = .true.
    do i=num_words,2,-1
      if (det%v(i) .ne. 0_ik)  return
    enddo
    if (det%v(1) .gt. itgr)  return
    gt_int = .false.

  end function gt_int

  logical function gt_ik(det,itgr)
  ! A Holmes, 22 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer(ik),intent(in) :: itgr
    integer :: i

    gt_ik = .true.
    do i=num_words,2,-1
      if (det%v(i) .ne. 0_ik)  return
    enddo
    if (det%v(1) .gt. itgr)  return
    gt_ik = .false.

  end function gt_ik

  logical function le_vec(det1,det2)
  ! Returns true if det1 is less than or equal to det2, false otherwise
  ! A Holmes, 2 Apr 2013.
    type(ik_vec),intent(in) :: det1,det2
    integer :: i
 
    do i=num_words,1,-1
      if (det1%v(i).lt.det2%v(i)) then
        le_vec = .true.; return
      elseif (det1%v(i).gt.det2%v(i)) then
        le_vec = .false.; return
      endif
    enddo
    le_vec = .true. ! dets are equal

  end function le_vec

  logical function le_int(det,itgr)
  ! A Holmes, 22 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer,intent(in) :: itgr
    integer :: i

    le_int = .false.
    do i=num_words,2,-1
      if (det%v(i) .ne. 0_ik)  return
    enddo
    if (det%v(1) .gt. itgr)  return
    le_int = .true.

  end function le_int

  logical function le_ik(det,itgr)
  ! A Holmes, 22 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer(ik),intent(in) :: itgr
    integer :: i

    le_ik = .false.
    do i=num_words,2,-1
      if (det%v(i) .ne. 0_ik)  return
    enddo
    if (det%v(1) .gt. itgr)  return
    le_ik = .true.

  end function le_ik

  logical function ge_vec(det1,det2)
  ! Returns true if det1 is greater than or equal to det2, false otherwise
  ! A Holmes, 2 Apr 2013.
    type(ik_vec),intent(in) :: det1,det2
    integer :: i
 
    do i=num_words,1,-1
      if (det1%v(i).gt.det2%v(i)) then
        ge_vec = .true.; return
      elseif (det1%v(i).lt.det2%v(i)) then
        ge_vec = .false.; return
      endif
    enddo
    ge_vec = .true. ! dets are equal

  end function ge_vec

  logical function ge_int(det,itgr)
  ! A Holmes, 22 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer,intent(in) :: itgr
    integer :: i

    ge_int = .true.
    if (det%v(1) .ge. itgr)  return
    do i=2,num_words
      if (det%v(i) .ne. 0_ik)  return
    enddo
    ge_int = .false.

  end function ge_int

  logical function ge_ik(det,itgr)
  ! A Holmes, 22 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer(ik),intent(in) :: itgr
    integer :: i

    ge_ik = .true.
    if (det%v(1) .ge. itgr)  return
    do i=2,num_words
      if (det%v(i) .ne. 0_ik)  return
    enddo
    ge_ik = .false.

  end function ge_ik

  function det_plus_det(det1,det2) result (ans)
  ! A Holmes, 30 Apr 2013.
    type(ik_vec),intent(in) :: det1,det2
    type(ik_vec) :: ans
    integer :: i

    ans = det1
    do i=1,num_words
      if (maskr(bits_per_word,ik)-ans%v(i).ge.det2%v(i)) then
        ans%v(i) = ans%v(i) + det2%v(i)
      else
        ans%v(i) = ans%v(i) + det2%v(i)
        if (btest(ans%v(i),i*(bits_per_word+1)-1)) then
          ans%v(i+1) = ans%v(i+1)+1_ik
        endif
      endif
    enddo

  end function det_plus_det

  function det_plus_int(det,itgr) result (ans)
  ! A Holmes, 23 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer,intent(in) :: itgr
    type(ik_vec) :: ans
    integer :: i
    integer :: remainder

    ans = det
    remainder = itgr
    do i=1,num_words
      ans%v(i) = det%v(i) + remainder
      if (maskr(bits_per_word,ik)-remainder.ge.det%v(i)) then
        return
      else
        remainder = 1
      endif
    enddo
    
  end function det_plus_int

  function int_plus_det(itgr,det) result (ans)
  ! A Holmes, 23 Apr 2013.
    integer,intent(in) :: itgr
    type(ik_vec),intent(in) :: det
    type(ik_vec) :: ans
    integer :: i
    integer :: remainder

    ans = det
    remainder = itgr
    do i=1,num_words
      ans%v(i) = det%v(i) + remainder
      if (maskr(bits_per_word,ik)-remainder.ge.det%v(i)) then
        return
      else
        remainder = 1
      endif
    enddo
   
  end function int_plus_det

  function det_plus_ik(det,itgr) result (ans)
  ! A Holmes, 23 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer(ik),intent(in) :: itgr
    type(ik_vec) :: ans
    integer :: i
    integer(ik) :: remainder

    ans = det
    remainder = itgr
    do i=1,num_words
      ans%v(i) = det%v(i) + remainder
      if (maskr(bits_per_word,ik)-remainder.ge.det%v(i)) then
        return
      else
        remainder = 1_ik
      endif
    enddo

  end function det_plus_ik

  function ik_plus_det(itgr,det) result (ans)
  ! A Holmes, 23 Apr 2013.
    integer(ik),intent(in) :: itgr
    type(ik_vec),intent(in) :: det
    type(ik_vec) :: ans
    integer :: i
    integer(ik) :: remainder

    ans = det
    remainder = itgr
    do i=1,num_words
      ans%v(i) = det%v(i) + remainder
      if (huge(ik)-remainder.ge.det%v(i)) then
        return
      else
        remainder = 1_ik
      endif
    enddo

  end function ik_plus_det

  function negative_det(det) result (ans)
  ! A Holmes, 23 Apr 2013.
  ! Uses the fact that negative integers are represented
  ! in binary as two's complement, i.e., -b = not(b) + 1
    type(ik_vec),intent(in) :: det
    type(ik_vec) :: ans
    integer :: i

    do i=1,num_words
      ans%v(i) = not(det%v(i))
    enddo
    do i=1,num_words
      ans%v(i)=ibclr(ans%v(i),bits_per_word)
    enddo
    ans = ans + 1_ik

  end function negative_det

  function det_minus_det(det1,det2) result (ans)
  ! A Holmes, 23 Apr 2013.
    type(ik_vec),intent(in) :: det1,det2
    type(ik_vec) :: ans
    integer :: i

    ans = det1
    do i=1,num_words
      if (ans%v(i)>=det2%v(i)) then
        ans%v(i) = ans%v(i) - det2%v(i)
      else
        ans%v(i) = ans%v(i) - det2%v(i)
        ans%v(i+1) = ans%v(i+1) - 1_ik
      endif
    enddo

  end function det_minus_det

  function det_minus_int(det,itgr) result (ans)
  ! A Holmes, 23 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer,intent(in) :: itgr
    type(ik_vec) :: ans
    integer :: remainder
    integer :: i

    ans = det
    remainder = itgr
    do i=1,num_words
      if (det%v(i).ge.remainder) then
        ans%v(i) = det%v(i) - remainder
        return
      else
        ans%v(i) = det%v(i)-remainder
        remainder = 1
      endif
    enddo
    
  end function det_minus_int

  function det_minus_ik(det,itgr) result (ans)
  ! A Holmes, 23 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer(ik),intent(in) :: itgr
    type(ik_vec) :: ans
    integer(ik) :: remainder
    integer :: i

    ans = det
    remainder = itgr
    do i=1,num_words
      if (det%v(i).ge.remainder) then
        ans%v(i) = det%v(i) - remainder
        return
      else
        ans%v(i) = det%v(i)-remainder
        remainder = 1_ik
      endif
    enddo
    
  end function det_minus_ik

  function det_times_rk(det,x) result (ans)
  ! Good to a few parts in 10^16.
  ! A Holmes, 23 Apr 2013.
    type(ik_vec),intent(in) :: det
    real(rk),intent(in) :: x
    real(rk) :: ans
    integer :: i

    ans = 0._rk
    do i=1,num_words*bits_per_word
      if (btest(det,i-1))  ans = ans + 2._rk**(i-1)
    enddo
    ans = ans * x

  end function det_times_rk

  function int_times_det(i,det) result (ans)
  ! Returns an integer(ik), the idea being that this will only
  ! ever be called in the hubbard density matrix code.
  ! A Holmes, 24 Apr 2013.
    integer,intent(in) :: i
    type(ik_vec),intent(in) :: det
    integer(ik) :: ans

    ans = i*det%v(1)

  end function int_times_det

  function ik_times_det(i,det) result (ans)
  ! Returns an integer(ik), the idea being that this will only
  ! ever be called in the hubbard density matrix code.
  ! A Holmes, 24 Apr 2013.
    integer(ik),intent(in) :: i
    type(ik_vec),intent(in) :: det
    integer(ik) :: ans

    ans = i*det%v(1)

  end function ik_times_det

  function det_div_det(det1,det2) result (Q)
  ! Integer division for type(ik_vec).
  ! Assumes that det1,det2 > 0.
  ! A Holmes, 23 Apr 2013.
    type(ik_vec),intent(in) :: det1,det2
    type(ik_vec) :: Q,R
    integer :: i

    if (det2 == 0) then
      write (6,*) "STOP Attempted division by zero"; call flush(6)
      stop "Attempted division by zero"
    endif

    Q = 0
    R = 0
    do i = num_words*bits_per_word-1,0,-1
      R = ishft(R,1)
      if (btest(det1,i))  R = ibset(R,0)
      if (R >= det2) then
        R = R - det2
        Q = ibset(Q,i)
      endif
    enddo

  end function det_div_det

  function abs_vec(det) result (a)
  ! Absolute value for type(ik_vec).
  ! This is only used because we now have to sort
  ! by absolute value of det label, since we want
  ! to use norb = number of bits allowed by signed
  ! integer type.
  ! A Holmes 29 Jan 2015.
    type(ik_vec),intent(in) :: det
    type(ik_vec) :: a
    integer :: i

    if (det%v(num_words) <= -1_ik) then
      do i=1,num_words
        a%v(i) = -det%v(i)
      enddo
    else
      a = det
    endif

  end function abs_vec

  logical function btest_vec(det,i)
  ! Same as btest, but for vectors of ik integers
  ! A Holmes, 2 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer,intent(in) :: i
    integer :: j,k

    k=i
    do j=1,num_words
      if (k<bits_per_word) then
        btest_vec = btest(det%v(j),k)
        return
      else
        k = k - bits_per_word
      endif
    enddo
    write (6,*) "STOP btest_vec not working"; call flush(6)
    stop "btest_vec not working"

  end function btest_vec

  function ibset_vec(det,i) result (newdet)
  ! Same as ibset, but for vectors of ik integers
  ! A Holmes, 2 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer,intent(in) :: i
    type(ik_vec) :: newdet
    integer :: j,k

    newdet = det
    k=i
    do j=1,num_words
      if (k<bits_per_word) then
        newdet%v(j) = ibset(det%v(j),k)
        return
      else
        k = k - bits_per_word
      endif
    enddo
    write (6,*) "i=",i
    write (6,*) "STOP ibset_vec not working"; call flush(6)
    stop "ibset_vec not working"

  end function ibset_vec

  function ibclr_vec(det,i) result (newdet)
  ! Same as ibclr, but for vectors of ik integers
  ! A Holmes, 2 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer,intent(in) :: i
    type(ik_vec) :: newdet
    integer :: j,k

    newdet = det
    k=i
    do j=1,num_words
      if (k<bits_per_word) then
        newdet%v(j) = ibclr(det%v(j),k)
        return
      else
        k = k - bits_per_word
      endif
    enddo
    write (6,*) "STOP ibclr_vec not working"; call flush(6)
    stop "ibclr_vec not working"

  end function ibclr_vec

  integer function trailz_vec(det)
  ! Returns number of trailing zero bits of a vector of ik integers
  ! A Holmes, 2 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer :: i

    trailz_vec = 0
    do i=1,num_words
      trailz_vec = trailz_vec + trailz(det%v(i))
      if (det%v(i).ne.0_ik)  return
      trailz_vec = trailz_vec - 1 ! Don't count unused bit.
    enddo

  end function trailz_vec

  function maskr_vec(i) result (det)
  ! A Holmes, 23 Apr 2013.
    integer,intent(in) :: i
    type(ik_vec) :: det
    integer :: j,k

    det = 0
    j = i/bits_per_word
    do k=1,j
      det%v(k) = maskr(bits_per_word,ik)
    enddo
    det%v(j+1) = maskr(i-j*bits_per_word,ik)

  end function maskr_vec

  function iand_vec(det1,det2)
  ! Returns bit-wise 'and' of two vectors of ik integers
  ! A Holmes, 2 Apr 2013.
    type(ik_vec),intent(in) :: det1,det2
    type(ik_vec) :: iand_vec
    integer :: i

    do i=1,num_words
      iand_vec%v(i) = iand(det1%v(i),det2%v(i))
    enddo

  end function iand_vec

  function ior_vec(det1,det2)
  ! Returns bit-wise 'or' of two vectors of ik integers
  ! A Holmes, 2 Apr 2013.
    type(ik_vec),intent(in) :: det1,det2
    type(ik_vec) :: ior_vec
    integer :: i

    do i=1,num_words
      ior_vec%v(i) = ior(det1%v(i),det2%v(i))
    enddo

  end function ior_vec

  function ieor_vec(det1,det2)
  ! Returns bit-wise 'or' of two vectors of ik integers
  ! A Holmes, 2 Apr 2013.
    type(ik_vec),intent(in) :: det1,det2
    type(ik_vec) :: ieor_vec
    integer :: i

    do i=1,num_words
      ieor_vec%v(i) = ieor(det1%v(i),det2%v(i))
    enddo

  end function ieor_vec

  function not_vec(det)
  ! A Holmes 23 Apr 2013.
    type(ik_vec),intent(in) :: det
    type(ik_vec) :: not_vec
    integer :: i

    do i=1,num_words
      not_vec%v(i) = not(det%v(i))
    enddo

  end function not_vec

  recursive function ishft_vec(det,i) result (answer)
  ! Same as ishft, but for vectors of ik integers, and subroutine instead of function
  ! A Holmes, 2 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer,intent(in) :: i
    type(ik_vec) :: answer
    integer :: j,k

    if (i==0) then
      answer = det
    elseif (i==bits_per_word) then
      do j=2,num_words
        answer%v(j) = det%v(j-1)
      enddo
      answer%v(1) = 0_ik
    elseif (i==-bits_per_word) then
      do j=2,num_words
        answer%v(j-1) = det%v(j)
      enddo
      answer%v(num_words) = 0_ik
    elseif (i>bits_per_word) then
      answer = ishft_vec(ishft_vec(det,i-bits_per_word),bits_per_word)
    elseif (i<-bits_per_word) then
      answer = ishft_vec(ishft_vec(det,i+bits_per_word),-bits_per_word)
    elseif (i>0) then
      do j=num_words,1,-1
        answer%v(j) = ishft(det%v(j),i)
        if (j>1) then
          do k=1,i
            if (btest(det%v(j-1),bits_per_word-k)) then
              answer%v(j) = ibset(answer%v(j),i-k)
            endif
          enddo
        endif
      enddo
    else ! i<0
      do j=1,num_words
        answer%v(j) = ishft(det%v(j),i)
        if (j<num_words) then
          do k=1,abs(i)
            if (btest(det%v(j+1),k-1)) then
              answer%v(j) = ibset(answer%v(j),bits_per_word-abs(i)+k-1)
            endif
          enddo
        endif
      enddo
    endif

    do j=1,num_words
      if (btest(answer%v(j),bits_per_word)) then
        answer%v(j) = ibclr(answer%v(j),bits_per_word) ! This bit is unused anyway.
      endif
    enddo

  end function ishft_vec

  function popcnt_vec(det)
  ! A Holmes, 23 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer :: popcnt_vec
    integer :: i

    popcnt_vec = popcnt(det%v(1))
    do i=2,num_words
      popcnt_vec = popcnt_vec + popcnt(det%v(i))
    enddo

  end function popcnt_vec

  function poppar_vec(det)
  ! A Holmes, 23 Apr 2013.
    type(ik_vec),intent(in) :: det
    integer :: poppar_vec
    integer :: i

    poppar_vec = poppar(det%v(1))
    do i=2,num_words
      poppar_vec = abs(poppar_vec - poppar(det%v(i)))
    enddo

  end function poppar_vec

  function max_vec(det1,det2) result (ans)
  ! A Holmes, 24 Apr 2013.
    type(ik_vec),intent(in) :: det1,det2
    type(ik_vec) :: ans

    if (det1>det2) then
      ans = det1
    else
      ans = det2
    endif

  end function max_vec

  function maxval_vec(det) result (ans)
  ! A Holmes, 24 Apr 2013.
    type(ik_vec),intent(in) :: det(:)
    type(ik_vec) :: ans
    integer :: i

    ans = det(1)
    do i=2,size(det)
      if (det(i)>ans)  ans = det(i)
    enddo

  end function maxval_vec

  ! Test functions:

! function i16b_to_ik_vec(from) result(to)
! ! A Holmes, 26 Apr 2013.
!   use types, only : i16b
!   integer(i16b),intent(in) :: from
!   type(ik_vec) :: to
!   integer :: i

!   to = 0_ik
!   do i=1,128
!     if (btest(from,i-1))  to = ibset(to,i-1)
!   enddo

! end function i16b_to_ik_vec

! function ik_vec_to_i16b(from) result(to)
! ! A Holmes, 26 Apr 2013.
!   use types, only : i16b
!   type(ik_vec),intent(in) :: from
!   integer(i16b) :: to
!   integer :: i

!   to = 0_ik
!   do i=1,128
!     if (btest(from,i-1))  to = ibset(to,i-1)
!   enddo

! end function ik_vec_to_i16b

! logical function disp_int(itgr)
! ! A Holmes, 26 Apr 2013.
!   use types, only : i16b
!   integer(i16b),intent(in) :: itgr
!   integer :: i,j,v1(63),v2(63),v(63) !v(128)

!   v1(:) = 0
!   v2(:) = 0
!   do i=1,bits_per_word
!     v(i)=mod(i,10)
!     if (btest(itgr,i-1))  v1(i)=1
!   enddo
!   do i=bits_per_word+1,2*bits_per_word
!     if (btest(itgr,i-1))  v2(i-bits_per_word)=1
!   enddo
!   write(6,'(/,63i1,'' '',63i1)') (v(i),i=1,63),(v(j),j=1,63)
!   write(6,'(/,63i1,'' '',63i1)') (v1(i),i=1,63),(v2(j),j=1,63)
!   disp_int = .true.

! end function disp_int

! logical function disp_vec(det)
! ! A Holmes, 26 Apr 2013.
!   use types, only : i16b
!   type(ik_vec),intent(in) :: det
!   integer :: i,j,v1(63),v2(63),v(63)

!   v1(:) = 0
!   v2(:) = 0
!   do i=1,bits_per_word
!     v(i)=mod(i,10)
!     if (btest(det,i-1))  v1(i)=1
!   enddo
!   do i=bits_per_word+1,2*bits_per_word
!     if (btest(det,i-1))  v2(i-bits_per_word)=1
!   enddo
!   write(6,'(/,63i1,'' '',63i1)') (v(i),i=1,63),(v(j),j=1,63)
!   write(6,'(/,63i1,'' '',63i1)') (v1(i),i=1,63),(v2(j),j=1,63)
!   disp_vec = .true.

! end function disp_vec

! logical function compare_int_to_vec(itgr,det)
! ! A Holmes, 29 Apr 2013.
!   use types, only : i16b
!   integer(i16b),intent(in) :: itgr
!   type(ik_vec),intent(in) :: det
!   integer :: i

!   compare_int_to_vec = .false.
!   do i=1,num_words*bits_per_word
!     if (btest(itgr,i-1).neqv.btest(det,i-1))  return
!   enddo
!   compare_int_to_vec = .true.

! end function compare_int_to_vec

! logical function compare_ik_to_vec(itgr,det)
! ! A Holmes, 29 Apr 2013.
!   integer(ik),intent(in) :: itgr
!   type(ik_vec),intent(in) :: det
!   integer :: i

!   compare_ik_to_vec = .false.
!   do i=1,num_words*bits_per_word
!     if (btest(itgr,i-1).neqv.btest(det,i-1))  return
!   enddo
!   compare_ik_to_vec = .true.

! end function compare_ik_to_vec

! logical function compare_int_to_int(itgr1,itgr2)
! ! A Holmes, 30 Apr 2013.
!   integer,intent(in) :: itgr1,itgr2

!   compare_int_to_int = (itgr1==itgr2)

! end function compare_int_to_int

! logical function compare_rk_to_rk(rk1,rk2)
! ! A Holmes, 30 Apr 2013.
!   real(rk),intent(in) :: rk1,rk2

!   compare_rk_to_rk = ((rk1==rk2).or.(2._rk*(abs(rk1-rk2)/(rk1+rk2))<(3._rk*10._rk**-16)))

! end function compare_rk_to_rk

! logical function compare_log_to_log(log1,log2)
! ! A Holmes, 30 Apr 2013.
!   logical,intent(in) :: log1,log2

!   compare_log_to_log = (log1.eqv.log2)

! end function compare_log_to_log

#endif

end module overload
