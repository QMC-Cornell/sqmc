! Define the module for the key type.
! Override the hash_value and == operator interface.
module ik_hash_module
  use types, only: ik,i16b
  implicit none

  interface hash_value
    module procedure djb_hash2
  end interface

  contains

    function djb_hash(det) result(hash)

      implicit none
      integer(kind=ik), intent(in) :: det
      integer(i16b) :: hash
      integer :: i,j
      integer(kind=ik) :: test_nk
      integer :: test_int
      integer, parameter :: n_bits_nk = bit_size(test_nk)
      integer, parameter :: n_bits = bit_size(test_int)
      INTEGER(i16b) PRIME ; PARAMETER (PRIME = 309485009821345068724781371_i16b)
      INTEGER(i16b) OFFSET ; PARAMETER (OFFSET = 144066263297769815596495629667062367629_i16b)
      hash = PRIME !91847591247101_ik

      do j = 0, n_bits_nk -1, 8
        test_nk = IEOR(det,X'5555555555555555') !X'555555555....' = 01010101 01010101 01010101 01010101 ...
        test_nk = test_nk * prime
        hash = (ishft(hash,5) + hash) + ishft(test_nk,-j)
      enddo


    end function djb_hash

    function djb_hash2(det) result(hash)
      implicit none
      integer(kind=ik), intent(in) :: det
      integer(i16b) :: hash
      INTEGER(i16b) PRIME ; PARAMETER (PRIME = 309485009821345068724781371_i16b)

      hash = det * prime

    end function djb_hash2


    function hash_value_ik(ints) result(hash)
      integer(ik), intent(in) :: ints
      integer :: hash
      integer :: i

      hash = ints
    end function

end module ik_hash_module

! Define the macros needed by fhash and include fhash.f90
#define KEY_TYPE integer(ik)
#define VALUE_TYPE type(int_vec)
#define KEY_USE use ik_hash_module
#define VALUE_USE use types, only: int_vec
#define SHORTNAME ik_int_list
#include "fhash.f90"
