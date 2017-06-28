module mpi_routines
  use types
  use common_run, only : diag_elem_info
  use common_walk, only: walk_dets_up, walk_dets_dn, walk_wt, matrix_elements, imp_distance
  !You have to use ISO_C_BINDING for the CPTR stuff
  use, intrinsic :: ISO_C_BINDING

implicit none
#ifdef MPI
INCLUDE 'mpif.h'
#endif
public :: master_core, master_core_node
public :: cluster_init, cluster_finalize
public :: mpi_bsend, mpi_bsend_diffroot
public :: mpi_allred, mpi_gath, mpi_agathv, mpi_gathv, mpi_red_max, mpi_red_min, mpi_red, mpi_stop
public :: mpi_barr, mpi_barr_in_node,t_sndlist, mpi_scattv, mpi_redscatt_real_dparray
public :: mpi_sendwalks, mpi_sendnewwalks, init_snd_buffers, mpi_merge_sort2
public :: snd_table, snd_displs, mpi_push_nwalk, snd_cnt
public :: get_owner, init_hash_owners, get_det_owner, mpi_snd
public :: shmem_allocate,t_det,MPI_DET_TYPE,MPI_DIAG_ELEM_INFO_TYPE, MPI_DET_MIX_DIAG_TYPE,conv_128_to_64,conv_64_to_128,mpi_allgatherv_new_dets
public :: shmem_deallocate,det_map,det_map_l,shmem_reallocate,mpi_distribute_remote_det_map
public :: mpi_bsend_between_nodes
private

type t_sndlist
    integer :: to,wlkid
end type t_sndlist

type t_walk
    real(rk)     :: wt
    integer(i8b) :: det_u(2*num_words),det_d(2*num_words) ! 2 i8b's instead of 1 i16b because MPI can only handle size i8b
    integer(i1b) :: imp_ini(2) !contains both imp_distance and initiator
    real(rk) :: e_num,e_den,diag_elems
end type t_walk

!MJO - Global information for the determinants - shared per node, only master_core_node should edit
type det_map
#ifdef NUM_ORBITALS_GT_127
   type(ik_vec),dimension(:),contiguous,pointer :: global_dets_up(:),global_dets_dn(:)
#else
   integer(ik),dimension(:),contiguous,pointer :: global_dets_up,global_dets_dn !MJO For shared memory global_dets - not implemented now
   !integer(ik),allocatable :: global_dets_up(:),global_dets_dn(:)
#endif
   !real(rk), allocatable   :: remote_wts(:,:)
   real(rk), dimension(:,:),contiguous,pointer :: remote_wts
 end type det_map


!MJO - Local information for the determinants - takes global list and slices out local part
type det_map_l
   integer(i8b) :: ndets_global,ndets,n_it,ndets_global_old
   real(rk)     :: average_connections
   integer(i8b), allocatable :: n_det_remote(:,:)
#ifdef NUM_ORBITALS_GT_127
   type(ik_vec) :: hf_up,hf_dn
#else
   integer(ik) :: hf_up,hf_dn
#endif
   integer (i8b), allocatable :: local_indices(:)
end type det_map_l

type t_det_mix_diag
    integer(i8b) :: det_up(2*num_words),det_dn(2*num_words) ! 2 i8b's instead of 1 i16b because MPI can only handle size i8b
    real(rk) :: e_mix_num,e_mix_den
    type(diag_elem_info) :: diag_elems_info
end type t_det_mix_diag

type t_det_mix_diag_term
    integer(i8b) :: det_up(2*num_words),det_dn(2*num_words) ! 2 i8b's instead of 1 i16b because MPI can only handle size i8b
    real(rk) :: e_mix_num,e_mix_den,term1_big,term2_big
    type(diag_elem_info) :: diag_elems_info
end type t_det_mix_diag_term

type t_det
    integer(i8b) :: det_up(2*num_words),det_dn(2*num_words) ! 2 i8b's instead of 1 i16b because MPI can only handle size i8b
end type t_det

!whoaminode is the rank within a node. master_core_node defines who is the master within a node
!sharedComm is the communicator within a node
integer,public :: whoami,mpierr,ncores,whoaminode,sharedComm,ndets_global,between_nodes_comm
integer,public :: ncores_on_node=1,nnodes=1
!declare shmemWindow here so it can be freed later
!A "window" is MPI's name for the memory that is open to other cores to see
integer, PARAMETER :: max_shmem_windows = 50 !Arbitrary maximum number of shared arrays MJO FIXME? Maybe make this dynamic
integer :: shmemWindow(max_shmem_windows) ! integrals + q,J,prob for each of same and opposite spin double excitations in heatbath
integer :: win_num
logical :: master_core,master_core_node
integer(i4b) :: nwbuff,loc_mwalk
integer,allocatable :: snd_table(:),snd_displs(:),snd_dsp(:),snd_cnt(:),rcv_cnt(:),recv_displs(:),list_cowner(:)
type(t_sndlist),allocatable :: snd_list(:)
integer istat
integer,parameter :: mpi_send_limit = 15000

!***buffers for send process [WARNING: put this in common_walk]
  real(rk), allocatable ::  snd_matrix_elements(:)

  type(t_walk), allocatable :: snd_buff(:),recv_buff(:)

  integer :: MPI_WALK_TYPE, MPI_DET_TYPE, MPI_DIAG_ELEM_INFO_TYPE, MPI_DET_MIX_DIAG_TYPE, MPI_DET_MIX_DIAG_TERM_TYPE !** needed for communications

  integer :: snd_count,recv_count


  interface mpi_merge_sort2
     module procedure mpi_merge_sort2_nowts, mpi_merge_sort2_num_denom_diag_elems_info, mpi_merge_sort2_num_denom_diag_elems_info_term,&
          &mpi_merge_sort2_nowts_replace, mpi_merge_sort2_num_denom_diag_elems_info_replace, mpi_merge_sort2_num_denom_diag_elems_info_term_replace
  end interface mpi_merge_sort2

  interface shmem_allocate
    module procedure shmem_allocate_rk, shmem_allocate_real, shmem_allocate_int, shmem_allocate_ik, shmem_allocate_ik_vec, shmem_allocate2_rk, shmem_allocate2_rs_absH, shmem_allocate_rs_absH
    module procedure shmem_allocate_logical
  end interface shmem_allocate

  interface shmem_deallocate
    module procedure shmem_deallocate_rk, shmem_deallocate_real, shmem_deallocate_int, shmem_deallocate_ik, shmem_deallocate2_rk, shmem_deallocate_ik_vec, shmem_deallocate_rs_absH
  end interface shmem_deallocate

  interface shmem_reallocate
    module procedure shmem_reallocate_rk, shmem_reallocate_ik, shmem_reallocate_ik_vec,shmem_reallocate2_rk, shmem_reallocate_rs_absH
  end interface shmem_reallocate

interface mpi_bsend
  module procedure mpi_bsend_int, mpi_bsend_real_sp,mpi_bsend_real_dp, mpi_bsend_complex,mpi_bsend_string,mpi_bsend_iarray, &!
                   mpi_bsend_logical,mpi_bsend_int64,mpi_bsend_int128,mpi_bsend_irand_seed,mpi_bsend_real_dparray,mpi_bsend_i64arr,mpi_bsend_i128arr,mpi_bsend_ikvec,mpi_bsend_ikvec_arr
end interface

interface mpi_bsend_between_nodes
  module procedure mpi_bsend_between_nodes_dp_array,mpi_bsend_between_nodes_i128arr,mpi_bsend_between_nodes_ikvec_arr, mpi_bsend_between_nodes_tdet_array
end interface mpi_bsend_between_nodes

interface mpi_snd
  module procedure mpi_snd_int,mpi_snd_real_sp,mpi_snd_real_dp,mpi_snd_real_dparray,mpi_snd_iarray,mpi_snd_i64array,&
                   mpi_snd_int64,mpi_snd_i128arr,mpi_snd_ikvec,mpi_snd_ikvec_array
end interface

interface mpi_bsend_diffroot
  module procedure mpi_bsend_diffroot_int, mpi_bsend_diffroot_real_sp,mpi_bsend_diffroot_real_dp,mpi_bsend_diffroot_logical,&
                    mpi_bsend_diffroot_int64,mpi_bsend_diffroot_ik_vec,mpi_bsend_diffroot_int128,mpi_bsend_diffroot_int8,&
                    mpi_bsend_diffroot_sndlist,mpi_bsend_diffroot_iarray
end interface

interface mpi_gath
  module procedure mpi_gath_int, mpi_gath_real_sp, mpi_gath_real_dp, mpi_gath_int8 !, mpi_gath_ik_vec, mpi_gath_i128
end interface

interface mpi_gathv
  module procedure mpi_gathv_int, mpi_gathv_real_sp, mpi_gathv_real_dp, mpi_gathv_int8, mpi_gathv_ik_vec,mpi_gathv_i128
end interface

interface mpi_scattv
  module procedure mpi_scattv_int, mpi_scattv_real_sp, mpi_scattv_real_dp, mpi_scattv_int8, mpi_scattv_ik_vec,mpi_scattv_i128,mpi_scattv_iwalk
end interface

interface mpi_agathv
  module procedure mpi_agathv_int, mpi_agathv_real_sp, mpi_agathv_real_dp, mpi_agathv_int8, mpi_agathv_ik_vec, &
                   mpi_agathv_i128
end interface

interface mpi_allred
   module procedure mpi_allred_int,mpi_allred_i8b, mpi_allred_real_sp,mpi_allred_real_dp,mpi_allred_real_dparray,mpi_allred_iarray,mpi_allred_real_dparray2!,mpi_bsend_complex,mpi_bsend_string, &
                 !  mpi_bsend_logical,mpi_bsend_int64,mpi_bsend_irand_seed!,mpi_bsend_3string
end interface

interface mpi_red_max
  module procedure mpi_red_max_real_sp, mpi_red_max_real_dp, mpi_red_max_real_dparray
end interface

interface mpi_red_min
  module procedure mpi_red_min_real_sp, mpi_red_min_real_dp, mpi_red_min_real_dparray
end interface

interface mpi_red
  module procedure mpi_red_real_sp,mpi_red_real_dp,mpi_red_real_dparray
end interface

interface conv_128_to_64
  module procedure conv_128_to_64,conv_064_to_64
end interface

interface conv_64_to_128
  module procedure conv_64_to_128,conv_64_to_064
end interface


!---- HASH TABLE part
integer :: RandomHash1(0:255) !,RandomHash2(0:255)
!---- HASH TABLE part
!--- for 128-bit integers send
integer*8,allocatable  :: i128_high(:),i128_low(:),i128_buff(:)
!--- for 128-bit integers send

contains

subroutine init_snd_buffers(mwalk)
    integer(i8b), intent(in) :: mwalk
    integer :: icore
#ifdef MPI
    loc_mwalk=int(mwalk,i4b)
    allocate(snd_list(mwalk),stat=istat)
    if(istat.ne.0) stop 'failed to allocate snd_list'
    allocate(list_cowner(mwalk),stat=istat)
    if(istat.ne.0) stop 'failed to allocate list_cowner'
    allocate(snd_matrix_elements(mwalk),stat=istat)
    if(istat.ne.0) stop 'failed to allocate snd_matrix_elements'

    allocate(snd_buff(mwalk),stat=istat)
    if(istat.ne.0) stop 'failed to allocate snd_buff'
    allocate(recv_buff(mwalk),stat=istat)
    if(istat.ne.0) stop 'failed to allocate recv_buff'

    nwbuff = int(mwalk/ncores,i4b)

    do icore=0,(ncores-1)
      snd_dsp(icore+1)=icore*nwbuff
    enddo

#endif

end subroutine init_snd_buffers
!------------------------------------------------------

!---- HASH TABLE part [based on George Booth's code]
subroutine init_hash_owners()
    implicit none
    integer :: i,j,map
    logical :: found
    real(rk) rannyu

#ifdef MPI
    if(ncores>1) then
      if(master_core) then
        do i=0,255
            RandomHash1(i) = i
        enddo

        !Find mapping integer function from range 0 -> 255 to unique integers in range 0 -> 255*20,000
        do i=0,255
            found=.false.
            do while(.not.found)
                map=int(255*rannyu()*20000)

                do j=0,i-1
                    if(RandomHash1(j).eq.map) exit
                enddo
                if(j.eq.i) found=.true.
            enddo
            RandomHash1(i)=map
        enddo
      endif
      call mpi_bsend(RandomHash1)
    endif
#endif
end subroutine
!------------------------------------------------------

!Return a hash between the range 0 -> 'range-1' from the bit-representation 'det'
subroutine hash(det,range,hashx,RandomHash)
    use types, only : ik
    implicit none
    integer(kind=ik), intent(in) :: det(2*num_words)
    integer, intent(in) :: range,RandomHash(0:255)
    integer, intent(out) :: hashx
    integer(kind=ik) :: acc
    integer(kind=ik) :: test_nk
    integer :: test_int,i,j
    integer :: val
    integer, parameter :: n_bits_nk = bit_size(test_nk)
    integer, parameter :: n_bits = bit_size(test_int)

!    acc = 0!91847591247101_ik
!    do i=1,2*num_words    !run over integers defining the bit-string
!        do j = 0, n_bits_nk -1, 8
!            !val = int(iand(ishft(det(i),-j), int(255,ik)),kind(test_int))
!            val = int(iand(ishft(det(i)*1099511628211_ik,-j), int(255,ik)),kind(test_int)) ! This change from Matt Otten has better load balancing
!            !val is now a number between 0 -> 255
!            !Ensure that RandomHash has a mapping for the 0th element too
!            !1099511628211_ik = ibset(0_ik,27)
!            acc = (1099511628211_ik * acc) + (RandomHash(val) * (i + j))
!        enddo
!    enddo
    ! hashx = int(abs(mod(acc,int(range,ik))),kind(test_int))

    ! call FNV128(det,32,acc)
    ! hashx = int(abs(mod(acc,int(range,ik))),kind(test_int))

    acc = djb_hash(det) !MJO Temporary other hash function
    hashx = int(abs(mod(acc,int(range,ik))),kind(test_int))

end subroutine hash
!------------------------------------------------------

!Return a hash between the range 0 -> 'range-1' from the bit-representation 'det'
subroutine hash_1(det,range,hashx,RandomHash)
    use types, only : ik
    implicit none
    integer(kind=ik), intent(in) :: det(num_words)
    integer, intent(in) :: range,RandomHash(0:255)
    integer, intent(out) :: hashx
    integer(kind=ik) :: acc
    integer(kind=ik) :: test_nk
    integer :: test_int,i,j
    integer :: val
    integer, parameter :: n_bits_nk = bit_size(test_nk)
    integer, parameter :: n_bits = bit_size(test_int)

!    acc = 0!91847591247101_ik
!    do i=1,2*num_words    !run over integers defining the bit-string
!        do j = 0, n_bits_nk -1, 8
!            !val = int(iand(ishft(det(i),-j), int(255,ik)),kind(test_int))
!            val = int(iand(ishft(det(i)*1099511628211_ik,-j), int(255,ik)),kind(test_int)) ! This change from Matt Otten has better load balancing
!            !val is now a number between 0 -> 255
!            !Ensure that RandomHash has a mapping for the 0th element too
!            !1099511628211_ik = ibset(0_ik,27)
!            acc = (1099511628211_ik * acc) + (RandomHash(val) * (i + j))
!        enddo
!    enddo
    ! hashx = int(abs(mod(acc,int(range,ik))),kind(test_int))

    ! call FNV128(det,32,acc)
    ! hashx = int(abs(mod(acc,int(range,ik))),kind(test_int))

    acc = djb_hash_1(det) !MJO Temporary other hash function
    hashx = int(abs(mod(acc,int(range,ik))),kind(test_int))

end subroutine hash_1
!------------------------------------------------------
function djb_hash_1(det) result(hash)

  implicit none
  integer(kind=ik), intent(in) :: det(num_words)
  integer(i16b) :: hash
  integer :: i,j
  integer(kind=ik) :: test_nk
  integer :: test_int
  integer, parameter :: n_bits_nk = bit_size(test_nk)
  integer, parameter :: n_bits = bit_size(test_int)
  INTEGER(i16b) PRIME ; PARAMETER (PRIME = 309485009821345068724781371_i16b)
  INTEGER(i16b) OFFSET ; PARAMETER (OFFSET = 144066263297769815596495629667062367629_i16b)
  hash = PRIME !91847591247101_ik

  do i=1,num_words
    do j = 0, n_bits_nk -1, 8
      test_nk = IEOR(det(i),X'5555555555555555') !X'555555555....' = 01010101 01010101 01010101 01010101 ...
      test_nk = test_nk * prime
      !test_nk = det(i)
      !hash = IEOR((ishft(hash,5) + hash),ishft(test_nk,-j))
      hash = (ishft(hash,5) + hash) + ishft(test_nk,-j)
    enddo
  enddo

end function DJB_hash_1


function djb_hash(det) result(hash)

    implicit none
    integer(kind=ik), intent(in) :: det(2*num_words)
    integer(i16b) :: hash
    integer :: i,j
    integer(kind=ik) :: test_nk,tmp
    integer :: test_int
    integer, parameter :: n_bits_nk = bit_size(test_nk)
    integer, parameter :: n_bits = bit_size(test_int)
    INTEGER(i16b) PRIME ; PARAMETER (PRIME = 309485009821345068724781371_i16b)
    INTEGER(i16b) OFFSET ; PARAMETER (OFFSET = 144066263297769815596495629667062367629_i16b)
    hash = PRIME !91847591247101_ik
    do i=1,2*num_words
      tmp = det(i)
      if (i.eq.2) tmp = tmp + offset
      do j = 0, n_bits_nk -1, 8
        test_nk = IEOR(tmp,X'5555555555555555') !X'555555555....' = 01010101 01010101 01010101 01010101 ...
        test_nk = test_nk * prime
        !test_nk = det(i)
        !hash = IEOR((ishft(hash,5) + hash),ishft(test_nk,-j))
        hash = (ishft(hash,5) + hash) + ishft(test_nk,-j)
      enddo
    enddo

end function DJB_hash

!------------------------------------------------------

SUBROUTINE FNV128 (BUFFER, LENGTH, HASH)
  IMPLICIT NONE
  INTEGER(i16b) HASH
  INTEGER LENGTH
  INTEGER (ik) BUFFER,B
  DIMENSION BUFFER(LENGTH)

  INTEGER(i16b) PRIME ; PARAMETER (PRIME = 309485009821345068724781371_i16b)
  INTEGER(i16b) OFFSET ; PARAMETER (OFFSET = 144066263297769815596495629667062367629_i16b)
  INTEGER I, J, K, L


  ! *#######################################################################
  ! *                begin
  ! *#######################################################################
  ! *          FNV-1a hash each octet in the buffer
  HASH = OFFSET
  DO J = 1, LENGTH
    DO L = 1,16 !16 because 16 byte input
      B = BUFFER(J)
      K = 0
      DO I = 0, 7           ! copy each bit from B to K
        IF (BTEST(B, I+(L-1)*16)) K = IBSET(K, I)
      ENDDO

      !*          xor the bottom with the current octet
      HASH = IEOR(HASH, K)

      !*          multiply by the 32 bit FNV magic prime mod 2^32
      HASH = HASH * PRIME
!      HASH = IAND(HASH, X'FFFFFFFF')      ! discard > 32 bits
    ENDDO
  ENDDO
END SUBROUTINE FNV128

!------------------------------------------------------
function get_det_owner(det_up,det_dn) result(coreid)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: det_up,det_dn
#else
    integer(kind=ik), intent(in) :: det_up,det_dn
#endif
    integer :: coreid
    integer(kind=ik) :: det(2*num_words)
#ifdef MPI
    if(ncores==1) then
      coreid=0
    else
#ifdef NUM_ORBITALS_GT_127
      det(1:num_words)=det_up%v(1:num_words)
      det(num_words+1:2*num_words)=det_dn%v(1:num_words)
#else
      det(1)=det_up
      det(2)=det_dn
      !det(2)=not(det_up)
#endif
      call hash(det,ncores,coreid,RandomHash1)
      !call hash_1(det,ncores,coreid,RandomHash1)
    endif
#else
    coreid = 0
#endif
end function get_det_owner
!------------------------------------------------------

function get_owner(windex,wtot) result(coreid)
integer, intent(in) :: windex,wtot
integer :: coreid
integer :: cscale,ic

    cscale=wtot/ncores
    coreid=-1
    do ic=0,(ncores-2)
        if((windex .ge. (ic*cscale+1)).and.(windex .le. ((ic+1)*cscale))) then
            coreid=ic
            return
        endif
    enddo
    if(coreid < 0) coreid=ncores-1 !*just for now that the last one has different upper-bound

!coreid=MOD(windex,2)
end function get_owner
!------------------------------------------------------
!---- End of HASH TABLE part --------------------------

subroutine init_mpitype_walker

  integer :: block_len(0:6)
  integer :: oldtypes(0:6),offsets(0:6)
  integer :: i8extent,r8extent,i1extent

#ifdef MPI
  call MPI_TYPE_EXTENT(MPI_INTEGER8, i8extent, mpierr)
  call MPI_TYPE_EXTENT(MPI_REAL8, r8extent, mpierr)
  call MPI_TYPE_EXTENT(MPI_INTEGER1, i1extent, mpierr)

  oldtypes(0) = MPI_REAL8    ! wt
  oldtypes(1) = MPI_INTEGER8 ! det_up
  oldtypes(2) = MPI_INTEGER8 ! det_dn
  oldtypes(3) = MPI_INTEGER1 ! combination of imp_distance and initiator (?)
  oldtypes(4) = MPI_REAL8    ! e_num
  oldtypes(5) = MPI_REAL8    ! e_den
  oldtypes(6) = MPI_REAL8    ! diagonal Hamiltonian element

  block_len(0) = 1
#ifdef NUM_ORBITALS_GT_127
  block_len(1) = 2*num_words
  block_len(2) = 2*num_words
#else
  block_len(1) = 2
  block_len(2) = 2
#endif
  block_len(3) = 2
  block_len(4) = 1
  block_len(5) = 1
  block_len(6) = 1

  offsets(0) = 0
  offsets(1) = r8extent*block_len(0)
  offsets(2) = offsets(1) + i8extent*block_len(1)
  offsets(3) = offsets(2) + i8extent*block_len(2)
  offsets(4) = offsets(3) + i1extent*block_len(3)
  offsets(5) = offsets(4) + r8extent*block_len(4)
  offsets(6) = offsets(5) + r8extent*block_len(5)

  call MPI_TYPE_STRUCT(7,block_len,offsets,oldtypes,MPI_WALK_TYPE,mpierr)
  call MPI_TYPE_COMMIT(MPI_WALK_TYPE,mpierr)
#endif

end subroutine init_mpitype_walker

!------------------------------------------------------
subroutine init_mpitype_det_mix_diag

  integer :: block_len(5)
  integer :: oldtypes(5),offsets(5)
  integer :: i8extent,r8extent,diag_elem_info_extent

#ifdef MPI
  call MPI_TYPE_EXTENT(MPI_INTEGER8, i8extent, mpierr)
  call MPI_TYPE_EXTENT(MPI_REAL8, r8extent, mpierr)
  call MPI_TYPE_EXTENT(MPI_DIAG_ELEM_INFO_TYPE, diag_elem_info_extent, mpierr)


  oldtypes(1) = MPI_INTEGER8 ! det_up
  oldtypes(2) = MPI_INTEGER8 ! det_dn
  oldtypes(3) = MPI_REAL8    ! e_num
  oldtypes(4) = MPI_REAL8    ! e_den
  oldtypes(5) = MPI_DIAG_ELEM_INFO_TYPE    ! diagonal Hamiltonian element

#ifdef NUM_ORBITALS_GT_127
  block_len(1) = 2*num_words
  block_len(2) = 2*num_words
#else
  block_len(1) = 2
  block_len(2) = 2
#endif
  block_len(3) = 1
  block_len(4) = 1
  block_len(5) = 1

  offsets(1) = 0
  offsets(2) = offsets(1) + i8extent*block_len(1)
  offsets(3) = offsets(2) + i8extent*block_len(2)
  offsets(4) = offsets(3) + r8extent*block_len(3)
  offsets(5) = offsets(4) + r8extent*block_len(4)

  call MPI_TYPE_STRUCT(5,block_len,offsets,oldtypes,MPI_DET_MIX_DIAG_TYPE,mpierr)
  call MPI_TYPE_COMMIT(MPI_DET_MIX_DIAG_TYPE,mpierr)
#endif

end subroutine init_mpitype_det_mix_diag

!------------------------------------------------------
subroutine init_mpitype_det_mix_term_diag

  integer :: block_len(7)
  integer :: oldtypes(7),offsets(7)
  integer :: i8extent,r8extent,diag_elem_info_extent

#ifdef MPI
  call MPI_TYPE_EXTENT(MPI_INTEGER8, i8extent, mpierr)
  call MPI_TYPE_EXTENT(MPI_REAL8, r8extent, mpierr)
  call MPI_TYPE_EXTENT(MPI_DIAG_ELEM_INFO_TYPE, diag_elem_info_extent, mpierr)


  oldtypes(1) = MPI_INTEGER8 ! det_up
  oldtypes(2) = MPI_INTEGER8 ! det_dn
  oldtypes(3) = MPI_REAL8    ! e_num
  oldtypes(4) = MPI_REAL8    ! e_den
  oldtypes(5) = MPI_REAL8    ! e_num
  oldtypes(6) = MPI_REAL8    ! e_num
  oldtypes(7) = MPI_DIAG_ELEM_INFO_TYPE    ! diagonal Hamiltonian element

#ifdef NUM_ORBITALS_GT_127
  block_len(1) = 2*num_words
  block_len(2) = 2*num_words
#else
  block_len(1) = 2
  block_len(2) = 2
#endif
  block_len(3) = 1
  block_len(4) = 1
  block_len(5) = 1
  block_len(6) = 1
  block_len(7) = 1


  offsets(1) = 0
  offsets(2) = offsets(1) + i8extent*block_len(1)
  offsets(3) = offsets(2) + i8extent*block_len(2)
  offsets(4) = offsets(3) + r8extent*block_len(3)
  offsets(5) = offsets(4) + r8extent*block_len(4)
  offsets(6) = offsets(5) + r8extent*block_len(5)
  offsets(7) = offsets(6) + r8extent*block_len(6)

  call MPI_TYPE_STRUCT(7,block_len,offsets,oldtypes,MPI_DET_MIX_DIAG_TERM_TYPE,mpierr)
  call MPI_TYPE_COMMIT(MPI_DET_MIX_DIAG_TERM_TYPE,mpierr)
#endif

end subroutine init_mpitype_det_mix_term_diag

subroutine init_mpitype_diag_elem_info

  integer :: block_len(5)
  integer :: oldtypes(5),offsets(5)
  integer :: iextent,r8extent

#ifdef MPI
  call MPI_TYPE_EXTENT(MPI_INTEGER, iextent, mpierr)
  call MPI_TYPE_EXTENT(MPI_REAL8, r8extent, mpierr)

  oldtypes(1) = MPI_REAL8 ! old_diag_elem
  oldtypes(2) = MPI_INTEGER ! p
  oldtypes(3) = MPI_INTEGER ! q
  oldtypes(4) = MPI_INTEGER ! r
  oldtypes(5) = MPI_INTEGER ! s

  block_len(1) = 1
  block_len(2) = 1
  block_len(3) = 1
  block_len(4) = 1
  block_len(5) = 1

  offsets(1) = 0
  offsets(2) = offsets(1) + r8extent*block_len(1)
  offsets(3) = offsets(2) + iextent*block_len(2)
  offsets(4) = offsets(3) + iextent*block_len(3)
  offsets(5) = offsets(4) + iextent*block_len(4)

  call MPI_TYPE_STRUCT(5,block_len,offsets,oldtypes,MPI_DIAG_ELEM_INFO_TYPE,mpierr)
  call MPI_TYPE_COMMIT(MPI_DIAG_ELEM_INFO_TYPE,mpierr)
#endif

end subroutine init_mpitype_diag_elem_info

!------------------------------------------------------
subroutine init_mpitype_det

  integer :: block_len(2)
  integer :: oldtypes(2),offsets(2)
  integer :: i8extent

#ifdef MPI
  call MPI_TYPE_EXTENT(MPI_INTEGER8, i8extent, mpierr)

  oldtypes(1) = MPI_INTEGER8 ! det_up
  oldtypes(2) = MPI_INTEGER8 ! det_dn

#ifdef NUM_ORBITALS_GT_127
  block_len(1) = 2*num_words
  block_len(2) = 2*num_words
#else
  block_len(1) = 2
  block_len(2) = 2
#endif
  offsets(1) = 0
  offsets(2) = offsets(1) + i8extent*block_len(1)

  call MPI_TYPE_STRUCT(2,block_len,offsets,oldtypes,MPI_DET_TYPE,mpierr)
  call MPI_TYPE_COMMIT(MPI_DET_TYPE,mpierr)
#endif

end subroutine init_mpitype_det


!------------------------------------------------------


!*** these routines shall be moved elsewhere
subroutine conv_128_to_64(i128,i64_1,i64_2)
integer*16, intent(in) :: i128
integer*8, intent(out) :: i64_1,i64_2
! integer :: ibit

integer*8 :: arr(2)

arr = TRANSFER(i128, arr)
i64_1=arr(1)
i64_2=arr(2)
!----
!i64_1=0
!i64_2=0
!do ibit = 0,63
!    if(btest(i128,ibit)) i64_1=ibset(i64_1,ibit)
!enddo
!do ibit = 64,127
!    if(btest(i128,ibit)) i64_2=ibset(i64_2,ibit-64)
!enddo
end subroutine conv_128_to_64
!------------------------------------------------------

subroutine conv_64_to_128(i128,i64_1,i64_2)
integer*16, intent(out) :: i128
integer*8, intent(in) :: i64_1,i64_2
integer :: ibit

! integer*16 :: lowbits

!i128=0
!lowbits=0
!i128 = i64_2
!i128=ISHFT(i128,63)
!lowbits = i64_1
!i128=IOR(i128,lowbits)
!---
i128=0
do ibit = 0,63
    if(btest(i64_1,ibit)) i128=ibset(i128,ibit)
enddo
do ibit = 0,63
    if(btest(i64_2,ibit)) i128=ibset(i128,ibit+64)
enddo
end subroutine conv_64_to_128
!------------------------------------------------------

! These are same as above, but for when ik = i8b
subroutine conv_064_to_64(i128,i64_1,i64_2)
integer*8, intent(in) :: i128
integer*8, intent(out) :: i64_1,i64_2
! integer :: ibit

integer*8 :: arr(2)

arr = TRANSFER(i128, arr)
i64_1=arr(1)
i64_2=arr(2)
!----
!i64_1=0
!i64_2=0
!do ibit = 0,63
!    if(btest(i128,ibit)) i64_1=ibset(i64_1,ibit)
!enddo
!do ibit = 64,127
!    if(btest(i128,ibit)) i64_2=ibset(i64_2,ibit-64)
!enddo
end subroutine conv_064_to_64
!------------------------------------------------------

subroutine conv_64_to_064(i128,i64_1,i64_2)
integer*8, intent(out) :: i128
integer*8, intent(in) :: i64_1,i64_2
integer :: ibit

! integer*16 :: lowbits

!i128=0
!lowbits=0
!i128 = i64_2
!i128=ISHFT(i128,63)
!lowbits = i64_1
!i128=IOR(i128,lowbits)
!---
i128=0
do ibit = 0,63
    if(btest(i64_1,ibit)) i128=ibset(i128,ibit)
enddo
do ibit = 0,63
    if(btest(i64_2,ibit)) i128=ibset(i128,ibit+64)
enddo
end subroutine conv_64_to_064
!------------------------------------------------------

subroutine cluster_init

  character(len=16) filename

  integer            :: narg, iarg
  character(len=256) :: command_line_arguments
  character(len=256) :: executable_name
  character(len=256) :: argument
  character(len=256) :: input_file_name = ''
  integer iostat

#ifdef MPI
  call MPI_INIT(mpierr)

  if(mpierr /= 0) then
    write(6,*) "MPI init failed!! errorcode:",mpierr
    stop "MPI init failed!! errorcode:"
  endif

  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, ncores, mpierr )
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, whoami, mpierr )

  !Set up a communicator for each node
  CALL MPI_COMM_SPLIT_TYPE( MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,sharedComm,mpierr )
  !Figure out the rank for each processor, within the node
  CALL MPI_COMM_RANK( sharedComm,whoaminode,mpierr )
  CALL MPI_COMM_SIZE( sharedComm, ncores_on_node, mpierr )
  !MJO create a communicator for all of the master_core_node
  !    The second argument is the 'color'. All cores who pass the same integer
  !    will be put on the same communicator. This means all whoaminode==0 are
  !    are on the same comm, between_nodes_comm. The third argument selects ranks,
  !    but it is unimportant for our purposes.
  if (whoaminode.eq.0) then
    CALL MPI_COMM_SPLIT( MPI_COMM_WORLD,whoaminode,whoami,between_nodes_comm,mpierr )
    CALL MPI_COMM_SIZE ( between_nodes_comm,nnodes,mpierr)
  else
    CALL MPI_COMM_SPLIT( MPI_COMM_WORLD,MPI_UNDEFINED,whoami,between_nodes_comm,mpierr )
  endif
  call mpi_bsend(nnodes) ! Send the nnodes from master core to everyone
!**setup sending process
  allocate(snd_table(ncores*ncores),stat=istat)
  if(istat.ne.0) stop 'failed to allocate snd_table'
  allocate(snd_displs(ncores*ncores),stat=istat)
  if(istat.ne.0) stop 'failed to allocate snd_displs'
  allocate(snd_dsp(ncores),stat=istat)
  if(istat.ne.0) stop 'failed to allocate snd_dsp'
  allocate(rcv_cnt(ncores),stat=istat)
  if(istat.ne.0) stop 'failed to allocate rcv_cnt'

  allocate(recv_displs(ncores),stat=istat)
  if(istat.ne.0) stop 'failed to allocate recv_displs'

  call init_mpitype_walker !** to define MPI_WALK_TYPE
  call init_mpitype_det !** To define MPI_DET_TYPE
  call init_mpitype_diag_elem_info!** To define MPI_DIAG_ELEM_INFO_TYPE
  call init_mpitype_det_mix_diag !** To define MPI_DET_MIX_DIAG_TYPE
  call init_mpitype_det_mix_term_diag !** To define MPI_DET_MIX_DIAG_TERM_TYPE

  master_core=.false.
  if(whoami==0) master_core=.true.

  !Set master_core_node to true if you are the mastercore for the node
  master_core_node=.false.
  if(whoaminode==0) master_core_node=.true.
!** added redirection of stdout for all but master_core [09/12/13]
  if(.not.master_core) then
    close(6)
! Uncomment the following lines and comment the /dev/null line if you want output from slave processes
  ! if(whoami.le.9) then
  !   write(filename,'(''slave.'',i1)') whoami
  ! elseif(whoami.le.99) then
  !   write(filename,'(''slave.'',i2)') whoami
  ! elseif(whoami.le.999) then
  !   write(filename,'(''slave.'',i3)') whoami
  ! elseif(whoami.le.9999) then
  !   write(filename,'(''slave.'',i4)') whoami
  ! elseif(whoami.le.99999) then
  !   write(filename,'(''slave.'',i4)') whoami
  ! endif
  ! open(6,file=filename)
    open(6,file='/dev/null')
  endif

  win_num = 0 !MJO set counter for number or mpi windows (shared memory arrays)


  write(6,'(''MPI is ENABLED [nc='',i5,'']'')') ncores
#else
  whoami = 0
  master_core = .true.
  whoaminode = 0
  master_core_node = .true.
  ncores = 1
  write(6,*) "MPI is DISABLED!"
#endif
  !MJO - Input file setup, for both MPI and non MPI, is here in cluster init
  !      Maybe it should be elsewhere?
  if (master_core) then
     !MJO get input file name
     !    modified from champ/qmc/champ.f90

     ! number of arguments on command line
     narg = command_argument_count()

     ! get command line arguments
     call get_command (command_line_arguments)

     ! get executable name
     call get_command_argument (0, executable_name)

     write(6,'(2a)') 'Executable: ',trim(executable_name)
     write(6,'(2a)') 'Command line arguments: ',trim(command_line_arguments)

     iarg = 0
     do
        iarg = iarg + 1
        if (iarg > narg) exit
        call get_command_argument (iarg, argument)
        if (trim(argument) == '-i' .or. trim(argument) == '-input') then
           iarg = iarg + 1
           call get_command_argument (iarg, argument)
           if (iarg > narg) then
              write(6,*) 'value for option "-i{input}" missing.'
              stop
           endif
        endif
        input_file_name  = trim(argument)
     enddo

#ifdef MPI
     ! check input file
     if (trim(input_file_name) == '') then
        write(6,*)'input file name must be given in the command line by "-i{nput} filename" when using MPI'
        stop
     endif
#endif
     if (trim(input_file_name).ne.'') then
       !write(6,'(3a)') 'Opening input file >',trim(input_file_name),'<.'
        open(5, file=trim(input_file_name), status='old', iostat=iostat)
        if (iostat /= 0) then
           write(6,*) 'error on opening file >',trim(input_file_name),'<'
        endif
     endif
  endif

  allocate(snd_cnt(ncores),stat=istat)
  if(istat.ne.0) stop 'failed to allocate snd_cnt'
end subroutine cluster_init
!------------------------------------------------------

subroutine cluster_finalize
  integer :: i
#ifdef MPI
  !free the window
  do i=1,win_num
     call MPI_WIN_FREE(shmemWindow(i),mpierr)
  enddo
  call MPI_FINALIZE(mpierr)
  if(mpierr /= 0) then
    write(6,*) "MPI finalize failed!! errorcode:",mpierr
    stop "MPI finalize failed!! errorcode:"
  endif

!  deallocate(snd_table,snd_displs,recv_displs)
#endif
end subroutine cluster_finalize
!------------------------------------------------------
subroutine mpi_barr()
#ifdef MPI
  call MPI_Barrier(MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_barr

!------------------------------------------------------
subroutine mpi_barr_in_node()
#ifdef MPI
  call MPI_Barrier(sharedComm,mpierr)
#endif
end subroutine mpi_barr_in_node!--------------------------------------------------------------------------
!--- MPI_ALLTOALLV ---------------------------------------------------------

subroutine mpi_alltoallv_iwalk(sndbuf,sndcounts,sdispls,recvbuf,recvcounts,rdispls)
type(t_walk), intent(in) :: sndbuf(:)
integer, intent(in) :: sndcounts(:)
integer, intent(in) :: sdispls(:),rdispls(:)
type(t_walk), intent(inout) :: recvbuf(:)
integer, intent(inout) :: recvcounts(:)
#ifdef MPI
call MPI_ALLTOALLV(sndbuf,sndcounts,sdispls,MPI_WALK_TYPE,recvbuf,recvcounts,rdispls,MPI_WALK_TYPE,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_alltoallv_iwalk
!--------------------------------------------------------------------------

!--- MPI_AGATH ---------------------------------------------------------
subroutine mpi_agath_int(sndbuf,sndcount,recvbuf,recvcount)
integer, intent(in) :: sndbuf(:)
integer, intent(in) :: sndcount
integer, intent(inout) :: recvbuf(:)
integer, intent(inout) :: recvcount
#ifdef MPI
call MPI_ALLGATHER(sndbuf,sndcount,MPI_INTEGER,recvbuf,recvcount,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_agath_int
!--------------------------------------------------------------------------

!--- MPI_AGATHV ---------------------------------------------------------
subroutine mpi_agathv_int(sndbuf,sndcount,recvbuf,recvcounts,displs)
    integer, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount
    integer, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcounts(:)
    integer, intent(inout) :: displs(:)
#ifdef MPI
    call MPI_ALLGATHERV(sndbuf,sndcount,MPI_INTEGER,recvbuf,recvcounts,displs,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_agathv_int
!------------------------------------------------------

subroutine mpi_agathv_real_sp(sndbuf,sndcount,recvbuf,recvcounts,displs)
    real, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount
    real, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcounts(:)
    integer, intent(inout) :: displs(:)
#ifdef MPI
    call MPI_ALLGATHERV(sndbuf,sndcount,MPI_REAL,recvbuf,recvcounts,displs,MPI_REAL,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_agathv_real_sp
!------------------------------------------------------

subroutine mpi_agathv_real_dp(sndbuf,sndcount,recvbuf,recvcounts,displs)
    real*8, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount
    real*8, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcounts(:)
    integer, intent(inout) :: displs(:)
#ifdef MPI
    call MPI_ALLGATHERV(sndbuf,sndcount,MPI_REAL8,recvbuf,recvcounts,displs,MPI_REAL8,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_agathv_real_dp
!------------------------------------------------------

subroutine mpi_agathv_int8(sndbuf,sndcount,recvbuf,recvcounts,displs)
    integer*1, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount
    integer*1, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcounts(:)
    integer, intent(inout) :: displs(:)
#ifdef MPI
    call MPI_ALLGATHERV(sndbuf,sndcount,MPI_INTEGER1,recvbuf,recvcounts,displs,MPI_INTEGER1,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_agathv_int8
!------------------------------------------------------

subroutine mpi_agathv_ik_vec(sndbuf,sndcount,recvbuf,recvcounts,displs)
    type(ik_vec), intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount
    type(ik_vec), intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcounts(:)
    integer, intent(inout) :: displs(:)

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii,tot_snd
#ifdef MPI
    tot_snd=sum(recvcounts(1:ncores))
    allocate(high(tot_snd))
    allocate(low(tot_snd))
    allocate(buff(tot_snd))

    do ii=1,sndcount
        call conv_128_to_64(sndbuf(ii)%v(1),low(ii),high(ii))
    enddo
    call MPI_ALLGATHERV(low,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,MPI_COMM_WORLD,mpierr)
    low=buff !***now buff contains the full vector
    call MPI_ALLGATHERV(high,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,MPI_COMM_WORLD,mpierr)
    high=buff
    do ii=1,tot_snd
        call conv_64_to_128(recvbuf(ii)%v(1),low(ii),high(ii))
    enddo

    do ii=1,sndcount
        call conv_128_to_64(sndbuf(ii)%v(2),low(ii),high(ii))
    enddo
    call MPI_ALLGATHERV(low,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,MPI_COMM_WORLD,mpierr)
    low=buff !***now buff contains the full vector
    call MPI_ALLGATHERV(high,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,MPI_COMM_WORLD,mpierr)
    high=buff
    do ii=1,tot_snd
        call conv_64_to_128(recvbuf(ii)%v(2),low(ii),high(ii))
    enddo

    deallocate(high,low,buff)
#endif
end subroutine mpi_agathv_ik_vec
!------------------------------------------------------

subroutine mpi_agathv_i128(sndbuf,sndcount,recvbuf,recvcounts,displs)
    integer*16, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount
    integer*16, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcounts(:)
    integer, intent(inout) :: displs(:)

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii,tot_snd
#ifdef MPI
    tot_snd=sum(recvcounts(1:ncores))
    allocate(high(tot_snd))
    allocate(low(tot_snd))
    allocate(buff(tot_snd))

    do ii=1,sndcount
        call conv_128_to_64(sndbuf(ii),low(ii),high(ii))
    enddo
    call MPI_ALLGATHERV(low,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,MPI_COMM_WORLD,mpierr)
    low=buff !***now buff contains the full vector
    call MPI_ALLGATHERV(high,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,MPI_COMM_WORLD,mpierr)
    high=buff
    do ii=1,tot_snd
        call conv_64_to_128(recvbuf(ii),low(ii),high(ii))
    enddo

    deallocate(high,low,buff)
#endif
end subroutine mpi_agathv_i128
!--------------------------------------------------------------------------

!--- MPI_GATH ---------------------------------------------------------
subroutine mpi_gath_int(sndbuf,sndcount,recvbuf,recvcount,root)
    integer, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    integer, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcount
#ifdef MPI
    call MPI_GATHER(sndbuf,sndcount,MPI_INTEGER,recvbuf,recvcount,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
#else
    recvbuf(1:recvcount)=sndbuf(1:recvcount)
#endif
end subroutine mpi_gath_int
!------------------------------------------------------

subroutine mpi_gath_real_sp(sndbuf,sndcount,recvbuf,recvcount,root)
    real, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    real, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcount
#ifdef MPI
    call MPI_GATHER(sndbuf,sndcount,MPI_REAL,recvbuf,recvcount,MPI_REAL,root,MPI_COMM_WORLD,mpierr)
#else
    recvbuf(1:recvcount)=sndbuf(1:recvcount)
#endif
end subroutine mpi_gath_real_sp
!------------------------------------------------------

subroutine mpi_gath_real_dp(sndbuf,sndcount,recvbuf,recvcount,root)
    real*8, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    real*8, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcount
#ifdef MPI
    call MPI_GATHER(sndbuf,sndcount,MPI_REAL8,recvbuf,recvcount,MPI_REAL8,root,MPI_COMM_WORLD,mpierr)
#else
    recvbuf(1:recvcount)=sndbuf(1:recvcount)
#endif
end subroutine mpi_gath_real_dp
!------------------------------------------------------

subroutine mpi_gath_int8(sndbuf,sndcount,recvbuf,recvcount,root)
    integer*1, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    integer*1, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcount
#ifdef MPI
    call MPI_GATHER(sndbuf,sndcount,MPI_INTEGER1,recvbuf,recvcount,MPI_INTEGER1,root,MPI_COMM_WORLD,mpierr)
#else
    recvbuf(1:recvcount)=sndbuf(1:recvcount)
#endif
end subroutine mpi_gath_int8
!------------------------------------------------------

!subroutine mpi_gath_ik_vec(sndbuf,sndcount,recvbuf,recvcount,root)
!    type(ik_vec), intent(in) :: sndbuf(:)
!    integer, intent(in) :: sndcount,root
!    type(ik_vec), intent(out) :: recvbuf(:)
!    integer, intent(in) :: recvcount
!
!    integer*8,allocatable  :: high(:),low(:),buff(:)
!    integer :: ii,tot_snd
!#ifdef MPI
!    tot_snd=sum(recvcount(1:ncores))
!    allocate(high(tot_snd))
!    allocate(low(tot_snd))
!    allocate(buff(tot_snd))
!
!    do ii=1,sndcount
!        call conv_128_to_64(sndbuf(ii)%v(1),low(ii),high(ii))
!    enddo
!    call MPI_GATHER(low,sndcount,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
!    low=buff !***now buff contains the full vector
!    call MPI_GATHER(high,sndcount,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
!    high=buff
!    do ii=1,tot_snd
!        call conv_64_to_128(recvbuf(ii)%v(1),low(ii),high(ii))
!    enddo
!
!    do ii=1,sndcount
!        call conv_128_to_64(sndbuf(ii)%v(2),low(ii),high(ii))
!    enddo
!    call MPI_GATHER(low,sndcount,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
!    low=buff !***now buff contains the full vector
!    call MPI_GATHER(high,sndcount,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
!    high=buff
!    do ii=1,tot_snd
!        call conv_64_to_128(recvbuf(ii)%v(2),low(ii),high(ii))
!    enddo
!    deallocate(high,low,buff)
!#endif
!end subroutine mpi_gath_ik_vec
!!------------------------------------------------------
!
!subroutine mpi_gath_i128(sndbuf,sndcount,recvbuf,recvcount,root)
!    integer*16, intent(in) :: sndbuf(:)
!    integer, intent(in) :: sndcount,root
!    integer*16, intent(out) :: recvbuf(:)
!    integer, intent(in) :: recvcount
!
!    integer*8,allocatable  :: high(:),low(:),buff(:)
!    integer :: ii,tot_snd
!#ifdef MPI
!    tot_snd=sum(recvcount(1:ncores))
!    allocate(high(tot_snd))
!    allocate(low(tot_snd))
!    allocate(buff(tot_snd))
!
!    do ii=1,sndcount
!        call conv_128_to_64(sndbuf(ii),low(ii),high(ii))
!    enddo
!    call MPI_GATHER(low,sndcount,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
!    low=buff !***now buff contains the full vector
!    call MPI_GATHER(high,sndcount,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
!    high=buff
!    do ii=1,tot_snd
!        call conv_64_to_128(recvbuf(ii),low(ii),high(ii))
!    enddo
!
!    deallocate(high,low,buff)
!#endif
!end subroutine mpi_gath_i128
!!--------------------------------------------------------------------------

!--- MPI_GATHV ---------------------------------------------------------
subroutine mpi_gathv_int(sndbuf,sndcount,recvbuf,recvcounts,displs,root)
    integer, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    integer, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcounts(:)
    integer, intent(in) :: displs(:)
#ifdef MPI
    call MPI_GATHERV(sndbuf,sndcount,MPI_INTEGER,recvbuf,recvcounts,displs,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_gathv_int
!------------------------------------------------------

subroutine mpi_gathv_real_sp(sndbuf,sndcount,recvbuf,recvcounts,displs,root)
    real, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    real, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcounts(:)
    integer, intent(in) :: displs(:)
#ifdef MPI
    call MPI_GATHERV(sndbuf,sndcount,MPI_REAL,recvbuf,recvcounts,displs,MPI_REAL,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_gathv_real_sp
!------------------------------------------------------

subroutine mpi_gathv_real_dp(sndbuf,sndcount,recvbuf,recvcounts,displs,root)
    real*8, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    real*8, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcounts(:)
    integer, intent(in) :: displs(:)
#ifdef MPI
    call MPI_GATHERV(sndbuf,sndcount,MPI_REAL8,recvbuf,recvcounts,displs,MPI_REAL8,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_gathv_real_dp
!------------------------------------------------------

subroutine mpi_gathv_int8(sndbuf,sndcount,recvbuf,recvcounts,displs,root)
    integer*1, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    integer*1, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcounts(:)
    integer, intent(in) :: displs(:)
#ifdef MPI
    call MPI_GATHERV(sndbuf,sndcount,MPI_INTEGER1,recvbuf,recvcounts,displs,MPI_INTEGER1,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_gathv_int8
!------------------------------------------------------

subroutine mpi_gathv_ik_vec(sndbuf,sndcount,recvbuf,recvcounts,displs,root)
    type(ik_vec), intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    type(ik_vec), intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcounts(:)
    integer, intent(in) :: displs(:)

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii,tot_snd
#ifdef MPI
    tot_snd=sum(recvcounts(1:ncores))
    allocate(high(tot_snd))
    allocate(low(tot_snd))
    allocate(buff(tot_snd))

    do ii=1,sndcount
        call conv_128_to_64(sndbuf(ii)%v(1),low(ii),high(ii))
    enddo
    call MPI_GATHERV(low,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    low=buff !***now buff contains the full vector
    call MPI_GATHERV(high,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    high=buff
    do ii=1,tot_snd
        call conv_64_to_128(recvbuf(ii)%v(1),low(ii),high(ii))
    enddo

    do ii=1,sndcount
        call conv_128_to_64(sndbuf(ii)%v(2),low(ii),high(ii))
    enddo
    call MPI_GATHERV(low,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    low=buff !***now buff contains the full vector
    call MPI_GATHERV(high,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    high=buff
    do ii=1,tot_snd
        call conv_64_to_128(recvbuf(ii)%v(2),low(ii),high(ii))
    enddo
    deallocate(high,low,buff)
#endif
end subroutine mpi_gathv_ik_vec
!------------------------------------------------------

subroutine mpi_gathv_i128(sndbuf,sndcount,recvbuf,recvcounts,displs,root)
    integer*16, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcount,root
    integer*16, intent(out) :: recvbuf(:)
    integer, intent(in) :: recvcounts(:)
    integer, intent(in) :: displs(:)

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii,tot_snd
#ifdef MPI
    tot_snd=sum(recvcounts(1:ncores))
    allocate(high(tot_snd))
    allocate(low(tot_snd))
    allocate(buff(tot_snd))

    do ii=1,sndcount
        call conv_128_to_64(sndbuf(ii),low(ii),high(ii))
    enddo
    call MPI_GATHERV(low,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    low=buff !***now buff contains the full vector
    call MPI_GATHERV(high,sndcount,MPI_INTEGER8,buff,recvcounts,displs,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    high=buff
    do ii=1,tot_snd
        call conv_64_to_128(recvbuf(ii),low(ii),high(ii))
    enddo

    deallocate(high,low,buff)
#endif
end subroutine mpi_gathv_i128
!--------------------------------------------------------------------------

!--- MPI_SCATTV ---------------------------------------------------------
subroutine mpi_scattv_int(sndbuf,sndcounts,displs,recvbuf,recvcount,root)
    integer, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcounts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: root
    integer, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcount
#ifdef MPI
    call MPI_SCATTERV(sndbuf,sndcounts,displs,MPI_INTEGER,recvbuf,recvcount,MPI_INTEGER,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_scattv_int
!------------------------------------------------------

subroutine mpi_scattv_real_dp(sndbuf,sndcounts,displs,recvbuf,recvcount,root)
    real*8, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcounts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: root
    real*8, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcount
#ifdef MPI
    call MPI_SCATTERV(sndbuf,sndcounts,displs,MPI_REAL8,recvbuf,recvcount,MPI_REAL8,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_scattv_real_dp
!------------------------------------------------------

subroutine mpi_scattv_real_sp(sndbuf,sndcounts,displs,recvbuf,recvcount,root)
    real, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcounts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: root
    real, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcount
#ifdef MPI
    call MPI_SCATTERV(sndbuf,sndcounts,displs,MPI_REAL,recvbuf,recvcount,MPI_REAL,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_scattv_real_sp
!------------------------------------------------------

subroutine mpi_scattv_int8(sndbuf,sndcounts,displs,recvbuf,recvcount,root)
    integer*1, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcounts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: root
    integer*1, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcount
#ifdef MPI
    call MPI_SCATTERV(sndbuf,sndcounts,displs,MPI_INTEGER1,recvbuf,recvcount,MPI_INTEGER1,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_scattv_int8
!------------------------------------------------------

subroutine mpi_scattv_iwalk(sndbuf,sndcounts,displs,recvbuf,recvcount,root)
    type(t_walk), intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcounts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: root
    type(t_walk), intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcount
#ifdef MPI
    call MPI_SCATTERV(sndbuf,sndcounts,displs,MPI_WALK_TYPE,recvbuf,recvcount,MPI_WALK_TYPE,root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_scattv_iwalk
!------------------------------------------------------

subroutine mpi_scattv_ik_vec(sndbuf,sndcounts,displs,recvbuf,recvcount,root)
    type(ik_vec), intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcounts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: root
    type(ik_vec), intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcount

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii,tot_snd
#ifdef MPI
    tot_snd=sum(sndcounts(1:ncores))
    allocate(high(tot_snd))
    allocate(low(tot_snd))
    allocate(buff(tot_snd))

    do ii=1,tot_snd
        call conv_128_to_64(sndbuf(ii)%v(1),low(ii),high(ii))
    enddo
    call MPI_SCATTERV(low,sndcounts,displs,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    low=buff !***now buff contains the full vector
    call MPI_SCATTERV(high,sndcounts,displs,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    high=buff
    do ii=1,recvcount
        call conv_64_to_128(recvbuf(ii)%v(1),low(ii),high(ii))
    enddo

    do ii=1,tot_snd
        call conv_128_to_64(sndbuf(ii)%v(2),low(ii),high(ii))
    enddo
    call MPI_SCATTERV(low,sndcounts,displs,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    low=buff !***now buff contains the full vector
    call MPI_SCATTERV(high,sndcounts,displs,MPI_INTEGER8,buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    high=buff
    do ii=1,recvcount
        call conv_64_to_128(recvbuf(ii)%v(2),low(ii),high(ii))
    enddo
    deallocate(high,low,buff)
#endif
end subroutine mpi_scattv_ik_vec
!------------------------------------------------------

subroutine mpi_scattv_i128(sndbuf,sndcounts,displs,recvbuf,recvcount,root)
    integer*16, intent(in) :: sndbuf(:)
    integer, intent(in) :: sndcounts(:)
    integer, intent(in) :: displs(:)
    integer, intent(in) :: root
    integer*16, intent(inout) :: recvbuf(:)
    integer, intent(inout) :: recvcount

    integer :: ii,tot_snd
#ifdef MPI
    tot_snd=sum(sndcounts(1:ncores))

    if(.not.allocated(i128_low)) then
      allocate(i128_low(loc_mwalk),stat=istat)
      if(istat.ne.0) stop 'failed to allocate i128_low'
    endif
    if(.not.allocated(i128_high)) then
      allocate(i128_high(loc_mwalk),stat=istat)
      if(istat.ne.0) stop 'failed to allocate i128_high'
    endif
    if(.not.allocated(i128_buff)) then
      allocate(i128_buff(loc_mwalk),stat=istat)
      if(istat.ne.0) stop 'failed to allocate i128_buff'
    endif
    if(whoami==root) then
        do ii=1,tot_snd
            call conv_128_to_64(sndbuf(ii),i128_low(ii),i128_high(ii))
        enddo
    endif
    call MPI_SCATTERV(i128_low,sndcounts,displs,MPI_INTEGER8,i128_buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    i128_low(1:recvcount)=i128_buff(1:recvcount)
    call MPI_SCATTERV(i128_high,sndcounts,displs,MPI_INTEGER8,i128_buff,recvcount,MPI_INTEGER8,root,MPI_COMM_WORLD,mpierr)
    do ii=1,recvcount
        call conv_64_to_128(recvbuf(ii),i128_low(ii),i128_buff(ii))
    enddo
#endif
end subroutine mpi_scattv_i128
!--------------------------------------------------------------------------
!--- MPI_REDUCE for max --------------------------------------------------------
subroutine mpi_red_max_real_sp(rspvar_in,rspvar_out,root)
real, intent(in) :: rspvar_in
real, intent(out) :: rspvar_out
integer, intent(in) :: root
#ifdef MPI
  call MPI_REDUCE(rspvar_in,rspvar_out,1,MPI_REAL,MPI_MAX,root,MPI_COMM_WORLD,mpierr)
#else
  rspvar_out=rspvar_in
#endif
end subroutine mpi_red_max_real_sp
!------------------------------------------------------

subroutine mpi_red_max_real_dp(rdpvar_in,rdpvar_out,root)
real*8, intent(in) :: rdpvar_in
real*8, intent(out) :: rdpvar_out
integer, intent(in) :: root
#ifdef MPI
  call MPI_REDUCE(rdpvar_in,rdpvar_out,1,MPI_REAL8,MPI_MAX,root,MPI_COMM_WORLD,mpierr)
#else
  rdpvar_out=rdpvar_in
#endif
end subroutine mpi_red_max_real_dp
!------------------------------------------------------

subroutine mpi_red_max_real_dparray(rdparr_in,rdparr_out,root)
real*8, intent(in) :: rdparr_in(:)
real*8, intent(out) :: rdparr_out(:)
integer, intent(in) :: root
#ifdef MPI
  call MPI_REDUCE(rdparr_in,rdparr_out,size(rdparr_in),MPI_REAL8,MPI_MAX,root,MPI_COMM_WORLD,mpierr)
#else
  rdparr_out=rdparr_in
#endif
end subroutine mpi_red_max_real_dparray
!--------------------------------------------------------------------------

!--- MPI_REDUCE for min --------------------------------------------------------
subroutine mpi_red_min_real_sp(rspvar_in,rspvar_out,root)
real, intent(in) :: rspvar_in
real, intent(out) :: rspvar_out
integer, intent(in) :: root
#ifdef MPI
  call MPI_REDUCE(rspvar_in,rspvar_out,1,MPI_REAL,MPI_MIN,root,MPI_COMM_WORLD,mpierr)
#else
  rspvar_out=rspvar_in
#endif
end subroutine mpi_red_min_real_sp
!------------------------------------------------------

subroutine mpi_red_min_real_dp(rdpvar_in,rdpvar_out,root)
real*8, intent(in) :: rdpvar_in
real*8, intent(out) :: rdpvar_out
integer, intent(in) :: root
#ifdef MPI
  call MPI_REDUCE(rdpvar_in,rdpvar_out,1,MPI_REAL8,MPI_MIN,root,MPI_COMM_WORLD,mpierr)
#else
  rdpvar_out=rdpvar_in
#endif
end subroutine mpi_red_min_real_dp
!------------------------------------------------------

subroutine mpi_red_min_real_dparray(rdparr_in,rdparr_out,root)
real*8, intent(in) :: rdparr_in(:)
real*8, intent(out) :: rdparr_out(:)
integer, intent(in) :: root
#ifdef MPI
  call MPI_REDUCE(rdparr_in,rdparr_out,size(rdparr_in),MPI_REAL8,MPI_MIN,root,MPI_COMM_WORLD,mpierr)
#else
  rdparr_out=rdparr_in
#endif
end subroutine mpi_red_min_real_dparray
!--------------------------------------------------------------------------

!--- MPI_REDUCE for sum --------------------------------------------------------
subroutine mpi_red_real_sp(rspvar,root)
real, intent(inout) :: rspvar
integer, intent(in) :: root
real :: rspsum
#ifdef MPI
  call MPI_REDUCE(rspvar,rspsum,1,MPI_REAL,MPI_SUM,root,MPI_COMM_WORLD,mpierr)
  if(whoami==root) rspvar=rspsum
#endif
end subroutine mpi_red_real_sp
!------------------------------------------------------

subroutine mpi_red_real_dp(rdpvar,root)
real*8, intent(inout) :: rdpvar
integer, intent(in) :: root
real*8 :: rdpsum
#ifdef MPI
  call MPI_REDUCE(rdpvar,rdpsum,1,MPI_REAL8,MPI_SUM,root,MPI_COMM_WORLD,mpierr)
  if(whoami==root) rdpvar=rdpsum
#endif
end subroutine mpi_red_real_dp
!------------------------------------------------------

subroutine mpi_red_real_dparray(rdparr,root)
real*8, intent(inout) :: rdparr(:)
integer, intent(in) :: root
real*8 :: rdpsum(size(rdparr))
#ifdef MPI
  call MPI_REDUCE(rdparr,rdpsum,size(rdparr),MPI_REAL8,MPI_SUM,root,MPI_COMM_WORLD,mpierr)
  if(whoami==root) rdparr=rdpsum
#endif
end subroutine mpi_red_real_dparray
!--------------------------------------------------------------------------

!--- MPI_REDUCE_SCATTER for sum ---------------------------------------------------
subroutine mpi_redscatt_real_dparray(rdparr,counts)
real*8, intent(inout) :: rdparr(:)
integer, intent(in) :: counts(:)
real*8  :: rdpsum(counts(whoami+1))
#ifdef MPI
  call MPI_REDUCE_SCATTER(rdparr,rdpsum,counts, MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
  rdparr(1:counts(whoami+1))=rdpsum(1:counts(whoami+1))
#endif
end subroutine mpi_redscatt_real_dparray
!--------------------------------------------------------------------------

!--- MPI_ALLREDUCE --------------------------------------------------------
subroutine mpi_allred_int(ivar)
integer, intent(inout) :: ivar
integer :: isum
#ifdef MPI
  call MPI_ALLREDUCE(ivar,isum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
  ivar=isum
#endif
end subroutine mpi_allred_int
!------------------------------------------------------

subroutine mpi_allred_i8b(ivar)
  integer (i8b), intent(inout) :: ivar
  integer (i8b):: isum
#ifdef MPI
  call MPI_ALLREDUCE(ivar,isum,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,mpierr)
  ivar=isum
#endif
end subroutine mpi_allred_i8b
!------------------------------------------------------

subroutine mpi_allred_real_sp(rspvar)
real, intent(inout) :: rspvar
real :: rspsum
#ifdef MPI
  call MPI_ALLREDUCE(rspvar,rspsum,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,mpierr)
  rspvar=rspsum
#endif
end subroutine mpi_allred_real_sp
!------------------------------------------------------

subroutine mpi_allred_real_dp(rdpvar)
real*8, intent(inout) :: rdpvar
real*8 :: rdpsum
#ifdef MPI
  call MPI_ALLREDUCE(rdpvar,rdpsum,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
  rdpvar=rdpsum
#endif
end subroutine mpi_allred_real_dp
!------------------------------------------------------

subroutine mpi_allred_real_dparray(rdparr)
real*8, intent(inout) :: rdparr(:)
real*8 :: rdpsum(size(rdparr))
#ifdef MPI
  call MPI_ALLREDUCE(rdparr,rdpsum,size(rdparr),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
  rdparr=rdpsum
#endif
end subroutine mpi_allred_real_dparray
!------------------------------------------------------

subroutine mpi_allred_real_dparray2(rdparr)
  real*8, intent(inout) :: rdparr(:,:)
  real*8 :: rdpsum(size(rdparr,1),size(rdparr,2))
  integer :: total_size
  total_size = size(rdparr,1)*size(rdparr,2)
#ifdef MPI
  call MPI_ALLREDUCE(rdparr,rdpsum,total_size,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
  rdparr=rdpsum
#endif
end subroutine mpi_allred_real_dparray2
!------------------------------------------------------

subroutine mpi_allred_iarray(iarr)
integer, intent(inout) :: iarr(:)
integer :: isum(size(iarr))
#ifdef MPI
  call MPI_ALLREDUCE(iarr,isum,size(iarr),MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpierr)
  iarr=isum
#endif
end subroutine mpi_allred_iarray
!--------------------------------------------------------------------------

!--- MPI_BSEND_DIFFROOT ------------------------------------------------------
subroutine mpi_bsend_diffroot_int(ivar,id_root)
integer, intent(inout) :: ivar
integer, intent(in) :: id_root
#ifdef MPI
  call MPI_Bcast(ivar,1,MPI_INTEGER,id_root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_diffroot_int
!------------------------------------------------------

subroutine mpi_bsend_diffroot_iarray(iarr,id_root)
integer,intent(inout) :: iarr(:)
integer, intent(in) :: id_root
#ifdef MPI
call MPI_Bcast(iarr,size(iarr),MPI_INTEGER,id_root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_diffroot_iarray
!------------------------------------------------------

subroutine mpi_bsend_diffroot_real_sp(rspvar,id_root)
real, intent(inout) :: rspvar
integer, intent(in) :: id_root
#ifdef MPI
  call MPI_Bcast(rspvar,1,MPI_REAL,id_root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_diffroot_real_sp
!------------------------------------------------------

subroutine mpi_bsend_diffroot_real_dp(rdpvar,id_root)
real*8, intent(inout) :: rdpvar
integer, intent(in) :: id_root
#ifdef MPI
  call MPI_Bcast(rdpvar,1,MPI_REAL8,id_root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_diffroot_real_dp
!------------------------------------------------------

subroutine mpi_bsend_diffroot_logical(lvar,id_root)
logical, intent(inout) :: lvar
integer, intent(in) :: id_root
#ifdef MPI
  call MPI_Bcast(lvar,1,MPI_LOGICAL,id_root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_diffroot_logical
!------------------------------------------------------

subroutine mpi_bsend_diffroot_int64(ivar64, id_root)
integer*8, intent(inout) :: ivar64
integer, intent(in) :: id_root
#ifdef MPI
  call MPI_Bcast(ivar64,1,MPI_INTEGER8,id_root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_diffroot_int64
!------------------------------------------------------

subroutine mpi_bsend_diffroot_sndlist(isndlist, id_root)
type(t_sndlist), intent(inout) :: isndlist
integer, intent(in) :: id_root
#ifdef MPI
  call MPI_Bcast(isndlist%to,1,MPI_INTEGER,id_root,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(isndlist%wlkid,1,MPI_INTEGER,id_root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_diffroot_sndlist
!------------------------------------------------------

subroutine mpi_bsend_diffroot_int128(ivar128, id_root)
integer*16, intent(inout) :: ivar128
integer, intent(in) :: id_root
integer*8  :: high,low
#ifdef MPI

  call conv_128_to_64(ivar128,low,high)
  call MPI_Bcast(low,1,MPI_INTEGER8,id_root,MPI_COMM_WORLD,mpierr)
  call MPI_Bcast(high,1,MPI_INTEGER8,id_root,MPI_COMM_WORLD,mpierr)
  call conv_64_to_128(ivar128,low,high)

#endif
end subroutine mpi_bsend_diffroot_int128
!------------------------------------------------------

subroutine mpi_bsend_diffroot_int8(ivar8, id_root)
integer*1, intent(inout) :: ivar8
integer, intent(in) :: id_root
#ifdef MPI
  call MPI_Bcast(ivar8,1,MPI_INTEGER1,id_root,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_diffroot_int8
!------------------------------------------------------

subroutine mpi_bsend_diffroot_ik_vec(ik_vecvar,id_root)
type(ik_vec), intent(inout) :: ik_vecvar
integer, intent(in) :: id_root
integer*8  :: high,low
integer :: i
#ifdef MPI

  if (ik==selected_int_kind(38)) then
    do i=1,num_words
      call conv_128_to_64(ik_vecvar%v(i),low,high)
      call MPI_Bcast(low,1,MPI_INTEGER8,id_root,MPI_COMM_WORLD,mpierr)
      call MPI_Bcast(high,1,MPI_INTEGER8,id_root,MPI_COMM_WORLD,mpierr)
      call conv_64_to_128(ik_vecvar%v(i),low,high)
    enddo
  else
    do i=1,num_words
      call MPI_Bcast(ik_vecvar%v(i),1,MPI_INTEGER8,id_root,MPI_COMM_WORLD,mpierr)
    enddo
  endif

#endif
end subroutine mpi_bsend_diffroot_ik_vec
!--------------------------------------------------------------------------

!--- MPI_BSEND ------------------------------------------------------------
subroutine mpi_bsend_int(ivar)
integer, intent(inout) :: ivar
#ifdef MPI
  call MPI_Bcast(ivar,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_int
!------------------------------------------------------

subroutine mpi_bsend_real_sp(rspvar)
real, intent(inout) :: rspvar
#ifdef MPI
  call MPI_Bcast(rspvar,1,MPI_REAL,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_real_sp
!------------------------------------------------------

subroutine mpi_bsend_real_dp(rdpvar)
real*8, intent(inout) :: rdpvar
#ifdef MPI
  call MPI_Bcast(rdpvar,1,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_real_dp
!------------------------------------------------------

subroutine mpi_bsend_real_dparray(rdparr)
real*8, intent(inout) :: rdparr(:)
#ifdef MPI
  call MPI_Bcast(rdparr,size(rdparr),MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_real_dparray
!------------------------------------------------------

subroutine mpi_bsend_complex(cvar)
complex, intent(inout) :: cvar
#ifdef MPI
  call MPI_Bcast(cvar,1,MPI_COMPLEX,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_complex
!------------------------------------------------------

subroutine mpi_bsend_string(svar)
character, intent(inout) :: svar*16
#ifdef MPI
 call MPI_Bcast(svar,len(svar),MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_string
!------------------------------------------------------

subroutine mpi_bsend_iarray(iarr)
integer,intent(inout) :: iarr(:)
#ifdef MPI
call MPI_Bcast(iarr,size(iarr),MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_iarray
!------------------------------------------------------

!ad-hoc for rand seed, should be substitute with one for matrices
subroutine mpi_bsend_irand_seed(irand)
integer,intent(inout),dimension(4,2) :: irand
#ifdef MPI
call MPI_Bcast(irand(:,1),4,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
call MPI_Bcast(irand(:,2),4,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_irand_seed
!------------------------------------------------------

subroutine mpi_bsend_logical(lvar)
logical, intent(inout) :: lvar
#ifdef MPI
  call MPI_Bcast(lvar,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_logical
!------------------------------------------------------

subroutine mpi_bsend_int64(ivar64)
integer*8, intent(inout) :: ivar64
#ifdef MPI
  call MPI_Bcast(ivar64,1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_int64
!------------------------------------------------------

subroutine mpi_bsend_i64arr(i64arr)
integer*8, intent(inout) :: i64arr(:)
#ifdef MPI
    call MPI_Bcast(i64arr,size(i64arr),MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
#endif
end subroutine mpi_bsend_i64arr
!------------------------------------------------------

subroutine mpi_bsend_int128(i128)
integer*16, intent(inout) :: i128
#ifdef MPI
    integer*8  :: high,low
    if(master_core) then
      call conv_128_to_64(i128,low,high)
    endif
    call MPI_Bcast(low,1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(high,1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
    call conv_64_to_128(i128,low,high)
#endif
end subroutine mpi_bsend_int128
!------------------------------------------------------

subroutine mpi_bsend_i128arr(i128arr)
integer*16, intent(inout) :: i128arr(:)

    integer*8,allocatable  :: high(:),low(:)
    integer :: ii
#ifdef MPI
    allocate(high(size(i128arr)))
    allocate(low(size(i128arr)))
    if(master_core) then
        do ii=1,size(i128arr)
            call conv_128_to_64(i128arr(ii),low(ii),high(ii))
        enddo
    endif
    call MPI_Bcast(low,size(low),MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(high,size(high),MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
    do ii=1,size(i128arr)
        call conv_64_to_128(i128arr(ii),low(ii),high(ii))
    enddo
#endif
end subroutine mpi_bsend_i128arr
!------------------------------------------------------

subroutine mpi_bsend_ikvec(ikv)
    type(ik_vec), intent(inout) :: ikv
    integer*8,allocatable  :: high(:),low(:)
    integer :: ii
#ifdef MPI
  if (ik==selected_int_kind(38)) then
    allocate(high(num_words))
    allocate(low(num_words))
    if(master_core) then
        do ii=1,num_words
            call conv_128_to_64(ikv%v(ii),low(ii),high(ii))
        enddo
    endif
    call MPI_Bcast(low,size(low),MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(high,size(high),MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
    do ii=1,num_words
        call conv_64_to_128(ikv%v(ii),low(ii),high(ii))
    enddo
    deallocate(high,low)
  else
    do ii=1,num_words
      call MPI_Bcast(ikv%v(ii),1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
    enddo
  endif

#endif
end subroutine mpi_bsend_ikvec
!------------------------------------------------------

subroutine mpi_bsend_ikvec_arr(ikvarr)
  type(ik_vec), intent(inout) :: ikvarr(:)
  integer*8,allocatable  :: high(:),low(:)
  integer :: ii,jj,kk
#ifdef MPI
  if (ik==selected_int_kind(38)) then
    allocate(high(num_words*size(ikvarr)))
    allocate(low(num_words*size(ikvarr)))
    if(master_core) then
        kk=1
        do jj=1,size(ikvarr)
            do ii=1,num_words
                call conv_128_to_64(ikvarr(jj)%v(ii),low(kk),high(kk))
                kk=kk+1
            enddo
        enddo
    endif
    call MPI_Bcast(low,size(low),MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
    call MPI_Bcast(high,size(high),MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
    kk=1
    do jj=1,size(ikvarr)
        do ii=1,num_words
            call conv_64_to_128(ikvarr(jj)%v(ii),low(kk),high(kk))
            kk=kk+1
        enddo
    enddo
    deallocate(high,low)
  else
    do jj=1,size(ikvarr)
      do ii=1,num_words
        call MPI_Bcast(ikvarr(jj)%v(ii),1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpierr)
      enddo
    enddo
  endif

#endif
end subroutine mpi_bsend_ikvec_arr
!------------------------------------------------------

subroutine mpi_bsend_between_nodes_ikvec_arr(ikvarr)
  type(ik_vec), intent(inout) :: ikvarr(:)
  integer*8,allocatable  :: high(:),low(:)
  integer :: ii,jj,kk
#ifdef MPI
  if (ik==selected_int_kind(38)) then
    allocate(high(num_words*size(ikvarr)))
    allocate(low(num_words*size(ikvarr)))
    if(master_core) then
        kk=1
        do jj=1,size(ikvarr)
            do ii=1,num_words
                call conv_128_to_64(ikvarr(jj)%v(ii),low(kk),high(kk))
                kk=kk+1
            enddo
        enddo
    endif
    call MPI_Bcast(low,size(low),MPI_INTEGER8,0,between_nodes_comm,mpierr)
    call MPI_Bcast(high,size(high),MPI_INTEGER8,0,between_nodes_comm,mpierr)
    kk=1
    do jj=1,size(ikvarr)
        do ii=1,num_words
            call conv_64_to_128(ikvarr(jj)%v(ii),low(kk),high(kk))
            kk=kk+1
        enddo
    enddo
    deallocate(high,low)
  else
    do jj=1,size(ikvarr)
      do ii=1,num_words
        call MPI_Bcast(ikvarr(jj)%v(ii),1,MPI_INTEGER8,0,between_nodes_comm,mpierr)
      enddo
    enddo
  endif

#endif
end subroutine mpi_bsend_between_nodes_ikvec_arr
!------------------------------------------------------

subroutine mpi_bsend_between_nodes_i128arr(i128arr)
  integer*16, intent(inout) :: i128arr(:)

  integer*8,allocatable  :: high(:),low(:),buff(:)
  integer :: ii
#ifdef MPI
  allocate(high(size(i128arr)))
  allocate(low(size(i128arr)))
  allocate(buff(size(i128arr)))
  if(master_core) then
    do ii=1,size(i128arr)
      call conv_128_to_64(i128arr(ii),low(ii),high(ii))
    enddo
  endif
  call MPI_Bcast(low,size(low),MPI_INTEGER8,0,between_nodes_comm,mpierr)
  call MPI_Bcast(high,size(high),MPI_INTEGER8,0,between_nodes_comm,mpierr)
  do ii=1,size(i128arr)
    call conv_64_to_128(i128arr(ii),low(ii),high(ii))
  enddo
  deallocate(high,low,buff)
#endif
end subroutine mpi_bsend_between_nodes_i128arr
!------------------------------------------------------

subroutine mpi_bsend_between_nodes_dp_array(rdparr)
  real*8, intent(inout) :: rdparr(:)
#ifdef MPI
  call MPI_Bcast(rdparr,size(rdparr),MPI_REAL8,0,between_nodes_comm,mpierr)
#endif
end subroutine mpi_bsend_between_nodes_dp_array
!------------------------------------------------------

subroutine mpi_bsend_between_nodes_tdet_array(tdetarr)
  type(t_det), intent(inout) :: tdetarr(:)
#ifdef MPI
  call MPI_Bcast(tdetarr,size(tdetarr),MPI_DET_TYPE,0,between_nodes_comm,mpierr)
#endif
end subroutine mpi_bsend_between_nodes_tdet_array
!------------------------------------------------------

!-------------------------------------------------------------------------

!--- MPI_SND ------------------------------------------------------------
!-- assuming master_core is the sender!!
!-------------------------------------------------------------------------
subroutine mpi_snd_int(ivar,target_core)
integer, intent(inout) :: ivar
integer, intent(in)    :: target_core
#ifdef MPI
integer :: cstatus(MPI_STATUS_SIZE)
  if(master_core) call mpi_send(ivar,1,MPI_INTEGER,target_core,69,MPI_COMM_WORLD,mpierr)
  if(whoami==target_core) call mpi_recv(ivar,1,MPI_INTEGER,0,69,MPI_COMM_WORLD,cstatus,mpierr)
#endif
end subroutine mpi_snd_int
!------------------------------------------------------

subroutine mpi_snd_real_sp(rspvar,target_core)
real, intent(inout) :: rspvar
integer, intent(in)    :: target_core
#ifdef MPI
integer :: cstatus(MPI_STATUS_SIZE)
  if(master_core) call mpi_send(rspvar,1,MPI_REAL,target_core,69,MPI_COMM_WORLD,mpierr)
  if(whoami==target_core) call mpi_recv(rspvar,1,MPI_REAL,0,69,MPI_COMM_WORLD,cstatus,mpierr)
#endif
end subroutine mpi_snd_real_sp
!------------------------------------------------------

subroutine mpi_snd_real_dp(rdpvar,target_core)
real*8, intent(inout) :: rdpvar
integer, intent(in)    :: target_core
#ifdef MPI
integer :: cstatus(MPI_STATUS_SIZE)
  if(master_core) call mpi_send(rdpvar,1,MPI_REAL8,target_core,69,MPI_COMM_WORLD,mpierr)
  if(whoami==target_core) call mpi_recv(rdpvar,1,MPI_REAL8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
#endif
end subroutine mpi_snd_real_dp
!------------------------------------------------------

subroutine mpi_snd_real_dparray(rdparr,target_core)
real*8, intent(inout) :: rdparr(:)
integer, intent(in)    :: target_core
#ifdef MPI
integer :: cstatus(MPI_STATUS_SIZE)
integer :: total_size,istart,iend,i
  total_size=size(rdparr)
  if (total_size<=mpi_send_limit) then
    if(master_core) call mpi_send(rdparr,total_size,MPI_REAL8,target_core,69,MPI_COMM_WORLD,mpierr)
    if(whoami==target_core) call mpi_recv(rdparr,total_size,MPI_REAL8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
  else
    do i=1,(total_size-1)/mpi_send_limit+1
      istart=(i-1)*mpi_send_limit+1
      iend=min(total_size,i*mpi_send_limit)
     !write (6,*) "Sending connections",istart,"through",iend; call flush(6)
      if(master_core) call mpi_send(rdparr(istart:iend),iend-istart+1,MPI_REAL8,target_core,69,MPI_COMM_WORLD,mpierr)
      if(whoami==target_core) call mpi_recv(rdparr(istart:iend),iend-istart+1,MPI_REAL8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
    enddo
  endif
#endif
end subroutine mpi_snd_real_dparray
!------------------------------------------------------

subroutine mpi_snd_iarray(iarr,target_core)
integer,intent(inout) :: iarr(:)
integer, intent(in)    :: target_core
#ifdef MPI
integer :: cstatus(MPI_STATUS_SIZE)
integer :: total_size,istart,iend,i
  total_size=size(iarr)
  if (total_size<=mpi_send_limit) then
    if(master_core) call mpi_send(iarr,total_size,MPI_REAL8,target_core,69,MPI_COMM_WORLD,mpierr)
    if(whoami==target_core) call mpi_recv(iarr,total_size,MPI_REAL8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
  else
    do i=1,(total_size-1)/mpi_send_limit+1
      istart=(i-1)*mpi_send_limit+1
      iend=min(total_size,i*mpi_send_limit)
     !write (6,*) "Sending connections",istart,"through",iend; call flush(6)
      if(master_core) call mpi_send(iarr(istart:iend),iend-istart+1,MPI_REAL8,target_core,69,MPI_COMM_WORLD,mpierr)
      if(whoami==target_core) call mpi_recv(iarr(istart:iend),iend-istart+1,MPI_REAL8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
    enddo
  endif
#endif
end subroutine mpi_snd_iarray
!------------------------------------------------------

subroutine mpi_snd_i64array(i64arr,target_core)
integer*8,intent(inout) :: i64arr(:)
integer, intent(in)    :: target_core
#ifdef MPI
integer :: cstatus(MPI_STATUS_SIZE)
integer :: total_size,istart,iend,i
  total_size=size(i64arr)
  if (total_size<=mpi_send_limit) then
    if(master_core) call mpi_send(i64arr,total_size,MPI_REAL8,target_core,69,MPI_COMM_WORLD,mpierr)
    if(whoami==target_core) call mpi_recv(i64arr,total_size,MPI_REAL8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
  else
    do i=1,(total_size-1)/mpi_send_limit+1
      istart=(i-1)*mpi_send_limit+1
      iend=min(total_size,i*mpi_send_limit)
     !write (6,*) "Sending connections",istart,"through",iend; call flush(6)
      if(master_core) call mpi_send(i64arr(istart:iend),iend-istart+1,MPI_REAL8,target_core,69,MPI_COMM_WORLD,mpierr)
      if(whoami==target_core) call mpi_recv(i64arr(istart:iend),iend-istart+1,MPI_REAL8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
    enddo
  endif
#endif
end subroutine mpi_snd_i64array
!------------------------------------------------------

subroutine mpi_snd_int64(int64,target_core)
integer*8,intent(inout) :: int64
integer, intent(in)    :: target_core
#ifdef MPI
integer :: cstatus(MPI_STATUS_SIZE)
  if(master_core) call mpi_send(int64,1,MPI_INTEGER8,target_core,69,MPI_COMM_WORLD,mpierr)
  if(whoami==target_core) call mpi_recv(int64,1,MPI_INTEGER8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
#endif
end subroutine mpi_snd_int64

subroutine mpi_snd_i128arr(i128arr,target_core)
    integer*16, intent(inout) :: i128arr(:)
    integer, intent(in) :: target_core

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii
#ifdef MPI
integer :: cstatus(MPI_STATUS_SIZE)
integer :: total_size,istart,iend,i
  total_size=size(i128arr)
  if(master_core .or. (whoami==target_core)) then
    if (total_size<=mpi_send_limit) then
      allocate(high(total_size))
      allocate(low(total_size))
      allocate(buff(total_size))
      if(master_core) then
          do ii=1,total_size
              call conv_128_to_64(i128arr(ii),low(ii),high(ii))
          enddo
          call mpi_send(low,total_size,MPI_INTEGER8,target_core,69,MPI_COMM_WORLD,mpierr)
          call mpi_send(high,total_size,MPI_INTEGER8,target_core,79,MPI_COMM_WORLD,mpierr)
      endif
      if(whoami == target_core) then
          call mpi_recv(low,total_size,MPI_INTEGER8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
          call mpi_recv(high,total_size,MPI_INTEGER8,0,79,MPI_COMM_WORLD,cstatus,mpierr)
          do ii=1,total_size
              call conv_64_to_128(i128arr(ii),low(ii),high(ii))
          enddo
      endif
    else
      allocate(high(mpi_send_limit))
      allocate(low(mpi_send_limit))
      allocate(buff(mpi_send_limit))
      do i=1,(total_size-1)/mpi_send_limit+1
        istart=(i-1)*mpi_send_limit+1
        iend=min(total_size,i*mpi_send_limit)
       !write (6,*) "Sending connections",istart,"through",iend; call flush(6)
        if(master_core) then
            do ii=istart,iend
                call conv_128_to_64(i128arr(ii),low(ii-istart+1),high(ii-istart+1))
            enddo
            call mpi_send(low(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,target_core,69,MPI_COMM_WORLD,mpierr)
            call mpi_send(high(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,target_core,79,MPI_COMM_WORLD,mpierr)
        endif
        if(whoami == target_core) then
            call mpi_recv(low(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
            call mpi_recv(high(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,0,79,MPI_COMM_WORLD,cstatus,mpierr)
            do ii=istart,iend
                call conv_64_to_128(i128arr(ii),low(ii-istart+1),high(ii-istart+1))
            enddo
        endif
      enddo
    endif
    deallocate(high,low,buff)
  endif
#endif
end subroutine mpi_snd_i128arr
!------------------------------------------------------

subroutine mpi_snd_ikvec(ikv,target_core)
    type(ik_vec), intent(inout) :: ikv
    integer, intent(in) :: target_core

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii
#ifdef MPI
  integer :: cstatus(MPI_STATUS_SIZE)
  integer :: total_size,istart,iend,i
  total_size=num_words
  if(master_core .or. (whoami==target_core)) then
    if (total_size<=mpi_send_limit) then
      allocate(high(total_size))
      allocate(low(total_size))
      allocate(buff(total_size))
      if(master_core) then
          do ii=1,total_size
              call conv_128_to_64(ikv%v(ii),low(ii),high(ii))
          enddo
          call mpi_send(low,total_size,MPI_INTEGER8,target_core,69,MPI_COMM_WORLD,mpierr)
          call mpi_send(high,total_size,MPI_INTEGER8,target_core,79,MPI_COMM_WORLD,mpierr)
      endif
      if(whoami == target_core) then
          call mpi_recv(low,total_size,MPI_INTEGER8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
          call mpi_recv(high,total_size,MPI_INTEGER8,0,79,MPI_COMM_WORLD,cstatus,mpierr)
          do ii=1,total_size
              call conv_64_to_128(ikv%v(ii),low(ii),high(ii))
          enddo
      endif
    else
      allocate(high(mpi_send_limit))
      allocate(low(mpi_send_limit))
      allocate(buff(mpi_send_limit))
      do i=1,(total_size-1)/mpi_send_limit+1
        istart=(i-1)*mpi_send_limit+1
        iend=min(total_size,i*mpi_send_limit)
       !write (6,*) "Sending connections",istart,"through",iend; call flush(6)
        if(master_core) then
            do ii=istart,iend
                call conv_128_to_64(ikv%v(ii),low(ii-istart+1),high(ii-istart+1))
            enddo
            call mpi_send(low(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,target_core,69,MPI_COMM_WORLD,mpierr)
            call mpi_send(high(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,target_core,79,MPI_COMM_WORLD,mpierr)
        endif
        if(whoami == target_core) then
            call mpi_recv(low(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
            call mpi_recv(high(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,0,79,MPI_COMM_WORLD,cstatus,mpierr)
            do ii=istart,iend
                call conv_64_to_128(ikv%v(ii),low(ii-istart+1),high(ii-istart+1))
            enddo
        endif
      enddo
    endif
    deallocate(high,low,buff)
  endif
#endif
end subroutine mpi_snd_ikvec
!------------------------------------------------------

subroutine mpi_snd_ikvec_array(ikvarr,target_core)
    type(ik_vec), intent(inout) :: ikvarr(:)
    integer, intent(in) :: target_core

    integer*8,allocatable  :: high(:),low(:),buff(:)
    integer :: ii,jj,kk
#ifdef MPI
  integer :: cstatus(MPI_STATUS_SIZE)
  integer :: total_size,istart,iend,i
  total_size=num_words*size(ikvarr)
  if(master_core .or. (whoami==target_core)) then
    if (total_size<=mpi_send_limit) then
      allocate(high(total_size))
      allocate(low(total_size))
      allocate(buff(total_size))
      if(master_core) then
          kk=1
          do jj=1,size(ikvarr)
              do ii=1,num_words
                  call conv_128_to_64(ikvarr(jj)%v(ii),low(kk),high(kk))
                  kk=kk+1
              enddo
          enddo
          call mpi_send(low,total_size,MPI_INTEGER8,target_core,69,MPI_COMM_WORLD,mpierr)
          call mpi_send(high,total_size,MPI_INTEGER8,target_core,79,MPI_COMM_WORLD,mpierr)
      endif
      if(whoami == target_core) then
          call mpi_recv(low,total_size,MPI_INTEGER8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
          call mpi_recv(high,total_size,MPI_INTEGER8,0,79,MPI_COMM_WORLD,cstatus,mpierr)
          kk=1
          do jj=1,size(ikvarr)
              do ii=1,num_words
                  call conv_64_to_128(ikvarr(jj)%v(ii),low(kk),high(kk))
                  kk=kk+1
              enddo
          enddo
      endif
    else
      allocate(high(mpi_send_limit))
      allocate(low(mpi_send_limit))
      allocate(buff(mpi_send_limit))
      do i=1,(total_size-1)/mpi_send_limit+1
        istart=(i-1)*mpi_send_limit+1
        iend=min(total_size,i*mpi_send_limit)
       !write (6,*) "Sending connections",istart,"through",iend; call flush(6)
        if(master_core) then
            do ii=istart,iend
                call conv_128_to_64(ikvarr(int((ii-1)/num_words)+1)%v(mod(ii-1,num_words)+1),low(ii-istart+1),high(ii-istart+1))
            enddo
            call mpi_send(low(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,target_core,69,MPI_COMM_WORLD,mpierr)
            call mpi_send(high(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,target_core,79,MPI_COMM_WORLD,mpierr)
        endif
        if(whoami == target_core) then
            call mpi_recv(low(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,0,69,MPI_COMM_WORLD,cstatus,mpierr)
            call mpi_recv(high(1:iend-istart+1),iend-istart+1,MPI_INTEGER8,0,79,MPI_COMM_WORLD,cstatus,mpierr)
            do ii=istart,iend
                call conv_64_to_128(ikvarr(int((ii-1)/num_words)+1)%v(mod(ii-1,num_words)+1),low(ii-istart+1),high(ii-istart+1))
            enddo
        endif
      enddo
    endif
    deallocate(high,low,buff)
  endif
#endif
end subroutine mpi_snd_ikvec_array
!------------------------------------------------------------------------------

!**WARNING: may break consistency if we change the kind of some of the walker arrays, should work with a generic object walker [TODO]
!my_cnt: number of walks to be sent
!my_nwalk: current population on core, used as offset for addresses [buf_snd=buf_snd(my_nwalk+1:my_nwalk+my_cnt)]
!nwalk: total population, used for current get_owner [CHANGE WHEN USING HASH_TABLE]
!my_nwalk_1 : first element of my new population [shall be set to zero when everyone works with its own walks, just use my_nwalk instead]
! Only used once at beginning of run

subroutine mpi_sendwalks(isnd_start,my_nwalk,my_cnt,my_nwalk_1,initiator,e_num,e_den,diag_elems)
  integer,intent(in) :: isnd_start,my_nwalk_1,my_cnt
  integer(i1b), intent(inout) :: initiator(:)
  real(rk),intent(inout) :: e_num(:),e_den(:),diag_elems(:)

  integer,intent(inout) :: my_nwalk

  integer :: isnd,iown,tot_snd,icore,my_ic,my_own,iw
#ifdef NUM_ORBITALS_GT_127
  integer :: islice
#endif

#ifdef MPI
  snd_table=0
  do isnd=1,my_cnt
    iown=get_det_owner(walk_dets_up(isnd_start+isnd),walk_dets_dn(isnd_start+isnd))
    snd_list(isnd)%to=iown
    snd_list(isnd)%wlkid=isnd_start+isnd
    snd_table(whoami*ncores + iown + 1) = snd_table(whoami*ncores + iown + 1) + 1
  enddo
    call mpi_allred(snd_table)

  recv_displs=0
  do icore=1,(ncores-1)
    recv_displs(icore+1)=snd_table((icore-1)*ncores+whoami+1)
  enddo
  do icore=2,ncores
    recv_displs(icore)=recv_displs(icore)+recv_displs(icore-1)
  enddo
  recv_displs(1)=0

  recv_displs=recv_displs+my_nwalk_1-1

  tot_snd=sum(snd_table(whoami*ncores + 1:(whoami+1)*ncores))

  snd_displs=0
  do icore=2,ncores
    snd_displs(whoami*ncores + icore)=sum(snd_table(whoami*ncores + 1:whoami*ncores + icore-1))
  enddo

  if(tot_snd > 0) then
    do icore=0,(ncores-1)
        iown=0
        my_ic=whoami*ncores + icore + 1
        do isnd=1,(tot_snd+my_own)
          if(snd_list(isnd)%to == icore) then
            iown=iown+1
            list_cowner(iown)=snd_list(isnd)%wlkid
          endif
          if(iown.ge.snd_table(my_ic)) exit
        enddo
        if(iown /= snd_table(my_ic)) stop "something's WRONG!!!"

        snd_buff(snd_displs(my_ic)+1:snd_displs(my_ic)+snd_table(my_ic))%wt=walk_wt(list_cowner(1:snd_table(my_ic)))
        snd_buff(snd_displs(my_ic)+1:snd_displs(my_ic)+snd_table(my_ic))%imp_ini(1) = imp_distance(list_cowner(1:snd_table(my_ic)))
        snd_buff(snd_displs(my_ic)+1:snd_displs(my_ic)+snd_table(my_ic))%imp_ini(2) = initiator(list_cowner(1:snd_table(my_ic)))
        snd_buff(snd_displs(my_ic)+1:snd_displs(my_ic)+snd_table(my_ic))%e_num = e_num(list_cowner(1:snd_table(my_ic)))
        snd_buff(snd_displs(my_ic)+1:snd_displs(my_ic)+snd_table(my_ic))%e_den = e_den(list_cowner(1:snd_table(my_ic)))
        snd_buff(snd_displs(my_ic)+1:snd_displs(my_ic)+snd_table(my_ic))%diag_elems = diag_elems(list_cowner(1:snd_table(my_ic)))
        do iw=1,snd_table(my_ic)
#ifdef NUM_ORBITALS_GT_127
            do islice=1,num_words
                call conv_128_to_64(walk_dets_up(list_cowner(iw))%v(islice),snd_buff(snd_displs(my_ic)+iw)%det_u(2*(islice-1)+1),snd_buff(snd_displs(my_ic)+iw)%det_u(2*(islice-1)+2))
                call conv_128_to_64(walk_dets_dn(list_cowner(iw))%v(islice),snd_buff(snd_displs(my_ic)+iw)%det_d(2*(islice-1)+1),snd_buff(snd_displs(my_ic)+iw)%det_d(2*(islice-1)+2))
            enddo
#else
            call conv_128_to_64(walk_dets_up(list_cowner(iw)),snd_buff(snd_displs(my_ic)+iw)%det_u(1),snd_buff(snd_displs(my_ic)+iw)%det_u(2))
            call conv_128_to_64(walk_dets_dn(list_cowner(iw)),snd_buff(snd_displs(my_ic)+iw)%det_d(1),snd_buff(snd_displs(my_ic)+iw)%det_d(2))
#endif
        enddo
        snd_matrix_elements(snd_displs(my_ic)+1:snd_displs(my_ic)+snd_table(my_ic)) = matrix_elements(list_cowner(1:snd_table(my_ic)))
    enddo
  endif

  do icore=0,(ncores-1)
    if(sum(snd_table(icore*ncores + 1:(icore+1)*ncores)) > 0) then !**do the following only if icore has something for the others
      call mpi_scattv(snd_buff(1:tot_snd),snd_table(icore*ncores + 1:(icore+1)*ncores),snd_displs(icore*ncores + 1:(icore+1)*ncores),recv_buff(recv_displs(icore+1)+1:recv_displs(icore+1)+snd_table(icore*ncores+ whoami +1)),snd_table(icore*ncores+ whoami +1),icore)
      call mpi_scattv(snd_matrix_elements(1:tot_snd),snd_table(icore*ncores + 1:(icore+1)*ncores),snd_displs(icore*ncores + 1:(icore+1)*ncores),matrix_elements(recv_displs(icore+1)+1:recv_displs(icore+1)+snd_table(icore*ncores+ whoami +1)),snd_table(icore*ncores+ whoami +1),icore)
    endif
  enddo

  do icore=0,(ncores-1)
    my_nwalk = my_nwalk +snd_table(icore*ncores+whoami+1)
  enddo

  walk_wt(recv_displs(1)+1:my_nwalk)=recv_buff(recv_displs(1)+1:my_nwalk)%wt
  do iw=recv_displs(1)+1,my_nwalk
#ifdef NUM_ORBITALS_GT_127
    do islice=1,num_words
      call conv_64_to_128(walk_dets_up(iw)%v(islice),recv_buff(iw)%det_u(2*(islice-1)+1),recv_buff(iw)%det_u(2*(islice-1)+2))
      call conv_64_to_128(walk_dets_dn(iw)%v(islice),recv_buff(iw)%det_d(2*(islice-1)+1),recv_buff(iw)%det_d(2*(islice-1)+2))
    enddo
#else
    call conv_64_to_128(walk_dets_up(iw),recv_buff(iw)%det_u(1),recv_buff(iw)%det_u(2))
    call conv_64_to_128(walk_dets_dn(iw),recv_buff(iw)%det_d(1),recv_buff(iw)%det_d(2))
#endif
  enddo
  imp_distance(recv_displs(1)+1:my_nwalk)=recv_buff(recv_displs(1)+1:my_nwalk)%imp_ini(1)
  initiator(recv_displs(1)+1:my_nwalk)=recv_buff(recv_displs(1)+1:my_nwalk)%imp_ini(2)
  e_num(recv_displs(1)+1:my_nwalk)=recv_buff(recv_displs(1)+1:my_nwalk)%e_num
  e_den(recv_displs(1)+1:my_nwalk)=recv_buff(recv_displs(1)+1:my_nwalk)%e_den
  diag_elems(recv_displs(1)+1:my_nwalk)=recv_buff(recv_displs(1)+1:my_nwalk)%diag_elems

  snd_table=0
#endif
end subroutine mpi_sendwalks
!------------------------------------------------------------------------------

subroutine mpi_push_nwalk(iw,iinitiator,e_num,e_den,diag_elem)
integer,intent(in) :: iw
integer(i1b), intent(in) :: iinitiator
real(rk),intent(in) :: e_num,e_den,diag_elem

integer :: iown,sdispls
#ifdef NUM_ORBITALS_GT_127
integer :: islice
#endif

#ifdef MPI
iown=get_det_owner(walk_dets_up(iw),walk_dets_dn(iw))
snd_cnt(iown+1) = snd_cnt(iown+1) + 1

if(snd_cnt(iown+1) > nwbuff) stop "Need to enlarge dimension of SND_BUFFs!"

sdispls=nwbuff*iown+snd_cnt(iown+1)

snd_buff(sdispls)%wt=walk_wt(iw)
snd_buff(sdispls)%imp_ini(1)=imp_distance(iw)
snd_buff(sdispls)%imp_ini(2)=iinitiator
snd_buff(sdispls)%e_num=e_num
snd_buff(sdispls)%e_den=e_den
snd_buff(sdispls)%diag_elems=diag_elem

#ifdef NUM_ORBITALS_GT_127
do islice=1,num_words
    call conv_128_to_64(walk_dets_up(iw)%v(islice),snd_buff(sdispls)%det_u(2*(islice-1)+1),snd_buff(sdispls)%det_u(2*(islice-1)+2))
    call conv_128_to_64(walk_dets_dn(iw)%v(islice),snd_buff(sdispls)%det_d(2*(islice-1)+1),snd_buff(sdispls)%det_d(2*(islice-1)+2))
enddo
#else
call conv_128_to_64(walk_dets_up(iw),snd_buff(sdispls)%det_u(1),snd_buff(sdispls)%det_u(2))
call conv_128_to_64(walk_dets_dn(iw),snd_buff(sdispls)%det_d(1),snd_buff(sdispls)%det_d(2))
#endif

#endif
end subroutine mpi_push_nwalk
!------------------------------------------------------------------------------

!**improved version to use when we know that the walkers are newly generated:
!** - don't sends matrix_elements but set them to 1e51_rk
!** - makes just one scatterv with a derived type
subroutine mpi_sendnewwalks(my_nwalk,initiator,e_num,e_den,diag_elems)
    integer(i1b), intent(inout) :: initiator(:)
    integer,intent(inout) :: my_nwalk
    real(rk),intent(inout) :: e_num(:),e_den(:),diag_elems(:)

    integer :: tot_snd,icore,my_own,iw
#ifdef NUM_ORBITALS_GT_127
    integer :: islice
#endif
#ifdef MPI
    call mpi_agath_int(snd_cnt,ncores,snd_table,ncores)

    recv_displs=0
    do icore=1,(ncores-1)
        recv_displs(icore+1)=snd_table((icore-1)*ncores+whoami+1)
    enddo
    do icore=2,ncores
        recv_displs(icore)=recv_displs(icore)+recv_displs(icore-1)
    enddo
    recv_displs(1)=0

    recv_displs=recv_displs+my_nwalk

    my_own=snd_table(whoami*ncores+whoami+1)
    snd_table(whoami*ncores+whoami+1) = 0 !**don't send back my own walkers

    tot_snd=sum(snd_table(whoami*ncores + 1:(whoami+1)*ncores))

    do icore=0,(ncores-1)
      rcv_cnt(1+icore)=snd_table(icore*ncores+ whoami +1)
    enddo

    call mpi_alltoallv_iwalk(snd_buff(1:tot_snd),snd_table(whoami*ncores + 1:(whoami+1)*ncores),snd_dsp,recv_buff,rcv_cnt,recv_displs)

    my_nwalk = my_nwalk + sum(rcv_cnt(1:ncores))

    my_nwalk = my_nwalk+my_own

    if(.not.master_core) then
      walk_wt(recv_displs(1)+1:recv_displs(whoami+1))=recv_buff(recv_displs(1)+1:recv_displs(whoami+1))%wt
      do iw=recv_displs(1)+1,recv_displs(whoami+1)
#ifdef NUM_ORBITALS_GT_127
        do islice=1,num_words
          call conv_64_to_128(walk_dets_up(iw)%v(islice),recv_buff(iw)%det_u(2*(islice-1)+1),recv_buff(iw)%det_u(2*(islice-1)+2))
          call conv_64_to_128(walk_dets_dn(iw)%v(islice),recv_buff(iw)%det_d(2*(islice-1)+1),recv_buff(iw)%det_d(2*(islice-1)+2))
        enddo
#else
        call conv_64_to_128(walk_dets_up(iw),recv_buff(iw)%det_u(1),recv_buff(iw)%det_u(2))
        call conv_64_to_128(walk_dets_dn(iw),recv_buff(iw)%det_d(1),recv_buff(iw)%det_d(2))
#endif
      enddo
      imp_distance(recv_displs(1)+1:recv_displs(whoami+1))=recv_buff(recv_displs(1)+1:recv_displs(whoami+1))%imp_ini(1)
      initiator(recv_displs(1)+1:recv_displs(whoami+1))=recv_buff(recv_displs(1)+1:recv_displs(whoami+1))%imp_ini(2)
      e_num(recv_displs(1)+1:recv_displs(whoami+1))=recv_buff(recv_displs(1)+1:recv_displs(whoami+1))%e_num
      e_den(recv_displs(1)+1:recv_displs(whoami+1))=recv_buff(recv_displs(1)+1:recv_displs(whoami+1))%e_den
      diag_elems(recv_displs(1)+1:recv_displs(whoami+1))=recv_buff(recv_displs(1)+1:recv_displs(whoami+1))%diag_elems
    endif

    walk_wt(recv_displs(whoami+1)+1:recv_displs(whoami+1)+my_own)=snd_buff(snd_dsp(whoami+1)+1:snd_dsp(whoami+1)+my_own)%wt
    do iw=1,my_own
#ifdef NUM_ORBITALS_GT_127
        do islice=1,num_words
          call conv_64_to_128(walk_dets_up(recv_displs(whoami+1)+iw)%v(islice),snd_buff(snd_dsp(whoami+1)+iw)%det_u(2*(islice-1)+1),snd_buff(snd_dsp(whoami+1)+iw)%det_u(2*(islice-1)+2))
          call conv_64_to_128(walk_dets_dn(recv_displs(whoami+1)+iw)%v(islice),snd_buff(snd_dsp(whoami+1)+iw)%det_d(2*(islice-1)+1),snd_buff(snd_dsp(whoami+1)+iw)%det_d(2*(islice-1)+2))
        enddo
#else
        call conv_64_to_128(walk_dets_up(recv_displs(whoami+1)+iw),snd_buff(snd_dsp(whoami+1)+iw)%det_u(1),snd_buff(snd_dsp(whoami+1)+iw)%det_u(2))
        call conv_64_to_128(walk_dets_dn(recv_displs(whoami+1)+iw),snd_buff(snd_dsp(whoami+1)+iw)%det_d(1),snd_buff(snd_dsp(whoami+1)+iw)%det_d(2))
#endif
    enddo
    imp_distance(recv_displs(whoami+1)+1:recv_displs(whoami+1)+my_own)=snd_buff(snd_dsp(whoami+1)+1:snd_dsp(whoami+1)+my_own)%imp_ini(1)
    initiator(recv_displs(whoami+1)+1:recv_displs(whoami+1)+my_own)=snd_buff(snd_dsp(whoami+1)+1:snd_dsp(whoami+1)+my_own)%imp_ini(2)
    e_num(recv_displs(whoami+1)+1:recv_displs(whoami+1)+my_own)=snd_buff(snd_dsp(whoami+1)+1:snd_dsp(whoami+1)+my_own)%e_num
    e_den(recv_displs(whoami+1)+1:recv_displs(whoami+1)+my_own)=snd_buff(snd_dsp(whoami+1)+1:snd_dsp(whoami+1)+my_own)%e_den
    diag_elems(recv_displs(whoami+1)+1:recv_displs(whoami+1)+my_own)=snd_buff(snd_dsp(whoami+1)+1:snd_dsp(whoami+1)+my_own)%diag_elems

    if(whoami+2 .le. ncores) then
      walk_wt(recv_displs(whoami+2)+1:my_nwalk)=recv_buff(recv_displs(whoami+2)+1:my_nwalk)%wt
      do iw=recv_displs(whoami+2)+1,my_nwalk
#ifdef NUM_ORBITALS_GT_127
        do islice=1,num_words
          call conv_64_to_128(walk_dets_up(iw)%v(islice),recv_buff(iw)%det_u(2*(islice-1)+1),recv_buff(iw)%det_u(2*(islice-1)+2))
          call conv_64_to_128(walk_dets_dn(iw)%v(islice),recv_buff(iw)%det_d(2*(islice-1)+1),recv_buff(iw)%det_d(2*(islice-1)+2))
        enddo
#else
        call conv_64_to_128(walk_dets_up(iw),recv_buff(iw)%det_u(1),recv_buff(iw)%det_u(2))
        call conv_64_to_128(walk_dets_dn(iw),recv_buff(iw)%det_d(1),recv_buff(iw)%det_d(2))
#endif
      enddo
      imp_distance(recv_displs(whoami+2)+1:my_nwalk)=recv_buff(recv_displs(whoami+2)+1:my_nwalk)%imp_ini(1)
      initiator(recv_displs(whoami+2)+1:my_nwalk)=recv_buff(recv_displs(whoami+2)+1:my_nwalk)%imp_ini(2)
      e_num(recv_displs(whoami+2)+1:my_nwalk)=recv_buff(recv_displs(whoami+2)+1:my_nwalk)%e_num
      e_den(recv_displs(whoami+2)+1:my_nwalk)=recv_buff(recv_displs(whoami+2)+1:my_nwalk)%e_den
      diag_elems(recv_displs(whoami+2)+1:my_nwalk)=recv_buff(recv_displs(whoami+2)+1:my_nwalk)%diag_elems
    endif

    matrix_elements(recv_displs(1)+1:my_nwalk) = 1e51_rk
#else
    my_nwalk = my_nwalk + sum(snd_cnt)
#endif
end subroutine mpi_sendnewwalks

!------------------------------------------------------------------------------------------------------------
subroutine mpi_merge_sort2_nowts(old_dets_up,old_dets_dn,n_det,dets_up,dets_dn,n_allocate,iorder,temp_i16_up,temp_i16_dn,temp_i_2)
! ===========================================================================================================
!
!
! -----------------------------------------------------------------------------------------------------------
  use tools, only : merge_sort2_up_dn,merge_original_with_spawned3 ! Note: this merge_original_with_spawned3 is very different from the routine in do_walk.f90!
  implicit none
  integer, allocatable, intent(inout) :: iorder(:)
  integer, intent(inout) :: n_det
  integer(i8b), intent(inout) :: n_allocate
  integer, allocatable, intent(inout) :: temp_i_2(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), allocatable, intent(inout) :: old_dets_up(:), old_dets_dn(:)
  type(ik_vec), allocatable, intent(inout) :: dets_up(:), dets_dn(:)
  type(ik_vec), allocatable, intent(inout) :: temp_i16_up(:), temp_i16_dn(:)
  type(ik_vec),allocatable :: send_buf_core_up(:,:),send_buf_core_dn(:,:)
#else
  integer(ik), allocatable, intent(inout) :: old_dets_up(:), old_dets_dn(:)
  integer(ik), allocatable, intent(inout) :: dets_up(:), dets_dn(:)
  integer(ik), allocatable, intent(inout) :: temp_i16_up(:), temp_i16_dn(:)
  integer(ik),allocatable :: send_buf_core_up(:,:),send_buf_core_dn(:,:)
#endif
  !Local variables
  integer(i8b) :: n_allocate_new
  integer :: core_id,i,j,k,l,send_size,recv_size,mpierr,n_det_remote(0:ncores-1)
  integer :: send_buf_cnt(0:ncores-1),recv_buf_cnt(0:ncores-1),send_buf_disp(0:ncores-1),recv_buf_disp(0:ncores-1)

  type(t_det),allocatable :: send_buf(:),recv_buf(:)

  !MJO remote merge/sort
  !MJO - pack buffers, send, recv, do another merge / sort
#ifdef MPI
  !Loop through list and figure out array sizes
  send_buf_cnt = 0
  do i=1,n_det
    core_id = get_det_owner(old_dets_up(i),old_dets_dn(i))
    send_buf_cnt(core_id) = send_buf_cnt(core_id) + 1
  enddo

  !Allocate arrays - We don't have to reallocate this way, but we have to loop over n_det twice
  allocate(send_buf_core_up(maxval(send_buf_cnt),0:ncores-1))
  allocate(send_buf_core_dn(maxval(send_buf_cnt),0:ncores-1))

  send_buf_cnt = 0
  do i=1,n_det
    core_id = get_det_owner(old_dets_up(i),old_dets_dn(i))
    !Fill up send buffers - we'll arrange them later
    send_buf_cnt(core_id) = send_buf_cnt(core_id) + 1
    send_buf_core_up(send_buf_cnt(core_id),core_id) = old_dets_up(i)
    send_buf_core_dn(send_buf_cnt(core_id),core_id) = old_dets_dn(i)
  enddo
  call flush(6)
  !Send over send_buf_cnts
  call MPI_alltoall(send_buf_cnt,1,MPI_INTEGER,recv_buf_cnt,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)

  !calculate displacements
  send_buf_disp(0) = 0
  recv_buf_disp(0) = 0
  do i=1,ncores-1
    send_buf_disp(i) = send_buf_cnt(i-1)+send_buf_disp(i-1)
    recv_buf_disp(i) = recv_buf_cnt(i-1)+recv_buf_disp(i-1)
  enddo

  send_size = sum(send_buf_cnt)
  recv_size = sum(recv_buf_cnt)

  !allocate send and receive bufs
  allocate(send_buf(0:send_size))
  allocate(recv_buf(0:recv_size))


  k=0
  do j=0,ncores-1
    do i=1,send_buf_cnt(j)
#ifdef NUM_ORBITALS_GT_127
      do l=1,num_words
        call conv_128_to_64(send_buf_core_up(i,j)%v(l),send_buf(k)%det_up(2*l-1),send_buf(k)%det_up(2*l))
        call conv_128_to_64(send_buf_core_dn(i,j)%v(l),send_buf(k)%det_dn(2*l-1),send_buf(k)%det_dn(2*l))
      enddo
#else
      call conv_128_to_64(send_buf_core_up(i,j),send_buf(k)%det_up(1),send_buf(k)%det_up(2))
      call conv_128_to_64(send_buf_core_dn(i,j),send_buf(k)%det_dn(1),send_buf(k)%det_dn(2))
#endif
      k = k+1
    enddo
  enddo

  !FIXME MJO - need to support norb.gt.127
  call MPI_alltoallv(send_buf,send_buf_cnt,send_buf_disp,MPI_DET_TYPE,&
       recv_buf,recv_buf_cnt,recv_buf_disp,MPI_DET_TYPE,&
       MPI_COMM_WORLD,mpierr)

  !Reallocate arrays if needed
  if (recv_size > n_allocate) then
    n_allocate_new = n_det+recv_size
    write(6,'(''1mpi: reallocating dets_up etc. arrays of size'',i10,'' ='',es11.4)') n_allocate_new, real(n_allocate_new)
    call flush(6)
    if(allocated(old_dets_up)) deallocate(old_dets_up,old_dets_dn)
    if(allocated(dets_up)) deallocate(dets_up,dets_dn)

    if(allocated(iorder)) deallocate(iorder)
    if(allocated(temp_i16_up)) deallocate(temp_i16_up)
    if(allocated(temp_i16_dn)) deallocate(temp_i16_dn)
    if(allocated(temp_i_2)) deallocate(temp_i_2)

    allocate(iorder(n_allocate_new))
    allocate(temp_i16_up((n_allocate_new+1)/2))
    allocate(temp_i16_dn((n_allocate_new+1)/2))
    allocate(temp_i_2((n_allocate_new+1)/2))

    allocate(old_dets_up(n_allocate_new))
    allocate(old_dets_dn(n_allocate_new))
    allocate(dets_up(n_allocate_new))
    allocate(dets_dn(n_allocate_new))
    write(6,'(''1mpi: after reallocating dets_up etc. arrays of size'',i10,'' ='',es11.4)') n_allocate_new, real(n_allocate_new)
    call flush(6)
  endif

  ! !Unpack buffers - includes local part
  do i=1,recv_size
#ifdef NUM_ORBITALS_GT_127
    do l=1,num_words
      call conv_64_to_128(dets_up(i)%v(l),recv_buf(i-1)%det_up(2*l-1),recv_buf(i-1)%det_up(2*l))
      call conv_64_to_128(dets_dn(i)%v(l),recv_buf(i-1)%det_dn(2*l-1),recv_buf(i-1)%det_dn(2*l))
    enddo
#else
    call conv_64_to_128(dets_up(i),recv_buf(i-1)%det_up(1),recv_buf(i-1)%det_up(2))
    call conv_64_to_128(dets_dn(i),recv_buf(i-1)%det_dn(1),recv_buf(i-1)%det_dn(2))
#endif
  enddo


  !merge sort the determinants from other cores, only if we received new determinants

  call merge_sort2_up_dn(old_dets_up(1:recv_size),old_dets_dn(1:recv_size), iorder, recv_size, temp_i16_up, temp_i16_dn, temp_i_2)

  n_det = recv_size !Set our new number of determinants

  call merge_original_with_spawned3(old_dets_up(1:recv_size),old_dets_dn(1:recv_size),n_det,dets_up,dets_dn)

  deallocate(send_buf)
  deallocate(recv_buf)
  deallocate(send_buf_core_dn)
  deallocate(send_buf_core_up)

  call MPI_ALLGATHER(n_det,1,MPI_INTEGER,n_det_remote,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
  ndets_global = sum(n_det_remote)
#endif

end subroutine mpi_merge_sort2_nowts

subroutine mpi_merge_sort2_nowts_replace(n_det,dets_up,dets_dn,n_allocate,iorder,temp_i16_up,temp_i16_dn,temp_i_2)
! ===========================================================================================================
!
!
! -----------------------------------------------------------------------------------------------------------
  use tools, only : merge_sort2_up_dn,merge_original_with_spawned3 ! Note: this merge_original_with_spawned3 is very different from the routine in do_walk.f90!
  implicit none
  integer, allocatable, intent(inout) :: iorder(:)
  integer, intent(inout) :: n_det
  integer(i8b), intent(inout) :: n_allocate
  integer, allocatable, intent(inout) :: temp_i_2(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), allocatable, intent(inout) :: dets_up(:), dets_dn(:)
  type(ik_vec), allocatable, intent(inout) :: temp_i16_up(:), temp_i16_dn(:)
#else
  integer(ik), allocatable, intent(inout) :: dets_up(:), dets_dn(:)
  integer(ik), allocatable, intent(inout) :: temp_i16_up(:), temp_i16_dn(:)
#endif

  !Local variables
  integer(i8b) :: n_allocate_new
  integer :: core_id,i,j,k,l,send_size,recv_size,mpierr,n_det_remote(0:ncores-1),index
  integer :: send_buf_cnt(0:ncores-1),recv_buf_cnt(0:ncores-1),send_buf_disp(0:ncores-1),recv_buf_disp(0:ncores-1)
  integer, allocatable :: send_buf_core_ind(:,:)
  type(t_det),allocatable :: send_buf(:),recv_buf(:)
  !MJO remote merge/sort
  !MJO - pack buffers, send, recv, do another merge / sort
#ifdef MPI
  !Loop through list and figure out array sizes
  send_buf_cnt = 0
  do i=1,n_det
    core_id = get_det_owner(dets_up(i),dets_dn(i))
    send_buf_cnt(core_id) = send_buf_cnt(core_id) + 1
  enddo
  !Allocate arrays - We don't have to reallocate this way, but we have to loop over n_det twice
  allocate(send_buf_core_ind(maxval(send_buf_cnt),0:ncores-1))

  send_buf_cnt = 0
  do i=1,n_det
    core_id = get_det_owner(dets_up(i),dets_dn(i))
    !Fill up send buffers - we'll arrange them later
    send_buf_cnt(core_id) = send_buf_cnt(core_id) + 1
    send_buf_core_ind(send_buf_cnt(core_id),core_id) = i
  enddo
  call flush(6)
  !Send over send_buf_cnts
  call MPI_alltoall(send_buf_cnt,1,MPI_INTEGER,recv_buf_cnt,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)

  !calculate displacements
  send_buf_disp(0) = 0
  recv_buf_disp(0) = 0
  do i=1,ncores-1
    send_buf_disp(i) = send_buf_cnt(i-1)+send_buf_disp(i-1)
    recv_buf_disp(i) = recv_buf_cnt(i-1)+recv_buf_disp(i-1)
  enddo

  send_size = sum(send_buf_cnt)
  recv_size = sum(recv_buf_cnt)

  !allocate send and receive bufs
  allocate(send_buf(0:send_size))
  allocate(recv_buf(0:recv_size))

  k=0
  do j=0,ncores-1
    do i=1,send_buf_cnt(j)
      index = send_buf_core_ind(i,j)
#ifdef NUM_ORBITALS_GT_127
      do l=1,num_words
        call conv_128_to_64(dets_up(index)%v(l),send_buf(k)%det_up(2*l-1),send_buf(k)%det_up(2*l))
        call conv_128_to_64(dets_dn(index)%v(l),send_buf(k)%det_dn(2*l-1),send_buf(k)%det_dn(2*l))
      enddo
#else
      call conv_128_to_64(dets_up(index),send_buf(k)%det_up(1),send_buf(k)%det_up(2))
      call conv_128_to_64(dets_dn(index),send_buf(k)%det_dn(1),send_buf(k)%det_dn(2))
#endif
      k = k+1
    enddo
  enddo

  call MPI_alltoallv(send_buf,send_buf_cnt,send_buf_disp,MPI_DET_TYPE,&
       recv_buf,recv_buf_cnt,recv_buf_disp,MPI_DET_TYPE,&
       MPI_COMM_WORLD,mpierr)

  !Reallocate arrays if needed
  if (recv_size > n_allocate) then
    n_allocate_new = n_det+recv_size
    write(6,'(''2mpi: reallocating dets_up etc. arrays of size'',i10,'' ='',es11.4)') n_allocate_new, real(n_allocate_new); call flush(6)
    if(allocated(dets_up)) deallocate(dets_up,dets_dn)

    if(allocated(iorder)) deallocate(iorder)
    if(allocated(temp_i16_up)) deallocate(temp_i16_up)
    if(allocated(temp_i16_dn)) deallocate(temp_i16_dn)
    if(allocated(temp_i_2)) deallocate(temp_i_2)

    allocate(iorder(n_allocate_new))
    allocate(temp_i16_up((n_allocate_new+1)/2))
    allocate(temp_i16_dn((n_allocate_new+1)/2))
    allocate(temp_i_2((n_allocate_new+1)/2))

    allocate(dets_up(n_allocate_new))
    allocate(dets_dn(n_allocate_new))
  endif

  ! !Unpack buffers - includes local part
  do i=1,recv_size
#ifdef NUM_ORBITALS_GT_127
    do l=1,num_words
      call conv_64_to_128(dets_up(i)%v(l),recv_buf(i-1)%det_up(2*l-1),recv_buf(i-1)%det_up(2*l))
      call conv_64_to_128(dets_dn(i)%v(l),recv_buf(i-1)%det_dn(2*l-1),recv_buf(i-1)%det_dn(2*l))
    enddo
#else
    call conv_64_to_128(dets_up(i),recv_buf(i-1)%det_up(1),recv_buf(i-1)%det_up(2))
    call conv_64_to_128(dets_dn(i),recv_buf(i-1)%det_dn(1),recv_buf(i-1)%det_dn(2))
#endif
  enddo


  !merge sort the determinants from other cores, only if we received new determinants
  call merge_sort2_up_dn(dets_up(1:recv_size),dets_dn(1:recv_size), iorder, recv_size, temp_i16_up, temp_i16_dn, temp_i_2)

  n_det = recv_size !Set our new number of determinants

  call merge_original_with_spawned3(n_det,dets_up(1:recv_size),dets_dn(1:recv_size))

  deallocate(send_buf)
  deallocate(recv_buf)

  call MPI_ALLGATHER(n_det,1,MPI_INTEGER,n_det_remote,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
  ndets_global = sum(n_det_remote)
#endif

end subroutine mpi_merge_sort2_nowts_replace


subroutine mpi_allgatherv_new_dets(dets_up,dets_dn,n_det,n_allocate,remote_det_map,local_det_map)
! =================================================================================================
!
!
! --------------------------------------------------------------------------------------------------
  implicit none
  integer, allocatable :: iorder(:)
  integer, intent(inout) :: n_det
  integer(i8b), intent(inout) :: n_allocate
  integer, allocatable :: temp_i_2(:)
  type(det_map), intent(inout) :: remote_det_map
  type(det_map_l), intent(inout) :: local_det_map
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), allocatable, intent(inout) :: dets_up(:), dets_dn(:)
  type(ik_vec), allocatable :: tmp_dets_up(:),tmp_dets_dn(:)
#else
  integer(ik), allocatable, intent(inout) :: dets_up(:), dets_dn(:)
  integer(ik), allocatable :: tmp_dets_up(:),tmp_dets_dn(:)
#endif
  !Local variables
  integer(i8b) :: n_allocate_new,new_det_index
  integer(i8b),allocatable :: local_indices_tmp(:)
  integer :: core_id,i,j,k,l,send_size,recv_size,mpierr,n_it,i_it,i_core
  integer :: offset,ndets_global_old,ndets_old,ndet_tmp
  integer :: send_buf_cnt(0:ncores-1),recv_buf_cnt(0:ncores-1),recv_buf_disp(0:ncores-1)
  type(t_det),allocatable :: send_buf(:),recv_buf(:)
  integer   :: n_det_remote(0:ncores-1),ncores_before_offset
  !MJO - pack buffers, send, recv
#ifdef MPI
  !Send over send_buf_cnts
  call MPI_allgather(n_det,1,MPI_INTEGER,recv_buf_cnt,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
  n_det_remote = recv_buf_cnt
  !calculate displacements
  recv_buf_disp(0) = 0
  do i=1,ncores-1
    recv_buf_disp(i) = recv_buf_cnt(i-1)+recv_buf_disp(i-1)
  enddo
  ! do i=0,ncores-1
  !   write(6,*) 'recv',recv_buf_cnt(i)
  ! enddo
  send_size = n_det
  recv_size = sum(recv_buf_cnt)
  write(6,'(''Newly spawned determinants (global)'',i10,'' ='',es11.4)') recv_size, real(recv_size); call flush(6)
  !allocate send and receive bufs
  allocate(send_buf(send_size))

  ! Only master_core_node will recv all of the determinants
  if (master_core_node) then
    allocate(recv_buf(recv_size))
  endif

  do i=1,n_det
#ifdef NUM_ORBITALS_GT_127
    do j=1,num_words
      call conv_128_to_64(dets_up(i)%v(j),send_buf(i)%det_up(2*j-1),send_buf(i)%det_up(2*j))
      call conv_128_to_64(dets_dn(i)%v(j),send_buf(i)%det_dn(2*j-1),send_buf(i)%det_dn(2*j))
    enddo
#else
    call conv_128_to_64(dets_up(i),send_buf(i)%det_up(1),send_buf(i)%det_up(2))
    call conv_128_to_64(dets_dn(i),send_buf(i)%det_dn(1),send_buf(i)%det_dn(2))
#endif
  enddo

  !Gather all dets to the master_core
  call MPI_gatherv(send_buf,n_det,MPI_DET_TYPE,&
    recv_buf,recv_buf_cnt,recv_buf_disp,MPI_DET_TYPE,&
    0,MPI_COMM_WORLD,mpierr)

  !Broadcast recv_buf to all master_core_node s
  !and unpack the buffers
  if(master_core_node) then
    call mpi_bsend_between_nodes(recv_buf)


    !Reallocate arrays if needed
    if (recv_size > n_allocate) then
      n_allocate_new = n_det+recv_size
      write(6,'(''3mpi: reallocating dets_up etc. arrays of size'',i10,'' ='',es11.4)') n_allocate_new, real(n_allocate_new); call flush(6)
      if(allocated(dets_up)) deallocate(dets_up,dets_dn)
      if(allocated(dets_up)) deallocate(dets_up,dets_dn)
      allocate(dets_up(n_allocate_new))
      allocate(dets_dn(n_allocate_new))
    endif

    ! !Unpack buffers - includes local part
    do i=1,recv_size
#ifdef NUM_ORBITALS_GT_127
      do j=1,num_words
        call conv_64_to_128(dets_up(i)%v(j),recv_buf(i)%det_up(2*j-1),recv_buf(i)%det_up(2*j))
        call conv_64_to_128(dets_dn(i)%v(j),recv_buf(i)%det_dn(2*j-1),recv_buf(i)%det_dn(2*j))
      enddo
#else
      call conv_64_to_128(dets_up(i),recv_buf(i)%det_up(1),recv_buf(i)%det_up(2))
      call conv_64_to_128(dets_dn(i),recv_buf(i)%det_dn(1),recv_buf(i)%det_dn(2))
#endif
    enddo
    deallocate(recv_buf)
  endif
  call mpi_barr_in_node()

  deallocate(send_buf)


  ndets_global_old = local_det_map%ndets_global

  local_det_map%ndets_global_old = ndets_global_old


  ndets_old        = local_det_map%ndets

  ! increment ndets_global, ndets
  ndet_tmp = sum(n_det_remote)
  local_det_map%ndets_global = local_det_map%ndets_global + ndet_tmp
  local_det_map%ndets = local_det_map%ndets + n_det_remote(whoami)

  !MJO - Resize, append the new determinants to the old global list
  call shmem_reallocate(remote_det_map%global_dets_up,int(local_det_map%ndets_global,i8b))
  call shmem_reallocate(remote_det_map%global_dets_dn,int(local_det_map%ndets_global,i8b))

  ! allocate(tmp_dets_up(ndets_global_old))
  ! allocate(tmp_dets_dn(ndets_global_old))
  ! tmp_dets_up(1:ndets_global_old) = remote_det_map%global_dets_up(1:ndets_global_old)
  ! tmp_dets_dn(1:ndets_global_old) = remote_det_map%global_dets_dn(1:ndets_global_old)

  ! if(allocated(remote_det_map%global_dets_up)) deallocate(remote_det_map%global_dets_up)
  ! if(allocated(remote_det_map%global_dets_dn)) deallocate(remote_det_map%global_dets_dn)
  ! allocate(remote_det_map%global_dets_up(local_det_map%ndets_global))
  ! allocate(remote_det_map%global_dets_dn(local_det_map%ndets_global))
  ! remote_det_map%global_dets_up(1_i8b:ndets_global_old) = tmp_dets_up(1:ndets_global_old)
  ! remote_det_map%global_dets_dn(1_i8b:ndets_global_old) = tmp_dets_dn(1:ndets_global_old)

  ! deallocate(tmp_dets_up)
  ! deallocate(tmp_dets_dn)

  new_det_index = ndets_global_old
  local_det_map%n_it = local_det_map%n_it + 1
  !Only master_core_node writes to shared array global_dets
  if (master_core_node) then
    do j=1,maxval(n_det_remote)
      do i_core = 0,ncores-1
        if (j<=n_det_remote(i_core)) then
          !MJO Append a determinant
          !    To get our interleaved structure, we loop through dets_up,
          !    adding one det per core at a time. To get the right det,
          !    we have to use offsets in dets_up and dets_dn
          new_det_index = new_det_index + 1
          remote_det_map%global_dets_up(new_det_index) = dets_up(recv_buf_disp(i_core)+j)
          remote_det_map%global_dets_dn(new_det_index) = dets_dn(recv_buf_disp(i_core)+j)
        endif
      enddo
    enddo
  endif
  call mpi_barr_in_node()
  !MJO Append n_det_remote_new to the full n_det_remote
  n_it = local_det_map%n_it
  do i_core=0,ncores-1
     local_det_map%n_det_remote(n_it,i_core) = n_det_remote(i_core)
  enddo

  !MJO append local_indices list
  !    now resize, append the indices to the end of the local list

  allocate(local_indices_tmp(ndets_old))
  local_indices_tmp(1:ndets_old) = local_det_map%local_indices(1:ndets_old)
  if(allocated(local_det_map%local_indices)) deallocate(local_det_map%local_indices)
  allocate(local_det_map%local_indices(local_det_map%ndets))

  local_det_map%local_indices(1_i8b:ndets_old) = local_indices_tmp(1:ndets_old)
  deallocate(local_indices_tmp)

  ncores_before_offset = 0
  do i=1,n_det_remote(whoami)
    !MJO This global index choice gives perfect interleaving of the determinants
    !    This gives much better load balancing than straight tiling
    do i_core = 0,whoami-1
      !Add the interleaved determinants from previous cores
      if (n_det_remote(i_core)>=i) then
        ncores_before_offset = ncores_before_offset + 1
      endif
    enddo
    !Add this cores determinant
    local_det_map%local_indices(ndets_old+i) = ndets_global_old + ncores_before_offset + i
    do i_core = whoami+1,ncores-1
      !Add the interleaved dets from ccores after this one
      if (n_det_remote(i_core)>=i) then
        ncores_before_offset = ncores_before_offset+1
      endif
    enddo
  enddo
  

  ! MJO - Below does global indices for tiled data - we now use interleaved
  ! l = 1
  ! do i_core=0,ncores-1
  !   do i_it = 1,n_it
  !     offset = 0
  !     ! sum the offset contributions from the tiles from all cores' previous iterations
  !     do j=0,ncores-1
  !       do k=1,i_it-1
  !         offset = offset + local_det_map%n_det_remote(k,j)
  !       enddo
  !     enddo
  !     ! now sum the offset contributions from this iteration's previous cores
  !     do j=0,i_core-1
  !       offset = offset + local_det_map%n_det_remote(i_it,j)
  !     enddo
  !     ! Now that we have the offsets, set the global_indices appropriately
  !     do j=1,local_det_map%n_det_remote(i_it,i_core)
  !       remote_det_map%global_indices(l) = offset+j
  !       l=l+1
  !     enddo
  !   enddo
  ! enddo
#endif
end subroutine mpi_allgatherv_new_dets

! ===========================================================================================================
! Take in the global list from remote_det_map, distribute it,
! then fill up the local_dets by hashing
! ===========================================================================================================
subroutine mpi_distribute_remote_det_map(remote_det_map,local_det_map)
  type(det_map), intent(inout) :: remote_det_map
  type(det_map_l), intent(inout) :: local_det_map
  integer, allocatable :: local_indices(:)
  integer ndets_global,det_owner,j,ndets,i,l

  ndets_global = local_det_map%ndets_global_old
  !allocate send and receive bufs
  write(6,*) 'ndets_global,ncores',ndets_global,ncores; call flush(6)
  allocate(local_indices(10*ndets_global/ncores)) !Allocate temporary space to store local indices

  !fixme MJO - need to support norb.gt.127

  !Send the full list to everyone
  if(master_core_node) then
    call mpi_bsend_between_nodes(remote_det_map%global_dets_up)
    call mpi_bsend_between_nodes(remote_det_map%global_dets_dn)
  endif
  call mpi_barr_in_node()

  j=0
  !Loop through the global list to find all of this cores local values
  do i=1,ndets_global
    det_owner = get_det_owner(remote_det_map%global_dets_up(i),remote_det_map%global_dets_dn(i))
    if (det_owner.eq.whoami) then
      j = j+1
      !This determinant belongs to us, store its index in our temporary local_indices array
      local_indices(j) = i
    endif
  enddo

  !  ndets = j-1
  ndets = j

  local_det_map%ndets = ndets
  allocate(local_det_map%local_indices(ndets))
  local_det_map%local_indices(1:ndets) = local_indices(1:ndets)

end subroutine mpi_distribute_remote_det_map

!------------------------------------------------------------------------------------------------------------
subroutine mpi_merge_sort2_num_denom_diag_elems_info(old_dets_up,old_dets_dn,n_det,dets_up,dets_dn,n_allocate,old_e_mix_num,old_e_mix_den,old_diag_elems_info,e_mix_num,e_mix_den,diag_elems_info,iorder,temp_i16_up,temp_i16_dn,temp_i_2)
! ===========================================================================================================
!
!
! -----------------------------------------------------------------------------------------------------------
  use tools, only : merge_sort2_up_dn,merge_original_with_spawned3 ! Note: this merge_original_with_spawned3 is very different from the routine in do_walk.f90!
  use common_run, only : diag_elem_info

  implicit none
  integer, allocatable, intent(inout) :: iorder(:)
  integer, intent(inout) :: n_det
  real(rk),allocatable,intent(inout) :: old_e_mix_num(:)
  real(rk),allocatable,intent(inout) :: old_e_mix_den(:)
  type(diag_elem_info),allocatable,intent(inout) :: old_diag_elems_info(:)
  real(rk),allocatable,intent(inout) :: e_mix_num(:)
  real(rk),allocatable,intent(inout) :: e_mix_den(:)
  type(diag_elem_info),allocatable,intent(inout) :: diag_elems_info(:)

  integer(i8b), intent(inout) :: n_allocate
  integer, allocatable, intent(inout) :: temp_i_2(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), allocatable, intent(inout) :: old_dets_up(:), old_dets_dn(:)
  type(ik_vec), allocatable, intent(inout) :: dets_up(:), dets_dn(:)
  type(ik_vec), allocatable, intent(inout) :: temp_i16_up(:), temp_i16_dn(:)
  type(ik_vec),allocatable :: send_buf_core_up(:,:),send_buf_core_dn(:,:)
#else
  integer(ik), allocatable, intent(inout) :: old_dets_up(:), old_dets_dn(:)
  integer(ik), allocatable, intent(inout) :: dets_up(:), dets_dn(:)
  integer(ik), allocatable, intent(inout) :: temp_i16_up(:), temp_i16_dn(:)
  integer(ik),allocatable :: send_buf_core_up(:,:),send_buf_core_dn(:,:)
#endif
  !Local variables
  integer(i8b) :: n_allocate_new
  integer :: core_id,i,j,k,l,send_size,recv_size,mpierr,n_det_remote(0:ncores-1)
  integer :: send_buf_cnt(0:ncores-1),recv_buf_cnt(0:ncores-1),send_buf_disp(0:ncores-1),recv_buf_disp(0:ncores-1)

  real(rk),allocatable :: send_buf_core_e_mix_num(:,:),send_buf_core_e_mix_den(:,:)
  type(diag_elem_info), allocatable :: send_buf_core_diag_elems_info(:,:)
  type(t_det_mix_diag),allocatable :: send_buf(:),recv_buf(:)

  !MJO remote merge/sort
  !MJO - pack buffers, send, recv, do another merge / sort
#ifdef MPI
  !Loop through list and figure out array sizes
  send_buf_cnt = 0
  do i=1,n_det
    core_id = get_det_owner(old_dets_up(i),old_dets_dn(i))
    send_buf_cnt(core_id) = send_buf_cnt(core_id) + 1
  enddo

  !Allocate arrays - We don't have to reallocate this way, but we have to loop over n_det twice
  allocate(send_buf_core_up(maxval(send_buf_cnt),0:ncores-1))
  allocate(send_buf_core_dn(maxval(send_buf_cnt),0:ncores-1))
  allocate(send_buf_core_e_mix_num(maxval(send_buf_cnt),0:ncores-1))
  allocate(send_buf_core_e_mix_den(maxval(send_buf_cnt),0:ncores-1))
  allocate(send_buf_core_diag_elems_info(maxval(send_buf_cnt),0:ncores-1))

  send_buf_cnt = 0
  do i=1,n_det
    core_id = get_det_owner(old_dets_up(i),old_dets_dn(i))
    !Fill up send buffers - we'll arrange them later
    send_buf_cnt(core_id)                            = send_buf_cnt(core_id) + 1
    send_buf_core_up(send_buf_cnt(core_id),core_id)  = old_dets_up(i)
    send_buf_core_dn(send_buf_cnt(core_id),core_id)  = old_dets_dn(i)
    send_buf_core_e_mix_num(send_buf_cnt(core_id),core_id) = old_e_mix_num(i)
    send_buf_core_e_mix_den(send_buf_cnt(core_id),core_id) = old_e_mix_den(i)
    send_buf_core_diag_elems_info(send_buf_cnt(core_id),core_id) = old_diag_elems_info(i)
  enddo
  call flush(6)
  !Send over send_buf_cnts
  call MPI_alltoall(send_buf_cnt,1,MPI_INTEGER,recv_buf_cnt,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)

  !calculate displacements
  send_buf_disp(0) = 0
  recv_buf_disp(0) = 0
  do i=1,ncores-1
    send_buf_disp(i) = send_buf_cnt(i-1)+send_buf_disp(i-1)
    recv_buf_disp(i) = recv_buf_cnt(i-1)+recv_buf_disp(i-1)
  enddo

  send_size = sum(send_buf_cnt)
  recv_size = sum(recv_buf_cnt)

  !allocate send and receive bufs
  allocate(send_buf(0:send_size))
  allocate(recv_buf(0:recv_size))


  k=0
  do j=0,ncores-1
    do i=1,send_buf_cnt(j)
#ifdef NUM_ORBITALS_GT_127
      do l=1,num_words
        call conv_128_to_64(send_buf_core_up(i,j)%v(l),send_buf(k)%det_up(2*l-1),send_buf(k)%det_up(2*l))
        call conv_128_to_64(send_buf_core_dn(i,j)%v(l),send_buf(k)%det_dn(2*l-1),send_buf(k)%det_dn(2*l))
      enddo
#else
      call conv_128_to_64(send_buf_core_up(i,j),send_buf(k)%det_up(1),send_buf(k)%det_up(2))
      call conv_128_to_64(send_buf_core_dn(i,j),send_buf(k)%det_dn(1),send_buf(k)%det_dn(2))
#endif
      send_buf(k)%e_mix_num = send_buf_core_e_mix_num(i,j)
      send_buf(k)%e_mix_den = send_buf_core_e_mix_den(i,j)
      send_buf(k)%diag_elems_info = send_buf_core_diag_elems_info(i,j)
      k = k+1
    enddo
  enddo

  !FIXME MJO - need to support norb.gt.127
  call MPI_alltoallv(send_buf,send_buf_cnt,send_buf_disp,MPI_DET_MIX_DIAG_TYPE,&
       recv_buf,recv_buf_cnt,recv_buf_disp,MPI_DET_MIX_DIAG_TYPE,&
       MPI_COMM_WORLD,mpierr)

  !Reallocate arrays if needed
  if (recv_size > n_allocate) then
    n_allocate_new = n_det+recv_size
    write(6,'(''4mpi: reallocating dets_up etc. arrays of size'',i10,'' ='',es11.4)') n_allocate_new, real(n_allocate_new)
    if(allocated(old_dets_up)) deallocate(old_dets_up,old_dets_dn)
    if(allocated(dets_up)) deallocate(dets_up,dets_dn)
    if(allocated(old_e_mix_num)) deallocate(old_e_mix_num,old_e_mix_den)
    if(allocated(e_mix_num)) deallocate(e_mix_num,e_mix_den)
    if(allocated(old_diag_elems_info)) deallocate(old_diag_elems_info)
    if(allocated(diag_elems_info)) deallocate(diag_elems_info)

    if(allocated(iorder)) deallocate(iorder)
    if(allocated(temp_i16_up)) deallocate(temp_i16_up)
    if(allocated(temp_i16_dn)) deallocate(temp_i16_dn)
    if(allocated(temp_i_2)) deallocate(temp_i_2)

    allocate(iorder(n_allocate_new))
    allocate(temp_i16_up((n_allocate_new+1)/2))
    allocate(temp_i16_dn((n_allocate_new+1)/2))
    allocate(temp_i_2((n_allocate_new+1)/2))

    allocate(old_dets_up(n_allocate_new))
    allocate(old_dets_dn(n_allocate_new))
    allocate(dets_up(n_allocate_new))
    allocate(dets_dn(n_allocate_new))
    allocate(old_e_mix_num(n_allocate_new))
    allocate(old_e_mix_den(n_allocate_new))
    allocate(e_mix_num(n_allocate_new)) !necessary? Doesn't look like it?
    allocate(e_mix_den(n_allocate_new)) !necessary? Doesn't look like it?
    allocate(old_diag_elems_info(n_allocate_new))
    allocate(diag_elems_info(n_allocate_new))
  endif

  ! !Unpack buffers - includes local part
  do i=1,recv_size
#ifdef NUM_ORBITALS_GT_127
    do l=1,num_words
      call conv_64_to_128(dets_up(i)%v(l),recv_buf(i-1)%det_up(2*l-1),recv_buf(i-1)%det_up(2*l))
      call conv_64_to_128(dets_dn(i)%v(l),recv_buf(i-1)%det_dn(2*l-1),recv_buf(i-1)%det_dn(2*l))
    enddo
#else
    call conv_64_to_128(dets_up(i),recv_buf(i-1)%det_up(1),recv_buf(i-1)%det_up(2))
    call conv_64_to_128(dets_dn(i),recv_buf(i-1)%det_dn(1),recv_buf(i-1)%det_dn(2))
#endif
    old_e_mix_num(i) = recv_buf(i-1)%e_mix_num
    old_e_mix_den(i) = recv_buf(i-1)%e_mix_den
    old_diag_elems_info(i) = recv_buf(i-1)%diag_elems_info
  enddo

  !merge sort the determinants from other cores, only if we received new determinants
  do j=1,recv_size
    iorder(j)=j
  enddo

  call merge_sort2_up_dn(old_dets_up(1:recv_size),old_dets_dn(1:recv_size), iorder, recv_size, temp_i16_up, temp_i16_dn, temp_i_2)

  n_det = recv_size !Set our new number of determinants

  old_e_mix_num(1:n_det) = old_e_mix_num(iorder(1:n_det))
  if(allocated(e_mix_num)) deallocate(e_mix_num)

  old_e_mix_den(1:n_det) = old_e_mix_den(iorder(1:n_det))
  if(allocated(e_mix_den)) deallocate(e_mix_den)

  old_diag_elems_info(1:n_det) = old_diag_elems_info(iorder(1:n_det))
  call merge_original_with_spawned3(old_dets_up,old_dets_dn,n_det,dets_up,dets_dn,old_e_mix_num,old_e_mix_den,old_diag_elems_info,e_mix_num,e_mix_den,diag_elems_info)

  deallocate(send_buf)
  deallocate(recv_buf)
  deallocate(send_buf_core_dn)
  deallocate(send_buf_core_up)

  call MPI_ALLGATHER(n_det,1,MPI_INTEGER,n_det_remote,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
  ndets_global = sum(n_det_remote)
#endif

end subroutine mpi_merge_sort2_num_denom_diag_elems_info


subroutine mpi_merge_sort2_num_denom_diag_elems_info_replace(n_det,dets_up,dets_dn,n_allocate,e_mix_num,e_mix_den,diag_elems_info,iorder,temp_i16_up,temp_i16_dn,temp_i_2)
! ===========================================================================================================
!
!
! -----------------------------------------------------------------------------------------------------------
  use tools, only : merge_sort2_up_dn,merge_original_with_spawned3 ! Note: this merge_original_with_spawned3 is very different from the routine in do_walk.f90!
  use common_run, only : diag_elem_info

  implicit none
  integer, allocatable, intent(inout) :: iorder(:)
  integer, intent(inout) :: n_det
  real(rk),allocatable,intent(inout) :: e_mix_num(:)
  real(rk),allocatable,intent(inout) :: e_mix_den(:)
  type(diag_elem_info),allocatable,intent(inout) :: diag_elems_info(:)

  integer(i8b), intent(inout) :: n_allocate
  integer, allocatable, intent(inout) :: temp_i_2(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), allocatable, intent(inout) :: dets_up(:), dets_dn(:)
  type(ik_vec), allocatable, intent(inout) :: temp_i16_up(:), temp_i16_dn(:)
#else
  integer(ik), allocatable, intent(inout) :: dets_up(:), dets_dn(:)
  integer(ik), allocatable, intent(inout) :: temp_i16_up(:), temp_i16_dn(:)
#endif
  !Local variables
  integer(i8b) :: n_allocate_new
  integer :: core_id,i,j,k,l,send_size,recv_size,mpierr,n_det_remote(0:ncores-1),index
  integer :: send_buf_cnt(0:ncores-1),recv_buf_cnt(0:ncores-1),send_buf_disp(0:ncores-1),recv_buf_disp(0:ncores-1)
  integer, allocatable :: send_buf_core_ind(:,:)
  type(t_det_mix_diag),allocatable :: send_buf(:),recv_buf(:)

  !MJO remote merge/sort
  !MJO - pack buffers, send, recv, do another merge / sort
#ifdef MPI
  !Loop through list and figure out array sizes
  send_buf_cnt = 0
  do i=1,n_det
    core_id = get_det_owner(dets_up(i),dets_dn(i))
    send_buf_cnt(core_id) = send_buf_cnt(core_id) + 1
  enddo

  !Allocate arrays - We don't have to reallocate this way, but we have to loop over n_det twice
  allocate(send_buf_core_ind(maxval(send_buf_cnt),0:ncores-1))

  send_buf_cnt = 0
  do i=1,n_det
    core_id = get_det_owner(dets_up(i),dets_dn(i))
    !Fill up send buffers - we'll arrange them later
    send_buf_cnt(core_id)                            = send_buf_cnt(core_id) + 1
    !Store index so we can fill the buffers properly later
    send_buf_core_ind(send_buf_cnt(core_id),core_id)  = i
  enddo


  !Send over send_buf_cnts
  call MPI_alltoall(send_buf_cnt,1,MPI_INTEGER,recv_buf_cnt,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
  !calculate displacements
  send_buf_disp(0) = 0
  recv_buf_disp(0) = 0
  do i=1,ncores-1
    send_buf_disp(i) = send_buf_cnt(i-1)+send_buf_disp(i-1)
    recv_buf_disp(i) = recv_buf_cnt(i-1)+recv_buf_disp(i-1)
  enddo

  send_size = sum(send_buf_cnt)
  recv_size = sum(recv_buf_cnt)


  !allocate send and receive bufs
  allocate(send_buf(0:send_size))
  allocate(recv_buf(0:recv_size))


  k=0
  do j=0,ncores-1
    do i=1,send_buf_cnt(j)
      index = send_buf_core_ind(i,j)
#ifdef NUM_ORBITALS_GT_127
      do l=1,num_words
        call conv_128_to_64(dets_up(index)%v(l),send_buf(k)%det_up(2*l-1),send_buf(k)%det_up(2*l))
        call conv_128_to_64(dets_dn(index)%v(l),send_buf(k)%det_dn(2*l-1),send_buf(k)%det_dn(2*l))
      enddo
#else
      call conv_128_to_64(dets_up(index),send_buf(k)%det_up(1),send_buf(k)%det_up(2))
      call conv_128_to_64(dets_dn(index),send_buf(k)%det_dn(1),send_buf(k)%det_dn(2))
#endif
      send_buf(k)%e_mix_num = e_mix_num(index)
      send_buf(k)%e_mix_den = e_mix_den(index)
      send_buf(k)%diag_elems_info = diag_elems_info(index)
      k = k+1
    enddo
  enddo

  !FIXME MJO - need to support norb.gt.127
  call MPI_alltoallv(send_buf,send_buf_cnt,send_buf_disp,MPI_DET_MIX_DIAG_TYPE,&
       recv_buf,recv_buf_cnt,recv_buf_disp,MPI_DET_MIX_DIAG_TYPE,&
       MPI_COMM_WORLD,mpierr)

  !Reallocate arrays if needed
  if (recv_size > n_allocate) then
    n_allocate_new = n_det+recv_size
    write(6,'(''5mpi: reallocating dets_up etc. arrays of size'',i10,'' ='',es11.4)') n_allocate_new, real(n_allocate_new)
    if(allocated(dets_up)) deallocate(dets_up,dets_dn)
    if(allocated(e_mix_num)) deallocate(e_mix_num,e_mix_den)
    if(allocated(diag_elems_info)) deallocate(diag_elems_info)

    if(allocated(iorder)) deallocate(iorder)
    if(allocated(temp_i16_up)) deallocate(temp_i16_up)
    if(allocated(temp_i16_dn)) deallocate(temp_i16_dn)
    if(allocated(temp_i_2)) deallocate(temp_i_2)

    allocate(iorder(n_allocate_new))
    allocate(temp_i16_up((n_allocate_new+1)/2))
    allocate(temp_i16_dn((n_allocate_new+1)/2))
    allocate(temp_i_2((n_allocate_new+1)/2))

    allocate(dets_up(n_allocate_new))
    allocate(dets_dn(n_allocate_new))
    allocate(e_mix_num(n_allocate_new))
    allocate(e_mix_den(n_allocate_new))
    allocate(diag_elems_info(n_allocate_new))
  endif

  ! !Unpack buffers - includes local part
  do i=1,recv_size
#ifdef NUM_ORBITALS_GT_127
    do l=1,num_words
      call conv_64_to_128(dets_up(i)%v(l),recv_buf(i-1)%det_up(2*l-1),recv_buf(i-1)%det_up(2*l))
      call conv_64_to_128(dets_dn(i)%v(l),recv_buf(i-1)%det_dn(2*l-1),recv_buf(i-1)%det_dn(2*l))
    enddo
#else
    call conv_64_to_128(dets_up(i),recv_buf(i-1)%det_up(1),recv_buf(i-1)%det_up(2))
    call conv_64_to_128(dets_dn(i),recv_buf(i-1)%det_dn(1),recv_buf(i-1)%det_dn(2))
#endif
    e_mix_num(i) = recv_buf(i-1)%e_mix_num
    e_mix_den(i) = recv_buf(i-1)%e_mix_den
    diag_elems_info(i) = recv_buf(i-1)%diag_elems_info
  enddo

  !merge sort the determinants from other cores, only if we received new determinants
  do j=1,recv_size
    iorder(j)=j
  enddo

  call merge_sort2_up_dn(dets_up(1:recv_size),dets_dn(1:recv_size), iorder, recv_size, temp_i16_up, temp_i16_dn, temp_i_2)

  n_det = recv_size !Set our new number of determinants

  e_mix_num(1:n_det) = e_mix_num(iorder(1:n_det))
  !if(allocated(e_mix_num)) deallocate(e_mix_num)

  e_mix_den(1:n_det) = e_mix_den(iorder(1:n_det))
  !if(allocated(e_mix_den)) deallocate(e_mix_den)

  diag_elems_info(1:n_det) = diag_elems_info(iorder(1:n_det))

  call merge_original_with_spawned3(n_det, dets_up, dets_dn, e_mix_num, e_mix_den, diag_elems_info)

  deallocate(send_buf)
  deallocate(recv_buf)

  call MPI_ALLGATHER(n_det,1,MPI_INTEGER,n_det_remote,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
  ndets_global = sum(n_det_remote)
#endif

end subroutine mpi_merge_sort2_num_denom_diag_elems_info_replace


subroutine mpi_merge_sort2_num_denom_diag_elems_info_term(old_dets_up,old_dets_dn,n_det,dets_up,dets_dn,n_allocate,old_e_mix_num,old_e_mix_den,old_diag_elems_info,old_term1_big,old_term2_big,e_mix_num,e_mix_den,diag_elems_info,term1_big,term2_big,iorder,temp_i16_up,temp_i16_dn,temp_i_2)
! ===========================================================================================================
!
!
! -----------------------------------------------------------------------------------------------------------
  use tools, only : merge_sort2_up_dn,merge_original_with_spawned3 ! Note: this merge_original_with_spawned3 is very different from the routine in do_walk.f90!
  use common_run, only : diag_elem_info

  implicit none
  integer, allocatable, intent(inout) :: iorder(:)
  integer, intent(inout) :: n_det
  real(rk),allocatable,intent(inout) :: old_e_mix_num(:),old_term1_big(:)
  real(rk),allocatable,intent(inout) :: old_e_mix_den(:),old_term2_big(:)
  type(diag_elem_info),allocatable,intent(inout) :: old_diag_elems_info(:)
  real(rk),allocatable,intent(inout) :: e_mix_num(:),term1_big(:)
  real(rk),allocatable,intent(inout) :: e_mix_den(:),term2_big(:)
  type(diag_elem_info),allocatable,intent(inout) :: diag_elems_info(:)

  integer(i8b), intent(inout) :: n_allocate
  integer, allocatable, intent(inout) :: temp_i_2(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), allocatable, intent(inout) :: old_dets_up(:), old_dets_dn(:)
  type(ik_vec), allocatable, intent(inout) :: dets_up(:), dets_dn(:)
  type(ik_vec), allocatable, intent(inout) :: temp_i16_up(:), temp_i16_dn(:)
  type(ik_vec), allocatable :: send_buf_core_up(:,:),send_buf_core_dn(:,:)
#else
  integer(ik), allocatable, intent(inout) :: old_dets_up(:), old_dets_dn(:)
  integer(ik), allocatable, intent(inout) :: dets_up(:), dets_dn(:)
  integer(ik), allocatable, intent(inout) :: temp_i16_up(:), temp_i16_dn(:)
  integer(ik), allocatable :: send_buf_core_up(:,:),send_buf_core_dn(:,:)
#endif
  !Local variables
  integer(i8b) :: n_allocate_new
  integer :: core_id,i,j,k,l,send_size,recv_size,mpierr,n_det_remote(0:ncores-1)
  integer :: send_buf_cnt(0:ncores-1),recv_buf_cnt(0:ncores-1),send_buf_disp(0:ncores-1),recv_buf_disp(0:ncores-1)
  real(rk),allocatable :: send_buf_core_e_mix_num(:,:),send_buf_core_e_mix_den(:,:)
  real(rk),allocatable :: send_buf_core_term1_big(:,:),send_buf_core_term2_big(:,:)
  type(diag_elem_info), allocatable :: send_buf_core_diag_elems_info(:,:)
  type(t_det_mix_diag_term),allocatable :: send_buf(:),recv_buf(:)

  !MJO remote merge/sort
  !MJO - pack buffers, send, recv, do another merge / sort
#ifdef MPI
  !Loop through list and figure out array sizes
  send_buf_cnt = 0
  do i=1,n_det
    core_id = get_det_owner(old_dets_up(i),old_dets_dn(i))
    send_buf_cnt(core_id) = send_buf_cnt(core_id) + 1
  enddo

  !Allocate arrays - We don't have to reallocate this way, but we have to loop over n_det twice
  allocate(send_buf_core_up(maxval(send_buf_cnt),0:ncores-1))
  allocate(send_buf_core_dn(maxval(send_buf_cnt),0:ncores-1))
  allocate(send_buf_core_e_mix_num(maxval(send_buf_cnt),0:ncores-1))
  allocate(send_buf_core_e_mix_den(maxval(send_buf_cnt),0:ncores-1))
  allocate(send_buf_core_term1_big(maxval(send_buf_cnt),0:ncores-1))
  allocate(send_buf_core_term2_big(maxval(send_buf_cnt),0:ncores-1))
  allocate(send_buf_core_diag_elems_info(maxval(send_buf_cnt),0:ncores-1))

  send_buf_cnt = 0
  do i=1,n_det
    core_id = get_det_owner(old_dets_up(i),old_dets_dn(i))
    !Fill up send buffers - we'll arrange them later
    send_buf_cnt(core_id)                            = send_buf_cnt(core_id) + 1
    send_buf_core_up(send_buf_cnt(core_id),core_id)  = old_dets_up(i)
    send_buf_core_dn(send_buf_cnt(core_id),core_id)  = old_dets_dn(i)
    send_buf_core_e_mix_num(send_buf_cnt(core_id),core_id) = old_e_mix_num(i)
    send_buf_core_e_mix_den(send_buf_cnt(core_id),core_id) = old_e_mix_den(i)
    send_buf_core_diag_elems_info(send_buf_cnt(core_id),core_id) = old_diag_elems_info(i)
    send_buf_core_term1_big(send_buf_cnt(core_id),core_id) = old_term1_big(i)
    send_buf_core_term2_big(send_buf_cnt(core_id),core_id) = old_term2_big(i)
  enddo
  call flush(6)
  !Send over send_buf_cnts

  call MPI_alltoall(send_buf_cnt,1,MPI_INTEGER,recv_buf_cnt,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)

  !calculate displacements
  send_buf_disp(0) = 0
  recv_buf_disp(0) = 0
  do i=1,ncores-1
    send_buf_disp(i) = send_buf_cnt(i-1)+send_buf_disp(i-1)
    recv_buf_disp(i) = recv_buf_cnt(i-1)+recv_buf_disp(i-1)
  enddo

  send_size = sum(send_buf_cnt)
  recv_size = sum(recv_buf_cnt)


  !allocate send and receive bufs
  allocate(send_buf(0:send_size))
  allocate(recv_buf(0:recv_size))


  k=0
  do j=0,ncores-1
    do i=1,send_buf_cnt(j)
#ifdef NUM_ORBITALS_GT_127
      do l=1,num_words
        call conv_128_to_64(send_buf_core_up(i,j)%v(l),send_buf(k)%det_up(2*l-1),send_buf(k)%det_up(2*l))
        call conv_128_to_64(send_buf_core_dn(i,j)%v(l),send_buf(k)%det_dn(2*l-1),send_buf(k)%det_dn(2*l))
      enddo
#else
      call conv_128_to_64(send_buf_core_up(i,j),send_buf(k)%det_up(1),send_buf(k)%det_up(2))
      call conv_128_to_64(send_buf_core_dn(i,j),send_buf(k)%det_dn(1),send_buf(k)%det_dn(2))
#endif
      send_buf(k)%e_mix_num = send_buf_core_e_mix_num(i,j)
      send_buf(k)%e_mix_den = send_buf_core_e_mix_den(i,j)
      send_buf(k)%term1_big = send_buf_core_term1_big(i,j)
      send_buf(k)%term2_big = send_buf_core_term2_big(i,j)
      send_buf(k)%diag_elems_info = send_buf_core_diag_elems_info(i,j)
      k = k+1
    enddo
  enddo

  !FIXME MJO - need to support norb.gt.127
  call MPI_alltoallv(send_buf,send_buf_cnt,send_buf_disp,MPI_DET_MIX_DIAG_TERM_TYPE,&
       recv_buf,recv_buf_cnt,recv_buf_disp,MPI_DET_MIX_DIAG_TERM_TYPE,&
       MPI_COMM_WORLD,mpierr)

  !Reallocate arrays if needed
  if (recv_size > n_allocate) then
    n_allocate_new = n_det+recv_size
    write(6,'(''6mpi: reallocating dets_up etc. arrays of size'',i10,'' ='',es11.4)') n_allocate_new, real(n_allocate_new)
    if(allocated(old_dets_up)) deallocate(old_dets_up,old_dets_dn)
    if(allocated(dets_up)) deallocate(dets_up,dets_dn)
    if(allocated(old_e_mix_num)) deallocate(old_e_mix_num,old_e_mix_den)
    if(allocated(old_term1_big)) deallocate(old_term1_big,old_term2_big)
    if(allocated(e_mix_num)) deallocate(e_mix_num,e_mix_den)
    if(allocated(term1_big)) deallocate(term1_big,term2_big)
    if(allocated(old_diag_elems_info)) deallocate(old_diag_elems_info)
    if(allocated(diag_elems_info)) deallocate(diag_elems_info)

    if(allocated(iorder)) deallocate(iorder)
    if(allocated(temp_i16_up)) deallocate(temp_i16_up)
    if(allocated(temp_i16_dn)) deallocate(temp_i16_dn)
    if(allocated(temp_i_2)) deallocate(temp_i_2)

    allocate(iorder(n_allocate_new))
    allocate(temp_i16_up((n_allocate_new+1)/2))
    allocate(temp_i16_dn((n_allocate_new+1)/2))
    allocate(temp_i_2((n_allocate_new+1)/2))

    allocate(old_dets_up(n_allocate_new))
    allocate(old_dets_dn(n_allocate_new))
    allocate(dets_up(n_allocate_new))
    allocate(dets_dn(n_allocate_new))
    allocate(old_e_mix_num(n_allocate_new))
    allocate(old_e_mix_den(n_allocate_new))
    allocate(old_term1_big(n_allocate_new))
    allocate(old_term2_big(n_allocate_new))
    allocate(e_mix_num(n_allocate_new)) !necessary? Doesn't look like it?
    allocate(e_mix_den(n_allocate_new)) !necessary? Doesn't look like it?
    allocate(term1_big(n_allocate_new))
    allocate(term2_big(n_allocate_new))
    allocate(old_diag_elems_info(n_allocate_new))
    allocate(diag_elems_info(n_allocate_new))
  endif

  ! !Unpack buffers - includes local part
  do i=1,recv_size
#ifdef NUM_ORBITALS_GT_127
    do l=1,num_words
      call conv_64_to_128(dets_up(i)%v(l),recv_buf(i-1)%det_up(2*l-1),recv_buf(i-1)%det_up(2*l))
      call conv_64_to_128(dets_dn(i)%v(l),recv_buf(i-1)%det_dn(2*l-1),recv_buf(i-1)%det_dn(2*l))
    enddo
#else
    call conv_64_to_128(dets_up(i),recv_buf(i-1)%det_up(1),recv_buf(i-1)%det_up(2))
    call conv_64_to_128(dets_dn(i),recv_buf(i-1)%det_dn(1),recv_buf(i-1)%det_dn(2))
#endif
    old_e_mix_num(i) = recv_buf(i-1)%e_mix_num
    old_e_mix_den(i) = recv_buf(i-1)%e_mix_den
    old_term1_big(i) = recv_buf(i-1)%term1_big
    old_term2_big(i) = recv_buf(i-1)%term2_big
    old_diag_elems_info(i) = recv_buf(i-1)%diag_elems_info
  enddo


  !merge sort the determinants from other cores, only if we received new determinants
  do j=1,recv_size
    iorder(j)=j
  enddo

  call merge_sort2_up_dn(old_dets_up(1:recv_size),old_dets_dn(1:recv_size), iorder, recv_size, temp_i16_up, temp_i16_dn, temp_i_2)

  n_det = recv_size !Set our new number of determinants

  old_e_mix_num(1:n_det) = old_e_mix_num(iorder(1:n_det))
  if(allocated(e_mix_num)) deallocate(e_mix_num)

  old_e_mix_den(1:n_det) = old_e_mix_den(iorder(1:n_det))
  if(allocated(e_mix_den)) deallocate(e_mix_den)

  old_term1_big(1:n_det) = old_term1_big(iorder(1:n_det))
  if(allocated(term1_big)) deallocate(term1_big)

  old_term2_big(1:n_det) = old_term2_big(iorder(1:n_det))
  if(allocated(term2_big)) deallocate(term2_big)


  old_diag_elems_info(1:n_det) = old_diag_elems_info(iorder(1:n_det))
  call merge_original_with_spawned3(old_dets_up,old_dets_dn,n_det,dets_up,dets_dn,old_e_mix_num,old_e_mix_den,old_diag_elems_info,old_term1_big,old_term2_big,e_mix_num,e_mix_den,diag_elems_info,term1_big,term2_big)

  deallocate(send_buf)
  deallocate(recv_buf)
  deallocate(send_buf_core_dn)
  deallocate(send_buf_core_up)

  call MPI_ALLGATHER(n_det,1,MPI_INTEGER,n_det_remote,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
  ndets_global = sum(n_det_remote)
#endif

end subroutine mpi_merge_sort2_num_denom_diag_elems_info_term

subroutine test_mpi_merge_sort2_term()
  use common_run, only : diag_elem_info
  use overload
  implicit none
  integer, allocatable :: iorder(:)
  real(rk),allocatable :: e_mix_num(:),term1_big(:)
  real(rk),allocatable :: e_mix_den(:),term2_big(:)
  type(diag_elem_info),allocatable :: diag_elems_info(:)

  integer, allocatable :: temp_i_2(:)

#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), allocatable :: dets_up(:), dets_dn(:)
  type(ik_vec), allocatable :: temp_i16_up(:), temp_i16_dn(:)
#else
  integer(ik), allocatable :: dets_up(:), dets_dn(:)
  integer(ik), allocatable :: temp_i16_up(:), temp_i16_dn(:)
#endif
  integer ndet,j,i,number_excitations_up,number_excitations_dn,l
  integer(i8b) n_allocate
  ndet = 50
  n_allocate = 50
  allocate(dets_up(ndet))
  allocate(dets_dn(ndet))
  allocate(e_mix_num(ndet))
  allocate(e_mix_den(ndet))
  allocate(term1_big(ndet))
  allocate(term2_big(ndet))
  allocate(diag_elems_info(ndet))
  allocate(temp_i_2((ndet+1)/2))
  allocate(temp_i16_up((ndet+1)/2))
  allocate(temp_i16_dn((ndet+1)/2))
  allocate(iorder(ndet))

#if DNUM_ORB_GT_127
  do j=1,ndet
    do l=1,num_words
      dets_up(j)%v(l) = 0
      dets_dn(j)%v(l) = 0
    enddo
    do i=1,10
      dets_up(j) = ibset(dets_up(j),j+mod(j,i)*9+whoami)
      dets_dn(j) = ibset(dets_dn(j),j+mod(j,i)*10+whoami)
    enddo
  enddo

  do j=1,ndet
    number_excitations_up = popcnt(dets_up(j))
    number_excitations_dn = popcnt(dets_dn(j))
    write(6,*) 'num_excitations_before',number_excitations_up,number_excitations_dn
    !write(6,'(4B128)') dets_up%v(1),dets_dn%v(2), dets_up%v(1),dets_up%v(2)
  enddo


  call mpi_merge_sort2_num_denom_diag_elems_info_term_replace(ndet,dets_up,dets_dn,n_allocate,e_mix_num,e_mix_den,diag_elems_info,term1_big,term2_big,iorder,temp_i16_up,temp_i16_dn,temp_i_2)

  write(6,*) 'n_after',ndet
  do j=1,ndet
    number_excitations_up = popcnt(dets_up(j))
    number_excitations_dn = popcnt(dets_dn(j))
    write(6,*) 'num_excitations_after',number_excitations_up,number_excitations_dn
    !write(6,'(4B128)') dets_up%v(1),dets_dn%v(2), dets_up%v(1),dets_up%v(2)
  enddo
  call flush(6)
  call mpi_barr()
  stop
#endif
end subroutine test_mpi_merge_sort2_term




subroutine mpi_merge_sort2_num_denom_diag_elems_info_term_replace(n_det,dets_up,dets_dn,n_allocate,e_mix_num,e_mix_den,diag_elems_info,term1_big,term2_big,iorder,temp_i16_up,temp_i16_dn,temp_i_2)
! ===========================================================================================================
!
! -----------------------------------------------------------------------------------------------------------
  use tools, only : merge_sort2_up_dn,merge_original_with_spawned3 ! Note: this merge_original_with_spawned3 is very different from the routine in do_walk.f90!
  use common_run, only : diag_elem_info

  implicit none
  integer, allocatable, intent(inout) :: iorder(:)
  integer, intent(inout) :: n_det
  real(rk),allocatable,intent(inout) :: e_mix_num(:),term1_big(:)
  real(rk),allocatable,intent(inout) :: e_mix_den(:),term2_big(:)
  type(diag_elem_info),allocatable,intent(inout) :: diag_elems_info(:)

  integer(i8b), intent(inout) :: n_allocate
  integer, allocatable, intent(inout) :: temp_i_2(:)
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec), allocatable, intent(inout) :: dets_up(:), dets_dn(:)
  type(ik_vec), allocatable, intent(inout) :: temp_i16_up(:), temp_i16_dn(:)
#else
  integer(ik), allocatable, intent(inout) :: dets_up(:), dets_dn(:)
  integer(ik), allocatable, intent(inout) :: temp_i16_up(:), temp_i16_dn(:)
#endif
  !Local variables
  integer(i8b) :: n_allocate_new
  integer :: core_id,i,j,k,l,send_size,recv_size,mpierr,n_det_remote(0:ncores-1),index
  integer :: send_buf_cnt(0:ncores-1),recv_buf_cnt(0:ncores-1),send_buf_disp(0:ncores-1),recv_buf_disp(0:ncores-1)
  integer, allocatable  :: send_buf_core_ind(:,:)
  integer(i8b), allocatable :: high_send(:),low_send(:),high_recv(:),low_recv(:)
  real(rk), allocatable :: tmp_rk_send1(:), tmp_rk_send2(:), tmp_rk_recv1(:), tmp_rk_recv2(:)
  type(diag_elem_info), allocatable :: diag_elems_send(:),diag_elems_recv(:)
  !MJO remote merge/sort
  !MJO - pack buffers, send, recv, do another merge / sort
#ifdef MPI
  !Loop through list and figure out array sizes
  send_buf_cnt = 0
  do i=1,n_det
    core_id = get_det_owner(dets_up(i),dets_dn(i))
    send_buf_cnt(core_id) = send_buf_cnt(core_id) + 1
  enddo

  !Deallocate our sorting arrays initially; they will be reallocated later
  if(allocated(iorder)) deallocate(iorder)
  if(allocated(temp_i16_up)) deallocate(temp_i16_up)
  if(allocated(temp_i16_dn)) deallocate(temp_i16_dn)
  if(allocated(temp_i_2)) deallocate(temp_i_2)

  !Allocate arrays - We don't have to reallocate this way, but we have to loop over n_det twice
  allocate(send_buf_core_ind(maxval(send_buf_cnt),0:ncores-1))

  send_buf_cnt = 0
  do i=1,n_det
    core_id = get_det_owner(dets_up(i),dets_dn(i))
    !Fill up send buffers - we'll arrange them later
    send_buf_cnt(core_id)                            = send_buf_cnt(core_id) + 1
    !Store index so we can fill the buffers properly later
    send_buf_core_ind(send_buf_cnt(core_id),core_id) = i
  enddo

  !Send over send_buf_cnts

  call MPI_alltoall(send_buf_cnt,1,MPI_INTEGER,recv_buf_cnt,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)

  !Increase the size by num_words, to support norb_gt_127
  send_buf_cnt = send_buf_cnt*num_words
  recv_buf_cnt = recv_buf_cnt*num_words
  !calculate displacements
  send_buf_disp(0) = 0
  recv_buf_disp(0) = 0
  do i=1,ncores-1
    send_buf_disp(i) = send_buf_cnt(i-1)+send_buf_disp(i-1)
    recv_buf_disp(i) = recv_buf_cnt(i-1)+recv_buf_disp(i-1)
  enddo
  send_size = sum(send_buf_cnt)
  recv_size = sum(recv_buf_cnt)
  write(6,'(''7mpi: reallocating dets_up etc. arrays of size'',i10,'' ='',es11.4)') recv_size, real(recv_size); call flush(6)

  !allocate send and receive bufs for det_up
  allocate(high_send(0:send_size))
  allocate(low_send(0:send_size))
  allocate(high_recv(0:recv_size))
  allocate(low_recv(0:recv_size))

  !Fill det_up
  k=0
  do j=0,ncores-1
    do i=1,send_buf_cnt(j)/num_words !We want to loop through determinants, not their words
      index = send_buf_core_ind(i,j)
#ifdef NUM_ORBITALS_GT_127
      do l=1,num_words
        call conv_128_to_64(dets_up(index)%v(l),low_send(k),high_send(k))
        k=k+1
      enddo
#else
      call conv_128_to_64(dets_up(index),low_send(k),high_send(k))
      k = k+1
#endif
    enddo
  enddo


  call MPI_alltoallv(high_send,send_buf_cnt,send_buf_disp,MPI_INTEGER8,&
    high_recv,recv_buf_cnt,recv_buf_disp,MPI_INTEGER8,&
    MPI_COMM_WORLD,mpierr)

  call MPI_alltoallv(low_send,send_buf_cnt,send_buf_disp,MPI_INTEGER8,&
    low_recv,recv_buf_cnt,recv_buf_disp,MPI_INTEGER8,&
    MPI_COMM_WORLD,mpierr)

  !Reallocate arrays
  n_allocate_new = recv_size/num_words !MJO divide since we received num_words slices per determinant

  if(allocated(dets_up)) deallocate(dets_up)
  allocate(dets_up(n_allocate_new))

  ! Unpack buffers - includes local part
  do i=1,recv_size/num_words !Divide by num_words because that is how many dets we actually have
#ifdef NUM_ORBITALS_GT_127
    do l=1,num_words
      ! 2*(i-1) + l-1 gives us the num_words contiguously packed parts of a single det
      call conv_64_to_128(dets_up(i)%v(l),low_recv(2*(i-1)+l-1),high_recv(2*(i-1)+l-1))
    enddo
#else
    call conv_64_to_128(dets_up(i),low_recv(i-1),high_recv(i-1))
#endif
  enddo

  !Fill det_dn
  k=0
  do j=0,ncores-1
    do i=1,send_buf_cnt(j)/num_words !We want to loop through determinants
      index = send_buf_core_ind(i,j)
#ifdef NUM_ORBITALS_GT_127
      do l=1,num_words
        call conv_128_to_64(dets_dn(index)%v(l),low_send(k),high_send(k))
        k=k+1
      enddo
#else
      call conv_128_to_64(dets_dn(index),low_send(k),high_send(k))
      k = k+1
#endif
    enddo
  enddo


  call MPI_alltoallv(high_send,send_buf_cnt,send_buf_disp,MPI_INTEGER8,&
    high_recv,recv_buf_cnt,recv_buf_disp,MPI_INTEGER8,&
    MPI_COMM_WORLD,mpierr)
  call MPI_alltoallv(low_send,send_buf_cnt,send_buf_disp,MPI_INTEGER8,&
    low_recv,recv_buf_cnt,recv_buf_disp,MPI_INTEGER8,&
    MPI_COMM_WORLD,mpierr)

  n_allocate_new = recv_size/num_words
  if(allocated(dets_dn)) deallocate(dets_dn)
  allocate(dets_dn(n_allocate_new))


  ! Unpack buffers - includes local part
  do i=1,recv_size/num_words
#ifdef NUM_ORBITALS_GT_127
    do l=1,num_words
      call conv_64_to_128(dets_dn(i)%v(l),low_recv(2*(i-1)+l-1),high_recv(2*(i-1)+l-1))
    enddo
#else
    call conv_64_to_128(dets_dn(i),low_recv(i-1),high_recv(i-1))
#endif
  enddo

  if(allocated(high_send)) deallocate(high_send,low_send,high_recv,low_recv)

  send_buf_cnt = send_buf_cnt/num_words
  recv_buf_cnt = recv_buf_cnt/num_words
  !recalculate displacements
  send_buf_disp(0) = 0
  recv_buf_disp(0) = 0
  do i=1,ncores-1
    send_buf_disp(i) = send_buf_cnt(i-1)+send_buf_disp(i-1)
    recv_buf_disp(i) = recv_buf_cnt(i-1)+recv_buf_disp(i-1)
  enddo

  !Recalculate send/recv size without the extra num_words
  send_size = sum(send_buf_cnt)
  recv_size = sum(recv_buf_cnt)

  ! Allocate space for pairs of reals
  ! We do two at a time because the memory cost for 1 det
  ! is 16 bytes; two rk is the same. This cleans up
  ! the coding a bit
  allocate(tmp_rk_send1(send_size))
  allocate(tmp_rk_send2(send_size))
  allocate(tmp_rk_recv1(recv_size))
  allocate(tmp_rk_recv2(recv_size))

  ! Fill e_mix_num, e_mix_den
  k=1
  do j=0,ncores-1
    do i=1,send_buf_cnt(j)
      index = send_buf_core_ind(i,j)
      tmp_rk_send1(k) = e_mix_num(index)
      tmp_rk_send2(k) = e_mix_den(index)
      k = k+1
    enddo
  enddo

  call MPI_alltoallv(tmp_rk_send1,send_buf_cnt,send_buf_disp,MPI_REAL8,&
    tmp_rk_recv1,recv_buf_cnt,recv_buf_disp,MPI_REAL8,&
    MPI_COMM_WORLD,mpierr)
  call MPI_alltoallv(tmp_rk_send2,send_buf_cnt,send_buf_disp,MPI_REAL8,&
    tmp_rk_recv2,recv_buf_cnt,recv_buf_disp,MPI_REAL8,&
    MPI_COMM_WORLD,mpierr)


  n_allocate_new = recv_size
  if(allocated(e_mix_num)) deallocate(e_mix_num,e_mix_den)
  allocate(e_mix_num(n_allocate_new))
  allocate(e_mix_den(n_allocate_new))

  ! Unpack buffers - includes local part
  do i=1,recv_size
    e_mix_num(i) = tmp_rk_recv1(i)
    e_mix_den(i) = tmp_rk_recv2(i)
  enddo

  ! Fill term1_big, term2_big
  k=1
  do j=0,ncores-1
    do i=1,send_buf_cnt(j)
      index = send_buf_core_ind(i,j)
      tmp_rk_send1(k) = term1_big(index)
      tmp_rk_send2(k) = term2_big(index)
      k = k+1
    enddo
  enddo


  call MPI_alltoallv(tmp_rk_send1,send_buf_cnt,send_buf_disp,MPI_REAL8,&
    tmp_rk_recv1,recv_buf_cnt,recv_buf_disp,MPI_REAL8,&
    MPI_COMM_WORLD,mpierr)
  call MPI_alltoallv(tmp_rk_send2,send_buf_cnt,send_buf_disp,MPI_REAL8,&
    tmp_rk_recv2,recv_buf_cnt,recv_buf_disp,MPI_REAL8,&
    MPI_COMM_WORLD,mpierr)

  !Reallocate arrays if needed
  n_allocate_new = recv_size
  if(allocated(term1_big)) deallocate(term1_big,term2_big)
  allocate(term1_big(n_allocate_new))
  allocate(term2_big(n_allocate_new))

  ! Unpack buffers - includes local part
  do i=1,recv_size
    term1_big(i) = tmp_rk_recv1(i)
    term2_big(i) = tmp_rk_recv2(i)
  enddo

  if(allocated(tmp_rk_send1)) deallocate(tmp_rk_send1,tmp_rk_send2,tmp_rk_recv1,tmp_rk_recv2)

  ! Fill diag_elems
  allocate(diag_elems_send(send_size))
  allocate(diag_elems_recv(recv_size))

  k=1
  do j=0,ncores-1
    do i=1,send_buf_cnt(j)
      index = send_buf_core_ind(i,j)
      diag_elems_send(k) = diag_elems_info(index)
      k = k+1
    enddo
  enddo
  call MPI_alltoallv(diag_elems_send,send_buf_cnt,send_buf_disp,MPI_DIAG_ELEM_INFO_TYPE,&
    diag_elems_recv,recv_buf_cnt,recv_buf_disp,MPI_DIAG_ELEM_INFO_TYPE,&
    MPI_COMM_WORLD,mpierr)



  n_allocate_new = recv_size
  if(allocated(diag_elems_info)) deallocate(diag_elems_info)
  allocate(diag_elems_info(n_allocate_new))

  !Also reallocate the auxiliary sorting arrays

  allocate(iorder(n_allocate_new))
  allocate(temp_i16_up((n_allocate_new+1)/2))
  allocate(temp_i16_dn((n_allocate_new+1)/2))
  allocate(temp_i_2((n_allocate_new+1)/2))

  ! Unpack buffers - includes local part
  do i=1,recv_size
    diag_elems_info(i) = diag_elems_recv(i)
  enddo

  if(allocated(diag_elems_send)) deallocate(diag_elems_send,diag_elems_recv)

  !merge sort the determinants from other cores, only if we received new determinants
  do j=1,recv_size
    iorder(j)=j
  enddo

  call merge_sort2_up_dn(dets_up(1:recv_size),dets_dn(1:recv_size), iorder, recv_size, temp_i16_up, temp_i16_dn, temp_i_2)

  n_det = recv_size !Set our new number of determinants

  e_mix_num(1:n_det) = e_mix_num(iorder(1:n_det))
  e_mix_den(1:n_det) = e_mix_den(iorder(1:n_det))

  term1_big(1:n_det) = term1_big(iorder(1:n_det))
  term2_big(1:n_det) = term2_big(iorder(1:n_det))

  diag_elems_info(1:n_det) = diag_elems_info(iorder(1:n_det))

  call merge_original_with_spawned3(n_det, dets_up, dets_dn, e_mix_num, e_mix_den, diag_elems_info, term1_big, term2_big)

  call MPI_ALLGATHER(n_det,1,MPI_INTEGER,n_det_remote,1,MPI_INTEGER,MPI_COMM_WORLD,mpierr)
  ndets_global = sum(n_det_remote)

#endif

end subroutine mpi_merge_sort2_num_denom_diag_elems_info_term_replace


!---------------------------------------------------------------

!------------Shared Memory Allocation---------------------------
!-- This routine uses MPI3 routines the allocate shared memory--
!-- MPI routines used:
!-- MPI_WIN_ALLOCATE_SHARED, MPI_WIN_SHARED_QUERY
!-- Additionally, C_F_POINTER (from ISO_C_BINDING) is used
!-- to set up the pointers needed.
!--
!--Use:
!-- call shmem_allocate(array,numberOfElements)
!-- in: numberOfElements: the number of elements of the array
!--                       to be allocated
!-- out: array:           pointer to the sharedMemory array
!--                       allocated
!--
!-- call shmem_deallocate(array)
!-- in: array:            pointer to memory to be deallocated
!--
!--------------------------------------------------------------
!--Added by: Matthew Otten, April 26, 2014
!--Edited:   Matthew Otten, January 24, 2017 - Added shmem_deallocate
!---------------------------------------------------------------
!---------------------------------------------------------------
!-------------How to Implement----------------------------------
!-- for any array for which you want shared, declare it as
!-- real(rk),dimension(:),pointer :: array
!-- and then call shmem_allocate(array,numberOfElements)
!--
!-- All shmem_allocated arrays are automatically freed at
!-- the end of the run. If you want to specifically free
!-- a shmem_allocated array during the run,
!-- call shmem_deallocate(array)
!--
!-- To change from arrays that are currently allocated:
!-- change real(rk),allocatable :: array(:) to
!-- real(rk),dimension(:),pointer :: array
!-- and change allocate(array(numberOfElements))
!-- to call shmem_allocate(array,numberOfElements)
!---------------------------------------------------------------

subroutine shmem_deallocate2_rk(array)
  real(rk),dimension(:,:),pointer,intent(inout) :: array
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  integer                                     :: i, this_win_num
  type(C_PTR)                                 :: array_cptr
#ifdef MPI
  !Find which window this array belongs too
  this_win_num = -1
  do i=1,win_num
     call MPI_WIN_SHARED_QUERY(shmemWindow(i),0,shmem_buf_size,shmem_buf_size,array_cptr,mpierr)
     ! C_ASSOCIATED checks if two c-like pointers point to the same array
     ! C_LOC gives the c-like pointer of a fortran array
     ! array(1) is because we want the pointer to the first element
     if (C_ASSOCIATED(array_cptr,C_LOC(array(1,1)))) then
        this_win_num = i
        exit
     endif
  enddo

  if (this_win_num.eq.-1) then
     ! If we didn't find a shmem_window for this array,
     ! print a warning and deallocate it normally
     write(6,*) "Warning! Tried to deallocate an array that was not allocated via shmem_allocate."
     write(6,*) "         Deallocating in the normal way, but this may not be what you intended."
     deallocate(array)
  else
     !Free this specific window
     call MPI_WIN_FREE(shmemWindow(this_win_num),mpierr)

     ! We need to 'pop' this window out of our window list
     do i=this_win_num,win_num-1
        shmemWindow(i) = shmemWindow(i+1)
     enddo
     shmemWindow(win_num) = 0 !0 out window so we don't reference it later
     win_num = win_num - 1 !Decrement the total number of shmem windows
  endif
#else
  !just deallocate the array normally
  deallocate(array)
#endif
end subroutine shmem_deallocate2_rk

subroutine shmem_deallocate_rk(array)
  real(rk),dimension(:),pointer,intent(inout) :: array
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  integer                                     :: i, this_win_num
  type(C_PTR)                                 :: array_cptr
#ifdef MPI
  !Find which window this array belongs too
  this_win_num = -1
  do i=1,win_num
    call MPI_WIN_SHARED_QUERY(shmemWindow(i),0,shmem_buf_size,shmem_buf_size,array_cptr,mpierr)
    ! C_ASSOCIATED checks if two c-like pointers point to the same array
    ! C_LOC gives the c-like pointer of a fortran array
    ! array(1) is because we want the pointer to the first element
    if (C_ASSOCIATED(array_cptr,C_LOC(array(1)))) then
      this_win_num = i
      exit
    endif
  enddo

  if (this_win_num.eq.-1) then
     ! If we didn't find a shmem_window for this array,
     ! print a warning and deallocate it normally
     write(6,*) "Warning! Tried to deallocate an array that was not allocated via shmem_allocate."
     write(6,*) "         Deallocating in the normal way, but this may not be what you intended."
     call flush(6)
     deallocate(array)
   else
     !Free this specific window
     call MPI_WIN_FREE(shmemWindow(this_win_num),mpierr)

     ! We need to 'pop' this window out of our window list
     do i=this_win_num,win_num-1
        shmemWindow(i) = shmemWindow(i+1)
     enddo
     shmemWindow(win_num) = 0 !0 out window so we don't reference it later
     win_num = win_num - 1 !Decrement the total number of shmem windows
  endif
#else
  !just deallocate the array normally
  deallocate(array)
#endif
end subroutine shmem_deallocate_rk

subroutine shmem_deallocate_int(array)
  integer,dimension(:),pointer,intent(inout) :: array
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  integer                                     :: i, this_win_num
  type(C_PTR)                                 :: array_cptr
#ifdef MPI
  !Find which window this array belongs too
  this_win_num = -1
  do i=1,win_num
     call MPI_WIN_SHARED_QUERY(shmemWindow(i),0,shmem_buf_size,shmem_buf_size,array_cptr,mpierr)
     ! C_ASSOCIATED checks if two c-like pointers point to the same array
     ! C_LOC gives the c-like pointer of a fortran array
     ! array(1) is because we want the pointer to the first element
     if (C_ASSOCIATED(array_cptr,C_LOC(array(1)))) then
        this_win_num = i
        exit
     endif
  enddo

  if (this_win_num.eq.-1) then
     ! If we didn't find a shmem_window for this array,
     ! print a warning and deallocate it normally
     write(6,*) "Warning! Tried to deallocate an array that was not allocated via shmem_allocate."
     write(6,*) "         Deallocating in the normal way, but this may not be what you intended."
     call flush(6)
     deallocate(array)
  else
     !Free this specific window
     call MPI_WIN_FREE(shmemWindow(this_win_num),mpierr)

     ! We need to 'pop' this window out of our window list
     do i=this_win_num,win_num-1
        shmemWindow(i) = shmemWindow(i+1)
     enddo
     shmemWindow(win_num) = 0 !0 out window so we don't reference it later
     win_num = win_num - 1 !Decrement the total number of shmem windows
  endif
#else
  !just deallocate the array normally
  deallocate(array)
#endif
end subroutine shmem_deallocate_int

subroutine shmem_deallocate_real(array)
  real,dimension(:),pointer,intent(inout) :: array
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  integer                                     :: i, this_win_num
  type(C_PTR)                                 :: array_cptr
#ifdef MPI
  !Find which window this array belongs too
  this_win_num = -1
  do i=1,win_num
     call MPI_WIN_SHARED_QUERY(shmemWindow(i),0,shmem_buf_size,shmem_buf_size,array_cptr,mpierr)
     ! C_ASSOCIATED checks if two c-like pointers point to the same array
     ! C_LOC gives the c-like pointer of a fortran array
     ! array(1) is because we want the pointer to the first element
     if (C_ASSOCIATED(array_cptr,C_LOC(array(1)))) then
        this_win_num = i
        exit
     endif
  enddo

  if (this_win_num.eq.-1) then
     ! If we didn't find a shmem_window for this array,
     ! print a warning and deallocate it normally
     write(6,*) "Warning! Tried to deallocate an array that was not allocated via shmem_allocate."
     write(6,*) "         Deallocating in the normal way, but this may not be what you intended."
     call flush(6)
     deallocate(array)
  else
     !Free this specific window
     call MPI_WIN_FREE(shmemWindow(this_win_num),mpierr)

     ! We need to 'pop' this window out of our window list
     do i=this_win_num,win_num-1
        shmemWindow(i) = shmemWindow(i+1)
     enddo
     shmemWindow(win_num) = 0 !0 out window so we don't reference it later
     win_num = win_num - 1 !Decrement the total number of shmem windows
  endif
#else
  !just deallocate the array normally
  deallocate(array)
#endif
end subroutine shmem_deallocate_real

subroutine shmem_deallocate_ik_vec(array)
  type(ik_vec),dimension(:),pointer,intent(inout) :: array
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  integer                                     :: i, this_win_num
  type(C_PTR)                                 :: array_cptr
#ifdef MPI
  !Find which window this array belongs too
  this_win_num = -1
  do i=1,win_num
     call MPI_WIN_SHARED_QUERY(shmemWindow(i),0,shmem_buf_size,shmem_buf_size,array_cptr,mpierr)
     ! C_ASSOCIATED checks if two c-like pointers point to the same array
     ! C_LOC gives the c-like pointer of a fortran array
     ! array(1) is because we want the pointer to the first element
     if (C_ASSOCIATED(array_cptr,C_LOC(array(1)))) then
        this_win_num = i
        exit
     endif
  enddo

  if (this_win_num.eq.-1) then
     ! If we didn't find a shmem_window for this array,
     ! print a warning and deallocate it normally
     write(6,*) "Warning! Tried to deallocate an array that was not allocated via shmem_allocate."
     write(6,*) "         Deallocating in the normal way, but this may not be what you intended."
     call flush(6)
     deallocate(array)
  else
     !Free this specific window
     call MPI_WIN_FREE(shmemWindow(this_win_num),mpierr)

     ! We need to 'pop' this window out of our window list
     do i=this_win_num,win_num-1
       shmemWindow(i) = shmemWindow(i+1)
     enddo
     shmemWindow(win_num) = 0 !0 out window so we don't reference it later
     win_num = win_num - 1 !Decrement the total number of shmem windows
  endif
#else
  !just deallocate the array normally
  deallocate(array)
#endif
end subroutine shmem_deallocate_ik_vec


subroutine shmem_deallocate_ik(array)
  integer(ik),dimension(:),pointer,intent(inout) :: array
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  integer                                     :: i, this_win_num
  type(C_PTR)                                 :: array_cptr
#ifdef MPI
  !Find which window this array belongs too
  this_win_num = -1
  do i=1,win_num
     call MPI_WIN_SHARED_QUERY(shmemWindow(i),0,shmem_buf_size,shmem_buf_size,array_cptr,mpierr)
     ! C_ASSOCIATED checks if two c-like pointers point to the same array
     ! C_LOC gives the c-like pointer of a fortran array
     ! array(1) is because we want the pointer to the first element
     if (C_ASSOCIATED(array_cptr,C_LOC(array(1)))) then
        this_win_num = i
        exit
     endif
  enddo

  if (this_win_num.eq.-1) then
     ! If we didn't find a shmem_window for this array,
     ! print a warning and deallocate it normally
     write(6,*) "Warning! Tried to deallocate an array that was not allocated via shmem_allocate."
     write(6,*) "         Deallocating in the normal way, but this may not be what you intended."
     call flush(6)
     deallocate(array)
  else
     !Free this specific window
     call MPI_WIN_FREE(shmemWindow(this_win_num),mpierr)

     ! We need to 'pop' this window out of our window list
     do i=this_win_num,win_num-1
       shmemWindow(i) = shmemWindow(i+1)
     enddo
     shmemWindow(win_num) = 0 !0 out window so we don't reference it later
     win_num = win_num - 1 !Decrement the total number of shmem windows
  endif
#else
  !just deallocate the array normally
  deallocate(array)
#endif
end subroutine shmem_deallocate_ik

subroutine shmem_deallocate_rs_absH(array)
  type(rs_absH),dimension(:),pointer,intent(inout) :: array
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  integer                                     :: i, this_win_num
  type(C_PTR)                                 :: array_cptr
#ifdef MPI
  !Find which window this array belongs too
  this_win_num = -1
  do i=1,win_num
     call MPI_WIN_SHARED_QUERY(shmemWindow(i),0,shmem_buf_size,shmem_buf_size,array_cptr,mpierr)
     ! C_ASSOCIATED checks if two c-like pointers point to the same array
     ! C_LOC gives the c-like pointer of a fortran array
     ! array(1) is because we want the pointer to the first element
     if (C_ASSOCIATED(array_cptr,C_LOC(array(1)))) then
        this_win_num = i
        exit
     endif
  enddo

  if (this_win_num.eq.-1) then
     ! If we didn't find a shmem_window for this array,
     ! print a warning and deallocate it normally
     write(6,*) "Warning! Tried to deallocate an array that was not allocated via shmem_allocate."
     write(6,*) "         Deallocating in the normal way, but this may not be what you intended."
     call flush(6)
     deallocate(array)
  else
     !Free this specific window
     call MPI_WIN_FREE(shmemWindow(this_win_num),mpierr)

     ! We need to 'pop' this window out of our window list
     do i=this_win_num,win_num-1
       shmemWindow(i) = shmemWindow(i+1)
     enddo
     shmemWindow(win_num) = 0 !0 out window so we don't reference it later
     win_num = win_num - 1 !Decrement the total number of shmem windows
  endif
#else
  !just deallocate the array normally
  deallocate(array)
#endif
end subroutine shmem_deallocate_rs_absH


subroutine shmem_allocate_general_setup(shmem_buf_size_input,array_local_cptr)

#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND),intent(in)   ::  shmem_buf_size_input
  integer(kind=MPI_ADDRESS_KIND)              ::  shmem_buf_size
#else
  integer,intent(in)   ::  shmem_buf_size_input
  integer              ::  shmem_buf_size
#endif
  type(C_PTR),intent(out)                     :: array_local_cptr
  integer                                        my_ncores_on_node,disp_unit

  !figure out how many cores on a node
#ifdef MPI
  call MPI_COMM_SIZE(sharedComm,my_ncores_on_node,mpierr)
#endif
  !set the buffer size; it should be the total size/ncores_on_node
  !1._rk is just chosen as an arbitrary number of the type we
  !want the sizeof
  if (whoaminode.eq.0) then
    shmem_buf_size = shmem_buf_size_input
  else
#ifdef MPI
    shmem_buf_size = 0_MPI_ADDRESS_KIND
#else
    shmem_buf_size = 0
#endif
  endif

  write (6,'(''Allocating shared memory of size'',i12,'' ='',es11.4,'' bytes'')') shmem_buf_size, real(shmem_buf_size) ; call flush(6)
  !Allocate shared memory. Arguments: bufferSize, elementSize, MPI_INFO
  !communicator,c_ptr for array,window,mpierr
  win_num = win_num + 1
  if (win_num > max_shmem_windows) then
     write(6,*) "ERROR! Tried to allocate too many shmem arrays!"
     write(6,*) "       Consider increasing max_shmem_windows in mpi_routines.f90"
     call mpi_stop("Error in shmem_allocate")
  endif

  disp_unit = 1
#ifdef MPI
  call MPI_WIN_ALLOCATE_SHARED(shmem_buf_size,disp_unit,MPI_INFO_NULL,sharedComm,array_local_cptr,shmemWindow(win_num),mpierr)
  !query the master processor (of the window [so, of the sharedComm, or of the node]) for
  !the location of the start of the memory.
  call MPI_WIN_SHARED_QUERY(shmemWindow(win_num),0,shmem_buf_size,disp_unit,array_local_cptr,mpierr)
#endif
end subroutine shmem_allocate_general_setup


! Reallocate a shmem array
subroutine shmem_reallocate_ik(array,numberOfElements)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  use types, only : i8b
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  integer                                        my_ncores_on_node,disp_unit,i
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_local_cptr,array_cptr
  !This array is a f90 pointer
  integer(ik),dimension(:),pointer,intent(inout) :: array
  integer(ik),dimension(:),pointer               :: tmp_array
  integer(i8b),intent(in)                          :: numberOfElements
  integer(i8b)                          :: numberOfElementsOld
  integer(ik) :: test
  integer :: info

 !if MPI is defined, use the MPI sharedmemory things
#ifdef MPI

  shmem_buf_size = numberOfElements*sizeof(test)
  call shmem_allocate_general_setup(shmem_buf_size,array_cptr)

  !Map the new array to a tmp_array, just for now
  call C_F_POINTER(array_cptr,tmp_array,[numberOfElements])
  !Now, copy the old data (array) into tmp_array

  if (master_core_node) then
     do i=1,size(array)
        tmp_array(i) = array(i)
     enddo
  endif
  call mpi_barr_in_node()

  !Deallocate array
  call shmem_deallocate(array)

  !Finally, map our point to array
  call C_F_POINTER(array_cptr,array,[numberOfElements])

#else
  !otherwise, just reallocate it like normal.
  numberOfElementsOld = size(array)
  if (master_core_node) then
     allocate(tmp_array(numberOfElementsOld))
     do i=1,numberOfElementsOld
        tmp_array(i) = array(i)
     enddo
  endif

  deallocate(array)
  allocate(array(numberOfElements))

  if (master_core_node) then
     do i=1,numberOfElements
        array(i) = tmp_array(i)
     enddo
     deallocate(tmp_array)
  endif
#endif
end subroutine shmem_reallocate_ik

! Reallocate a shmem array
subroutine shmem_reallocate_ik_vec(array,numberOfElements)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  use types, only : i8b
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  integer                                        my_ncores_on_node,disp_unit,i
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_local_cptr,array_cptr
  !This array is a f90 pointer
  type(ik_vec),dimension(:),pointer,intent(inout) :: array
  type(ik_vec),dimension(:),pointer               :: tmp_array
  integer(i8b),intent(in)                          :: numberOfElements
  integer(i8b)                          :: numberOfElementsOld
  type(ik_vec) :: test
  integer :: info

 !if MPI is defined, use the MPI sharedmemory things
#ifdef MPI

  shmem_buf_size = numberOfElements*sizeof(test)
  call shmem_allocate_general_setup(shmem_buf_size,array_cptr)

  !Map the new array to a tmp_array, just for now
  call C_F_POINTER(array_cptr,tmp_array,[numberOfElements])
  !Now, copy the old data (array) into tmp_array

  if (master_core_node) then
     do i=1,size(array)
        tmp_array(i) = array(i)
     enddo
  endif
  call mpi_barr_in_node()

  !Deallocate array
  call shmem_deallocate(array)

  !Finally, map our point to array
  call C_F_POINTER(array_cptr,array,[numberOfElements])

#else
  !otherwise, just reallocate it like normal.
  numberOfElementsOld = size(array)
  if (master_core_node) then
     allocate(tmp_array(numberOfElementsOld))
     do i=1,numberOfElementsOld
        tmp_array(i) = array(i)
     enddo
  endif

  deallocate(array)
  allocate(array(numberOfElements))

  if (master_core_node) then
     do i=1,numberOfElements
        array(i) = tmp_array(i)
     enddo
     deallocate(tmp_array)
  endif
#endif
end subroutine shmem_reallocate_ik_vec


! Reallocate a shmem array
subroutine shmem_reallocate_rk(array,numberOfElements)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  use types, only : i8b
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  integer                                        my_ncores_on_node,disp_unit,i
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_local_cptr,array_cptr
  !This array is a f90 pointer
  real(rk),dimension(:),pointer,intent(inout) :: array
  real(rk),dimension(:),pointer               :: tmp_array
  integer(i8b),intent(in)                          :: numberOfElements
  integer(i8b)                          :: numberOfElementsOld
  real(rk) :: test
  integer :: info

 !if MPI is defined, use the MPI sharedmemory things
#ifdef MPI

  shmem_buf_size = numberOfElements*sizeof(test)
  call shmem_allocate_general_setup(shmem_buf_size,array_cptr)

  !Map the new array to a tmp_array, just for now
  call C_F_POINTER(array_cptr,tmp_array,[numberOfElements])
  !Now, copy the old data (array) into tmp_array

  if (master_core_node) then
     do i=1,size(array)
        tmp_array(i) = array(i)
     enddo
  endif
  call mpi_barr_in_node()

  !Deallocate array
  call shmem_deallocate(array)

  !Finally, map our point to array
  call C_F_POINTER(array_cptr,array,[numberOfElements])

#else
  !otherwise, just reallocate it like normal.
  numberOfElementsOld = size(array)
  if (master_core_node) then
     allocate(tmp_array(numberOfElementsOld))
     do i=1,numberOfElementsOld
        tmp_array(i) = array(i)
     enddo
  endif

  deallocate(array)
  allocate(array(numberOfElements))

  if (master_core_node) then
     do i=1,numberOfElements
        array(i) = tmp_array(i)
     enddo
     deallocate(tmp_array)
  endif
#endif
end subroutine shmem_reallocate_rk

! Reallocate a shmem array

subroutine shmem_reallocate_rs_absH(array,numberOfElements)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  use types, only : i8b
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  integer                                        my_ncores_on_node,disp_unit,i
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_local_cptr,array_cptr
  !This array is a f90 pointer
  type(rs_absH),dimension(:),pointer,intent(inout) :: array
  type(rs_absH),dimension(:),pointer               :: tmp_array
  integer(i8b),intent(in)                          :: numberOfElements
  integer(i8b)                          :: numberOfElementsOld
  type(rs_absH) :: test
  integer :: info

 !if MPI is defined, use the MPI sharedmemory things
#ifdef MPI

  shmem_buf_size = numberOfElements*sizeof(test)
  call shmem_allocate_general_setup(shmem_buf_size,array_cptr)

  !Map the new array to a tmp_array, just for now
  call C_F_POINTER(array_cptr,tmp_array,[numberOfElements])
  !Now, copy the old data (array) into tmp_array

  if (master_core_node) then
     do i=1,min(size(tmp_array),size(array))
        tmp_array(i) = array(i)
     enddo
  endif
  call mpi_barr_in_node()

  !Deallocate array
  call shmem_deallocate(array)

  !Finally, map our point to array
  call C_F_POINTER(array_cptr,array,[numberOfElements])

#else
  !otherwise, just reallocate it like normal.
  numberOfElementsOld = size(array)
  if (master_core_node) then
     allocate(tmp_array(numberOfElementsOld))
     do i=1,numberOfElementsOld
        tmp_array(i) = array(i)
     enddo
  endif

  deallocate(array)
  allocate(array(numberOfElements))

  if (master_core_node) then
     do i=1,numberOfElements
        array(i) = tmp_array(i)
     enddo
     deallocate(tmp_array)
  endif
#endif
end subroutine shmem_reallocate_rs_absH

! Reallocate a shmem array

subroutine shmem_reallocate2_rk(array,numberOfElements1,numberOfElements2)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  use types, only : i8b
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  integer                                        my_ncores_on_node,disp_unit,i,j
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_local_cptr,array_cptr
  !This array is a f90 pointer
  real(rk),dimension(:,:),pointer,intent(inout) :: array
  real(rk),dimension(:,:),pointer               :: tmp_array
  integer(i8b),intent(in)                          :: numberOfElements1,numberOfElements2
  integer(i8b)                          :: numberOfElementsOld1,numberOfElementsOld2
  real(rk) :: test
  integer :: info
  integer, allocatable :: arrayshape(:)

 !if MPI is defined, use the MPI sharedmemory things
#ifdef MPI

  shmem_buf_size = numberOfElements1*numberOfElements2*sizeof(test)
  call shmem_allocate_general_setup(shmem_buf_size,array_cptr)

  allocate(arrayshape(2))
  arrayshape = (/ numberOfElements1,numberOfElements2 /)
  !Map the new array to a tmp_array, just for now
  call C_F_POINTER(array_cptr,tmp_array,arrayshape)
  !Now, copy the old data (array) into tmp_array

  if (master_core_node) then
    do j=1,size(array,2)
      do i=1,size(array,1)
        tmp_array(i,j) = array(i,j)
      enddo
    enddo
  endif
  call mpi_barr_in_node()

  !Deallocate array
  call shmem_deallocate(array)

  !Finally, map our point to array
  call C_F_POINTER(array_cptr,array,arrayshape)

#else
  !otherwise, just reallocate it like normal.
  numberOfElementsOld1 = size(array,1)
  numberOfElementsOld2 = size(array,2)
  if (master_core_node) then
    allocate(tmp_array(numberOfElementsOld1,numberOfElementsOld2))
    do j=1,numberOfElementsOld2
      do i=1,numberOfElementsOld1
        tmp_array(i,j) = array(i,j)
      enddo
    enddo
  endif

  deallocate(array)
  allocate(array(numberOfElements1,numberOfElements2))

  if (master_core_node) then
    do j=1,numberOfElementsOld2
      do i=1,numberOfElementsOld1
        array(i,j) = tmp_array(i,j)
      enddo
    enddo
    deallocate(tmp_array)
  endif
#endif
end subroutine shmem_reallocate2_rk


subroutine shmem_allocate_rk(array,numberOfElements)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  use types, only : i8b
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_cptr
  !This array is a f90 pointer
  real(rk),dimension(:),pointer,intent(inout) :: array
  integer(i8b),intent(in)                          :: numberOfElements
  real(rk) :: test

 !if MPI is defined, use the MPI sharedmemory things
#ifdef MPI

  shmem_buf_size = numberOfElements*sizeof(test)
  call shmem_allocate_general_setup(shmem_buf_size,array_cptr)

  !remap that c_ptr to the array
  call C_F_POINTER(array_cptr,array,[numberOfElements])
#else
  !otherwise, just allocate it like normal.
  allocate(array(numberOfElements))
#endif
end subroutine shmem_allocate_rk
!------------------------------------------------------------------------------

subroutine shmem_allocate_rs_absH(array,numberOfElements)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  use types, only : i8b, rs_absH
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_cptr
  !This array is a f90 pointer
  type(rs_absH),dimension(:),pointer,intent(inout) :: array
  !type(optional_integer),dimension(:,:),pointer,intent(inout) :: array
  integer(i8b),intent(in)                          :: numberOfElements
  type(rs_absH) :: test

 !if MPI is defined, use the MPI sharedmemory things
#ifdef MPI

  shmem_buf_size = numberOfElements*sizeof(test)
  call shmem_allocate_general_setup(shmem_buf_size,array_cptr)

  !remap that c_ptr to the array
  call C_F_POINTER(array_cptr,array,[numberOfElements])
#else
  !otherwise, just allocate it like normal.
  allocate(array(numberOfElements))
#endif
end subroutine shmem_allocate_rs_absH
!------------------------------------------------------------------------------

subroutine shmem_allocate_real(array,numberOfElements)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  ! Added type real to Matt's shmem_allocate routine
  ! A Holmes, 5 Apr 2015
  use types, only : i8b
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_cptr
  !This array is a f90 pointer
  real,dimension(:),pointer,intent(inout) :: array
  integer(i8b),intent(in)                          :: numberOfElements
  real :: test


 !if MPI is defined, use the MPI sharedmenmory things
#ifdef MPI

  shmem_buf_size = numberOfElements*sizeof(test)
  call shmem_allocate_general_setup(shmem_buf_size,array_cptr)

  !remap that c_ptr to the array
  call C_F_POINTER(array_cptr,array,[numberOfElements])
#else
  !otherwise, just allocate it like normal.
  allocate(array(numberOfElements))
#endif
end subroutine shmem_allocate_real
!------------------------------------------------------------------------------

subroutine shmem_allocate_int(array,numberOfElements)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  ! Added type int to Matt's shmem_allocate routine
  ! A Holmes, 5 Apr 2015
  use types, only : i8b
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_cptr
  !This array is a f90 pointer
  integer,dimension(:),pointer,intent(inout) :: array
  integer(i8b),intent(in)                          :: numberOfElements
  integer :: test

 !if MPI is defined, use the MPI sharedmenmory things
#ifdef MPI

  shmem_buf_size = numberOfElements*sizeof(test)
  call shmem_allocate_general_setup(shmem_buf_size,array_cptr)

  !remap that c_ptr to the array
  call C_F_POINTER(array_cptr,array,[numberOfElements])
#else
  !otherwise, just allocate it like normal.
  allocate(array(numberOfElements))
#endif
end subroutine shmem_allocate_int
!------------------------------------------------------------------------------
subroutine shmem_allocate_logical(array,numberOfElements)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  ! Added type int to Matt's shmem_allocate routine
  ! A Holmes, 5 Apr 2015
  use types, only : i8b
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_cptr
  !This array is a f90 pointer
  logical,dimension(:),pointer,intent(inout) :: array
  integer(i8b),intent(in)                          :: numberOfElements
  integer :: test

 !if MPI is defined, use the MPI sharedmenmory things
#ifdef MPI

  shmem_buf_size = numberOfElements*sizeof(test)
  call shmem_allocate_general_setup(shmem_buf_size,array_cptr)

  !remap that c_ptr to the array
  call C_F_POINTER(array_cptr,array,[numberOfElements])
#else
  !otherwise, just allocate it like normal.
  allocate(array(numberOfElements))
#endif
end subroutine shmem_allocate_logical
!------------------------------------------------------------------------------
subroutine shmem_allocate_ik(array,numberOfElements)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  ! Added type int to Matt's shmem_allocate routine
  ! A Holmes, 5 Apr 2015
  use types, only : i8b, ik
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_cptr
  !This array is a f90 pointer
  integer(ik),dimension(:),pointer,intent(inout) :: array
  integer(i8b),intent(in)                          :: numberOfElements
  integer(ik) :: test

 !if MPI is defined, use the MPI sharedmenmory things
#ifdef MPI
  shmem_buf_size = numberOfElements*sizeof(test)
  call shmem_allocate_general_setup(shmem_buf_size,array_cptr)

  !remap that c_ptr to the array
  call C_F_POINTER(array_cptr,array,[numberOfElements])
#else
  !otherwise, just allocate it like normal.
  allocate(array(numberOfElements))
#endif
end subroutine shmem_allocate_ik
!------------------------------------------------------------------------------

subroutine shmem_allocate_ik_vec(array,numberOfElements)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  ! Added type int to Matt's shmem_allocate routine
  ! A Holmes, 5 Apr 2015
  use types, only : i8b, ik
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_cptr
  !This array is a f90 pointer
  type(ik_vec),dimension(:),pointer,intent(inout) :: array
  integer(i8b),intent(in)                          :: numberOfElements
  type(ik_vec) :: test

 !if MPI is defined, use the MPI sharedmenmory things
#ifdef MPI
  shmem_buf_size = numberOfElements*sizeof(test)
  call shmem_allocate_general_setup(shmem_buf_size,array_cptr)
  !remap that c_ptr to the array
  call C_F_POINTER(array_cptr,array,[numberOfElements])
#else
  !otherwise, just allocate it like normal.
  allocate(array(numberOfElements))
#endif
end subroutine shmem_allocate_ik_vec
!------------------------------------------------------------------------------

subroutine shmem_allocate2_rk(array,numberOfElements1,numberOfElements2)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  use types, only : i8b
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_cptr
  !This array is a f90 pointer
  real(rk),dimension(:,:),pointer,intent(inout) :: array
  integer(i8b),intent(in)                          :: numberOfElements1,numberOfElements2
  real(rk) :: test
  integer, allocatable :: arrayshape(:)

 !if MPI is defined, use the MPI sharedmemory things
#ifdef MPI

  shmem_buf_size = numberOfElements1*numberOfElements2*sizeof(test)
  call shmem_allocate_general_setup(shmem_buf_size,array_cptr)

  !remap that c_ptr to the array
  allocate(arrayshape(2))
  arrayshape = (/ numberOfElements1,numberOfElements2 /)
  call C_F_POINTER(array_cptr,array,arrayshape)
#else
  !otherwise, just allocate it like normal.
  allocate(array(numberOfElements1,numberOfElements2))
#endif
end subroutine shmem_allocate2_rk
!------------------------------------------------------------------------------

subroutine shmem_allocate2_rs_absH(array,numberOfElements1,numberOfElements2)
  !MPI_ADDRESS_KIND is required for the MPI_WIN_ALLOCATE_SHARED
  use types, only : i8b, rs_absH
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND)                 shmem_buf_size
#endif
  !C_PTR are pointers like in the C language: they are just a memory location
  type(C_PTR)                                 :: array_cptr
  !This array is a f90 pointer
  type(rs_absH),dimension(:,:),pointer,intent(inout) :: array
  !type(optional_integer),dimension(:,:),pointer,intent(inout) :: array
  integer(i8b),intent(in)                          :: numberOfElements1,numberOfElements2
  type(rs_absH) :: test
  integer, allocatable :: arrayshape(:)

 !if MPI is defined, use the MPI sharedmemory things
#ifdef MPI

  shmem_buf_size = numberOfElements1*numberOfElements2*sizeof(test)
  call shmem_allocate_general_setup(shmem_buf_size,array_cptr)

  !remap that c_ptr to the array
  allocate(arrayshape(2))
  arrayshape = (/ numberOfElements1,numberOfElements2 /)
  call C_F_POINTER(array_cptr,array,arrayshape)
#else
  !otherwise, just allocate it like normal.
  allocate(array(numberOfElements1,numberOfElements2))
#endif
end subroutine shmem_allocate2_rs_absH
!------------------------------------------------------------------------------

subroutine mpi_stop(msg)
  USE ISO_FORTRAN_ENV, only : error_unit ! access computing environment
  character(len=*) msg
  integer mpi_err

  write(6,'(a)') msg ; call flush(6)
  write(error_unit,'(a)') msg ! write message to standard error
#ifdef MPI
  call MPI_ABORT(MPI_COMM_WORLD,1,mpi_err)
#else
  stop 'mpi_stop'
#endif

end subroutine mpi_stop

end module mpi_routines
