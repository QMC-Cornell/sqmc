module more_tools
  use types, only: ik, ik_vec, rk
#ifdef NUM_ORBITALS_GT_127
  use overload
#endif
  use tools, only: merge_sort_real
  implicit none
  save
  private
  public ::  signmod, constrained_dets, fermionic_phase, get_nbr,check_nbr, is_in_list, is_in_list_mod,                        &
           & print_real_array_with_index,                                                                                      &
           & print_walker,print_walker2, print_real_matrix, print_real_matrix_lower_triangular, print_complex_matrix, print_int_matrix, print_map_on_square,       &
           & sherman_morrison, print_bin,get_adj_list_12_cross,get_adj_list_44_square,get_nbr_from_adj_list,                   &
           & choose_entry_by_its_weight,choose_entry_by_its_weight_rank2,convert_integer_to_nary,                              &
           & real_symmetric_diagonalize, real_symmetric_diagonalize_ow_ham, real_general_diagonalize_ow_ham, real_sym_gen_eig,matrix_lanczos,                    &
           & sparse_matrix_multiply,fast_sparse_matrix_multiply,fast_sparse_matrix_multiply_upper_triangular,binary_search,binary_search_single,binary_search_list_and_update,linear_search_list_and_update,linear_search_list,           &
           & overlap_two_slater,gp_two_slater,create_kspace_sym_maps,create_rspace_sym_maps,generate_fourfold_k_configs,       &
           & generate_fourfold_k_configs_efficient,generate_fourfold_k_configs_efficient_given_locations,print_sym_configs,    &
           & get_rep_only,linear_search,sparse_fast_sparse_matrix_multiply_upper_triangular,reduce_matrix_indices,             &
           & fast_sparse_matrix_multiply_local_band,from_upper_triangular_to_band,fast_sparse_matrix_multiply_local_band_block,&
           & from_upper_triangular_to_band_shuffle,from_upper_triangular_to_band_shuffle_old,sample_discrete_distribution,gen_hist,add_to_hist,get_occ_orbs,setup_alias,sample_alias,&
           & cyrus_diagonalize_sparse, davidson_sparse, davidson_sparse_single, davidson_sparse_mpi, parpack_diagonalize, &
           & davidson_sparse_mpi2

  public :: binary_search_lbound, binary_search_rbound
!=====================================================================================================================
 interface matrix_lanczos
   module procedure matrix_lanczos_square, matrix_lanczos_sparse
 end interface matrix_lanczos
!=====================================================================================================================

!=====================================================================================================================
 interface convert_integer_to_nary
   module procedure convert_integer_to_nary_int4, convert_integer_to_nary_int16
 end interface convert_integer_to_nary
!=====================================================================================================================

!=====================================================================================================================
 interface print_sym_configs
#ifdef NUM_ORBITALS_GT_127
   module procedure print_sym_configs_scalar, print_sym_configs_vec
#else
   module procedure print_sym_configs_scalar
#endif
 end interface print_sym_configs
!=====================================================================================================================

!=====================================================================================================================
 interface setup_alias
   module procedure setup_alias_rk, setup_alias_real
 end interface setup_alias
!=====================================================================================================================

!=====================================================================================================================
 interface sample_alias
   module procedure sample_alias_rk, sample_alias_real
 end interface sample_alias
!=====================================================================================================================

!=====================================================================================================================
 interface binary_search
   module procedure binary_search_int, binary_search_i8b
 end interface binary_search
!=====================================================================================================================

!=====================================================================================================================
 interface get_occ_orbs
   module procedure get_occ_orbs_1det, get_occ_orbs_2dets
 end interface get_occ_orbs
!=====================================================================================================================

! Note: it is NOT a good idea to increase epsilon much, because the number of HCI iterations goes up. 1e-6 is too big.
 real(rk)                 :: epsilon=1.e-10 ! for Lanczos, Davidson

contains
  !==========================================================================
  function signmod(a,b)
  ! -----------------------------------------------------------------------
  ! Description : Returns the mod of a number ... i.e. 6 mod 8 = 6
  !               and consider                        -6 mod 8 = 2
  !               So signmod will return min(6,2) =2
  ! Author      : H.J. Changlani
  ! Comments    : Old code,Not needed by new main code -get_nbr used instead
  ! ------------------------------------------------------------------------
  implicit none
  integer,intent(in) :: a
  integer,intent(in) :: b
  integer               tmp1,tmp2
  integer               signmod

  if (b.eq. 0 ) then
      signmod=0
      return
  endif

  tmp1=modulo(a,b)
  tmp2=modulo(-a,b)

  signmod=min(tmp1,tmp2)
  return

  end function signmod
  !==========================================================================

  !=================================================================================================
  subroutine constrained_dets(n_alpha,num_configs,dets)
  ! ----------------------------------------------------------------------------------------------
  ! Description   : Gives us a neat way of generating all configs with a fixed
  !                 particle number
  !                 Only used in subroutines make_hubbard_matrix_2d,
  !                 make_hubbard_hamiltonian_k, lanczos_hubbardlanczos_hubbard,
  !                 arnoldi_hubbard_binary_search, lanczos_hubbard_binary_search,
  !                 and density_matrix_2by2 in hubbard.f90, i.e., only for small
  !                 Hilbert spaces. -> No need for a vector_ik version yet.
  ! Author        : F. Petruzielo's code used by H.J. Changlani (ref. Bit Twiddling Hacks)
  ! ----------------------------------------------------------------------------------------------

  implicit none
  integer,intent(in)      ::n_alpha
  integer,intent(in)      ::num_configs
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(out) ::dets(num_configs)
  type(ik_vec)               temp
#else
  integer(ik),intent(out) ::dets(num_configs)
  integer(ik)               temp
#endif
  integer                   i

  dets(1) = 2_ik ** n_alpha - 1
  do i = 2, num_configs
     temp = ior(dets(i-1), dets(i-1) - 1) + 1
     dets(i) = ior(temp,ishft(iand(temp, -temp) / iand(dets(i-1),-dets(i-1)),-1) - 1)
  enddo

  end subroutine constrained_dets


  !===========================================================================
  function fermionic_phase(config,site_1,site_2)
  ! ------------------------------------------------------------------------
  ! Description   : Compute Jordan Wigner phase by computing number of occupied sites (of same spin)
  !                 between 2 positions
  ! Author        : H.J. Changlani, 21 Nov 2010
  ! ------------------------------------------------------------------------


  implicit none
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in)  ::config
#else
  integer(ik),intent(in)  ::config
#endif
  integer,intent(in)      ::site_1,site_2
  integer                 ::i,num_betn
  integer                 ::fermionic_phase

  num_betn=0

  !print *,"min(site_1,site_2)",min(site_1,site_2)
  !print *,"max(site_1,site_2)",max(site_1,site_2)

  do i=min(site_1,site_2)+1,max(site_1,site_2)-1
      !print *,"Checking fermion phase"
      if (btest(config,i-1) .eqv. .true. ) then
          num_betn=num_betn+1
      endif
  enddo

  !fermionic_phase=-1**(mod(num_betn,2))
  if (mod(num_betn,2) .eq. 0 ) then
    fermionic_phase=1
  else
    fermionic_phase=-1
  endif

  !print *,"mod(num_betn,2)",mod(num_betn,2)
  !print *,"fermionic_phase=",fermionic_phase

  end function fermionic_phase

  !===========================================================================
  !===========================================================================
  subroutine check_nbr(l_x,l_y,pbc,site_1,site_2,nbr_types,num_types)
  ! ------------------------------------------------------------------------
  ! Description   : Checks if site_1 and site_2 are nbrs
  !                 If yes, what kind of neighbor is site_2 wrt site_1
  !                 1=LEFT 2=RIGHT 3=UP 4=DOWN
  !                 If not a neighbor nbr_type = 0
  !
  ! Author        : H.J. Changlani, 21 Nov 2010
  ! ----------------------------------------------------------------------------------------------
      implicit none

      ! Dummy variables
      integer, intent(in)     ::  l_x,l_y,site_1,site_2
      integer, intent(out)    ::  nbr_types(4)
      integer,intent(out)     ::  num_types
      logical, intent(in)     ::  pbc

      ! Local
      integer                 :: j,nbr
      logical                 :: allowed

      allowed=.false.
      nbr_types(:)=0
      num_types=0

      do j=1,4
          call get_nbr(l_x,l_y,pbc,site_1,j,nbr,allowed)
          if (allowed) then
              if (nbr .eq. site_2) then
                  num_types=num_types+1
                  nbr_types(num_types)=j
              endif
          endif
      enddo
      return

  end subroutine check_nbr

  !===========================================================================
  subroutine get_nbr(l_x,l_y,pbc,site,nbr_type,nbr,allowed)
  ! ------------------------------------------------------------------------
  ! Description   : Returns neighbor of a site on a square lattice and tells if its allowed or not
  !                 i.e. on the boundary a left or right or Up or down might
  !                 not be allowed!
  ! Note          : We use integers to indicate nbr type
  !                 1=LEFT 2=RIGHT 3=UP 4=DOWN
  !
  ! Author        : H.J. Changlani, 21 Nov 2010
  ! ----------------------------------------------------------------------------------------------
      implicit none

      ! Dummy variables
      integer, intent(in)     ::  l_x,l_y,site,nbr_type
      integer, intent(out)    ::  nbr
      logical, intent(in)     ::  pbc
      logical, intent(out)    ::  allowed

      ! Local variables
      integer                 ::  x_1,x_2,y_1,y_2

      ! Convert site to x,y notation
      y_1=((site-1)/(l_x))+1
      x_1=site-((y_1-1)*l_x)
      allowed=.true.

      ! Initialize neighbor

      y_2=y_1
      x_2=x_1

      !------------------------------------------------------!
      !                       LEFT                           !
      !------------------------------------------------------!

      if (nbr_type .eq. 1) then
         x_2=x_1-1
         y_2=y_1

         if (pbc .eqv. .false.) then         ! OPEN BOUNDARY CONDITIONS
             if ((x_2 .gt. 0)) then
                  allowed=.true.
             else
                  allowed=.false.
             endif
         else                                ! PERIODIC BOUNDARY CONDITIONS
              if (x_2 .eq. 0) then
                  x_2=l_x
              endif
              if (x_2 .eq. x_1) then
                  allowed=.false.
              endif
         endif
      endif

      !------------------------------------------------------!
      !                       RIGHT                          !
      !------------------------------------------------------!

      if (nbr_type .eq. 2) then
         x_2=x_1+1
         y_2=y_1

         if (pbc .eqv. .false.) then         ! OPEN BOUNDARY CONDITIONS
             if ((x_2 .le. l_x)) then
                  allowed=.true.
             else
                  allowed=.false.
             endif
         else                                ! PERIODIC BOUNDARY CONDITIONS
              if (x_2 .eq. l_x+1) then
                  x_2=1
              endif
              if (x_2 .eq. x_1) then
                  allowed=.false.
              endif
         endif

      endif

      !------------------------------------------------------!
      !                       UP                             !
      !------------------------------------------------------!

      if (nbr_type .eq. 3) then
          x_2=x_1
          y_2=y_1+1
          if (pbc .eqv. .false.) then         ! OPEN BOUNDARY CONDITIONS
             if ((y_2 .le. l_y)) then
                  allowed=.true.
             else
                  allowed=.false.
             endif
          else                                ! PERIODIC BOUNDARY CONDITIONS
              if (y_2 .eq. l_y+1) then
                  y_2=1
              endif
              if (y_2 .eq. y_1) then
                  allowed=.false.
              endif
         endif
      endif

      !------------------------------------------------------!
      !                       DOWN                           !
      !------------------------------------------------------!

      if (nbr_type .eq. 4) then
          x_2=x_1
          y_2=y_1-1
          if (pbc .eqv. .false.) then         ! OPEN BOUNDARY CONDITIONS
             if ((y_2 .gt. 0)) then
                  allowed=.true.
             else
                  allowed=.false.
             endif
          else                                ! PERIODIC BOUNDARY CONDITIONS
              if (y_2 .eq. 0) then
                  y_2=l_y
              endif
              if (y_2 .eq. y_1) then
                  allowed=.false.
              endif
         endif
      endif

      if (allowed .eqv. .true.) then
          nbr=((y_2-1)*(l_x))+x_2
      else
          nbr=-1
      endif

  end subroutine get_nbr
! -------------------------------------------------------------------------------------------------
  subroutine get_nbr_from_adj_list(adj_list,site,nbr_type,nbr,allowed)
      implicit none

      ! Dummy variables
      integer, intent(in)     ::  site,nbr_type
      integer, intent(out)    ::  nbr
      logical, intent(out)    ::  allowed
      integer,intent(in)      ::  adj_list(:,:)

      allowed=.true.

      if (site .gt. size(adj_list,1)) then
          nbr=0
          write (6,*) "Adj list error. It should never get here"
          stop
      else
          nbr=adj_list(site,nbr_type)
      endif

      if (nbr .eq. 0 ) allowed=.false.

   end subroutine get_nbr_from_adj_list

!---------------------------------------------------------------------------------------------------------
  subroutine get_adj_list_62_strip(pbc,adj_list)
  ! ----------------------------------------------------------------------------------------------------
  ! Description   : Adjacency list of a 6x2 strip

  !                  5---6---1---2----9---10
  !                  |   |   | x |    |    |
  !                  7---8---3---4----11---12

  ! Note          : We use integers to indicate nbr type
  !                 1=LEFT 2=RIGHT 3=UP 4=DOWN
  !
  ! Author        : H.J. Changlani, Oct 7 2011
  ! ----------------------------------------------------------------------------------------------------
      implicit none

      ! Dummy variables
      logical, intent(in)     ::  pbc
      integer,intent(out)     ::  adj_list(12,4)

      adj_list(1,1)=6
      adj_list(1,2)=2
      adj_list(1,3)=3
      adj_list(1,4)=3

      adj_list(2,1)=1
      adj_list(2,2)=9
      adj_list(2,3)=4
      adj_list(2,4)=4

      adj_list(3,1)=8
      adj_list(3,2)=4
      adj_list(3,3)=1
      adj_list(3,4)=1

      adj_list(4,1)=3
      adj_list(4,2)=11
      adj_list(4,3)=2
      adj_list(4,4)=2

      adj_list(5,1)=10
      adj_list(5,2)=6
      adj_list(5,3)=7
      adj_list(5,4)=7

      adj_list(6,1)=5
      adj_list(6,2)=1
      adj_list(6,3)=8
      adj_list(6,4)=8

      adj_list(7,1)=12
      adj_list(7,2)=8
      adj_list(7,3)=5
      adj_list(7,4)=5

      adj_list(8,1)=7
      adj_list(8,2)=3
      adj_list(8,3)=6
      adj_list(8,4)=6

      adj_list(9,1)=2
      adj_list(9,2)=10
      adj_list(9,3)=11
      adj_list(9,4)=11

      adj_list(10,1)=9
      adj_list(10,2)=5
      adj_list(10,3)=12
      adj_list(10,4)=12

      adj_list(11,1)=4
      adj_list(11,2)=12
      adj_list(11,3)=9
      adj_list(11,4)=9

      adj_list(12,1)=11
      adj_list(12,2)=7
      adj_list(12,3)=10
      adj_list(12,4)=10

      if (pbc .eqv. .false.) then
          ! Left
          adj_list(5,1)=0
          adj_list(7,1)=0

          ! Right
          adj_list(10,2)=0
          adj_list(12,2)=0

          ! Up
          adj_list(5,3)=0
          adj_list(6,3)=0
          adj_list(1,3)=0
          adj_list(2,3)=0
          adj_list(9,3)=0
          adj_list(10,3)=0

          ! Down
          adj_list(7,4)=0
          adj_list(8,4)=0
          adj_list(3,4)=0
          adj_list(4,4)=0
          adj_list(11,4)=0
          adj_list(12,4)=0

      endif

      end subroutine get_adj_list_62_strip

  subroutine get_adj_list_44_square(pbc,adj_list)
  ! ----------------------------------------------------------------------------------------------------
  ! Description   : Adjacency list of a 4x4 square lattice
  !                     1---2----5---6
  !                     | x |    |   |
  !                     3---4----7---8
  !                     |   |    |   |
  !                     9---10---13--14
  !                     |   |    |   |
  !                     11--12---15--16
  ! Note          : We use integers to indicate nbr type
  !                 1=LEFT 2=RIGHT 3=UP 4=DOWN
  !
  ! Author        : H.J. Changlani, Oct 7 2011
  ! ----------------------------------------------------------------------------------------------------
      implicit none

      ! Dummy variables
      logical, intent(in)     ::  pbc
      integer,intent(out)     ::  adj_list(16,4)

      adj_list(1,1)=6
      adj_list(1,2)=2
      adj_list(1,3)=11
      adj_list(1,4)=3

      adj_list(2,1)=1
      adj_list(2,2)=5
      adj_list(2,3)=12
      adj_list(2,4)=4

      adj_list(3,1)=8
      adj_list(3,2)=4
      adj_list(3,3)=1
      adj_list(3,4)=9

      adj_list(4,1)=3
      adj_list(4,2)=7
      adj_list(4,3)=2
      adj_list(4,4)=10

      adj_list(5,1)=2
      adj_list(5,2)=6
      adj_list(5,3)=15
      adj_list(5,4)=7

      adj_list(6,1)=5
      adj_list(6,2)=1
      adj_list(6,3)=16
      adj_list(6,4)=8

      adj_list(7,1)=4
      adj_list(7,2)=8
      adj_list(7,3)=5
      adj_list(7,4)=13

      adj_list(8,1)=7
      adj_list(8,2)=3
      adj_list(8,3)=6
      adj_list(8,4)=14

      adj_list(9,1)=14
      adj_list(9,2)=10
      adj_list(9,3)=3
      adj_list(9,4)=11

      adj_list(10,1)=9
      adj_list(10,2)=13
      adj_list(10,3)=4
      adj_list(10,4)=12

      adj_list(11,1)=16
      adj_list(11,2)=12
      adj_list(11,3)=9
      adj_list(11,4)=1

      adj_list(12,1)=11
      adj_list(12,2)=15
      adj_list(12,3)=10
      adj_list(12,4)=2

      adj_list(13,1)=10
      adj_list(13,2)=14
      adj_list(13,3)=7
      adj_list(13,4)=15

      adj_list(14,1)=13
      adj_list(14,2)=9
      adj_list(14,3)=8
      adj_list(14,4)=16

      adj_list(15,1)=12
      adj_list(15,2)=16
      adj_list(15,3)=13
      adj_list(15,4)=5

      adj_list(16,1)=15
      adj_list(16,2)=11
      adj_list(16,3)=14
      adj_list(16,4)=6

      if (pbc .eqv. .false.) then
          ! Left
          adj_list(1,1)=0
          adj_list(3,1)=0
          adj_list(9,1)=0
          adj_list(11,1)=0

          ! Right
          adj_list(6,2)=0
          adj_list(8,2)=0
          adj_list(14,2)=0
          adj_list(16,2)=0

          ! Up
          adj_list(1,3)=0
          adj_list(2,3)=0
          adj_list(5,3)=0
          adj_list(6,3)=0

          ! Down
          adj_list(11,4)=0
          adj_list(12,4)=0
          adj_list(15,4)=0
          adj_list(16,4)=0

      endif

      end subroutine get_adj_list_44_square

  !========================================================================================================
  subroutine get_adj_list_12_cross(pbc,adj_list)
  ! -----------------------------------------------------------------------------------------------------
  ! Description   : Returns neighbor of a site on the following lattice  and tells if its allowed or not
  !                         6----7
  !                         |    |
  !                     5---1----2---8
  !                     |   |    |   |
  !                     9---3----4---10
  !                         |    |
  !                         11---12
  ! Note          : We use integers to indicate nbr type
  !                 1=LEFT 2=RIGHT 3=UP 4=DOWN
  !
  ! Author        : H.J. Changlani, Oct 7 2011
  ! ---------------------------------------------------------------------------------------------------
      implicit none

      ! Dummy variables
      logical, intent(in)     ::  pbc
      integer,intent(out)     ::  adj_list(12,4)

      adj_list(1,1)=5
      adj_list(1,2)=2
      adj_list(1,3)=6
      adj_list(1,4)=3

      adj_list(2,1)=1
      adj_list(2,2)=8
      adj_list(2,3)=7
      adj_list(2,4)=4

      adj_list(3,1)=9
      adj_list(3,2)=4
      adj_list(3,3)=1
      adj_list(3,4)=11

      adj_list(4,1)=3
      adj_list(4,2)=10
      adj_list(4,3)=2
      adj_list(4,4)=12

      if (pbc) then
          adj_list(5,1)=8
      else
          adj_list(5,1)=0
      endif
      adj_list(5,2)=1
      if (pbc) then
          adj_list(5,3)=9
      else
          adj_list(5,3)=0
      endif
      adj_list(5,4)=9

      if (pbc) then
          adj_list(6,1)=7
      else
          adj_list(6,1)=0
      endif

      adj_list(6,2)=7
      if (pbc) then
         adj_list(6,3)=11
      else
         adj_list(6,3)=0
      endif
      adj_list(6,4)=1

      adj_list(7,1)=6
      if (pbc) then
          adj_list(7,2)=6
      else
          adj_list(7,2)=0
      endif
      if (pbc) then
          adj_list(7,3)=12
      else
          adj_list(7,3)=0
      endif
      adj_list(7,4)=2

      adj_list(8,1)=2
      if (pbc) then
          adj_list(8,2)=5
      else
          adj_list(8,2)=0
      endif
      if (pbc) then
          adj_list(8,3)=10
      else
          adj_list(8,3)=0
      endif
      adj_list(8,4)=10

      if (pbc) then
          adj_list(9,1)=10
      else
          adj_list(9,1)=0
      endif
      adj_list(9,2)=3
      adj_list(9,3)=5
      if (pbc) then
          adj_list(9,4)=5
      else
          adj_list(9,4)=0
      endif

      adj_list(10,1)=4
      if (pbc) then
          adj_list(10,2)=9
      else
          adj_list(10,2)=0
      endif
      adj_list(10,3)=8
      if (pbc) then
          adj_list(10,4)=8
      else
          adj_list(10,4)=0
      endif

      if (pbc) then
          adj_list(11,1)=12
      else
          adj_list(11,1)=0
      endif
      adj_list(11,2)=12
      adj_list(11,3)=3
      if (pbc) then
          adj_list(11,4)=6
      else
          adj_list(11,4)=0
      endif

      adj_list(12,1)=11
      if (pbc) then
          adj_list(12,2)=11
      else
          adj_list(12,2)=0
      endif
      adj_list(12,3)=4
      if (pbc) then
          adj_list(12,4)=7
      else
          adj_list(12,4)=0
      endif

      end subroutine get_adj_list_12_cross

  ! ===============================================================================================
  subroutine is_in_list(config,list_size,list,entry_present,location)
  ! ===============================================================================================
  ! ----------------------------------------------------------------------------------------------
  ! Description   : Checks if walker = config  present within
  !                 a list of dets
  !                 We return if the entry is present or not
  !                 and the location of the entry (default location = -1)
  ! Author        : H.J. Changlani
  ! ----------------------------------------------------------------------------------------------

  implicit none

  integer,intent(in)      ::list_size
  integer,intent(out)     ::location
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in)  ::config
  type(ik_vec),intent(in)  ::list(list_size)
#else
  integer(ik),intent(in)  ::config
  integer(ik),intent(in)  ::list(list_size)
#endif
  logical,intent(out)     ::entry_present
  integer                   i

  do i=1,list_size
      if (config .eq. list(i)) then
          entry_present=.true.
          location=i
          return
      endif
  enddo

  entry_present=.false.
  location=-1
  return

  end subroutine is_in_list
  ! ===============================================================================================

  subroutine is_in_list_mod(up_config, dn_config, up_list, dn_list, list_size,entry_present,location)
  ! ===============================================================================================
  ! ----------------------------------------------------------------------------------------------
  ! Description   : Checks if walker = up_config,down_config is present within
  !                 a list of walkers
  !                 We return if the entry is present or not
  !                 and the location of the entry (default location = -1)
  ! Author        : H.J. Changlani
  ! ----------------------------------------------------------------------------------------------

  implicit none

  integer,intent(in)      ::list_size
  integer,intent(out)     ::location
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in)  ::up_config
  type(ik_vec),intent(in)  ::dn_config
  type(ik_vec),intent(in)  ::up_list(list_size),dn_list(list_size)
#else
  integer(ik),intent(in)  ::up_config
  integer(ik),intent(in)  ::dn_config
  integer(ik),intent(in)  ::up_list(list_size),dn_list(list_size)
#endif
  logical,intent(out)     ::entry_present
  integer                   i

  do i=1,list_size
      if ((up_config .eq. up_list(i)) .and. (dn_config .eq. dn_list(i))) then
          entry_present=.true.
          location=i
          return
      endif
   enddo

  entry_present=.false.
  location=-1
  return

  end subroutine is_in_list_mod

  ! ===============================================================================================
  subroutine print_walker(l_x,l_y,up_config,dn_config)
  ! ===============================================================================================
  ! ----------------------------------------------------------------------------------------------
  ! Description   : Print walker in a b a b a b a b format
  !                 where a = site of alpha electron b = site of beta electron
  !                 Spin up would be                  a b =1 0
  !                 Spin down would be                a b =0 1
  !                 Double occupancy                  a b =1 1
  !                 Empty                             a b =0 0
  ! Author        : H.J. Changlani
  ! ----------------------------------------------------------------------------------------------


  implicit none
  integer,    intent(in)  ::l_x,l_y
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in)  ::up_config,dn_config
#else
  integer(ik),intent(in)  ::up_config,dn_config
#endif
  integer                 ::i
  integer                 ::config(2*l_x*l_y)

  do i=1,l_x*l_y
      if (btest(up_config,i-1) .eqv. .true.) then
          config(2*i-1)=1
      else
          config(2*i-1)=0
      endif

      if (btest(dn_config,i-1) .eqv. .true.) then
          config(2*i)=1
      else
          config(2*i)=0
      endif
  enddo

  do i=1,2*l_x*l_y
   write (*,'(i2)',advance="no") config(i)
  enddo
  print *,""

  end subroutine print_walker
  ! ===============================================================================================

  subroutine print_walker2(l_x,l_y,up_config,dn_config,full_config)
  ! ===============================================================================================
  ! ----------------------------------------------------------------------------------------------
  ! Description   : Print walker in a b a b a b a b format
  !                 where a = site of alpha electron b = site of beta electron
  !                 Spin up would be                  a b =1 0
  !                 Spin down would be                a b =0 1
  !                 Double occupancy                  a b =1 1
  !                 Empty                             a b =0 0
  ! Author        : H.J. Changlani
  ! ----------------------------------------------------------------------------------------------


  implicit none
  integer,    intent(in)  ::l_x,l_y
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in)  ::up_config,dn_config
  type(ik_vec), intent(inout):: full_config
#else
  integer(ik),intent(in)  ::up_config,dn_config
  integer(ik), intent(inout):: full_config
#endif
  integer                 ::i
  integer                 ::config(2*l_x*l_y)
  integer(ik)             :: base

  do i=1,l_x*l_y
      if (btest(up_config,i-1) .eqv. .true.) then
          config(i)=1
      else
          config(i)=0
      endif

      if (btest(dn_config,i-1) .eqv. .true.) then
          config(i+(l_x*l_y))=1
      else
          config(i+(l_x*l_y))=0
      endif
  enddo

  base=1_ik
  full_config=0_ik

  do i=2*l_x*l_y,1,-1
   full_config=full_config+(config(i)*base)
   base=base*2_ik
  enddo

  end subroutine print_walker2
  ! ===============================================================================================

  subroutine print_real_matrix(rows,cols,matrix,print_rows,print_cols)
  ! ===============================================================================================
  ! ----------------------------------------------------------------------------------------------
  ! Description   : Print matrix
  ! Author      : H.J. Changlani
  ! ----------------------------------------------------------------------------------------------
  implicit none
  integer,intent(in)  :: rows,cols
  real(rk),intent(in) :: matrix(rows,cols)
  integer,optional    :: print_rows,print_cols
  integer             :: i,j
  integer             :: rows_to_print,cols_to_print

  if (present (print_rows) ) then
      rows_to_print=print_rows
  else
      rows_to_print=rows
  endif

  if (present (print_cols)) then
      cols_to_print=print_cols
  else
      cols_to_print=cols
  endif

  !write(6,*)
  do i=1,rows_to_print
      write(6,'(1000f10.3)') (matrix(i,j),j=1,cols_to_print)
  enddo

  end subroutine print_real_matrix


  subroutine print_real_matrix_lower_triangular(rows,cols,matrix,print_rows,print_cols)
  ! ===============================================================================================
  ! ----------------------------------------------------------------------------------------------
  ! Description : Print lower triangular part of a matrix (assumed to be
  !               symmetric)
  ! Author      : A Holmes, 20 Dec 2016
  ! ----------------------------------------------------------------------------------------------
  implicit none
  integer,intent(in)  :: rows,cols
  real(rk),intent(in) :: matrix(rows,cols)
  integer,optional    :: print_rows,print_cols
  integer             :: i,j
  integer             :: rows_to_print,cols_to_print

  if (present (print_rows) ) then
      rows_to_print=print_rows
  else
      rows_to_print=rows
  endif

  if (present (print_cols)) then
      cols_to_print=print_cols
  else
      cols_to_print=cols
  endif

  do i=1,rows_to_print
           write(6,'(1000f13.6)') (matrix(i,j),j=1,min(i,cols_to_print))
  enddo

  end subroutine print_real_matrix_lower_triangular


  subroutine print_complex_matrix(rows,cols,matrix,print_rows,print_cols)
  ! ===============================================================================================
  ! ----------------------------------------------------------------------------------------------
  ! Description   : Print matrix
  ! Author      : H.J. Changlani
  ! ----------------------------------------------------------------------------------------------
  implicit none
  integer,intent(in)    :: rows,cols
  complex*16,intent(in) :: matrix(rows,cols)
  integer,optional      :: print_rows,print_cols
  integer               :: i,j
  integer               :: rows_to_print,cols_to_print

  if (present (print_rows) ) then
      rows_to_print=print_rows
  else
      rows_to_print=rows
  endif

  if (present (print_cols)) then
      cols_to_print=print_cols
  else
      cols_to_print=cols
  endif

  !write(6,*)
  write(6,*) "Real part"
  do i=1,rows_to_print
           write(6,'(1000f10.5)') (real(matrix(i,j)),j=1,cols_to_print)
  enddo

  write(6,*) "Imaginary part"
  do i=1,rows_to_print
           write(6,'(1000f10.5)') (aimag(matrix(i,j)),j=1,cols_to_print)
  enddo

  call flush(6)

  end subroutine print_complex_matrix


  subroutine print_real_array_with_index(nelts,matrix)
  ! ===============================================================================================
  ! ----------------------------------------------------------------------------------------------
  ! Description   : Print matrix
  ! Author      : H.J. Changlani
  ! ----------------------------------------------------------------------------------------------
  implicit none
  integer,intent(in)  :: nelts
  real(rk),intent(in) :: matrix(nelts)
  integer i

  !write(6,*)
  do i=1,nelts
           write(6,'(i5,f10.5)') i,matrix(i)
  enddo

  end subroutine print_real_array_with_index

! ===============================================================================================

  subroutine print_int_matrix(rows,cols,matrix)
  ! ===============================================================================================
  ! ----------------------------------------------------------------------------------------------
  ! Description   : Print matrix of integers
  ! Author      : A. Holmes, Jun 2011
  ! ----------------------------------------------------------------------------------------------
  implicit none
  integer,intent(in)  :: rows,cols
  integer,intent(in) :: matrix(rows,cols)
  integer i,j

  do i=1,rows
           write(6,'(100i5)') (matrix(i,j),j=1,cols)
  enddo
  write(6,*)

  end subroutine print_int_matrix


! ===============================================================================================

  subroutine sherman_morrison(n,rownum,oldrow,newrow,inv,newinv,ratio)
  ! ===============================================================================================
  ! ----------------------------------------------------------------------------------------------
  ! Description   : Sherman-Morrison formula for the case where only one row changes.
  !
  ! Created     : Adam Holmes, 27 May 2011
  ! ----------------------------------------------------------------------------------------------

  implicit none

  ! Dummy variables
  integer,intent(in) :: n,rownum
  real(rk),intent(out) :: ratio
  real(rk),intent(in) :: inv(n,n),oldrow(n),newrow(n)
  real(rk),intent(out) :: newinv(n,n)

  ! Local variables
  integer i,j
  real(rk),dimension(n) :: a,b
  real(rk),dimension(n,n) :: d

  a = inv(:,rownum) ! column vector is rownum'th column of inv
  b = matmul(newrow-oldrow,inv) ! row vector is change in row times inv

  ratio = dot_product(newrow,a) ! ratio of new det to old det

  do i=1,size(a)
    do j=1,size(b)
      d(i,j) = a(i)*b(j) ! matrix structure of change in inv is outer product of a and b
    enddo
  enddo

  newinv = inv - d/(1+dot_product(newrow-oldrow,a))

  end subroutine sherman_morrison


! ===============================================================================================

  subroutine print_bin(rows,cols,det)
  ! ===============================================================================================
  ! ----------------------------------------------------------------------------------------------
  ! Description   : Converts bitpacked determinant into binary and displays.
  !
  ! Created     : Adam Holmes, 26 Jun 2011
  ! ----------------------------------------------------------------------------------------------

  implicit none

  ! Dummy variables
  integer,intent(in) :: rows,cols
#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: det(:)
#else
  integer(ik),intent(in) :: det(:)
#endif

  ! Local variables
  integer :: i,j
  integer,allocatable :: disp(:,:)

  allocate(disp(rows,cols))

  disp(:,:) = 0

  !if (rows==1) then
  !  do i=1,cols
  !    if (btest(det,j-1))  disp(1,j)=1
  !  enddo
  !else
  do i=1,rows
    do j=1,cols
      if (btest(det(i),j-1))  disp(i,j)=1
    enddo
  enddo
  !endif

  call print_int_matrix(rows,cols,disp)

  deallocate(disp)
  end subroutine print_bin

! ===============================================================================================

! subroutine list_combinations(n,k,nchoosek,list)
! ! ===============================================================================================
! ! ----------------------------------------------------------------------------------------------
! ! Description   : Lists out the n_choose_k combinations of integers from 1 to n.
! !               This algorithm first finds all permutations, so could be improved.
! !               Not presently being used
! !
! ! Created     : Adam Holmes, 6 Dec 2011
! ! ----------------------------------------------------------------------------------------------

! use tools, only: n_choose_k
! use types, only: i4b

! implicit none

! integer,intent(in) :: n,k
! integer,intent(out) :: nchoosek
! integer,allocatable,intent(out) :: list(:,:)

! integer :: i,j,ind,kfact,indmod
! integer :: old_perm(n),new_perm(n)

! if (k==1) then
!   nchoosek = n
!   allocate(list(n,1))
!   do i=1,n
!   list(i,1) = i
!   enddo
!   return
! endif

! kfact = 1
! do i=1,n-k
!   kfact = kfact*i
! enddo

! nchoosek = int(n_choose_k(n,k),i4b)
! allocate(list(nchoosek,k))
! list(:,:) = 0

! do i=1,n
!   old_perm(i) = i
!   new_perm(i) = i
! enddo

! do i=1,k
!   list(1,i) = i
! enddo

! ind = 1
! indmod = 1

! do while (.true.)
!   i=n-1
!   do while (new_perm(i) .ge. new_perm(i+1))
!   i=i-1
!   enddo
!   j=n
!   do while (new_perm(j) .le. new_perm(i))
!   j=j-1
!   enddo
!   new_perm(i) = old_perm(j)
!   new_perm(j) = old_perm(i)
!   old_perm(i) = new_perm(i)
!   old_perm(j) = new_perm(j)
!   i=i+1
!   j=n
!   do while (i<j)
!   new_perm(i) = old_perm(j)
!   new_perm(j) = old_perm(i)
!   i=i+1
!   j=j-1
!   enddo
!   do i=1,k-1
!   if (new_perm(i)>new_perm(i+1))  goto 10
!   enddo
!   ind=ind+1
!   if (mod(ind,kfact)==1) then
!   indmod=indmod+1
!   do i=1,k
!     list(indmod,i) = new_perm(i)
!   enddo
!   if (indmod==nchoosek)  exit
!   endif
!0  do j=1,n
!   old_perm(j) = new_perm(j)
!   enddo
! enddo

! end subroutine list_combinations

!=====================================================================================================================
 subroutine choose_entry_by_its_weight(entries,weights,number_of_entries,output_entry,out_wt,success)
!=====================================================================================================================
 ! Description : Given a set of determinants (integers) and corresponding real positive weights for those entries
 !             this subroutine normalizes those weights and selects the entry with a probability proportional to its
 !             weight. I envision that such a code could be more generally used for heat bath algorithms.
 ! Created by  : Hitesh J. Changlani
 ! Date      : Feb 1, 2012

 implicit none
 ! Dummy
 integer,intent(in)    :: number_of_entries
 integer(ik),intent(in) :: entries(number_of_entries)
 real(rk),intent(inout)  :: weights(number_of_entries)
 integer(ik),intent(out) :: output_entry
 real(rk),intent(out)  :: out_wt
 logical,intent(out)   :: success

 ! Local
 real(rk)              :: rannyu,rnd
 real(rk)              :: cwt,total_weight
 integer               :: i,j

 total_weight=0._rk
 rnd=rannyu()
 i=1
 cwt=0._rk
 output_entry=0_ik
 out_wt=0._rk
 success=.true.

 ! Normalize weights
 do j=1,number_of_entries
  total_weight=total_weight+weights(j)
 enddo

 if (total_weight .le. 1.0e-7) then
  success=.false.
  return
 endif

 weights(:)=weights(:)/total_weight

 do
  if (rnd>cwt .and. rnd <cwt+weights(i)) then
      output_entry=entries(i)
      out_wt=weights(i)
      return
  else
      cwt=cwt+weights(i)
      i=i+1
  endif
 enddo

 if (output_entry .eq. 0_ik) then
  write (6,*) "Bug in choose_entry_by_its_weight"
  call flush(6)
  stop
 endif
 end subroutine choose_entry_by_its_weight

!=====================================================================================================================
 subroutine choose_entry_by_its_weight_rank2(weight_table,rows,cols,out_row,out_col,out_wt,success)
!=====================================================================================================================
 ! Description : Given a set of determinants (integers) and corresponding real positive weights for those entries
 !             this subroutine normalizes those weights and selects the entry with a probability proportional to its
 !             weight. I envision that such a code could be more generally used for heat bath algorithms.
 ! Created by  : Hitesh J. Changlani
 ! Date      : Feb 1, 2012

 implicit none
 ! Dummy
 integer,intent(in)    :: rows,cols
 real(rk),intent(in)   :: weight_table(rows,cols)
 integer(ik),intent(out) :: out_row,out_col
 real(rk),intent(out)  :: out_wt
 logical,intent(out)   :: success

 ! Local
 real(rk)              :: wt_table(rows,cols)
 real(rk)              :: rannyu,rnd
 real(rk)              :: cwt,total_weight,total_weight_inv
 integer               :: j,k

 total_weight=0._rk
 rnd=rannyu()
 cwt=0._rk
 out_row=0
 out_col=0
 out_wt=0._rk
 success=.true.

 !write (6,*) "Weight table" !write (6,*) weight_table !stop
 ! Normalize weights
 do j=1,rows
  do k=1,cols
      wt_table(j,k)=abs(weight_table(j,k))
      total_weight=total_weight+wt_table(j,k)
  enddo
 enddo

 if (total_weight .le. 1.0e-7) then
  !write (6,*) "Total weight can never be so small"
  !call print_real_matrix(rows,cols,weight_table)
  !call print_real_matrix(rows,cols,wt_table)
  !write (6,*) "rows=",rows,"cols=",cols
  !call flush(6)
  success=.false.
  return
  !stop
 endif
 total_weight_inv=1._rk/(total_weight)
 wt_table(:,:)=wt_table(:,:)*total_weight_inv
 !print *,"wt_table (after norm) ",wt_table
 !print *,"total_weight=",total_weight

 do j=1,rows
  do k=1,cols
      if ((rnd>cwt) .and. (rnd .le. cwt+wt_table(j,k))) then
          out_wt=wt_table(j,k)
          out_row=j
          out_col=k
          return
      else
          cwt=cwt+wt_table(j,k)
      endif
  enddo
 enddo

  if ( (out_row .eq. 0) .or. (out_col .eq. 0)) then
  print *,"out_row=",out_row
  print *,"out_col=",out_col
  print *,"out_wt=",out_wt
  stop
  endif

 end subroutine choose_entry_by_its_weight_rank2


!=====================================================================================================================
 subroutine convert_integer_to_nary_int4(num,base,num_bits,nary)
!=====================================================================================================================
 ! Description : Given a integer and base (eg. 2) and how many digits desired, compute (bi)nary number
 ! Created by  : Hitesh J. Changlani
 ! Date      : Feb 8, 2012

 implicit none

 ! Dummy variables
 integer,intent(in)      :: num
 integer,intent(in)      :: base
 integer,intent(in)      :: num_bits
 integer,intent(out)     :: nary(num_bits)

 ! Local variables
 integer                 :: num_copy
 integer                 :: ctr

 num_copy=num
 nary(:)=0
 if (base**num_bits<num) then
      write (6,*) "Input number greater than representable by num_bits, Expect errors"
 endif

 ! Local variables
 do ctr=1,num_bits
  nary(ctr)=mod(num_copy,base)
  num_copy=num_copy/base
 enddo

 end subroutine convert_integer_to_nary_int4

!=====================================================================================================================
 subroutine convert_integer_to_nary_int16(num,base,num_bits,nary)
!=====================================================================================================================
 ! Description : Given a integer and base (eg. 2) and how many digits desired, compute (bi)nary number
 ! Created by  : Hitesh J. Changlani
 ! Date      : Feb 8, 2012

 implicit none

 ! Dummy variables
 integer(ik),intent(in)  :: num
 integer,intent(in)      :: base
 integer,intent(in)      :: num_bits
 integer(ik),intent(out)   :: nary(num_bits)

 ! Local variables
 integer(ik)             :: num_copy
 integer                 :: ctr

 num_copy=num
 nary(:)=0
 if (base**num_bits<num) then
      write (6,*) "Input number greater than representable by num_bits, Expect errors"
 endif

 ! Local variables
 do ctr=1,num_bits
  nary(ctr)=mod(num_copy,base)
  num_copy=num_copy/base
 enddo

 end subroutine convert_integer_to_nary_int16

!=====================================================================================================================
 subroutine real_symmetric_diagonalize(n,ham,eigenvectors,eigenvalues)
 ! Created by : Hitesh Changlani, March 12,2012
!=====================================================================================================================

 implicit none

 integer,intent(in)   :: n
 real(rk),intent(in)  :: ham(n,n)
 real(rk),intent(out) :: eigenvectors(n,n),eigenvalues(n)

 integer            :: len_work,info
 real(rk),allocatable :: work(:)

 call my_second(1,'real_symmetric_diagonalize')

 len_work = 3 * n -1
 allocate(work(len_work))

 eigenvectors=ham

 call dsyev('V', 'U', n, eigenvectors, n, eigenvalues, work, len_work, info)
 deallocate(work)

 if (info /= 0) then
      write(6,*) "info = ", info
      stop 'Diagonalization Failed!'
 endif

 call my_second(2,'real_symmetric_diagonalize')

 end subroutine real_symmetric_diagonalize

!=====================================================================================================================
 subroutine real_symmetric_diagonalize_ow_ham(n,ham,eigenvalues)
 ! Created by : Hitesh Changlani, March 12,2012
!=====================================================================================================================

 implicit none

 integer,intent(in)   :: n
 real(rk),intent(inout) :: ham(n,n)
 real(rk),intent(out)   :: eigenvalues(n)

 integer              :: len_work,info
 real(rk),allocatable   :: work(:)

 call my_second(1,'real_symmetric_diagonalize_ow_ham')

 len_work = 3 * n -1
 allocate(work(len_work))

 call dsyev('V', 'U', n, ham, n, eigenvalues, work, len_work, info)
 deallocate(work)

 if (info /= 0) then
      write(6,*) "info = ", info
      stop 'Diagonalization Failed!'
 endif

 call my_second(2,'real_symmetric_diagonalize_ow_ham')

 end subroutine real_symmetric_diagonalize_ow_ham

!=====================================================================================================================
 subroutine real_general_diagonalize_ow_ham(n,ham,re_eigenvalues,im_eigenvalues)
 ! Created by : Hitesh Changlani, October 17,2012
!=====================================================================================================================

 implicit none

 integer,intent(in)   :: n
 real(rk),intent(inout) :: ham(n,n)
 real(rk),intent(out)   :: re_eigenvalues(n)
 real(rk)             :: im_eigenvalues(n)
 real(rk)             :: vl(n,n),vr(n,n)
 integer              :: i,len_work,info
 real(rk)             :: work(8*n)
 integer              :: iorder(n)
 real(rk), dimension((n+1)/2) :: temp_real
 integer,  dimension((n+1)/2) :: temp_i_2   ! temporary array actually neither in nor out

 call my_second(1,'real_general_diagonalize_ow_ham')

 len_work = 8* n

 call dgeev('N', 'V', n, ham, n, re_eigenvalues, im_eigenvalues,vl,n,vr,n, work, len_work, info)
 ! Only the right eigenvectors are computed
 ! However unlike the symmetric case, the eigenvalues are not sorted

 if (info /= 0) then
      write(6,*) "info = ", info
      stop 'Diagonalization Failed!'
 endif

 do i=1,n
  iorder(i)=i
 enddo

 !write (6,*) "Before sort"
 !call flush(6)
 !do i=1,n
 !   write(6,*) re_eigenvalues(i),im_eigenvalues(i)
 !   call flush(6)
 !enddo

 call merge_sort_real(re_eigenvalues, iorder, n, temp_real, temp_i_2)

 im_eigenvalues(1:n)=im_eigenvalues(iorder(1:n))
 do i=1,n
  ham(:,i)=vr(:,iorder(i))
 enddo

 !write (6,*) "After sort"
 !call flush(6)
 !do i=1,n
 !   write(6,*) re_eigenvalues(i),im_eigenvalues(i)
 !   call flush(6)
 !enddo

 call my_second(2,'real_general_diagonalize_ow_ham')

 end subroutine real_general_diagonalize_ow_ham
!=====================================================================================================================

!=====================================================================================================================
 subroutine matrix_lanczos_square(n,lowest_eigenvector,lowest_eigenvalue,ham,highest_eigenvalue,second_lowest_eigenvalue,initial_vector)
 ! Created by  : Hitesh Changlani, March 2012
 ! Modified by : A Holmes, 14 Feb 2013. Optional input initial_vector is the starting vector for Lanczos.
!=====================================================================================================================
 implicit none

 ! Dummy
 integer,intent(in)        :: n
 real(rk),intent(in)       :: ham(n,n)
 real(rk),intent(out)      :: lowest_eigenvector(:)
 real(rk),intent(out)      :: lowest_eigenvalue
 real(rk),intent(out),optional :: highest_eigenvalue,second_lowest_eigenvalue
 real(rk),intent(in),optional :: initial_vector(:)

 ! Local
 integer                  :: i,it
!real(rk)                 :: rannyu
 real(rk)                 :: energy_shift
 real(rk)                 :: norm,norm_inv
 real(rk),allocatable     :: w(:),v(:,:)
 real(rk),allocatable     :: alphas(:),betas(:)
 integer                  :: iterations
 real(rk)                 :: lowest_eigenvalue_prev
 logical                  :: converged=.false.
 integer                  :: len_work,info
 real(rk),allocatable     :: work(:),eigenvalues(:),tridiag(:,:)

  call my_second(1,'matrix_lanczos_square')

  iterations=50        ! User option
  iterations=min(n,iterations)
  len_work = 3*n-1
  allocate(work(len_work))
  allocate (v(n,iterations+1))
  allocate (w(n))
  allocate(alphas(iterations+1))
  allocate(betas(iterations+1))
  w(:)=0_rk

  if (present(initial_vector)) then
  v(:,1)=initial_vector(:)
  else
  v(:,1)=0
  v(1,1)=1
  endif

  energy_shift=0._rk
  betas(1)=0._rk
  allocate (tridiag(iterations,iterations))
  converged=.false.

  write(6,'(/,''Executing matrix_lanczos_square in more_tools.f90'')')
  write (6,*) "n,iterations=",n,iterations;call flush(6)
  if (n>1) then
    do it=1,iterations
       w(:)=matmul(ham,v(:,it))
       if (it .gt. 1) w(:)=w(:)-betas(it)*v(:,it-1)
       alphas(it)=dot_product(w,v(:,it))
       w(:)=w(:)-alphas(it)*v(:,it)
       norm=dot_product(w,w)
       if (norm<(1.e-12_rk))  converged=.true.
       betas(it+1)=norm**(0.5_rk)
       norm_inv=1._rk/betas(it+1)
       v(:,it+1)=w(:)*norm_inv
       w(:)=v(:,it+1)
       do i=1,it        ! Reorthogonalization
           norm=dot_product(v(:,it+1),v(:,i))
           !write (6,*) "q=",norm
           call flush(6)
           w(:)=w(:)-norm*v(:,i)
       enddo
       v(:,it+1)=w(:)
       w(:)=0._rk
       norm=dot_product(v(:,it+1),v(:,it+1))
       norm_inv=1._rk/(norm**(0.5_rk))
       v(:,it+1)=v(:,it+1)*norm_inv
      !if (allocated(tridiag)) deallocate(tridiag)
      !allocate(tridiag(it,it))
       tridiag(:,:)=0._rk

       if (allocated(eigenvalues)) deallocate(eigenvalues)
       allocate(eigenvalues(it))
       eigenvalues(:)=0._rk

       do i=1,it
           tridiag(i,i)=alphas(i)
           if (i<it) then
               tridiag(i,i+1)=betas(i+1)
               tridiag(i+1,i)=betas(i+1)
           endif
       enddo

       !diagonalize with lapack routine
       len_work = 3 * (it) -1
       if (allocated(work)) deallocate(work)
       allocate(work(len_work))

       !write (6,*) "it-npower",it-npower
       call dsyev('V', 'U', it, tridiag(1:it,1:it), it, eigenvalues, work, len_work, info)

       !deallocate(tridiag) !deallocate(work)
       lowest_eigenvalue=eigenvalues(1)
       if (present(highest_eigenvalue))  highest_eigenvalue=eigenvalues(it)
       if (present(second_lowest_eigenvalue).and.it>1)  second_lowest_eigenvalue=eigenvalues(2)
       !call print_real_matrix(size(eigenvalues),1,eigenvalues)
       if (it.gt.1 .and. abs(lowest_eigenvalue-lowest_eigenvalue_prev)<epsilon) then
           converged=.true.
           exit
          !write(6,*) "(Final) Eigenvalues"
          !call print_real_matrix(size(eigenvalues,1),1,eigenvalues)
          !call flush(6)
       else
           lowest_eigenvalue_prev=lowest_eigenvalue
           write(6,'(''Iteration, Eigenvalue='',i3,f15.9)') it, lowest_eigenvalue
          !call flush(6)
          !write(6,*) "Eigenvalues"
          !call print_real_matrix(size(eigenvalues,1),1,eigenvalues)
          !call flush(6)
       endif
       if (converged)  exit
    enddo

    it=min(it,iterations)
    write(6,'(''matrix_lanczos_square: iter, n, Lowest eigenvalue ='',i4,i10,f17.10)') it, n, lowest_eigenvalue ; call flush(6)

    v(:,1)=matmul(v(:,1:it),tridiag(1:it,1))

    if (allocated(eigenvalues)) deallocate(eigenvalues)
    if (allocated(tridiag))         deallocate(tridiag)
    if (allocated(work))        deallocate(work)

  else
    write (6,*) "Diagonalization attempted with n=1"
    lowest_eigenvalue = ham(1,1)
  endif

  lowest_eigenvector(:)=v(:,1)

  call my_second(2,'matrix_lanczos_square')

  end subroutine matrix_lanczos_square

!=====================================================================================================================
 subroutine matrix_lanczos_sparse(n,lowest_eigenvector,lowest_eigenvalue,matrix_indices,nelem_nonzero,matrix_values,highest_eigenvalue,second_lowest_eigenvalue,initial_vector)
 ! Created by  : Hitesh Changlani, March 2012
 ! Modified by : Adam Holmes, 5 April 2012
 !             Modified input to allow for sparse matrices as well.
 !             Also changed the initial guess to HF instead of random (should converge more quickly this way).
 ! Modified by : A Holmes, 14 Feb 2013. Optional input initial_vector is the starting vector for Lanczos.
!=====================================================================================================================
 use common_ham, only : hamiltonian_type
 use types, only : i8b

 implicit none

 ! Dummy
 integer,intent(in)        :: n
 integer(i8b),intent(in)   :: matrix_indices(:),nelem_nonzero(:)
 real(rk),intent(in)       :: matrix_values(:)
 real(rk),intent(out)      :: lowest_eigenvector(:)
 real(rk),intent(out)      :: lowest_eigenvalue
 real(rk),intent(out),optional :: highest_eigenvalue,second_lowest_eigenvalue
 real(rk),intent(in),optional :: initial_vector(:)

 ! Local
 integer                  :: i,it
!real(rk)                 :: rannyu
 real(rk)                 :: energy_shift
 real(rk)                 :: norm,norm_inv
 real(rk),allocatable     :: w(:),v(:,:)
 real(rk),allocatable     :: alphas(:),betas(:)
 integer                  :: iterations
 real(rk)                 :: lowest_eigenvalue_prev
 logical                  :: converged=.false.
 integer                  :: len_work,info
 real(rk),allocatable     :: work(:),eigenvalues(:),tridiag(:,:)

  call my_second(1,'matrix_lanczos_sparse')

  iterations=50        ! User option
  iterations=min(n,iterations)
  allocate (v(n,iterations+1))
  allocate (w(n))
  allocate(alphas(iterations+1))
  allocate(betas(iterations+1))
  w(:)=0._rk

  if (present(initial_vector)) then
  norm = 1._rk/sqrt(dot_product(initial_vector,initial_vector))
  v(:,1) = norm*initial_vector(:)
  else
  v(:,1)=0
  v(1,1)=1
  endif

  energy_shift=0._rk
  betas(1)=0._rk
  allocate (tridiag(iterations,iterations))

  allocate(eigenvalues(iterations))
  len_work = 3*iterations-1
  allocate(work(len_work))

  converged=.false.

  write(6,'(/,''Executing matrix_lanczos_sparse in more_tools.f90'')')
  if (n>1) then
    do it=1,iterations
      if (hamiltonian_type.eq.'hubbardk' .or. hamiltonian_type.eq.'chem' .or. hamiltonian_type.eq.'heg') then
         call fast_sparse_matrix_multiply_upper_triangular(n,matrix_indices,nelem_nonzero,matrix_values,v(:,it),w(:))
       else
         call fast_sparse_matrix_multiply(n,matrix_indices,nelem_nonzero,matrix_values,v(:,it),w(:))
       endif
       if (it .gt. 1) w(:)=w(:)-betas(it)*v(:,it-1)
       alphas(it)=dot_product(w,v(:,it))
       w(:)=w(:)-alphas(it)*v(:,it)
       norm=dot_product(w,w)
       if (norm<(1.e-12_rk))  converged=.true.
       betas(it+1)=norm**(0.5_rk)
       norm_inv=1._rk/betas(it+1)
       v(:,it+1)=w(:)*norm_inv
       w(:)=v(:,it+1)
       do i=1,it        ! Reorthogonalization
           norm=dot_product(v(:,it+1),v(:,i))
           call flush(6)
           w(:)=w(:)-norm*v(:,i)
       enddo
       v(:,it+1)=w(:)
       w(:)=0._rk
       norm=dot_product(v(:,it+1),v(:,it+1))
       norm_inv=1._rk/(norm**(0.5_rk))
       v(:,it+1)=v(:,it+1)*norm_inv
       tridiag(:,:)=0._rk

       eigenvalues(:)=0._rk

       do i=1,it
           tridiag(i,i)=alphas(i)
           if (i<it) then
               tridiag(i,i+1)=betas(i+1)
               tridiag(i+1,i)=betas(i+1)
           endif
       enddo

       !diagonalize with lapack routine
       len_work = 3*it-1

       call dsyev('V', 'U', it, tridiag(1:it,1:it), it, eigenvalues, work, len_work, info)

       !deallocate(tridiag) !deallocate(work)
       lowest_eigenvalue=eigenvalues(1)
       if (present(highest_eigenvalue))  highest_eigenvalue=eigenvalues(it)
       if (present(second_lowest_eigenvalue).and.it>1)  second_lowest_eigenvalue=eigenvalues(2)
       !call print_real_matrix(size(eigenvalues),1,eigenvalues)
       if (it.gt.1 .and. abs(lowest_eigenvalue-lowest_eigenvalue_prev)<epsilon) then
           converged=.true.
           exit
       else
           lowest_eigenvalue_prev=lowest_eigenvalue
           write(6,'(''Iteration, Eigenvalue='',i3,f16.9)') it, lowest_eigenvalue
           call flush(6)
       endif
       if (converged)  exit

    enddo

    it=min(it,iterations)
    write(6,'(''matrix_lanczos_sparse: iter, n, Lowest eigenvalue ='',i4,i10,f17.10)') it, n, lowest_eigenvalue ; call flush(6)

    v(:,1)=matmul(v(:,1:it),tridiag(1:it,1))

    if (allocated(eigenvalues)) deallocate(eigenvalues)
    if (allocated(tridiag))     deallocate(tridiag)
    if (allocated(work))        deallocate(work)

  else
    write (6,*) "Diagonalization attempted with n=1"
    lowest_eigenvalue = matrix_values(1)
  endif

  lowest_eigenvector(:)=v(:,1)

  call my_second(2,'matrix_lanczos_sparse')

  end subroutine matrix_lanczos_sparse

!=====================================================================================================================
 subroutine cyrus_diagonalize_sparse(n,lowest_eigenvector,lowest_eigenvalue,matrix_indices,nelem_nonzero,matrix_values,highest_eigenvalue,second_lowest_eigenvalue,initial_vector)
 ! Created by  : Cyrus Umrigar, November 2016
 ! Works very well for C_2 but not for Cr_2.  Best to use Davidson.
!=====================================================================================================================
 use common_ham, only : hamiltonian_type
 use types, only : i8b

 implicit none

 ! Dummy
 integer,intent(in)        :: n
 integer(i8b),intent(in)   :: matrix_indices(:),nelem_nonzero(:)
 real(rk),intent(in)       :: matrix_values(:)
 real(rk),intent(out)      :: lowest_eigenvector(:)
 real(rk),intent(out)      :: lowest_eigenvalue
 real(rk),intent(out),optional :: highest_eigenvalue,second_lowest_eigenvalue
 real(rk),intent(in),optional :: initial_vector(:)

 ! Local
 integer                  :: i, it, ind
 real(rk)                 :: norm, a
 real(rk),allocatable     :: w(:,:),v(:,:)
 integer                  :: iterations
 real(rk)                 :: lowest_eigenvalue_prev, PT_corr
 real(rk),allocatable     :: diag_elems(:)
 integer j

  call my_second(1,'cyrus_diagonalize_sparse')

  iterations=50        ! User option
  iterations=min(n,iterations)
  allocate (v(n,2))
  allocate (w(n,2))

! Store the diagonal elements of H.
! The highest_eigenvalue is an lower estimate since it is just the largest diagonal element
! The second_lowest_eigenvalue is a crude estimate
! The lowest_eigenvalue is the smallest diagonal element but will be corrected later.
! Here it is calculated just to get the crude estimate of the second_lowest eigenvalue
  allocate(diag_elems(n))
  lowest_eigenvalue=9.d99
  highest_eigenvalue=-9.d99
  ind = 1
  diag_elems(1) = matrix_values(1)
  do i=2,n
    ind = ind+nelem_nonzero(i-1)
    diag_elems(i) = matrix_values(ind)
    lowest_eigenvalue=min(lowest_eigenvalue,diag_elems(i))
    highest_eigenvalue=max(highest_eigenvalue,diag_elems(i))
  enddo
  second_lowest_eigenvalue=lowest_eigenvalue+(highest_eigenvalue-lowest_eigenvalue)/n
  write(6,'(''Estimated lowest, second_lowest and highest eigenvalue='',9f14.6)') lowest_eigenvalue, second_lowest_eigenvalue, highest_eigenvalue

! write(6,'(''diag_elems='',100f8.4)') diag_elems(1:n)

! For Cr2 the 1st determinant does not have the lowest diag. element, so if the input vector has a 1 on the 1st element change it to put a 1 on the one with the lowest diag elem.
  if (present(initial_vector) .and. abs(initial_vector(1)-1).gt.1.e-6) then ! If the initial vector is present and not one where the 1st elem is 1
    norm = 1._rk/sqrt(dot_product(initial_vector,initial_vector))
    v(:,2) = norm*initial_vector(:)
  else
    v(1:n,2)=0
    v(minloc(diag_elems(1:n)),2)=1
  endif
  write(6,'(''minloc='',9i9)') minloc(diag_elems(1:n))

  write(6,'(/,''Executing cyrus_diagonalize_sparse in more_tools.f90'')')
  write(6,'(''n='',i9)') n
  if(n>1) then
    do it=1,iterations
      if (hamiltonian_type.eq.'hubbardk' .or. hamiltonian_type.eq.'chem' .or. hamiltonian_type.eq.'heg') then
         call fast_sparse_matrix_multiply_upper_triangular(n,matrix_indices,nelem_nonzero,matrix_values,v(:,2),w(:,2))
       else
         call fast_sparse_matrix_multiply(n,matrix_indices,nelem_nonzero,matrix_values,v(:,2),w(:,2))
       endif
       lowest_eigenvalue=dot_product(v(1:n,2),w(1:n,2))/dot_product(v(1:n,2),v(1:n,2))
       write(6,'(''Iteration, Eigenvalue='',i3,f16.9)') it, lowest_eigenvalue ; call flush(6)
       write(6,'(''dot_prod='',9f14.6)') dot_product(v(1:n,2),w(1:n,2)),dot_product(v(1:n,2),v(1:n,2))

       PT_corr=0
       do i=1,n
         if(diag_elems(i)-lowest_eigenvalue.gt.epsilon) then
           PT_corr=PT_corr-(w(i,2)-lowest_eigenvalue*v(i,2))**2/(diag_elems(i)-lowest_eigenvalue)
         else
           write(6,'(''Warning: i, diag_elems(i)-lowest_eigenvalue'',i5,9es12.4)') i,diag_elems(i)-lowest_eigenvalue
         endif
       enddo
       write(6,'(''PT_corr, E_tot'',9f14.6)') PT_corr, lowest_eigenvalue+PT_corr

!      do i=1,n
!        w(i,2)=(diag_elems(i)*v(i,2)-w(i,2))/(diag_elems(i)-lowest_eigenvalue-PT_corr)
!      enddo
       w(1:n,2)=(diag_elems(1:n)*v(1:n,2)-w(1:n,2))/(diag_elems(1:n)-lowest_eigenvalue-PT_corr)
       w(:,2)=w(:,2)/norm2(w(:,2))
       if(it.eq.1) then
         v(1:n,1)=v(1:n,2)
         v(1:n,2)=0.9*v(1:n,2)+0.1*w(1:n,2)
       else
!        a=dot_product(residual(1:n,1),residual(1:n,1)-residual(1:n,2)) / dot_product(residual(1:n,1)-residual(1:n,2),residual(1:n,1)-residual(1:n,2))
         a=dot_product(w(1:n,1)-v(1:n,1), w(1:n,1)-v(1:n,1)-w(1:n,2)+v(1:n,2)) / dot_product(w(1:n,1)-v(1:n,1)-w(1:n,2)+v(1:n,2), w(1:n,1)-v(1:n,1)-w(1:n,2)+v(1:n,2))
         write(6,'(''a='',f8.4)') a
         if(a.lt.0.01) then
           a=0.01
           write(6,'(''a reset to 0.01'')')
         endif
         v(1:n,1)=v(1:n,2)
         v(1:n,2)=(1.-a)*w(1:n,1)+a*w(1:n,2)
       endif
       v(1:n,2)=v(1:n,2)/norm2(v(1:n,2))
       w(1:n,1)=w(1:n,2)
       if(it.gt.1 .and. abs(lowest_eigenvalue_prev-lowest_eigenvalue).lt.epsilon) then
         exit
       endif
       lowest_eigenvalue_prev=lowest_eigenvalue
    enddo

  else
    write (6,*) "Diagonalization attempted with n=1"
    lowest_eigenvalue = matrix_values(1)
  endif

  it=min(it,iterations)
  write(6,'(''cyrus_diagonalize: iter, n, Lowest eigenvalue ='',i4,i10,f17.10)') it, n, lowest_eigenvalue ; call flush(6)
  lowest_eigenvector(1:n)=v(1:n,2)

  deallocate(v)
  deallocate(w)

  call my_second(2,'cyrus_diagonalize_sparse')

  end subroutine cyrus_diagonalize_sparse

!=====================================================================================================================
 subroutine davidson_sparse(n,n_states,final_vector,lowest_eigenvalues,matrix_indices,nelem_nonzero,matrix_values,initial_vector)
 ! Diagonally pre-conditioned Davidson
 ! A Holmes, 17 Nov 2016
!=====================================================================================================================

 use common_ham, only : hamiltonian_type
 use types, only : i8b

 implicit none

 ! Dummy
 integer,intent(in)        :: n,n_states
 integer(i8b),intent(in)   :: matrix_indices(:),nelem_nonzero(:)
 real(rk),intent(in)       :: matrix_values(:)
 real(rk),intent(out)      :: final_vector(:,:)
 real(rk),intent(out)      :: lowest_eigenvalues(:)
 real(rk),intent(in),optional :: initial_vector(:,:)

 ! Local
 integer                  :: i,it
!real(rk)                 :: rannyu
 real(rk)                 :: energy_shift
 real(rk)                 :: norm,norm_inv
 real(rk),allocatable :: residual_norm(:)
 real(rk),allocatable     :: w(:,:),Hw(:,:),v(:,:),Hv(:,:)
 integer                  :: iterations
 real(rk),allocatable :: lowest_eigenvalues_prev(:)
 logical                  :: converged=.false.
 integer                  :: len_work,info
 real(rk),allocatable     :: work(:),eigenvalues(:),h_krylov(:,:),h_overwrite(:,:)
 real(rk),allocatable     :: diag_elems(:)
 integer :: ind,j
 integer :: niter
 integer :: n_diagonalize, it_circ

  ! v(:,1:n_states) are input vectors, v(:,i) for i>n_states are residuals
  ! w(:,1:n_states) are the best vectors so far

  call my_second(1,'davidson_sparse')

  iterations=50        ! User option
  iterations=min(n,iterations)
  allocate(lowest_eigenvalues_prev(n_states))
  allocate(v(n,n_states*iterations))
  allocate(Hv(n,n_states*iterations))
  allocate(w(n,n_states))
  allocate(Hw(n,n_states))
  allocate(residual_norm(n_states))

  if (present(initial_vector)) then
    do i=1,n_states
      norm = 1._rk/sqrt(dot_product(initial_vector(:,i),initial_vector(:,i)))
      v(:,i) = norm*initial_vector(:,i)
      if (i>1) then
        ! Orthogonalize
        do j=1,i-1
          norm=dot_product(v(:,i),v(:,j))
          v(:,i)=v(:,i)-norm*v(:,j)
        enddo
        ! Normalize
        norm=dot_product(v(:,i),v(:,i))
        norm_inv=1._rk/sqrt(norm)
        v(:,i)=v(:,i)*norm_inv
      endif
    enddo
  else
    ! Start with HF and random dets
    v(:,:)=0
    do i=1,n_states
      v(i,i)=1
    enddo
  endif

  energy_shift=0._rk
  allocate (h_krylov(n_states*iterations,n_states*iterations))
  allocate (h_overwrite(n_states*iterations,n_states*iterations))

  allocate(eigenvalues(n_states*iterations))
  len_work = 3*n_states*iterations-1
  allocate(work(len_work))

  converged=.false.

  ! w is the lowest energy vector so far
  if (n>1) then
    ! Get diagonal elements
    allocate(diag_elems(n))
    ind = 1
    diag_elems(1) = matrix_values(1)
    do i=2,n
      ind = ind+nelem_nonzero(i-1)
      diag_elems(i) = matrix_values(ind)
    enddo

    ! First iteration:
    do i=1,n_states
      if (hamiltonian_type.eq.'hubbardk' .or. hamiltonian_type.eq.'chem' .or. hamiltonian_type.eq.'heg') then
        call fast_sparse_matrix_multiply_upper_triangular(n,matrix_indices,nelem_nonzero,matrix_values,v(:,i),Hv(:,i))
      else
        call fast_sparse_matrix_multiply(n,matrix_indices,nelem_nonzero,matrix_values,v(:,i),Hv(:,i))
      endif
    enddo

    do i=1,n_states
      lowest_eigenvalues(i) = dot_product(v(:,i),Hv(:,i))
      h_krylov(i,i) = lowest_eigenvalues(i)
      do j=i+1,n_states
        h_krylov(i,j) = dot_product(v(:,i),Hv(:,j))
        h_krylov(j,i) = h_krylov(i,j)
      enddo
    enddo

    write(6,'(''Iteration, Eigenvalues='',i3,10f16.9)') 1, lowest_eigenvalues(:)

    w(:,:) = v(:,1:n_states)
    Hw(:,:) = Hv(:,1:n_states)

    residual_norm(:) = 1._rk ! so at least one iteration is done

    niter = min(n,n_states*(iterations)) !MJO Removed mysterious +1

    n_diagonalize = 1 ! just for printout

    do it=n_states+1,niter*10 !MJO TEMP for now, allow upto 10x as many iterations
       ! Compute residual for state i
       it_circ = mod(it-1,niter)+1
       if (it.gt.niter) then
         if (it_circ.eq.1) then
           ! MJO Put the current best eigenvectors into v[1:nstates]
           v(:,1:n_states) = w(:,1:n_states)
           Hv(:,1:n_states) = Hw(:,1:n_states)
           do i=1,n_states
             lowest_eigenvalues(i) = dot_product(v(:,i),Hv(:,i)) ! locally sized array
             h_krylov(i,i) = lowest_eigenvalues(i)
             ! Do local dot products
             do j=i+1,n_states
               h_krylov(i,j) = dot_product(v(:,i),Hv(:,j)) !locally sized array
             enddo
             ! Copy into other half of array
             do j=i+1,n_states
               h_krylov(j,i) = h_krylov(i,j)
             enddo
           enddo
           cycle
         endif
       endif

       i = mod(it_circ-1,n_states)+1
       v(:,it_circ) = (Hw(:,i) - lowest_eigenvalues(i)*w(:,i))/(lowest_eigenvalues(i) - diag_elems(:))
       do j=1,n
         if (abs(lowest_eigenvalues(i)-diag_elems(j))<1e-8_rk)  v(j,it_circ) = -1._rk  ! Since denominator could be 0
       enddo

       ! If residual small, converged
       residual_norm(i)=dot_product(v(:,it_circ),v(:,it_circ))
       if (sum(residual_norm)<(1.e-12_rk))  converged=.true.

       ! Orthogonalize
       do i=1,it_circ-1
         norm=dot_product(v(:,it_circ),v(:,i))
         v(:,it_circ)=v(:,it_circ)-norm*v(:,i)
       enddo

       ! Normalize
       norm=dot_product(v(:,it_circ),v(:,it_circ))
       norm_inv=1._rk/sqrt(norm)
       v(:,it_circ)=v(:,it_circ)*norm_inv

       ! Apply H once
       if (hamiltonian_type.eq.'hubbardk' .or. hamiltonian_type.eq.'chem' .or. hamiltonian_type.eq.'heg') then
         call fast_sparse_matrix_multiply_upper_triangular(n,matrix_indices,nelem_nonzero,matrix_values,v(:,it_circ),Hv(:,it_circ))
       else
         call fast_sparse_matrix_multiply(n,matrix_indices,nelem_nonzero,matrix_values,v(:,it_circ),Hv(:,it_circ))
       endif

       ! Construct Krylov matrix and diagonalize
       do i=1,it_circ
         h_krylov(i,it_circ) = dot_product(v(:,i),Hv(:,it_circ))
         h_krylov(it_circ,i) = h_krylov(i,it_circ)
       enddo

       ! Diagonalize with lapack routine after all new states added
       if (mod(it_circ,n_states)==0) then

         len_work = 3*it_circ-1
         h_overwrite(1:it_circ,1:it_circ) = h_krylov(1:it_circ,1:it_circ)
         call dsyev('V', 'U', it_circ, h_overwrite(1:it_circ,1:it_circ), it_circ, eigenvalues, work, len_work, info)

         lowest_eigenvalues(:)=eigenvalues(1:n_states)

         do i=1,n_states
           w(:,i)=matmul(v(:,1:it_circ),h_overwrite(1:it_circ,i))
           Hw(:,i)=matmul(Hv(:,1:it_circ),h_overwrite(1:it_circ,i))
         enddo

         if (it.gt.1 .and. maxval(abs(lowest_eigenvalues(:)-lowest_eigenvalues_prev(:)))<epsilon) then
             converged=.true.
             exit
         else
             lowest_eigenvalues_prev(:)=lowest_eigenvalues(:)
             n_diagonalize = n_diagonalize + 1
             write(6,'(''Iteration, Eigenvalues='',i3,10f16.9)') n_diagonalize, lowest_eigenvalues(:)
             call flush(6)
         endif
         if (converged)  exit

       endif ! Diagonalize

    enddo ! it

    write(6,'(''davidson_sparse: n, iter, Lowest eigenvalue ='',i10,i4,10f17.10)') n, it, lowest_eigenvalues(:)
    call flush(6)

    if (allocated(eigenvalues)) deallocate(eigenvalues)
    if (allocated(h_krylov))    deallocate(h_krylov)
    if (allocated(work))        deallocate(work)

  else
    write (6,*) "Diagonalization attempted with n=1"
    lowest_eigenvalues(1) = matrix_values(1)
  endif

  final_vector(:,:)=w(:,:)

  call my_second(2,'davidson_sparse')

  end subroutine davidson_sparse

!=====================================================================================================================
 subroutine davidson_sparse_mpi(remote_det_map,local_det_map,n_states,final_vector,lowest_eigenvalues,matrix_indices,nelem_nonzero,matrix_values,initial_vector)
 ! Diagonally pre-conditioned Davidson MPI version
 ! A Holmes, 17 Nov 2016 Base version
 ! M Otten,  5 Dec 2016 Added MPI support
!=====================================================================================================================

 use common_ham, only : hamiltonian_type
 use types, only : i8b
 use mpi_routines, only : whoami,ncores,mpierr,det_map,det_map_l,shmem_allocate,shmem_deallocate,master_core_node,between_nodes_comm,mpi_barr_in_node,mpi_barr
 implicit none
#ifdef MPI
 include 'mpif.h'
#endif


 ! Dummy
 integer,intent(in)        :: n_states
 integer(i8b),intent(in)   :: matrix_indices(:),nelem_nonzero(:)
 real(rk),intent(in)       :: matrix_values(:)
 real(rk),intent(out)      :: final_vector(:,:)
 real(rk),intent(out)      :: lowest_eigenvalues(:)
 real(rk),intent(in),optional :: initial_vector(:,:)
 type(det_map),intent(in)  :: remote_det_map
 type(det_map_l),intent(in)  :: local_det_map
 ! Local
 integer                  :: i,it
!real(rk)                 :: rannyu
 real(rk)                 :: energy_shift
 real(rk)                 :: norm,norm_inv
 real(rk),allocatable :: residual_norm(:)
 real(rk),allocatable     :: w(:,:),Hw(:,:),Hv_tmp(:)
 real(rk),dimension(:,:),pointer :: v,Hv
 integer                  :: iterations
 real(rk),allocatable :: lowest_eigenvalues_prev(:)
 logical                  :: converged=.false.
 integer                  :: len_work,info
 real(rk),allocatable     :: work(:),eigenvalues(:),h_krylov(:,:),h_overwrite(:,:)
 real(rk),allocatable     :: diag_elems(:)
 integer :: ind,j
 integer :: niter
 integer :: n_diagonalize

 !MJO Variables need for MPI calls
 integer  :: n_local,n_global
  ! v(:,1:n_states) are input vectors, v(:,i) for i>n_states are residuals
  ! w(:,1:n_states) are tlhe best vectors so far

  call my_second(1,'davidson_sparse_mpi')
  n_global = local_det_map%ndets_global
  n_local  = local_det_map%ndets

  iterations=50        ! User option
  iterations=min(n_global,iterations)
  allocate(lowest_eigenvalues_prev(n_states))
  !allocate(v(n_global,n_states*iterations))
  !allocate(Hv(n_global,n_states*iterations))
  call shmem_allocate(v,int(n_global,i8b),int(n_states*iterations,i8b))
  call shmem_allocate(Hv,int(n_global,i8b),int(n_states*iterations,i8b))
  allocate(Hv_tmp(n_global))
  allocate(w(n_global,n_states))
  allocate(Hw(n_global,n_states))
  allocate(residual_norm(n_states))

  if(master_core_node) then
     if (present(initial_vector)) then
        do i=1,n_states
           norm = 1._rk/sqrt(dot_product(initial_vector(:,i),initial_vector(:,i)))
           v(:,i) = norm*initial_vector(:,i)
           if (i>1) then
              ! Orthogonalize
              do j=1,i-1
                 norm=dot_product(v(:,i),v(:,j))
                 v(:,i)=v(:,i)-norm*v(:,j)
              enddo
              ! Normalize
              norm=dot_product(v(:,i),v(:,i))
              norm_inv=1._rk/sqrt(norm)
              v(:,i)=v(:,i)*norm_inv
           endif
        enddo
     else
        ! Start with HF and random dets
        v(:,:)=0
        do i=1,n_states
           v(i,i)=1
        enddo
     endif
  endif
  call mpi_barr_in_node()

  energy_shift=0._rk
  allocate (h_krylov(n_states*iterations,n_states*iterations))
  allocate (h_overwrite(n_states*iterations,n_states*iterations))

  allocate(eigenvalues(n_states*iterations))
  len_work = 3*n_states*iterations-1
  allocate(work(len_work))

  converged=.false.

  ! w is the lowest energy vector so far
  if (n_global>1) then
    ! Get diagonal elements
    allocate(diag_elems(n_global))
    diag_elems = 0
    diag_elems(local_det_map%local_indices(1)) = matrix_values(1)
    ind = 1
    do i=2,n_local
      ind = ind+nelem_nonzero(i-1)
      diag_elems(local_det_map%local_indices(i)) = matrix_values(ind)
    enddo

    if (master_core_node) then
       Hv = 0
    endif
    call mpi_barr_in_node()

    !call my_second(1,"mxm")
    ! First iteration:
    do i=1,n_states
       if (hamiltonian_type.eq.'hubbardk' .or. hamiltonian_type.eq.'chem') then
          !call fast_sparse_matrix_multiply_upper_triangular_mpi(n_local,matrix_indices,nelem_nonzero,matrix_values,v(:,i),Hv(:,i),local_det_map%local_indices)
          call fast_sparse_matrix_multiply_upper_triangular_mpi(n_local,matrix_indices,nelem_nonzero,matrix_values,v(:,i),Hv_tmp(:),local_det_map%local_indices)
       else
          call fast_sparse_matrix_multiply(n_local,matrix_indices,nelem_nonzero,matrix_values,v(:,i),Hv(:,i))
       endif
    enddo
    !MJO Reduce diag elems
#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,diag_elems,n_global,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
    !MJO collect all v - FIXME - do all n_states at once!
    do i=1,n_states

       ! MJO Reduce Hv_tmp
#ifdef MPI
       call MPI_ALLREDUCE(MPI_IN_PLACE,Hv_tmp(:),n_global,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
       !MJO Copy Hv_tmp into Hv
       if (master_core_node) then
          do j=1,n_global
             Hv(j,i) = Hv_tmp(j)
          enddo
       endif
    enddo
    call mpi_barr_in_node() !Make sure all communication has finished
    !call my_second(2,"mxm")

    do i=1,n_states
      lowest_eigenvalues(i) = dot_product(v(:,i),Hv(:,i))
      h_krylov(i,i) = lowest_eigenvalues(i)
      do j=i+1,n_states
        h_krylov(i,j) = dot_product(v(:,i),Hv(:,j))
        h_krylov(j,i) = h_krylov(i,j)
      enddo
    enddo

    write(6,'(''Iteration, Eigenvalues='',i3,10f16.9)') 1,lowest_eigenvalues(:)

    w(:,:) = v(:,1:n_states)
    Hw(:,:) = Hv(:,1:n_states)

    residual_norm(:) = 1._rk ! so at least one iteration is done

    niter = min(n_global,n_states*(iterations+1))

    n_diagonalize = 1 ! just for printout

    do it=n_states+1,niter
       ! Compute residual for state i
       i = mod(it-1,n_states)+1
       if (master_core_node) then
          v(:,it) = (Hw(:,i) - lowest_eigenvalues(i)*w(:,i))/(lowest_eigenvalues(i) - diag_elems(:))
          do j=1,n_global
             if (abs(lowest_eigenvalues(i)-diag_elems(j))<1e-8_rk)  v(j,it) = -1._rk  ! Since denominator could be 0
          enddo
       endif
       call mpi_barr_in_node()

       ! If residual small, converged
       residual_norm(i)=dot_product(v(:,it),v(:,it))
       if (sum(residual_norm)<(1.e-12_rk))  converged=.true.

       ! Orthogonalize
       if (master_core_node) then
          do i=1,it-1
             norm=dot_product(v(:,it),v(:,i))
             v(:,it)=v(:,it)-norm*v(:,i)
          enddo
          ! Normalize
          norm=dot_product(v(:,it),v(:,it))
          norm_inv=1._rk/sqrt(norm)
          v(:,it)=v(:,it)*norm_inv
       endif
       i = it ! Set i to the correct number for all cores
       call mpi_barr_in_node()

       ! Apply H once
!       call my_second(1,"mxm")

       if (hamiltonian_type.eq.'hubbardk' .or. hamiltonian_type.eq.'chem') then
          call fast_sparse_matrix_multiply_upper_triangular_mpi(n_local,matrix_indices,nelem_nonzero,matrix_values,v(:,i),&
               Hv_tmp(:),local_det_map%local_indices)
          !call fast_sparse_matrix_multiply_upper_triangular_mpi(n_local,matrix_indices,nelem_nonzero,matrix_values,v(:,i),&
          !     Hv(:,i),local_det_map%local_indices)
       else
         call fast_sparse_matrix_multiply(n_local,matrix_indices,nelem_nonzero,matrix_values,v(:,it),Hv(:,it))
       endif

       !MJO collect all v
#ifdef MPI
       call MPI_ALLREDUCE(MPI_IN_PLACE,Hv_tmp(:),local_det_map%ndets_global,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
       if (master_core_node) then
          do j=1,n_global
             Hv(j,i) = Hv_tmp(j)
          enddo
       endif
       call mpi_barr_in_node() !Make sure Hv is updated
!       call my_second(2,"mxm")


       ! Construct Krylov matrix and diagonalize
       do i=1,it
          h_krylov(i,it) = dot_product(v(:,i),Hv(:,it))
          h_krylov(it,i) = h_krylov(i,it)
       enddo

       ! Diagonalize with lapack routine after all new states added
       if (mod(it,n_states)==0) then

         len_work = 3*it-1
         h_overwrite(1:it,1:it) = h_krylov(1:it,1:it)
         call dsyev('V', 'U', it, h_overwrite(1:it,1:it), it, eigenvalues, work, len_work, info)

         lowest_eigenvalues(:)=eigenvalues(1:n_states)

         do i=1,n_states
           w(:,i)=matmul(v(:,1:it),h_overwrite(1:it,i))
           Hw(:,i)=matmul(Hv(:,1:it),h_overwrite(1:it,i))
         enddo

         if (it.gt.1 .and. maxval(abs(lowest_eigenvalues(:)-lowest_eigenvalues_prev(:)))<epsilon) then
             converged=.true.
             exit
         else
             lowest_eigenvalues_prev(:)=lowest_eigenvalues(:)
             n_diagonalize = n_diagonalize + 1
             write(6,'(''Iteration, Eigenvalues='',i3,10f16.9)') n_diagonalize, lowest_eigenvalues(:)
             call flush(6)
         endif
         if (converged)  exit

       endif ! Diagonalize

    enddo ! it

    write(6,'(''davidson_sparse: n, iter, Lowest eigenvalue ='',i10,i4,10f17.10)') n_global, it, lowest_eigenvalues(:)
    call flush(6)

    if (allocated(eigenvalues)) deallocate(eigenvalues)
    if (allocated(h_krylov))     deallocate(h_krylov)
    if (allocated(work))        deallocate(work)

  else
    write (6,*) "Diagonalization attempted with n=1"
    lowest_eigenvalues(1) = matrix_values(1)
  endif

  final_vector(:,:)=w(:,:)

  call shmem_deallocate(v)
  call shmem_deallocate(Hv)
  call my_second(2,'davidson_sparse_mpi')

end subroutine davidson_sparse_mpi

!=====================================================================================================================
subroutine davidson_sparse_mpi2(remote_det_map,local_det_map,n_states,final_vector,lowest_eigenvalues,matrix_indices,nelem_nonzero,matrix_values,initial_vector)
  ! Diagonally pre-conditioned Davidson MPI version
  ! A Holmes, 17 Nov 2016 Base version
  ! M Otten,  5 Dec 2016 Added MPI support
  !=====================================================================================================================

  use common_ham, only : hamiltonian_type
  use types, only : i8b
  use mpi_routines, only : whoami,ncores,mpierr,det_map,det_map_l,shmem_allocate,shmem_deallocate,master_core_node,between_nodes_comm,mpi_barr_in_node,mpi_barr,mpi_allred,shmem_reallocate
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif


  ! Dummy
  integer,intent(in)        :: n_states
  integer(i8b),intent(in)   :: matrix_indices(:),nelem_nonzero(:)
  real(rk),intent(in)       :: matrix_values(:)
  real(rk),intent(out)      :: final_vector(:,:)
  real(rk),intent(out)      :: lowest_eigenvalues(:)
  real(rk),intent(in)       :: initial_vector(:,:)
  type(det_map),intent(inout)  :: remote_det_map
  type(det_map_l),intent(in)  :: local_det_map
  ! Local
  integer                  :: i,it,it_circ
  !real(rk)                 :: rannyu
  real(rk)                 :: energy_shift
  real(rk)                 :: norm,norm_inv
  real(rk),allocatable :: residual_norm(:)
  real(rk),allocatable     :: w(:,:),Hw(:,:),Hv_tmp_global(:),v_tmp_global(:)
  real(rk),allocatable :: v(:,:),Hv(:,:)
  integer                  :: iterations,resid_norm_index
  real(rk),allocatable :: lowest_eigenvalues_prev(:)
  logical                  :: converged=.false.
  integer                  :: len_work,info
  real(rk),allocatable     :: work(:),eigenvalues(:),h_krylov(:,:),h_overwrite(:,:)
  real(rk),allocatable     :: diag_elems(:)
  integer :: ind,j,k
  integer :: niter
  integer :: n_diagonalize

  !MJO Variables need for MPI calls
  integer  :: n_local,n_global
  ! v(:,1:n_states) are input vectors, v(:,i) for i>n_states are residuals
  ! w(:,1:n_states) are the best vectors so far

  call my_second(1,'davidson_sparse_mpi2')
  n_global = local_det_map%ndets_global
  n_local  = local_det_map%ndets

  iterations=50       ! User option
  iterations=min(n_global,iterations)

  allocate(lowest_eigenvalues_prev(n_states))

  allocate(v(n_local,n_states*iterations))
  allocate(Hv(n_local,n_states*iterations))
  allocate(w(n_local,n_states))
  allocate(Hw(n_local,n_states))

  allocate(residual_norm(n_states))

  ! Tmp global storage for mxm
  allocate(Hv_tmp_global(n_global))
  allocate(v_tmp_global(n_global))


  do i=1,n_states
    norm = dot_product(initial_vector(:,i),initial_vector(:,i)) !initial vector is globally sized, slice out our part
    call mpi_allred(norm)
    norm = 1/sqrt(norm)
    do j=1,n_local
      v(j,i) = norm*initial_vector(j,i)
    enddo
    if (i>1) then
      ! Orthogonalize
      do j=1,i-1
        norm=dot_product(v(:,i),v(:,j)) !Locally sized array
        call mpi_allred(norm)
        v(:,i)=v(:,i)-norm*v(:,j) !Locally sized array
      enddo
      ! Normalize
      norm=dot_product(v(:,i),v(:,i)) !Locally sized array
      call mpi_allred(norm)
      norm_inv=1._rk/sqrt(norm)
      v(:,i)=v(:,i)*norm_inv !Locally sized array
    endif
  enddo

  energy_shift=0._rk
  allocate (h_krylov(n_states*iterations,n_states*iterations))
  allocate (h_overwrite(n_states*iterations,n_states*iterations))

  allocate(eigenvalues(n_states*iterations*10))
  len_work = 3*n_states*iterations-1
  allocate(work(len_work))

  converged=.false.

  ! w is the lowest energy vector so far
  if (n_global>1) then
    ! Get diagonal elements

    if (n_local.ge.1) then
      allocate(diag_elems(n_local))
      diag_elems(1) = matrix_values(1)
      ind = 1
      do i=2,n_local
        ind = ind+nelem_nonzero(i-1)
        diag_elems(i) = matrix_values(ind)
      enddo
    endif

    !call my_second(1,"mxm")
    ! First iteration:
    do i=1,n_states
      !Copy the vector into a tmp global vector and reduce
      v_tmp_global = 0
      do j=1,n_local
        v_tmp_global(local_det_map%local_indices(j)) =v(j,i)
      enddo
      call mpi_allred(v_tmp_global)

      if (hamiltonian_type.eq.'hubbardk' .or. hamiltonian_type.eq.'chem') then
        !call fast_sparse_matrix_multiply_upper_triangular_mpi(n_local,matrix_indices,nelem_nonzero,matrix_values,v(:,i),Hv(:,i),local_det_map%local_indices)
        call fast_sparse_matrix_multiply_upper_triangular_mpi(n_local,matrix_indices,nelem_nonzero,matrix_values,v_tmp_global(:),Hv_tmp_global(:),local_det_map%local_indices)
      else
        call fast_sparse_matrix_multiply(n_local,matrix_indices,nelem_nonzero,matrix_values,v(:,i),Hv(:,i))
      endif
      !MJO collect all v - FIXME - do all n_states at once!

      ! MJO Reduce Hv_tmp
      call mpi_allred(Hv_tmp_global)

      !MJO Copy Hv_tmp into Hv
      do j=1,n_local
        Hv(j,i) = Hv_tmp_global(local_det_map%local_indices(j))
      enddo
    enddo

    !call my_second(2,"mxm")

    !FIXME Do communication for all nstates at once!
    do i=1,n_states
      lowest_eigenvalues(i) = dot_product(v(:,i),Hv(:,i)) ! locally sized array
      call mpi_allred(lowest_eigenvalues(i))
      h_krylov(i,i) = lowest_eigenvalues(i)
      ! Do local dot products
      do j=i+1,n_states
        h_krylov(i,j) = dot_product(v(:,i),Hv(:,j)) !locally sized array
      enddo
      ! Collect all local dot products
      call mpi_allred(h_krylov(i,i+1:n_states))
#ifdef MPI
      !       call MPI_ALLREDUCE(MPI_IN_PLACE,h_krylov(i,i+1:n_states),(n_states-i),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
      ! Copy into other half of array
      do j=i+1,n_states
        h_krylov(j,i) = h_krylov(i,j)
      enddo
    enddo

    write(6,'(''Iteration, Eigenvalues='',i3,10f16.9)') 1,lowest_eigenvalues(:)

    w(:,:) = v(:,1:n_states) !Locally sized array
    Hw(:,:) = Hv(:,1:n_states) !Locally sized array

    residual_norm(:) = 1._rk ! so at least one iteration is done

    niter = min(n_global,n_states*(iterations)) !MJO Removed mysterious +1

    n_diagonalize = 1 ! just for printout

    do it=n_states+1,niter*10 !MJO TEMP for now, allow upto 10x as many iterations
      ! Compute residual for state i
      it_circ = mod(it-1,niter)+1

      if (it.gt.niter) then
        if (it_circ.eq.1) then
          ! MJO Put the current best eigenvectors into v[1:nstates]
          v(:,1:n_states) = w(:,1:n_states)
          Hv(:,1:n_states) = Hw(:,1:n_states)
          do i=1,n_states
            lowest_eigenvalues(i) = dot_product(v(:,i),Hv(:,i)) ! locally sized array
            call mpi_allred(lowest_eigenvalues(i))
            h_krylov(i,i) = lowest_eigenvalues(i)
            ! Do local dot products
            do j=i+1,n_states
              h_krylov(i,j) = dot_product(v(:,i),Hv(:,j)) !locally sized array
            enddo
            ! Collect all local dot products
            call mpi_allred(h_krylov(i,i+1:n_states))
#ifdef MPI
            !       call MPI_ALLREDUCE(MPI_IN_PLACE,h_krylov(i,i+1:n_states),(n_states-i),MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
          enddo
          ! Copy into other half of array
          do i=1,n_states
            do j=i+1,n_states
              h_krylov(j,i) = h_krylov(i,j)
            enddo
          enddo
          cycle
        elseif (it_circ<=n_states) then
          cycle
        endif
      endif
      i = mod(it_circ-1,n_states)+1

      !All are local arrays
      do j=1,n_local
        v(j,it_circ) = (Hw(j,i) - lowest_eigenvalues(i)*w(j,i))/(lowest_eigenvalues(i) - diag_elems(j))
        if (abs(lowest_eigenvalues(i)-diag_elems(j))<1e-8_rk)  v(j,it) = -1._rk  ! Since denominator could be 0
      enddo

      ! If residual small, converged
      residual_norm(i)=dot_product(v(:,it_circ),v(:,it_circ)) !Locally sized arrays
      resid_norm_index = i

      !MJO If we are larger than niter, we need to use the
      !    circular buffer, and we orthonormalize against
      !    every vector in the buffer except itself
      do i=1,it_circ-1
        norm=dot_product(v(:,it_circ),v(:,i)) !Locally sized array
        call mpi_allred(norm)
        v(:,it_circ)=v(:,it_circ)-norm*v(:,i) !Locally sized array
      enddo

      ! Normalize
      norm=dot_product(v(:,it_circ),v(:,it_circ)) !Locally sized array
      call mpi_allred(norm)
      call mpi_allred(residual_norm(resid_norm_index)) !Residual norm communication moved down to make one less sync point
      if (sum(residual_norm)<(1.e-12_rk))  converged=.true.

      norm_inv=1._rk/sqrt(norm)
      v(:,it_circ)=v(:,it_circ)*norm_inv !Locally sized array

      ! Apply H once
      !       call my_second(1,"mxm")

      if (hamiltonian_type.eq.'hubbardk' .or. hamiltonian_type.eq.'chem') then
        !Copy the vector into a tmp global vector and reduce
        v_tmp_global = 0
        do j=1,n_local
          v_tmp_global(local_det_map%local_indices(j)) = v(j,it_circ)
        enddo
        call mpi_allred(v_tmp_global)
        call fast_sparse_matrix_multiply_upper_triangular_mpi(n_local,matrix_indices,nelem_nonzero,matrix_values,&
          v_tmp_global(:),Hv_tmp_global(:),local_det_map%local_indices)
      else
        !MJO This is likely to not work
        call fast_sparse_matrix_multiply(n_local,matrix_indices,nelem_nonzero,matrix_values,v(:,it_circ),Hv(:,it_circ))
      endif
      !MJO collect all Hv
      call mpi_allred(Hv_tmp_global)

      !Store the result back in our local vector
      do j=1,n_local
        Hv(j,it_circ) = Hv_tmp_global(local_det_map%local_indices(j))
      enddo

      !       call my_second(2,"mxm")


      ! Construct Krylov mpatrix and diagonalize
      ! Do all local dot products first
      do i=1,it_circ
        h_krylov(i,it_circ) = dot_product(v(:,i),Hv(:,it_circ)) !locally sized arrays
      enddo
      call mpi_allred(h_krylov(1:it_circ,it_circ))
      do i=1,it_circ
        h_krylov(it_circ,i) = h_krylov(i,it_circ)
      enddo

      ! Diagonalize with lapack routine after all new states added
      if (mod(it_circ,n_states)==0) then

        len_work = 3*it_circ-1
        h_overwrite(1:it_circ,1:it_circ) = h_krylov(1:it_circ,1:it_circ)
        ! Everyone does the diagonalization redundantly, as it is very small
        call dsyev('V', 'U', it_circ, h_overwrite(1:it_circ,1:it_circ), it_circ, eigenvalues, work, len_work, info)

        lowest_eigenvalues(:)=eigenvalues(1:n_states)

        ! We only need the part of w that is local, and we have the part of v that is
        ! local, as well as all of h_overwrite.
        ! Therefore, this is a completely local operation - no communication needed
        do i=1,n_states
          w(:,i)=matmul(v(:,1:it_circ),h_overwrite(1:it_circ,i)) !Locally sized arrays
          Hw(:,i)=matmul(Hv(:,1:it_circ),h_overwrite(1:it_circ,i)) !Locally sized arrays
        enddo
        if (it.gt.1 .and. maxval(abs(lowest_eigenvalues(:)-lowest_eigenvalues_prev(:)))<epsilon) then
          converged=.true.
          exit
        else
          lowest_eigenvalues_prev(:)=lowest_eigenvalues(:)
          n_diagonalize = n_diagonalize + 1
          write(6,'(''Iteration, Eigenvalues='',i3,10f16.9)') n_diagonalize, lowest_eigenvalues(:)
          call flush(6)
        endif
        if (converged)  exit

      endif ! Diagonalize
!      if (it.gt.niter) stop
    enddo ! it

    write(6,'(''davidson_sparse: n, iter, Lowest eigenvalue ='',i10,i4,10f17.10)') n_global, it, lowest_eigenvalues(:)
    call flush(6)

  else
    write (6,*) "Diagonalization attempted with n=1"
    lowest_eigenvalues(1) = matrix_values(1)
  endif

  !Copy our data into the output vectors
  !FIXME Do all nstates at once!
  call shmem_reallocate(remote_det_map%remote_wts,int(n_global,i8b),int(n_states,i8b))
  if (master_core_node) then
    remote_det_map%remote_wts(:,:) = 0
  endif
  final_vector = 0
  do j=1,n_states
    ! Store the global weights in the remote_det_map
    Hv_tmp_global = 0
    do i=1,n_local
      final_vector(i,j) = w(i,j)
      ! Take local vector and put in a tmp global
      Hv_tmp_global(local_det_map%local_indices(i)) = w(i,j)
    enddo
    ! Reduce the weights
    call mpi_allred(Hv_tmp_global)
    if (master_core_node) then
      remote_det_map%remote_wts(:,j) = Hv_tmp_global(:)
    endif
    call mpi_barr_in_node()
  enddo

  call my_second(2,'davidson_sparse_mpi2')

end subroutine davidson_sparse_mpi2


!=====================================================================================================================
subroutine parpack_diagonalize(remote_det_map,local_det_map,n_states,final_vector,lowest_eigenvalues,matrix_indices,nelem_nonzero,matrix_values,initial_vector)
!=====================================================================================================================

  use common_ham, only : hamiltonian_type
  use types, only : i8b
  use mpi_routines, only : whoami,ncores,mpierr,det_map,det_map_l,shmem_allocate,shmem_deallocate,master_core_node,between_nodes_comm,mpi_barr_in_node,mpi_barr,master_core
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif

  ! Dummy
  integer,intent(in)        :: n_states
  integer(i8b),intent(in)   :: matrix_indices(:),nelem_nonzero(:)
  real(rk),intent(in)       :: matrix_values(:)
  real(rk),intent(out)      :: final_vector(:,:)
  real(rk),intent(out)      :: lowest_eigenvalues(:)
  real(rk),intent(in),optional :: initial_vector(:,:)
  type(det_map),intent(inout)  :: remote_det_map
  type(det_map_l),intent(in)  :: local_det_map
#ifdef PARPACK
  ! Local
  integer                  :: i
  !real(rk)                 :: rannyu
  real(rk),allocatable     :: resid(:),workl(:),lanczos_v(:,:),lanczos_d(:,:),workd(:)
  real(rk),allocatable     :: Hv_tmp(:),v_tmp(:)
  integer  :: n_local,n_global
  integer          ido, n, nev, ncv, lworkl, info, ierr, j, nx, nconv, maxitr, mode, ishfts
  real(rk) :: sigma,tol,norm,norm_inv
  character        bmat*1, which*2
  integer          iparam(11), ipntr(11)
  logical, allocatable :: select_evs(:)
  logical :: rvec


  ! v(:,1:n_states) are input vectors, v(:,i) for i>n_states are residuals
  ! w(:,1:n_states) are tlhe best vectors so far
  call my_second(1,'davidson_sparse_parpack')
  n_global = local_det_map%ndets_global
  n_local  = local_det_map%ndets
  allocate(Hv_tmp(n_global))
  allocate(v_tmp(n_global))
  if (present(initial_vector)) then
     !        do i=1,n_states
     norm = 1._rk/sqrt(dot_product(initial_vector(:,1),initial_vector(:,1)))
     v_tmp = norm*initial_vector(:,1)
  else
     ! Start with HF and random dets
     v_tmp = 0
     v_tmp(1) = 1
  endif

  ido = 0
  bmat = 'I' !Standard eigenvalue problem
  nev  = n_states   !number of eigenvalues
  which = 'LM' !Computer largest in magnitude eigenvalues
  tol   = epsilon  !Desired precision
  ncv   = 15   !Number of arnoldi vectors
  lworkl = ncv*(ncv+8) !Work space
  info = 1   !We will use an initial vector
  iparam(1) = 1 !use exact shifts
  iparam(3) = 300 !Max number of Arnoldi iterations
  iparam(7) = 1   !A*x=lambda*x, A is symmetric

  allocate(resid(n_local))
  allocate(workl(lworkl))
  allocate(lanczos_v(n_local,ncv))
  allocate(workd(3*n_local))
  do i=1,n_local
     resid(i) = v_tmp(local_det_map%local_indices(i))
  enddo

  do
     call pdsaupd ( MPI_COMM_WORLD, ido, bmat, n_local, which, nev, tol, resid,&
     &                 ncv, lanczos_v, n_local, iparam, ipntr, workd, workl,&
     &                 lworkl, info )

     !Do the mxm
     if (ido .eq. -1 .or. ido .eq. 1) then
        !Copy the vector into a tmp global vector and reduce
        v_tmp = 0
        do i=1,n_local
           v_tmp(local_det_map%local_indices(i)) = workd(ipntr(1)+i-1) 
        enddo
#ifdef MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE,v_tmp,n_global,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
        !Do the MXM
        Hv_tmp = 0
        if (hamiltonian_type.eq.'hubbardk' .or. hamiltonian_type.eq.'chem') then
           call fast_sparse_matrix_multiply_upper_triangular_mpi(n_local,matrix_indices,nelem_nonzero,matrix_values,v_tmp,&
                Hv_tmp,local_det_map%local_indices)
        else
!           call fast_sparse_matrix_multiply(n_local,matrix_indices,nelem_nonzero,matrix_values,v(:,it),Hv(:,it))
        endif
        !Gather the results
#ifdef MPI
        call MPI_ALLREDUCE(MPI_IN_PLACE,Hv_tmp,n_global,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
        !Store the result back in the work vector
        do i=1,n_local
           workd(ipntr(2)+i-1) = Hv_tmp(local_det_map%local_indices(i))
        enddo
     else
        !We are done or there was an error, exit
        exit
     endif
  enddo

  !Check for errors
  if (info.lt.0) then
     if(master_core) then
        write(6,*) "Error in pdsaupd!!",info; call flush(6)
        stop
     endif
  endif
  allocate(select_evs(ncv)) !Not used, but necessary to pass?
  rvec = .true.
  allocate(lanczos_d(ncv,2))
  sigma = 0 !not referenced
  !Get evalues and evecs
  call pdseupd ( MPI_COMM_WORLD, rvec, 'All', select_evs,&
  &        lanczos_d, lanczos_v, n_local, sigma,&
  &        bmat, n_local, which, nev, tol, resid, ncv, lanczos_v, n_local,&
  &        iparam, ipntr, workd, workl, lworkl, ierr )

  !Check for errors
  if (ierr.ne.0) then
     if(master_core) then
        write(6,*) "Error in pdseupd!!"
     endif
  endif

  nconv = iparam(5)

!   do j=1,nconv
!      !Computer the residual norm - just for checking?
!      !Copy the vector into a tmp global vector and reduce
!      v_tmp = 0
!      do i=1,n_local
!         v_tmp(local_det_map%local_indices(i)) = lanczos_v(i,j) !! Take local vector and put in a global
!      enddo
! #ifdef MPI
!      call MPI_ALLREDUCE(MPI_IN_PLACE,v_tmp,n_global,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
! #endif
!      !Do the MXM
!      if (hamiltonian_type.eq.'hubbardk' .or. hamiltonian_type.eq.'chem') then
!         call fast_sparse_matrix_multiply_upper_triangular_mpi(n_local,matrix_indices,nelem_nonzero,matrix_values,v_tmp,&
!              Hv_tmp,local_det_map%local_indices)
!      else
! !        call fast_sparse_matrix_multiply(n_local,matrix_indices,nelem_nonzero,matrix_values,v(:,it),Hv(:,it))
!      endif
!      !Gather the results
! #ifdef MPI
!      call MPI_ALLREDUCE(MPI_IN_PLACE,Hv_tmp,n_global,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
! #endif
!      call daxpy(n_local,-lanczos_d(j,1),lanczos_v(1,j),1,Hv_tmp(local_det_map%local_indices),1)
!      ! lanczos_d(j,2) = sqrt(dot_product(Hv_tmp(local_det_map%local_indices),Hv_tmp(local_det_map%local_indices)))
!      ! call MPI_ALL_REDUCE(MPI_IN_PLACE,lanczos_d(j,2),1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
!      lanczos_d(j,2) = pdnorm2(MPI_COMM_WORLD,n_local,Hv_tmp(local_det_map%local_indices),1)
!   enddo

!Copy our data into the output vectors
  do j=1,nconv
     lowest_eigenvalues(j) = lanczos_d(j,1)
     final_vector = 0
     do i=1,n_local
        final_vector(local_det_map%local_indices(i),j) = lanczos_v(i,j) !! Take local vector and put in a global
     enddo
#ifdef MPI
     call MPI_ALLREDUCE(MPI_IN_PLACE,final_vector(:,j),n_global,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
  enddo

  write(6,'(''parpack_diagonalize: n, Lowest eigenvalue ='',i10,&
       & 10f17.10)') n_global, lowest_eigenvalues(:)
  call flush(6)

  call my_second(2,'davidson_sparse_parpack')
#else
  write(6,*) "PARPACK not compiled, resorting to davidson_sparse_mpi2"
  call davidson_sparse_mpi2(remote_det_map,local_det_map,n_states,final_vector,lowest_eigenvalues,matrix_indices,nelem_nonzero,matrix_values,initial_vector)
#endif
end subroutine parpack_diagonalize


 subroutine davidson_sparse_single(n,lowest_eigenvector&
      &,lowest_eigenvalue,matrix_indices,nelem_nonzero,matrix_values&
      &,highest_eigenvalue,initial_vector)
 ! Diagonally pre-conditioned Davidson
 ! A Holmes, 17 Nov 2016
!=====================================================================================================================

 use common_ham, only : hamiltonian_type
 use types, only : i8b

 implicit none

 ! Dummy
 integer,intent(in)        :: n
 integer(i8b),intent(in)   :: matrix_indices(:),nelem_nonzero(:)
 real(rk),intent(in)       :: matrix_values(:)
 real(rk),intent(out)      :: lowest_eigenvector(:)
 real(rk),intent(out)      :: lowest_eigenvalue
 real(rk),intent(out),optional :: highest_eigenvalue
 real(rk),intent(in),optional :: initial_vector(:)

 ! Local
 integer                  :: i,it
 real(rk)                 :: energy_shift
 real(rk)                 :: norm,norm_inv
 real(rk),allocatable     :: w(:),Hw(:),v(:,:),Hv(:,:)
 integer                  :: iterations
 real(rk)                 :: lowest_eigenvalue_prev
 logical                  :: converged=.false.
 integer                  :: len_work,info
 real(rk),allocatable     :: work(:),eigenvalues(:),h_krylov(:,:)&
      &,h_overwrite(:,:)
 real(rk),allocatable     :: diag_elems(:)
 integer :: ind

  ! v(:,1) is input vector, v(:,i) for i>1 are residuals
  ! w(:) is best vector so far

  call my_second(1,'davidson_sparse_single')

  iterations=50        ! User option
  iterations=min(n,iterations)

  allocate(v(n,iterations))
  allocate(Hv(n,iterations))
  allocate(w(n))
  allocate(Hw(n))
  allocate(h_krylov(iterations,iterations))
  allocate(h_overwrite(iterations,iterations))
  allocate(eigenvalues(iterations))
  allocate(diag_elems(n))
  len_work = 3*iterations-1
  allocate(work(len_work))

  if (present(initial_vector)) then
  norm = 1._rk/sqrt(dot_product(initial_vector,initial_vector))
  v(:,1) = norm*initial_vector(:)
  else
  v(:,1)=0
  v(1,1)=1
  endif

  energy_shift=0._rk

  converged=.false.

  ! w is the lowest energy vector so far
  if (n>1) then
    ! Get diagonal elements
    ind = 1
    diag_elems(1) = matrix_values(1)
    highest_eigenvalue=diag_elems(1)
    do i=2,n
      ind = ind+nelem_nonzero(i-1)
      diag_elems(i) = matrix_values(ind)
      highest_eigenvalue=max(highest_eigenvalue,diag_elems(i))
    enddo

    ! First iteration:
      if (hamiltonian_type.eq.'hubbardk' .or. hamiltonian_type.eq.'chem' .or. hamiltonian_type.eq.'heg') then
      call fast_sparse_matrix_multiply_upper_triangular(n,matrix_indices,nelem_nonzero,matrix_values,v(:,1),Hv(:,1))
    else
      call fast_sparse_matrix_multiply(n,matrix_indices,nelem_nonzero&
           &,matrix_values,v(:,1),Hv(:,1))
    endif
    lowest_eigenvalue = dot_product(v(:,1),Hv(:,1))
    lowest_eigenvalue_prev = lowest_eigenvalue
    write(6,'(''Iteration, Eigenvalue='',i3,f16.9)') 1,&
         & lowest_eigenvalue
    w(:) = v(:,1)
    Hw(:) = Hv(:,1)
    h_krylov(1,1) = lowest_eigenvalue
    h_overwrite(1,1) = lowest_eigenvalue

    do it=2,iterations
       ! Compute residual
       v(:,it) = (Hw(:) - lowest_eigenvalue*w(:))/(lowest_eigenvalue &
            &- diag_elems(:))
       if (abs(lowest_eigenvalue-diag_elems(1))<1e-8_rk)  v(1,it) = &
            &-1._rk  ! Since denominator could be 0

       ! If residual small, converged
       norm=dot_product(v(:,it),v(:,it))
      !write (6,*) "Norm of residual=",norm
       if (norm<(1.e-12_rk))  converged=.true.
      !norm_inv=1._rk/sqrt(norm)
      !v(:,it)=v(:,it)*norm_inv

       ! Orthogonalize
       do i=1,it-1
         norm=dot_product(v(:,it),v(:,i))
         v(:,it)=v(:,it)-norm*v(:,i)
       enddo

       ! Normalize
       norm=dot_product(v(:,it),v(:,it))
       norm_inv=1._rk/sqrt(norm)
       v(:,it)=v(:,it)*norm_inv

       ! Apply H once
      if (hamiltonian_type.eq.'hubbardk' .or. hamiltonian_type.eq.'chem' .or. hamiltonian_type.eq.'heg') then
         call fast_sparse_matrix_multiply_upper_triangular(n,matrix_indices,nelem_nonzero,matrix_values,v(:,it),Hv(:,it))
       else
         call fast_sparse_matrix_multiply(n,matrix_indices&
              &,nelem_nonzero,matrix_values,v(:,it),Hv(:,it))
       endif

       ! Construct Krylov matrix and diagonalize
       do i=1,it
         h_krylov(i,it) = dot_product(v(:,i),Hv(:,it))
         h_krylov(it,i) = h_krylov(i,it)
       enddo
      !write (6,*) "Krylov matrix:"
      !call print_real_matrix(it,it,h_krylov(1:it,1:it))

       ! Diagonalize with lapack routine
       len_work = 3*it-1
       h_overwrite(1:it,1:it) = h_krylov(1:it,1:it)
       call dsyev('V', 'U', it, h_overwrite(1:it,1:it), it,&
            & eigenvalues, work, len_work, info)

       lowest_eigenvalue=eigenvalues(1)

       w(:)=matmul(v(:,1:it),h_overwrite(1:it,1))
       Hw(:)=matmul(Hv(:,1:it),h_overwrite(1:it,1))

       if (present(highest_eigenvalue))  highest_eigenvalue&
            &=max(highest_eigenvalue,eigenvalues(it))
       if (it.gt.1 .and. abs(lowest_eigenvalue&
            &-lowest_eigenvalue_prev)<epsilon) then
           converged=.true.
           exit
       else
           lowest_eigenvalue_prev=lowest_eigenvalue
           write(6,'(''Iteration, Eigenvalue='',i3,f16.9)') it,&
                & lowest_eigenvalue
           call flush(6)
       endif
       if (converged)  exit

    enddo

    it=min(it,iterations)
    write(6,'(''davidson_sparse: n, iter, Lowest eigenvalue ='',i10,i4,10f17.10)') n, it, lowest_eigenvalue ; call flush(6)

  else
    write (6,*) "Diagonalization attempted with n=1"
    lowest_eigenvalue = matrix_values(1)
  endif

  lowest_eigenvector(:)=w(:)

  deallocate (v,Hv,w,Hw,h_krylov,h_overwrite,eigenvalues,diag_elems&
       &,work)

  call my_second(2,'davidson_sparse_single')


  end subroutine davidson_sparse_single

!===================================================================
  !==================================================
  subroutine sparse_matrix_multiply(matrix_indices,matrix_values&
       &,vector,answer)
!===================================================================
    !==================================================
  ! Multiplies a sparse matrix by a (non-sparse) vector
  ! Returns a vector that is assumed to be the same size as the input
    !  vector
  ! Not presently being used
  ! A. Holmes, 24 Feb 2012

  integer,intent(in) :: matrix_indices(:,:)
  real(rk),intent(in) :: matrix_values(:)
  real(rk),intent(in) :: vector(:)
  real(rk),intent(out) :: answer(:)

  integer :: i

  answer(:) = 0._rk

  do i=1,size(matrix_values)
    answer(matrix_indices(i,1)) = answer(matrix_indices(i,1)) +&
         & matrix_values(i) * vector(matrix_indices(i,2))
  enddo

  end subroutine sparse_matrix_multiply


!===================================================================
  !==================================================
  subroutine fast_sparse_matrix_multiply(n,matrix_indices&
       &,nelem_nonzero,matrix_values,vector,answer)
!===================================================================
    !==================================================
  ! Multiplies a sparse matrix by a (non-sparse) vector
  ! Returns a vector that is assumed to be the same size as the input
    !  vector
  ! A. Holmes, 25 Feb 2012

  use types, only : i8b

  integer,intent(in) :: n ! size of vector
  integer(i8b),intent(in) :: matrix_indices(:) ! Vector of column
  !  numbers
  integer(i8b),intent(in) :: nelem_nonzero(:)
  real(rk),intent(in) :: matrix_values(:)
  real(rk),intent(in) :: vector(:)
  real(rk),intent(out) :: answer(:)

  integer :: i,j,k
  real(rk) :: x

  answer(:) = 0._rk

  k=0
  do i=1,n
    x=0
    do j=1,nelem_nonzero(i)
      k=k+1
      x = x + matrix_values(k)*vector(matrix_indices(k))
    enddo
    answer(i) = x
  enddo

  end subroutine fast_sparse_matrix_multiply

!===================================================================
  !==================================================
  subroutine from_upper_triangular_to_band(n,my_n,my_columns&
       &,matrix_indices,nelem_nonzero,matrix_values)!,my_rows)
!===================================================================
    !==================================================
  ! Takes as input an upper triangular matrix and returns the band
!  needed for local core
  ! A Roggero, Dec 2013

  integer, intent(in) :: n,my_n
  integer, intent(in) :: my_columns(my_n)
!  integer, optional,intent(out),allocatable :: my_rows(:)

  integer,allocatable ,intent(inout):: matrix_indices(:) ! Vector of column numbers
  integer,allocatable ,intent(inout):: nelem_nonzero(:)
  real(rk),allocatable,intent(inout):: matrix_values(:)

  real(rk) :: full_matrix(n,n)

  integer :: ii,jj,kk,indx
  integer :: istat,nnonzero

  !generate full matrix
  full_matrix(:,:)=0_rk
  kk = 0
  do ii=1,n
    do jj=1,nelem_nonzero(ii)
      kk = kk + 1
      indx = matrix_indices(kk)
      full_matrix(ii,indx) = matrix_values(kk)
      full_matrix(indx,ii) = matrix_values(kk)
    enddo
  enddo

 !deallocate(nelem_nonzero)
 !allocate(nelem_nonzero(n),stat=istat)
 !if(istat/=0) stop "failed to reallocate nelem_nonzero"
  !count number of nonzero elements needed
  nelem_nonzero(:) = 0
  do ii=1,n
    do jj=1,my_n
      if(full_matrix(ii,my_columns(jj)) /= 0_rk) nelem_nonzero(ii) = nelem_nonzero(ii) + 1
    enddo
  enddo

  nnonzero=sum(nelem_nonzero(1:n))
  deallocate(matrix_indices,matrix_values)
  allocate(matrix_indices(nnonzero),stat=istat)
  if(istat/=0) stop "failed to reallocate matrix_indices"
  allocate(matrix_values(nnonzero),stat=istat)
  if(istat/=0) stop "failed to reallocate matrix_values"

  !rebuild matrix as local band
  kk = 0
  do ii=1,n
    do jj=1,my_n
      if(full_matrix(ii,my_columns(jj)) /= 0_rk) then
        kk = kk + 1
        matrix_values(kk) = full_matrix(ii,my_columns(jj))
        matrix_indices(kk) = jj
      endif
    enddo
  enddo
  if(kk /= nnonzero) stop "problem in conversion from_upper_triangular_to_band" ! just for debug

end subroutine from_upper_triangular_to_band

!=====================================================================================================================
  subroutine from_upper_triangular_to_band_shuffle(iam,row_count,my_rows,shuffle_list,n,my_n,matrix_indices,nelem_nonzero,matrix_values)
!=====================================================================================================================
  ! Takes as input an upper triangular matrix and returns the band needed for local core
  ! Version of algorithm that only requires O(nonzero elements in deterministic projector) memory.
  ! A Holmes, 21 Mar 2014

  use types, only : i8b
  use generic_sort, only : quick_sort

  integer, intent(in) :: n,my_n,iam
  integer(i8b), intent(in) :: shuffle_list(:)
  integer(i8b), intent(inout) :: my_rows(:)
  integer(i8b), intent(inout) :: row_count(:)

  integer(i8b),allocatable ,intent(inout):: matrix_indices(:) ! Vector of column numbers
  integer(i8b),allocatable ,intent(inout):: nelem_nonzero(:)
  real(rk),allocatable,intent(inout):: matrix_values(:)

  integer :: ii
  integer(i8b) :: initial_stride
  integer(i8b) :: i,num_square_elems,num_tri_elems,j,row

  integer(i8b),allocatable :: square_indices(:),square_nelem_nonzero(:),tmp_indices(:),tmp_nelem_nonzero(:)
  real(rk),allocatable :: square_values(:),tmp_values(:)
  integer(i8b),allocatable :: nelem_nonzero_cum(:),shuf_nelem_nonzero_cum(:),shuffle_list_sparse(:),inverse_shuffle(:)

  initial_stride=sum(row_count(1:iam))
  do ii=1,my_n
    my_rows(ii)=initial_stride+ii
  enddo

  ! We need to store only (all rows, my_rows(1:my_n)) elements of the projector!

  ! Convert from upper triangular form to square form first
  allocate(square_nelem_nonzero(n))
  allocate(square_indices(2*size(matrix_indices)-n))
  allocate(square_values(2*size(matrix_values)-n))

  call from_upper_triangular_to_square(n,nelem_nonzero,matrix_indices,matrix_values, square_nelem_nonzero,square_indices,square_values,size(matrix_indices,1,i8b),num_square_elems)

  allocate(inverse_shuffle(size(shuffle_list)))
  do i=1,size(shuffle_list)
    inverse_shuffle(shuffle_list(i)) = i
  enddo

  ! Permute the rows and columns so that consecutive rows belong to the same core.
  allocate(nelem_nonzero_cum(n))
  allocate(shuf_nelem_nonzero_cum(n))
  nelem_nonzero_cum(1) = 0
  shuf_nelem_nonzero_cum(1) = 0
  do i=2,n
    nelem_nonzero_cum(i) = nelem_nonzero_cum(i-1) + square_nelem_nonzero(i-1)
    shuf_nelem_nonzero_cum(i) = shuf_nelem_nonzero_cum(i-1) + square_nelem_nonzero(shuffle_list(i-1)) !inverse_shuffle(i-1))
  enddo

  allocate(shuffle_list_sparse(num_square_elems))
  do i=1,n
    do j=1,square_nelem_nonzero(shuffle_list(i))
      shuffle_list_sparse(shuf_nelem_nonzero_cum(i)+j) = nelem_nonzero_cum(shuffle_list(i)) + j
    enddo
  enddo

  allocate(tmp_nelem_nonzero(n))
  allocate(tmp_indices(num_square_elems))
  allocate(tmp_values(num_square_elems))

  tmp_nelem_nonzero(1:n) = square_nelem_nonzero(shuffle_list(1:n))
  tmp_indices(1:num_square_elems) = square_indices(shuffle_list_sparse(1:num_square_elems))
  tmp_values(1:num_square_elems) = square_values(shuffle_list_sparse(1:num_square_elems))

  square_nelem_nonzero(:) = tmp_nelem_nonzero(:)
  square_indices(:) = tmp_indices(:)
  square_values(:) = tmp_values(:)

  ! Relabel the indices

  do i=1,num_square_elems
    tmp_indices(i) = inverse_shuffle(square_indices(i))
  enddo
  ! Now sort within each row
  do i=1,n
    ! sort shuf_nelem_nonzero_cum(i) + 1 through shuf_nelem_nonzero_cum(i) + square_nelem_nonzero(i)
    call quick_sort(shuf_nelem_nonzero_cum(i)+1,shuf_nelem_nonzero_cum(i)+square_nelem_nonzero(i),tmp_indices,tmp_values)
  enddo
  square_indices(:) = tmp_indices(:)
  square_values(:) = tmp_values(:)

  ! Count number of elems in band form
  num_tri_elems = 0
  row = 1
  do i=1,num_square_elems
    if (square_indices(i)>=my_rows(1).and.square_indices(i)<=my_rows(my_n))  num_tri_elems = num_tri_elems + 1
  enddo

  deallocate(matrix_indices,matrix_values)
  allocate(matrix_indices(num_tri_elems))
  allocate(matrix_values(num_tri_elems))

  call from_square_to_band(n, square_nelem_nonzero,square_indices,square_values,nelem_nonzero,matrix_indices,matrix_values,num_square_elems,num_tri_elems,my_rows(1),my_rows(my_n))

end subroutine from_upper_triangular_to_band_shuffle


!=====================================================================================================================
  subroutine from_upper_triangular_to_band_shuffle_old(iam,row_count,my_rows,shuffle_list,n,my_n,matrix_indices,nelem_nonzero,matrix_values)!,my_rows)
!=====================================================================================================================
  ! Takes as input an upper triangular matrix and returns the band needed for local core
  ! Not presently being used, but:
  ! I'm leaving this in here in case there is ever a bug in the above routine. Note that this routine requires O(n_dtm^2) memory
  ! A Holmes, 21 Mar 2014

  use types, only : i8b

  integer, intent(in) :: n,my_n,iam
  integer, intent(in) :: shuffle_list(:)
  integer(i8b), intent(inout) :: my_rows(:)
  integer(i8b), intent(inout) :: row_count(:)
!  integer, optional,intent(out),allocatable :: my_rows(:)

  integer,allocatable ,intent(inout):: matrix_indices(:) ! Vector of column numbers
  integer,allocatable ,intent(inout):: nelem_nonzero(:)
  real(rk),allocatable,intent(inout):: matrix_values(:)

  real(rk) :: full_matrix(n,n)

  integer :: ii,jj,kk,indx
  integer :: istat,nnonzero
  integer :: initial_stride

  !generate full matrix
  full_matrix(:,:)=0_rk
  kk = 0
  do ii=1,n
    do jj=1,nelem_nonzero(ii)
      kk = kk + 1
      indx = matrix_indices(kk)
      full_matrix(ii,indx) = matrix_values(kk)
      full_matrix(indx,ii) = matrix_values(kk)
    enddo
  enddo

  write (6,*) "Matrix="; call flush(6)
  do ii=1,n
    write (6,*) full_matrix(ii,1:n)
  enddo
  call flush(6)

  deallocate(nelem_nonzero)
  allocate(nelem_nonzero(n),stat=istat)
  if(istat/=0) stop "failed to reallocate nelem_nonzero"

  !shuffling rows to get them in correct order
  full_matrix(1:n,:)=full_matrix(shuffle_list(1:n),:)
  full_matrix(:,1:n)=full_matrix(:,shuffle_list(1:n))

  initial_stride=sum(row_count(1:iam))

  do ii=1,my_n
    my_rows(ii)=initial_stride+ii
  enddo
  !count number of nonzero elements needed
  nelem_nonzero(:) = 0
  do ii=1,n
    do jj=1,my_n
      if(full_matrix(ii,my_rows(jj)) /= 0_rk) nelem_nonzero(ii) = nelem_nonzero(ii) + 1
    enddo
  enddo

  nnonzero=sum(nelem_nonzero(1:n))
  deallocate(matrix_indices,matrix_values)
  allocate(matrix_indices(nnonzero),stat=istat)
  if(istat/=0) stop "failed to reallocate matrix_indices"
  allocate(matrix_values(nnonzero),stat=istat)
  if(istat/=0) stop "failed to reallocate matrix_values"

  !rebuild matrix as local band
  kk = 0
  do ii=1,n
    do jj=1,my_n
      if(full_matrix(ii,my_rows(jj)) /= 0_rk) then
        kk = kk + 1
        matrix_values(kk) = full_matrix(ii,my_rows(jj))
        matrix_indices(kk) = jj
      endif
    enddo
  enddo
  if(kk /= nnonzero) stop "problem in conversion from_upper_triangular_to_band" ! just for debug

end subroutine from_upper_triangular_to_band_shuffle_old

!=====================================================================================================================
  subroutine fast_sparse_matrix_multiply_local_band(n,matrix_indices,nelem_nonzero,matrix_values,vector,answer)
!=====================================================================================================================
  ! A Roggero, Dec 2013

  use types, only : i8b

  integer,intent(in) :: n
  integer(i8b),intent(in) :: matrix_indices(:) ! Vector of column numbers
  integer(i8b),intent(in) :: nelem_nonzero(:)
  real(rk),intent(in) :: matrix_values(:)
  real(rk),intent(in) :: vector(:)
  real(rk),intent(out) :: answer(:) !now answer has still dimension n --> reduce it later!

  integer :: i,j,k,m

  k=0
  do i=1,n
    answer(i)=0._rk
    do j=1,nelem_nonzero(i)
      k=k+1
      m = matrix_indices(k)
      answer(i) = answer(i) + matrix_values(k)*vector(m)
    enddo
  enddo

  end subroutine fast_sparse_matrix_multiply_local_band
!=====================================================================================================================!!=====================================================================================================================
  subroutine fast_sparse_matrix_multiply_local_band_block(block_offset,kstart,matrix_indices,nelem_nonzero,matrix_values,vector,answer)
!=====================================================================================================================
!just do one block of the dimension of answer by starting at the block_offset row, k_start is needed to find the correct entries in matrix_indices and matrix_values
  ! A Roggero, Dec 2013

  use types, only : i8b

  integer,intent(in) :: block_offset
  integer(i8b),intent(in) :: matrix_indices(:) ! Vector of column numbers
  integer(i8b),intent(in) :: nelem_nonzero(:)
  real(rk),intent(in) :: matrix_values(:)
  real(rk),intent(in) :: vector(:)
  real(rk),intent(out) :: answer(:) !now answer has still dimension n --> reduce it later!
  integer, intent(inout) :: kstart

  integer :: i,j,k,m,block_size

  block_size = size(answer)

  k=kstart!make sure at the beginning is set to 0
  do i=1,block_size
    answer(i)=0._rk
    do j=1,nelem_nonzero(i+block_offset)
      k=k+1
      m = matrix_indices(k)
      answer(i) = answer(i) + matrix_values(k)*vector(m)
    enddo
  enddo

  kstart=k !used for the next call to this routine

  end subroutine fast_sparse_matrix_multiply_local_band_block
!=====================================================================================================================!=====================================================================================================================
  subroutine fast_sparse_matrix_multiply_upper_triangular(n,matrix_indices,nelem_nonzero,matrix_values,vector,answer,diag_elems,den1)
!=====================================================================================================================
  ! Multiplies a sparse, symmetric matrix (where only the upper triangular part is stored) by a (non-sparse) vector
  ! Returns a vector that is assumed to be the same size as the input vector
  ! A Holmes, 14 May 2012
  ! Modified : A Holmes, 5 Jan 2013. Make matrix_indices and nelem_nonzero optional; if not present, assume that you're multiplying by only the first row and column. (for HF -> Psi_T)
  ! Also, new optional input diag_elems contains additional diagonal elements of the sparse matrix, and den1, which is the first psi trial denominator
  ! (since in the first row, we want E_num(1)/E_den(1),E_num(2),E_num(3),..., instead of E_num(1),E_num(2),E_num(3),...)

  use types, only : i8b

  integer,intent(in) :: n ! size of vector
  integer(i8b),optional,intent(in) :: matrix_indices(:) ! Vector of column numbers
  integer(i8b),optional,intent(in) :: nelem_nonzero(:)
  real(rk),intent(in) :: matrix_values(:)
  real(rk),intent(in) :: vector(:)
  real(rk),intent(out) :: answer(:)
  real(rk),optional,intent(in) :: diag_elems(:),den1

  integer(i8b) :: i,j,m
  integer(i8b) :: k

  answer(:) = 0._rk

  if (present(matrix_indices)) then ! general sparse upper triangular matrix
    k=0_i8b
    do i=1,n
      do j=1,nelem_nonzero(i)
        k=k+1_i8b
        m = matrix_indices(k)
        answer(i) = answer(i) + matrix_values(k)*vector(m)
        if (i.ne.m)  answer(m) = answer(m) + matrix_values(k)*vector(i)
      enddo
    enddo
  elseif (present(diag_elems)) then ! only first row and diagonal
    answer(1) = answer(1) + matrix_values(1)/den1*vector(1)
    do i=2,n
      answer(i) = answer(i) + matrix_values(i)*vector(1)
      answer(1) = answer(1) + matrix_values(i)*vector(i)
      answer(i) = answer(i) + diag_elems(i)*vector(i) ! only for states in psi_t_connected but outside deterministic space
    enddo
  else ! only first row, and zero the first element
    do i=2,n
      answer(i) = answer(i) + matrix_values(i)*vector(1)
      answer(1) = answer(1) + matrix_values(i)*vector(i)
    enddo
  endif

  end subroutine fast_sparse_matrix_multiply_upper_triangular
!=====================================================================================================================


  subroutine fast_sparse_matrix_multiply_upper_triangular_mpi(n,matrix_indices,nelem_nonzero,matrix_values,vector,answer,local_to_global)
!=====================================================================================================================
  ! Multiplies a sparse, symmetric matrix (where only the upper triangular part is stored) by a (non-sparse) vector
  ! Returns a vector that is assumed to be the same size as the input vector
  ! A Holmes, 14 May 2012
  ! Modified : A Holmes, 5 Jan 2013. Make matrix_indices and nelem_nonzero optional; if not present, assume that you're multiplying by only the first row and column. (for HF -> Psi_T)
  ! Also, new optional input diag_elems contains additional diagonal elements of the sparse matrix, and den1, which is the first psi trial denominator
  ! (since in the first row, we want E_num(1)/E_den(1),E_num(2),E_num(3),..., instead of E_num(1),E_num(2),E_num(3),...)

  use types, only : i8b

  integer,intent(in) :: n! size of vector, index map
  integer(i8b),intent(in) :: local_to_global(:) !index map
  integer(i8b),optional,intent(in) :: matrix_indices(:) ! Vector of column numbers
  integer(i8b),optional,intent(in) :: nelem_nonzero(:) 
  real(rk),intent(in) :: matrix_values(:)
  real(rk),intent(in) :: vector(:)
  real(rk),intent(out) :: answer(:)

  integer(i8b) :: i,j,m,i_global
  integer(i8b) :: k

  answer(:) = 0._rk

  if (present(matrix_indices)) then ! general sparse upper triangular matrix
    k=0_i8b
    do i=1,n
      i_global = local_to_global(i)
      do j=1,nelem_nonzero(i)
        k=k+1_i8b
        m = matrix_indices(k)
        answer(i_global) = answer(i_global) + matrix_values(k)*vector(m)
        if (i_global.ne.m)  answer(m) = answer(m) + matrix_values(k)*vector(i_global)
      enddo
    enddo

! MJO - I'm not sure that these lines are necessary for the MPI version
!  elseif (present(diag_elems)) then ! only first row and diagonal
!    answer(1) = answer(1) + matrix_values(1)/den1*vector(1)
!    do i=2,n
!      answer(i) = answer(i) + matrix_values(i)*vector(1)
!      answer(1) = answer(1) + matrix_values(i)*vector(i)
!      answer(i) = answer(i) + diag_elems(i)*vector(i) ! only for states in psi_t_connected but outside deterministic space
!    enddo
!  else ! only first row, and zero the first element
!    do i=2,n
!      answer(i) = answer(i) + matrix_values(i)*vector(1)
!      answer(1) = answer(1) + matrix_values(i)*vector(i)
!    enddo
  endif

  end subroutine fast_sparse_matrix_multiply_upper_triangular_mpi
!=====================================================================================================================

subroutine reduce_matrix_indices(n,locn,v_indices,matrix_indices,nelem_nonzero,matrix_values)
!=====================================================================================================================
! Reduces the number of indices in matrix_indices to match those present on local core (locn) and contained in v_indices, updates
! nelem_nonzero and reallocates matrix_values and matrix_indices accordingly
  ! A Roggero, Dec 2013

  integer,intent(in) :: n,locn,v_indices(:)
  integer,intent(inout),allocatable :: matrix_indices(:),nelem_nonzero(:)
  real(rk),intent(inout),allocatable:: matrix_values(:)

  integer :: i,j,k,m,l,nk,nonzero
  integer :: tmp_mindices(size(matrix_indices)),tmp_nelemnz(size(nelem_nonzero))
  real(rk) :: tmp_matrix_values(size(matrix_values))

  nonzero=locn
  tmp_nelemnz=0
  k=0; nk=0
  do i=1,n
    do j=1,nelem_nonzero(i)
      k=k+1
      m = matrix_indices(k)
      do l=1,nonzero
        if(m == v_indices(l)) then
          nk = nk + 1
          tmp_nelemnz(i) = tmp_nelemnz(i) + 1
          tmp_mindices(nk) = l
          tmp_matrix_values(nk) = matrix_values(k)
          exit
        endif
      enddo
    enddo
  enddo
  nelem_nonzero=tmp_nelemnz
!  deallocate(matrix_values,matrix_indices)
!  deallocate(matrix_values)
!  allocate(matrix_values(nk),stat=istat))
!  if(istat.ne.0) stop 'failed to reallocate matrix_values in reduce_matrix_indice'
!  allocate(matrix_indices(nk))
  matrix_indices(1:nk)=tmp_mindices(1:nk)
  matrix_values(1:nk)=tmp_matrix_values(1:nk)

end subroutine reduce_matrix_indices
!=====================================================================================================================
  subroutine sparse_fast_sparse_matrix_multiply_upper_triangular(v_indices,n,real_n,matrix_indices,nelem_nonzero,matrix_values,vector,answer,diag_elems,den1,v1)
!=====================================================================================================================
  ! Modified AR, 2/8/13: similar to previous one but for sparse vectors, first input v_indices contains non-zero indices of input vector
  !
  ! Multiplies a sparse, symmetric matrix (where only the upper triangular part is stored) by a (non-sparse) vector
  ! Returns a vector that is assumed to be the same size as the input vector
  ! A Holmes, 14 May 2012
  ! Modified : A Holmes, 5 Jan 2013. Make matrix_indices and nelem_nonzero optional; if not present, assume that you're multiplying by only the first row and column. (for HF -> Psi_T)
  ! Also, new optional input diag_elems contains additional diagonal elements of the sparse matrix, and den1, which is the first psi trial denominator
  ! (since in the first row, we want E_num(1)/E_den(1),E_num(2),E_num(3),..., instead of E_num(1),E_num(2),E_num(3),...)

  integer,intent(in) :: n ,real_n! size of out vector
  integer,optional,intent(in) :: matrix_indices(:) ! Vector of column numbers
  integer,optional,intent(in) :: nelem_nonzero(:)
  real(rk),intent(in) :: matrix_values(:)
  real(rk),intent(in) :: vector(:)
  real(rk),intent(out) :: answer(:)
  real(rk),optional,intent(in) :: diag_elems(:),den1,v1
  integer,intent(in) :: v_indices(:)

  integer :: i,j,k,m,l

  answer(:) = 0._rk

  if (present(matrix_indices)) then ! general sparse upper triangular matrix
    l=1
    k=0
    do i=1,n
      do j=1,nelem_nonzero(i)
        k=k+1
        m = matrix_indices(k)
        answer(i) = answer(i) + matrix_values(k)*vector(m)
        if ((i.ne.v_indices(m)).and.(i==v_indices(l))) then
          answer(v_indices(m)) = answer(v_indices(m)) + matrix_values(k)*vector(l)
          l = l + 1
        endif
      enddo
    enddo

  elseif (present(diag_elems)) then ! only first row and diagonal

    if(v_indices(1) == 1) answer(1)=answer(1)+matrix_values(1)/den1*v1
    do i=2,real_n
        answer(1+i) = answer(1+i) + matrix_values(i)*v1
    enddo
    if(v_indices(1) > 1) then
      answer(1) = answer(1) + matrix_values(1)*vector(1)
      answer(2) = answer(2) + diag_elems(1)*vector(1)+matrix_values(1)*v1 ! only for states in psi_t_connected but outside deterministic space
    endif

    do i=2,real_n
      answer(1) = answer(1) + matrix_values(i)*vector(i)
      answer(1+i) = answer(1+i) + diag_elems(i)*vector(i) ! only for states in psi_t_connected but outside deterministic space
    enddo

  else ! only first row, and zero the first element
    do i=2,real_n
      answer(1+i) = answer(1+i) + matrix_values(v_indices(i))*v1
    enddo
    if(v_indices(1) > 1) then
      answer(1) = answer(1) + matrix_values(v_indices(1))*vector(1)
      answer(2) = answer(2) + matrix_values(v_indices(1))*v1
    endif
    do i=2,real_n
      answer(1) = answer(1) + matrix_values(v_indices(i))*vector(i)
    enddo
  endif

  end subroutine sparse_fast_sparse_matrix_multiply_upper_triangular


  subroutine binary_search_lbound(item, arr, idx)
    implicit none
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: item, arr(:)
#else
    integer(ik), intent(in) :: item, arr(:)
#endif
    integer, intent(out) :: idx
    integer :: left, right, mid

    left = 0
    right = size(arr)
    do while(right - left > 1)
      mid = (left + right) / 2
      if (arr(mid) < item) then
        left = mid
      else
        right = mid
      endif
    enddo
    idx = right
  end subroutine binary_search_lbound

  subroutine binary_search_rbound(item, arr, idx)
    implicit none
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec), intent(in) :: item, arr(:)
#else
    integer(ik), intent(in) :: item, arr(:)
#endif
    integer, intent(out) :: idx
    integer :: left, right, mid

    left = 1
    right = size(arr) + 1
    do while(right - left > 1)
      mid = (left + right) / 2
      if (arr(mid) <= item) then
        left = mid
      else
        right = mid
      endif
    enddo
    idx = left
  end subroutine binary_search_rbound
!=====================================================================================================================
 subroutine binary_search_single(search_for_up,search_in_up,ind)
!=====================================================================================================================

    ! Searches for determinant search_for_up,dn in the list search_in_up,dn and returns its index in that list.
    ! If it isn't in the list, return 0.
    ! Currently only used in hubbard.f90
    ! A. Holmes, 13 Mar 2012

    implicit none

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: search_for_up  ! List of dets we are searching for (not necessarily in order)
    type(ik_vec),intent(in) :: search_in_up(:)! List of dets we are searching in (assumed to be in order)
#else
    integer(ik),intent(in) :: search_for_up  ! List of dets we are searching for (not necessarily in order)
    integer(ik),intent(in) :: search_in_up(:)! List of dets we are searching in (assumed to be in order)
#endif
    integer,intent(out) :: ind ! Index in search_in where search_for is located.

    integer :: lower_bound,upper_bound,midpoint,i

    lower_bound = 1
    upper_bound = size(search_in_up)
    midpoint = nint(.5*(lower_bound+upper_bound))

    if (search_for_up==search_in_up(lower_bound)) then
      ind = lower_bound
      return
    endif

    if (search_for_up==search_in_up(upper_bound)) then
      ind = upper_bound
      return
    endif

    do i=1,100 ! This is only here to prevent the program from running forever in case I made a mistake. The program should terminate on its own long before reaching i=100.
      if (search_for_up==search_in_up(midpoint)) then
        ind = midpoint
        return
      endif
      if (midpoint==upper_bound.or.midpoint==lower_bound) then
        ind = 0 ! Return 0 if we don't find it in the list.
        return
      endif
      if (abs(search_for_up)<abs(search_in_up(midpoint)).or.(search_for_up==search_in_up(midpoint))) then
        upper_bound=midpoint
        midpoint=nint(.5*(lower_bound+upper_bound))
        cycle
      endif
      if (abs(search_for_up)>abs(search_in_up(midpoint))) then
        lower_bound=midpoint
        midpoint=nint(.5*(lower_bound+upper_bound))
        cycle
      endif
    enddo

  end subroutine binary_search_single


!=====================================================================================================================
 subroutine binary_search_int(search_for_up,search_for_dn,search_in_up,search_in_dn,midpoint)
!=====================================================================================================================

    ! Searches for determinant search_for_up,dn in the list search_in_up,dn and returns its index in that list.
    ! If it isn't in the list, return 0.
    ! A. Holmes, 13 Mar 2012

    implicit none

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: search_for_up,search_for_dn ! List of dets we are searching for (not necessarily in order)
    type(ik_vec),intent(in) :: search_in_up(:),search_in_dn(:) ! List of dets we are searching in (assumed to be in order)
    type(ik_vec) :: temp_det_up, temp_det_dn, search_for_up_abs, search_for_dn_abs
#else
    integer(ik),intent(in) :: search_for_up,search_for_dn ! List of dets we are searching for (not necessarily in order)
    integer(ik),intent(in) :: search_in_up(:),search_in_dn(:) ! List of dets we are searching in (assumed to be in order)
    integer(ik) :: temp_det_up, temp_det_dn, search_for_up_abs, search_for_dn_abs
#endif
    integer,intent(out) :: midpoint ! Index in search_in where search_for is located.

    integer :: lower_bound,upper_bound

    lower_bound = 1
    upper_bound = size(search_in_up)

    search_for_up_abs = abs(search_for_up)
    search_for_dn_abs = abs(search_for_dn)

    do while (lower_bound <= upper_bound)
      midpoint = ishft(lower_bound + upper_bound, -1) ! Quick way of dividing by two.
      temp_det_up = search_in_up(midpoint)
      temp_det_dn = search_in_dn(midpoint)
      if (search_for_up_abs > abs(temp_det_up) .or. &
          & (search_for_up == temp_det_up .and. search_for_dn_abs > abs(temp_det_dn))) then
        lower_bound = midpoint + 1
      elseif (search_for_up == temp_det_up .and. search_for_dn == temp_det_dn) then
        return
      else !(search_for_up<search_in_up(midpoint).or.(search_for_up==search_in_up(midpoint).and.search_for_dn<search_in_dn(midpoint)))
        upper_bound=midpoint-1
      endif
    enddo

    midpoint=0
    return


  end subroutine binary_search_int


 subroutine binary_search_i8b(search_for_up,search_for_dn,search_in_up,search_in_dn,midpoint)
!=====================================================================================================================

    ! Searches for determinant search_for_up,dn in the list search_in_up,dn and returns its index in that list.
    ! If it isn't in the list, return 0.
    ! A. Holmes, 13 Mar 2012

    use types, only : i8b

    implicit none

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: search_for_up,search_for_dn ! List of dets we are searching for (not necessarily in order)
    type(ik_vec),intent(in) :: search_in_up(:),search_in_dn(:) ! List of dets we are searching in (assumed to be in order)
#else
    integer(ik),intent(in) :: search_for_up,search_for_dn ! List of dets we are searching for (not necessarily in order)
    integer(ik),intent(in) :: search_in_up(:),search_in_dn(:) ! List of dets we are searching in (assumed to be in order)
#endif
    integer(i8b),intent(out) :: midpoint ! Index in search_in where search_for is located.

    integer(i8b) :: lower_bound,upper_bound

    lower_bound = 1
    upper_bound = size(search_in_up)

    do while (lower_bound <= upper_bound)
      midpoint = (lower_bound+upper_bound)/2
      if (abs(search_for_up)>abs(search_in_up(midpoint)).or.(search_for_up==search_in_up(midpoint).and.abs(search_for_dn)>abs(search_in_dn(midpoint)))) then
        lower_bound=midpoint+1
      elseif (search_for_up==search_in_up(midpoint).and.search_for_dn==search_in_dn(midpoint)) then
        return
      else !(search_for_up<search_in_up(midpoint).or.(search_for_up==search_in_up(midpoint).and.search_for_dn<search_in_dn(midpoint)))
        upper_bound=midpoint-1
      endif
    enddo

    midpoint=0
    return


  end subroutine binary_search_i8b


!=====================================================================================================================
  subroutine binary_search_list_and_update(search_for_up,search_for_dn,walk_wt,search_in_up,search_in_dn,psi_t_connected_e_loc_num,psi_t_connected_e_loc_den,e_num_gen,e_den_gen,e_num2_cum,e_den2_cum,e_num_abs_cum,e_den_abs_cum,e_num_e_den_cum,e_num_walker,e_den_walker)
!=====================================================================================================================

    ! Searches for all of the elements in search_for_up,dn (assumed to be sorted) in search_in_up,dn, and updates all the running sums when it finds them.
    ! Scales as size(search_for_up)*log[size(search_in_up)]
    ! A. Holmes, 27 Mar 2012

    implicit none

    real(rk),intent(in) :: walk_wt(:)
#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: search_for_up(:),search_for_dn(:) ! List of dets we are searching for (assumed to be in order)
    type(ik_vec),intent(in) :: search_in_up(:),search_in_dn(:) ! List of dets we are searching in (assumed to be in order)
#else
    integer(ik),intent(in) :: search_for_up(:),search_for_dn(:) ! List of dets we are searching for (assumed to be in order)
    integer(ik),intent(in) :: search_in_up(:),search_in_dn(:) ! List of dets we are searching in (assumed to be in order)
#endif
    real(rk),intent(in) :: psi_t_connected_e_loc_num(:),psi_t_connected_e_loc_den(:)
    real(rk),intent(inout) :: e_num_walker(:),e_den_walker(:)
    real(rk),intent(inout) :: e_num_gen,e_den_gen,e_num2_cum,e_den2_cum,e_num_abs_cum,e_den_abs_cum,e_num_e_den_cum

    integer :: lower_bound,upper_bound,i,j
    real(rk) :: e_num,e_den

    lower_bound = 1
    upper_bound = size(search_in_up)

    do i=size(search_for_up),1,-1
      if (e_num_walker(i)>1.e50_rk) then
        call binary_search(search_for_up(i),search_for_dn(i),search_in_up(lower_bound:upper_bound),search_in_dn(lower_bound:upper_bound),j)
        if (j.eq.0) then
          e_num_walker(i) = 0._rk
          e_den_walker(i) = 0._rk
        else
          if (j.ne.upper_bound)  upper_bound=j-1
          e_num_walker(i) = psi_t_connected_e_loc_num(j)
          e_den_walker(i) = psi_t_connected_e_loc_den(j)
        endif
      endif

      e_num=e_num_walker(i)*walk_wt(i)
      e_den=e_den_walker(i)*walk_wt(i)

      if (e_num.ne.0._rk) then
        if (abs(e_den).lt.1.d-22)  e_den=abs(e_den) ! This is necessary because the sign of e_num depends on the sign of e_den, which can be either +0 or -0
        e_den_gen=e_den_gen+e_den
        e_den2_cum=e_den2_cum+e_den**2
        e_den_abs_cum=e_den_abs_cum+abs(e_den)
        e_num_gen=e_num_gen+e_num
        e_num2_cum=e_num2_cum+e_num**2
        e_num_abs_cum=e_num_abs_cum+e_num*sign(1._rk,e_den)
        e_num_e_den_cum=e_num_e_den_cum+e_num*e_den
        !write(6,'(''i,j,e_den,e_num,walk_wt(i)'',2i6,9es14.6)') i,j,e_den,e_num,walk_wt(i),e_num_walker(i),e_den_walker(i)
      endif
      if (lower_bound.gt.upper_bound) return
    enddo

  end subroutine binary_search_list_and_update


!=====================================================================================================================
 integer function sample_discrete_distribution(c_probs)
!=====================================================================================================================

    ! Samples from a discrete distribution specified by the vector of cumulative probabilities c_probs,
    ! assumed to be properly normalized.
    ! Uses binary search to locate the correct element, and returns its index.
    ! Scales as O(log n), so sample_alias (O(1)) is preferred.
    ! A Holmes, 21 Mar 2015

    implicit none

    real(rk),intent(in) :: c_probs(:)

    integer :: midpoint ! Index in search_in where search_for is located.
    integer :: lower_bound,upper_bound
    real(rk) :: r,rannyu

    r = rannyu()

    lower_bound = 1
    upper_bound = size(c_probs)

    do while (lower_bound < upper_bound)
      midpoint = (lower_bound+upper_bound)/2
      if (r<c_probs(midpoint)) then
        upper_bound=midpoint
      else
        lower_bound=midpoint+1
      endif
    enddo

    sample_discrete_distribution = lower_bound

  end function sample_discrete_distribution

!=====================================================================================================================
 subroutine real_sym_gen_eig(n,a,b,eigenvectors,eigenvalues)
!=====================================================================================================================
 ! Created by : Hitesh Changlani, March 19,2012
 ! Purpose    : Solve the generalized eigenproblem by calling LAPACK dsygv

 implicit none
 ! Dummy
 integer,intent(in)   :: n
 real(rk),intent(in)  :: a(n,n),b(n,n)
 real(rk),intent(out) :: eigenvectors(n,n),eigenvalues(n)
 ! Local
 integer              :: len_work,info
 real(rk),allocatable :: work(:)
 real(rk)             :: bcopy(n,n)

 len_work = 3 * n -1
 allocate(work(len_work))

 eigenvectors=a
 bcopy=b

 call dsygv(1,'V', 'U', n, eigenvectors, n, bcopy,n, eigenvalues, work, len_work, info)
 deallocate(work)

 if (info /= 0) then
        write(6,*) "info = ", info
        stop 'Generalized Eigenproblem Failed!'
 endif

 end subroutine real_sym_gen_eig

!====================================================================================================================
 subroutine overlap_two_slater(n_s,n_e,slater1,slater2,overlap)
!=====================================================================================================================
 ! Created by : Hitesh Changlani, March 19,2012
 ! Purpose    : Overlap between Mean Field solutions (Slater determinants)

 implicit none
 ! Dummy
 integer,intent(in)   :: n_s,n_e
 real(rk),intent(in)  :: slater1(n_s,n_e),slater2(n_s,n_e)
 real(rk),intent(out) :: overlap
 ! Local
 real(rk) :: s1ts2(n_e,n_e)

 s1ts2=matmul(transpose(slater1),slater2)  ! Refer Shiwei Zhang and general CPMC literature
 call determinant(s1ts2,n_e,overlap)

 end subroutine overlap_two_slater

!====================================================================================================================
 subroutine gp_two_slater(n_s,n_e,slater1,slater2,gp)
!=====================================================================================================================
 ! Created by : Hitesh Changlani, March 19,2012
 ! Purpose    : Solve the generalized eigenproblem by calling LAPACK dsygv

 implicit none
 ! Dummy
 integer,intent(in)   :: n_s,n_e
 real(rk),intent(in)  :: slater1(n_s,n_e),slater2(n_s,n_e)
 real(rk),intent(out) :: gp(n_s,n_s)
 ! Local
 real(rk) :: s1ts2(n_e,n_e),det

 s1ts2=matmul(transpose(slater1),slater2)               ! Refer Shiwei Zhang and general CPMC literature
 call matinv(s1ts2,n_e,det)
 gp=matmul(matmul(slater2,s1ts2),transpose(slater1))    ! G'_ij = <c_j^dagger c_i>

 end subroutine gp_two_slater

!=====================================================================================================================
 subroutine create_kspace_sym_maps(k_vectors,c4_map,reflection_map)
!=====================================================================================================================
 ! Created by             : Hitesh Changlani
 ! Date last edited       : April 3, 2012
 ! Purpose                : Construct in momentum space the C4 maps and the reflection map a new configuration with a new set of labels given by map. Such a routine is useful for spatial symmetries
 !                        : This is a one time computation done at the time of system_setup and so one need not worry too much about efficiency issues

 implicit none
 ! Dummy
 integer,intent(in)               :: k_vectors(:,:)
 integer,intent(out)              :: c4_map(:,:),reflection_map(:)

 ! Local
 integer                          :: i,l_x,nsites,row,col,orig_site,new_site
 integer,allocatable              :: inv_map(:)
 integer                          :: k_x,k_y,k_x_new,k_y_new

 nsites=size(k_vectors,2)
 l_x=int(sqrt(real(nsites)))
 allocate(inv_map(nsites))

 ! Create inverse map
 write (6,*) "Site (Energy order)    K_x           K_y   Location on original lattice"
 do i=1,nsites
    col=(k_vectors(1,i)+l_x)/2
    row=(k_vectors(2,i)+l_x)/2
    orig_site=l_x*(row-1)+col
    inv_map(orig_site)=i
    write (6,*) i,k_vectors(1,i),k_vectors(2,i),orig_site
 enddo

 ! C_1 map  (kx,ky) -> (ky,-kx)
 do i=1,nsites
    k_x_new=k_vectors(2,i)
    if (k_vectors(1,i) .eq. l_x) then
        k_y_new=k_vectors(1,i)
    else
        k_y_new=-k_vectors(1,i)
    endif
    col=(k_x_new+l_x)/2
    row=(k_y_new+l_x)/2
    orig_site=l_x*(row-1)+col
    new_site=inv_map(orig_site)
    c4_map(1,i)=new_site
 enddo

 ! C_2 map  (kx,ky) -> (-kx,-ky)
 do i=1,nsites
    if (k_vectors(1,i) .eq. l_x) then
        k_x_new=k_vectors(1,i)
    else
        k_x_new=-k_vectors(1,i)
    endif
    if (k_vectors(2,i) .eq. l_x) then
        k_y_new=k_vectors(2,i)
    else
        k_y_new=-k_vectors(2,i)
    endif
    col=(k_x_new+l_x)/2
    row=(k_y_new+l_x)/2
    orig_site=l_x*(row-1)+col
    new_site=inv_map(orig_site)
    c4_map(2,i)=new_site
 enddo

 ! C_3 map  (kx,ky) -> (-ky,kx)
 do i=1,nsites
    k_y_new=k_vectors(1,i)
    if (k_vectors(2,i) .eq. l_x) then
        k_x_new=k_vectors(2,i)
    else
        k_x_new=-k_vectors(2,i)
    endif
    col=(k_x_new+l_x)/2
    row=(k_y_new+l_x)/2
    orig_site=l_x*(row-1)+col
    new_site=inv_map(orig_site)
    c4_map(3,i)=new_site
 enddo

 ! Application of reflection about the diagonal (kx,ky)->(-k_y,-k_x)
 do i=1,nsites
    k_x=k_vectors(1,i)
    k_y=k_vectors(2,i)
    if (k_x .eq. l_x) then
         k_y_new=k_x
    else
        k_y_new=-k_x
    endif
    if (k_y .eq. l_x) then
         k_x_new=k_y
    else
        k_x_new=-k_y
    endif
    col=(k_x_new+l_x)/2
    row=(k_y_new+l_x)/2
    orig_site=l_x*(row-1)+col
    new_site=inv_map(orig_site)
    reflection_map(i)=new_site
 enddo

 end subroutine create_kspace_sym_maps

!=====================================================================================================================
 subroutine create_rspace_sym_maps(l_x,c4_map,reflection_map)
!=====================================================================================================================
 ! Created by             : Hitesh Changlani
 ! Date last edited       : April 3, 2012
 ! Purpose                : Construct in real space the C4 maps and the reflection map a new configuration with a new set of labels given by map. Such a routine is useful for spatial symmetries
 !                        : This is a one time computation done at the time of system_setup and so one need not worry too much about efficiency issues

 implicit none
 ! Dummy
 integer,intent(in)               :: l_x
 integer,intent(out)              :: c4_map(:,:),reflection_map(:)

 ! Local
 integer                          :: i,j,nsites
 integer,allocatable              :: c1mat(:,:),c2mat(:,:),c3mat(:,:)
 integer                          :: row,col,new_site
 logical                          :: found

 nsites=l_x*l_x
 allocate(c1mat(nsites,nsites))
 allocate(c2mat(nsites,nsites))
 allocate(c3mat(nsites,nsites))

 ! Consider a square 4x4 lattice

 ! 1  2  3  4
 ! 5  6  7  8
 ! 9  10 11 12
 ! 13 14 15 16

 ! Application of C^1 a 90 degree rigid rotation gives

 ! 13 9  5 1
 ! 14 10 6 2
 ! 15 11 7 3
 ! 16 12 8 4

 ! Create C_4 maps
 c1mat(:,:)=0
 c2mat(:,:)=0
 c3mat(:,:)=0

 ! C_1 map
 do i=1,nsites
    row=((i-1)/l_x)+1
    col=i-((row-1)*l_x)
    new_site=(l_x*(col-1))+(l_x-row+1)
    c4_map(1,new_site)=i
    c1mat(i,new_site)=1
 enddo

 c2mat=matmul(c1mat,c1mat)
 !c3mat=matmul(c1mat,c2mat)

 call print_int_matrix(nsites,nsites,c1mat)
 call print_int_matrix(nsites,nsites,c2mat)
 !call print_int_matrix(nsites,nsites,c3mat)

 ! Convert c2mat to a map - note they have only one 1 element in every row and column
 do i=1,nsites
   j=1
   found=.false.
   do
     if (c2mat(i,j) .eq. 1) then
        found=.true.
        exit
     endif
     j=j+1
   enddo
   c4_map(2,i)=j
 enddo

 ! Convert c3mat to a map - note they have only one 1 element in every row and column
! do i=1,nsites
!   j=1
!   found=.false.
!   do
!     if (c3mat(i,j) .eq. 1) then
!        found=.true.
!        exit
!     endif
!     j=j+1
!   enddo
!   c4_map(3,i)=j
! enddo

  do i=1,nsites
     row=((i-1)/l_x)+1
     col=i-((row-1)*l_x)
     new_site=(l_x*(l_x-col+1-1))+(row)
     c4_map(3,new_site)=i
!    c1mat(i,new_site)=1
 enddo


 ! Application of reflection around the y axis gives

 ! 4  3 | 2  1
 ! 8  7 | 6  5
 ! 12 11| 10 9
 ! 16 15| 14 13

 ! Create Reflection map
 do i=1,nsites
    row=((i-1)/l_x)+1
    col=i-((row-1)*l_x)
    reflection_map(l_x*(row-1)+(l_x-col+1))=i
 enddo

 end subroutine create_rspace_sym_maps

!=====================================================================================================================
 subroutine print_map_on_square(l_x,map)
!=====================================================================================================================
 ! Created by    : Hitesh Changlani
 ! Date created  : April 3 2012
 ! Purpose       : To print map on a square lattice. This is mainly to check the rotation+reflection+spin-inversion
 !                 symmetry adapted code for the Hubbard model

 implicit none
 integer,intent(in)  :: l_x
 integer,intent(in)  :: map(l_x*l_x)

 integer i,j

 do i=1,l_x
             write(6,'(100i5)') (map(((i-1)*l_x)+j),j=1,l_x)
 enddo
 call flush(6)

 end subroutine print_map_on_square

!============================================================================================================================
 subroutine relabel_efficient(map,up_locations,dn_locations,new_det_up,new_det_dn,new_up_locations,new_dn_locations,sign_map)
!============================================================================================================================
 ! Created by : Hitesh Changlani, April 3, 2012
 ! Purpose    : Construct a new configuration with a new set of labels given by map. Such a routine is useful for spatial symmetries
 use generic_sort, only : sort

 implicit none
 ! Dummy
 integer,intent(in)               :: map(:)
 integer,intent(in)               :: up_locations(:),dn_locations(:)
#ifdef NUM_ORBITALS_GT_127
 type(ik_vec),intent(out)          :: new_det_up,new_det_dn
#else
 integer(ik),intent(out)          :: new_det_up,new_det_dn
#endif
 real(rk),intent(out)             :: sign_map
 integer,intent(out)              :: new_up_locations(:),new_dn_locations(:)

 ! Local
 integer                          :: i,counter,nup_ops,ndn_ops

 new_det_up=0_ik
 new_det_dn=0_ik
 sign_map=1
 counter=0

 do i=1,size(up_locations)
        new_det_up=ibset(new_det_up,map(up_locations(i))-1)
        counter=counter+1
        new_up_locations(counter)=map(up_locations(i))
 enddo

 call sort(counter,new_up_locations,nup_ops)            ! This is a new order for the fermions
                                                        ! and has to be brought into the 1,2,3... form

 counter=0

 do i=1,size(dn_locations)
        new_det_dn=ibset(new_det_dn,map(dn_locations(i))-1)
        counter=counter+1
        new_dn_locations(counter)=map(dn_locations(i))
 enddo

 call sort(counter,new_dn_locations,ndn_ops)

 if (mod(nup_ops+ndn_ops,2) .ne. 0) sign_map=-1

 end subroutine relabel_efficient

!===================================================================
 subroutine relabel(map,det_up,det_dn,new_det_up,new_det_dn,sign_map)
!===================================================================
 ! Created by : Hitesh Changlani, April 3, 2012
 ! Purpose    : Construct a new configuration with a new set of labels given by map. Such a routine is useful for spatial symmetries
 use generic_sort, only : sort

 implicit none
 ! Dummy
 integer,intent(in)               :: map(:)
#ifdef NUM_ORBITALS_GT_127
 type(ik_vec),intent(in)           :: det_up,det_dn
 type(ik_vec),intent(out)          :: new_det_up,new_det_dn
#else
 integer(ik),intent(in)           :: det_up,det_dn
 integer(ik),intent(out)          :: new_det_up,new_det_dn
#endif
 real(rk),intent(out)             :: sign_map

 ! Local
 integer                          :: i,counter,nup_ops,ndn_ops
 integer                          :: labels(size(map,1))

 new_det_up=0_ik
 new_det_dn=0_ik
 sign_map=1
 counter=0

 do i=1,size(map)                                                      ! Can also utilize the number of up electrons to reduce this cost for very dilute lattices
      if (btest(det_up,i-1)) then
        new_det_up=ibset(new_det_up,map(i)-1)
        counter=counter+1
        labels(counter)=map(i)
      endif
 enddo

 call sort(counter,labels,nup_ops)                                    ! This is a new order for the fermions and has to be brought into the 1,2,3... form

 counter=0

 do i=1,size(map)                                                     ! Can also utilize the number of down electrons to reduce this cost for very dilute lattices
      if (btest(det_dn,i-1)) then
        counter=counter+1
        labels(counter)=map(i)
        new_det_dn=ibset(new_det_dn,map(i)-1)
      endif
 enddo

 call sort(counter,labels,ndn_ops)

 if (mod(nup_ops+ndn_ops,2) .ne. 0) sign_map=-1

 end subroutine relabel

!=======================================================================================================================================================
 subroutine print_sym_configs_scalar(c4_map,reflection_map,z,p,det_up_in,det_dn_in)
!=======================================================================================================================================================
 ! Created by : Hitesh Changlani, April 5,2012
 ! Purpose    : Generate all four fold rotated+ spatial inversion + time reversed configurations and print them
 use overload
 use types, only : num_words

 implicit none
 ! Dummy
 integer                  :: c4_map(:,:),reflection_map(:)     ! C_4 and Reflection Maps
 integer,intent(in)       :: z,p                               ! z and p are quantum numbers for time reversal (spin inversion) and parity (spatial reflection)
! type(ik_vec),intent(in)   :: det_up_in,det_dn_in                     ! Incoming representative determinant up and down
 integer(ik),intent(in)   :: det_up_in,det_dn_in                     ! Incoming representative determinant up and down
#ifdef NUM_ORBITALS_GT_127
 type(ik_vec)             :: det_up,det_dn
 type(ik_vec)              :: rep_up,rep_dn
 type(ik_vec)              :: new_dets_up(16),new_dets_dn(16)   ! Symmetry related determinants
#else
 integer(ik)             :: det_up,det_dn
 integer(ik)              :: rep_up,rep_dn
 integer(ik)              :: new_dets_up(16),new_dets_dn(16)   ! Symmetry related determinants
#endif

 ! Local
 real(rk)                 :: phases(16)                        ! Phase factors in the linear combination
 real(rk)                 :: phase_w_rep,norm                  ! Norm
 integer                  :: i,num_distinct
 character(len=2)         :: fmt

 det_up = det_up_in
 det_dn = det_dn_in

 call generate_fourfold_k_configs(c4_map,reflection_map,z,p,det_up,det_dn,new_dets_up,new_dets_dn,phases,rep_up,rep_dn,phase_w_rep,norm,num_distinct)

 write (fmt,'(i2)') 2*num_words+1
 do i=1,16
   write(6,'(' // trim(fmt) // 'i5)') new_dets_up(i),new_dets_dn(i),int(phases(i))
 enddo

 write (6,*) "norm=",norm
 write (6,*)

 end subroutine print_sym_configs_scalar

#ifdef NUM_ORBITALS_GT_127
!=======================================================================================================================================================
 subroutine print_sym_configs_vec(c4_map,reflection_map,z,p,det_up_in,det_dn_in)
!=======================================================================================================================================================
 ! Created by : Hitesh Changlani, April 5,2012
 ! Purpose    : Generate all four fold rotated+ spatial inversion + time reversed configurations and print them
 use overload
 use types, only : num_words

 implicit none
 ! Dummy
 integer                  :: c4_map(:,:),reflection_map(:)     ! C_4 and Reflection Maps
 integer,intent(in)       :: z,p                               ! z and p are quantum numbers for time reversal (spin inversion) and parity (spatial reflection)
 type(ik_vec),intent(in)   :: det_up_in,det_dn_in                     ! Incoming representative determinant up and down
#ifdef NUM_ORBITALS_GT_127
 type(ik_vec)              :: det_up,det_dn
 type(ik_vec)              :: rep_up,rep_dn
 type(ik_vec)              :: new_dets_up(16),new_dets_dn(16)   ! Symmetry related determinants
#else
 integer(ik)              :: det_up,det_dn
 integer(ik)              :: rep_up,rep_dn
 integer(ik)              :: new_dets_up(16),new_dets_dn(16)   ! Symmetry related determinants
#endif

 ! Local
 real(rk)                 :: phases(16)                        ! Phase factors in the linear combination
 real(rk)                 :: phase_w_rep,norm                  ! Norm
 integer                  :: i,num_distinct
 character(len=2)         :: fmt

 det_up = det_up_in
 det_dn = det_dn_in

 call generate_fourfold_k_configs(c4_map,reflection_map,z,p,det_up,det_dn,new_dets_up,new_dets_dn,phases,rep_up,rep_dn,phase_w_rep,norm,num_distinct)

 write (fmt,'(i2)') 2*num_words+1
 do i=1,16
   write(6,'(' // trim(fmt) // 'i5)') new_dets_up(i),new_dets_dn(i),int(phases(i))
 enddo

 write (6,*) "norm=",norm
 write (6,*)

 end subroutine print_sym_configs_vec
#endif

!=======================================================================================================================================================
 subroutine get_rep_only(c4_map,reflection_map,z,p,det_up,det_dn,rep_up,rep_dn)
!=======================================================================================================================================================
 ! Created by : Hitesh Changlani, April 5,2012
 ! Purpose    : Generate all four fold rotated+ spatial inversion + time reversed configurations and print them

 implicit none
 ! Dummy
 integer                  :: c4_map(:,:),reflection_map(:)     ! C_4 and Reflection Maps
 integer,intent(in)       :: z,p                               ! z and p are quantum numbers for time reversal (spin inversion) and parity (spatial reflection)
#ifdef NUM_ORBITALS_GT_127
 type(ik_vec),intent(in)   :: det_up,det_dn                     ! Incoming representative determinant up and down
 type(ik_vec),intent(out)  :: rep_up,rep_dn
 type(ik_vec)              :: new_dets_up(16),new_dets_dn(16)   ! Symmetry related determinants
#else
 integer(ik),intent(in)   :: det_up,det_dn                     ! Incoming representative determinant up and down
 integer(ik),intent(out)  :: rep_up,rep_dn
 integer(ik)              :: new_dets_up(16),new_dets_dn(16)   ! Symmetry related determinants
#endif

 ! Local
 real(rk)                 :: phases(16)                        ! Phase factors in the linear combination
 real(rk)                 :: phase_w_rep,norm                  ! Norm
 integer                  :: num_distinct

 call generate_fourfold_k_configs(c4_map,reflection_map,z,p,det_up,det_dn,new_dets_up,new_dets_dn,phases,rep_up,rep_dn,phase_w_rep,norm,num_distinct)

 end subroutine get_rep_only

!=================================================================================================================================================================================
 subroutine generate_fourfold_k_configs_efficient(c4_map,reflection_map,z,p,nup,ndn,det_up_in,det_dn_in,new_dets_up,new_dets_dn,phases,rep_up,rep_dn,phase_w_rep,norm,num_distinct)
!=================================================================================================================================================================================
 ! Created by : Hitesh Changlani, April 2,2012
 ! Purpose    : Generate all four fold rotated+ spatial inversion + time reversed configurations
 !              Also return the norm of the configuration as given by group theoretical considerations

 use overload
 implicit none
 ! Dummy
 integer,intent(in)       :: c4_map(:,:),reflection_map(:)              ! C_4 and Reflection Maps
 integer,intent(in)       :: nup,ndn
 integer,intent(in)       :: z,p                                        ! z and p are quantum numbers for time reversal (spin inversion) and parity (spatial reflection)
 real(rk),intent(out)     :: phases(16)                                 ! Phase factors in the linear combination
 real(rk),intent(out)     :: phase_w_rep,norm                           ! Norm
 integer,intent(out)      :: num_distinct                               ! Number of distinct configurations
#ifdef NUM_ORBITALS_GT_127
 type(ik_vec),intent(in)   :: det_up_in,det_dn_in                              ! Incoming representative determinant up and down
 type(ik_vec)             :: det_up,det_dn                              ! Incoming representative determinant up and down
 type(ik_vec),intent(out)  :: new_dets_up(16),new_dets_dn(16)            ! Symmetry related determinants
 type(ik_vec),intent(out)  :: rep_up,rep_dn                              ! Representatives
 type(ik_vec)              :: det
#else
 integer(ik),intent(in)   :: det_up_in,det_dn_in                              ! Incoming representative determinant up and down
 integer(ik)             :: det_up,det_dn                              ! Incoming representative determinant up and down
 integer(ik),intent(out)  :: new_dets_up(16),new_dets_dn(16)            ! Symmetry related determinants
 integer(ik),intent(out)  :: rep_up,rep_dn                              ! Representatives
 integer(ik)               :: det
#endif
 ! Local
 integer                  :: up_locations(16,nup),dn_locations(16,ndn)  ! Save locations of electrons
 integer                  :: i,counter

 det_up = det_up_in
 det_dn = det_dn_in

 ! Norm computation only true for Q=0 group which has G=A1 symmetry
 num_distinct=0
 phase_w_rep=1._rk

 rep_up=det_up
 rep_dn=det_dn

 ! We call lattice reflection as "parity", a convention used often

 ! Map i - C^i  need only i=1,2,3 i = 0 and 4 lead to identity transformation
 phases(1)=1
 new_dets_up(1)=det_up
 new_dets_dn(1)=det_dn
 norm=1._rk

! counter=0
! i=0
! do
!    i=i+1
!    if (btest(det_up,i-1)) then
!       counter=counter+1
!       up_locations(1,counter)=i
!    endif
!    if (counter .eq. nup) exit
! enddo

! counter=0
! i=0
! do
!    i=i+1
!    if (btest(det_dn,i-1)) then
!       counter=counter+1
!       dn_locations(1,counter)=i
!    endif
!    if (counter .eq. ndn) exit
! enddo

! New way of counting - taken from Adam
 counter=0
 i=0
 det = det_up
 do while (det.ne.0)
      i = trailz(det)+1
      counter=counter+1
      up_locations(1,counter)=i
      det = ibclr(det,i-1)
 enddo

 counter=0
 i=0
 det = det_dn
 do while (det.ne.0)
      i = trailz(det)+1
      counter=counter+1
      dn_locations(1,counter)=i
      det = ibclr(det,i-1)
 enddo
! End of edit

 counter=0
 i=0

 do i=1,3
    call relabel_efficient(c4_map(i,:),up_locations(1,:),dn_locations(1,:),new_dets_up(i+1),new_dets_dn(i+1),up_locations(i+1,:),dn_locations(i+1,:),phases(i+1))
    if ((new_dets_up(i+1) .eq. det_up) .and. (new_dets_dn(i+1) .eq. det_dn)) norm=norm+phases(i+1)
 enddo

 ! 4 fold + Spin inversion
 do i=5,8
    new_dets_up(i)=new_dets_dn(i-4)
    up_locations(i,:)=dn_locations(i-4,:)
    new_dets_dn(i)=new_dets_up(i-4)
    dn_locations(i,:)=up_locations(i-4,:)
    phases(i)=z*phases(i-4)
    if ((new_dets_up(i) .eq. det_up) .and. (new_dets_dn(i) .eq. det_dn)) norm=norm+phases(i)
 enddo

 ! Lattice reflection around y=-x axis
 do i=9,12
    call relabel_efficient(reflection_map(:),up_locations(i-8,:),dn_locations(i-8,:),new_dets_up(i),new_dets_dn(i),up_locations(i,:),dn_locations(i,:),phases(i))
    phases(i)=phases(i)*p*phases(i-8)
    if ((new_dets_up(i) .eq. det_up) .and. (new_dets_dn(i) .eq. det_dn)) norm=norm+phases(i)
 enddo

 ! 4 fold first followed by  Lattice reflection and then followed by spin inversion
 do i=13,16
    new_dets_up(i)=new_dets_dn(i-4)
    up_locations(i,:)=dn_locations(i-4,:)
    new_dets_dn(i)=new_dets_up(i-4)
    dn_locations(i,:)=up_locations(i-4,:)
    phases(i)=z*phases(i-4)
    if ((new_dets_up(i) .eq. det_up) .and. (new_dets_dn(i) .eq. det_dn)) norm=norm+phases(i)
 enddo

 do i=1,16
    if (new_dets_up(i) .lt. rep_up) then
        rep_up=new_dets_up(i)
        rep_dn=new_dets_dn(i)
        phase_w_rep=phases(i)
    elseif (new_dets_up(i) .eq. rep_up) then
        if (rep_dn .gt. new_dets_dn(i)) then
            phase_w_rep=phases(i)
            rep_dn=new_dets_dn(i)
        endif
    endif
 enddo

 if (norm .gt. 1.0e-10) then
    num_distinct=nint(16._rk/norm)
 endif

 ! Now evaluate norm
 norm=sqrt(abs(norm))

 end subroutine generate_fourfold_k_configs_efficient

!=================================================================================================================================================================================
 subroutine generate_fourfold_k_configs_efficient_given_locations(c4_map,reflection_map,z,p,nup,ndn,det_up,det_dn,given_up_locations,given_dn_locations,new_dets_up,new_dets_dn,phases,rep_up,rep_dn,phase_w_rep,norm,num_distinct)
!=================================================================================================================================================================================
 ! Created by : Hitesh Changlani, April 2,2012
 ! Purpose    : Generate all four fold rotated+ spatial inversion + time reversed configurations
 !              Also return the norm of the configuration as given by group theoretical considerations

 implicit none
 ! Dummy
 integer,intent(in)       :: c4_map(:,:),reflection_map(:)              ! C_4 and Reflection Maps
 integer,intent(in)       :: nup,ndn
 integer,intent(in)       :: z,p                                        ! z and p are quantum numbers for time reversal (spin inversion) and parity (spatial reflection)
#ifdef NUM_ORBITALS_GT_127
 type(ik_vec),intent(in)   :: det_up,det_dn                              ! Incoming representative determinant up and down
 type(ik_vec),intent(out)  :: new_dets_up(16),new_dets_dn(16)            ! Symmetry related determinants
 type(ik_vec),intent(out)  :: rep_up,rep_dn                              ! Representatives
#else
 integer(ik),intent(in)   :: det_up,det_dn                              ! Incoming representative determinant up and down
 integer(ik),intent(out)  :: new_dets_up(16),new_dets_dn(16)            ! Symmetry related determinants
 integer(ik),intent(out)  :: rep_up,rep_dn                              ! Representatives
#endif
 real(rk),intent(out)     :: phases(16)                                 ! Phase factors in the linear combination
 real(rk),intent(out)     :: phase_w_rep,norm                           ! Norm
 integer,intent(out)      :: num_distinct                               ! Number of distinct configurations
 integer,intent(in)       :: given_up_locations(:),given_dn_locations(:)

 ! Local
 integer                  :: i,counter
 integer                  :: up_locations(16,nup),dn_locations(16,ndn)  ! Save locations of electrons

 ! Norm computation only true for Q=0 group which has G=A1 symmetry
 num_distinct=0
 phase_w_rep=1._rk

 rep_up=det_up
 rep_dn=det_dn

 ! We call lattice reflection as "parity", a convention used often

 ! Map i - C^i  need only i=1,2,3 i = 0 and 4 lead to identity transformation
 phases(1)=1
 new_dets_up(1)=det_up
 new_dets_dn(1)=det_dn
 norm=1._rk

 up_locations(1,:)=given_up_locations(:)
 dn_locations(1,:)=given_dn_locations(:)

 counter=0
 i=0

 do i=1,3
    call relabel_efficient(c4_map(i,:),up_locations(1,:),dn_locations(1,:),new_dets_up(i+1),new_dets_dn(i+1),up_locations(i+1,:),dn_locations(i+1,:),phases(i+1))
    if ((new_dets_up(i+1) .eq. det_up) .and. (new_dets_dn(i+1) .eq. det_dn)) norm=norm+phases(i+1)
 enddo

 ! 4 fold + Spin inversion
 do i=5,8
    new_dets_up(i)=new_dets_dn(i-4)
    up_locations(i,:)=dn_locations(i-4,:)
    new_dets_dn(i)=new_dets_up(i-4)
    dn_locations(i,:)=up_locations(i-4,:)
    phases(i)=z*phases(i-4)
    if ((new_dets_up(i) .eq. det_up) .and. (new_dets_dn(i) .eq. det_dn)) norm=norm+phases(i)
 enddo

 ! Lattice reflection around y=-x axis
 do i=9,12
    call relabel_efficient(reflection_map(:),up_locations(i-8,:),dn_locations(i-8,:),new_dets_up(i),new_dets_dn(i),up_locations(i,:),dn_locations(i,:),phases(i))
    phases(i)=phases(i)*p*phases(i-8)
    if ((new_dets_up(i) .eq. det_up) .and. (new_dets_dn(i) .eq. det_dn)) norm=norm+phases(i)
 enddo

 ! 4 fold first followed by  Lattice reflection and then followed by spin inversion
 do i=13,16
    new_dets_up(i)=new_dets_dn(i-4)
    up_locations(i,:)=dn_locations(i-4,:)
    new_dets_dn(i)=new_dets_up(i-4)
    dn_locations(i,:)=up_locations(i-4,:)
    phases(i)=z*phases(i-4)
    if ((new_dets_up(i) .eq. det_up) .and. (new_dets_dn(i) .eq. det_dn)) norm=norm+phases(i)
 enddo

 do i=1,16
    if (new_dets_up(i) .lt. rep_up) then
        rep_up=new_dets_up(i)
        rep_dn=new_dets_dn(i)
        phase_w_rep=phases(i)
    elseif (new_dets_up(i) .eq. rep_up) then
        if (rep_dn .gt. new_dets_dn(i)) then
            phase_w_rep=phases(i)
            rep_dn=new_dets_dn(i)
        endif
    endif
 enddo

 if (norm .gt. 1.0e-10) then
    num_distinct=nint(16._rk/norm)
 endif

 ! Now evaluate norm
 norm=sqrt(abs(norm))

 end subroutine generate_fourfold_k_configs_efficient_given_locations

!==========================================================================================================================================================
 subroutine generate_fourfold_k_configs(c4_map,reflection_map,z,p,det_up,det_dn,new_dets_up,new_dets_dn,phases,rep_up,rep_dn,phase_w_rep,norm,num_distinct)
!==========================================================================================================================================================
 ! Created by : Hitesh Changlani, April 2,2012
 ! Purpose    : Generate all four fold rotated+ spatial inversion + time reversed configurations
 !              Also return the norm of the configuration as given by group theoretical considerations

 implicit none
 ! Dummy
 integer,intent(in)       :: c4_map(:,:),reflection_map(:)              ! C_4 and Reflection Maps
 integer,intent(in)       :: z,p                                        ! z and p are quantum numbers for time reversal (spin inversion) and parity (spatial reflection)
#ifdef NUM_ORBITALS_GT_127
 type(ik_vec),intent(in)   :: det_up,det_dn                              ! Incoming representative determinant up and down
 type(ik_vec),intent(out)  :: new_dets_up(16),new_dets_dn(16)            ! Symmetry related determinants
 type(ik_vec),intent(out)  :: rep_up,rep_dn                              ! Representatives
#else
 integer(ik),intent(in)   :: det_up,det_dn                              ! Incoming representative determinant up and down
 integer(ik),intent(out)  :: new_dets_up(16),new_dets_dn(16)            ! Symmetry related determinants
 integer(ik),intent(out)  :: rep_up,rep_dn                              ! Representatives
#endif
 real(rk),intent(out)     :: phases(16)                                 ! Phase factors in the linear combination
 real(rk),intent(out)     :: phase_w_rep,norm                           ! Norm
 integer,intent(out)      :: num_distinct                               ! Number of distinct configurations
 ! Local
 integer                  :: i

 ! Norm computation only true for Q=0 group which has G=A1 symmetry
 num_distinct=0
 phase_w_rep=1._rk

 rep_up=det_up
 rep_dn=det_dn

 ! We call lattice reflection as "parity", a convention used often

 ! Map i - C^i  need only i=1,2,3 i = 0 and 4 lead to identity transformation
 phases(1)=1
 new_dets_up(1)=det_up
 new_dets_dn(1)=det_dn
 norm=1._rk
 do i=1,3
    call relabel(c4_map(i,:),det_up,det_dn,new_dets_up(i+1),new_dets_dn(i+1),phases(i+1))
    if ((new_dets_up(i+1) .eq. det_up) .and. (new_dets_dn(i+1) .eq. det_dn)) norm=norm+phases(i+1)
 enddo

 ! 4 fold + Spin inversion
 do i=5,8
    new_dets_up(i)=new_dets_dn(i-4)
    new_dets_dn(i)=new_dets_up(i-4)
    phases(i)=z*phases(i-4)
    if ((new_dets_up(i) .eq. det_up) .and. (new_dets_dn(i) .eq. det_dn)) norm=norm+phases(i)
 enddo

 ! Lattice reflection around y=-x axis
 do i=9,12
    call relabel(reflection_map(:),new_dets_up(i-8),new_dets_dn(i-8),new_dets_up(i),new_dets_dn(i),phases(i))
    phases(i)=phases(i)*p*phases(i-8)
    if ((new_dets_up(i) .eq. det_up) .and. (new_dets_dn(i) .eq. det_dn)) norm=norm+phases(i)
 enddo

 ! 4 fold first followed by  Lattice reflection and then followed by spin inversion
 do i=13,16
    new_dets_up(i)=new_dets_dn(i-4)
    new_dets_dn(i)=new_dets_up(i-4)
    phases(i)=z*phases(i-4)
    if ((new_dets_up(i) .eq. det_up) .and. (new_dets_dn(i) .eq. det_dn)) norm=norm+phases(i)
 enddo

 do i=1,16
    if (new_dets_up(i) .lt. rep_up) then
        rep_up=new_dets_up(i)
        rep_dn=new_dets_dn(i)
        phase_w_rep=phases(i)
    elseif (new_dets_up(i) .eq. rep_up) then
        if (rep_dn .gt. new_dets_dn(i)) then
            phase_w_rep=phases(i)
            rep_dn=new_dets_dn(i)
        endif
    endif
 enddo

 if (norm .gt. 1.0e-10) then
    num_distinct=nint(16._rk/norm)
    !if (phases(1).lt.0) phases(:)=phases(:)*-1
 endif

 ! Now evaluate norm
 norm=sqrt(abs(norm))

 end subroutine generate_fourfold_k_configs

!====================================================================================================================================================================================================================
  subroutine linear_search_list_and_update(search_for_up,search_for_dn,walk_wt,search_in_up,search_in_dn,psi_t_connected_e_loc_num,psi_t_connected_e_loc_den,e_num_gen,e_den_gen,e_num2_cum,e_den2_cum,e_num_abs_cum,e_den_abs_cum,e_num_e_den_cum,e_num_walker,e_den_walker)
!====================================================================================================================================================================================================================

    ! Searches for all of the elements in search_for_up,dn (assumed to be sorted) in search_in_up,dn, and updates all the running sums when it finds them.
    ! Scales as size(search_for_up)+size(search_in_up)
    ! A. Holmes, 27 Mar 2012
    ! Modified : A Holmes, 14 Nov 2012. Added importance sampling with psi_g(psi_t_connected) = e_loc/e_trial.

    use common_psi_t, only : hf_to_psit

    implicit none

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: search_for_up(:),search_for_dn(:) ! List of dets (occupied ones) we are searching for (assumed to be in order)
    type(ik_vec),intent(in) :: search_in_up(:),search_in_dn(:) ! List of dets (with stored local energy pieces) we are searching in (assumed to be in order)
#else
    integer(ik),intent(in) :: search_for_up(:),search_for_dn(:) ! List of dets (occupied ones) we are searching for (assumed to be in order)
    integer(ik),intent(in) :: search_in_up(:),search_in_dn(:) ! List of dets (with stored local energy pieces) we are searching in (assumed to be in order)
#endif
    real(rk),intent(in) :: walk_wt(:)
    real(rk),intent(in) :: psi_t_connected_e_loc_num(:),psi_t_connected_e_loc_den(:)
    real(rk),intent(inout) :: e_num_walker(:),e_den_walker(:)
    real(rk),intent(inout) :: e_num_gen,e_den_gen,e_num2_cum,e_den2_cum,e_num_abs_cum,e_den_abs_cum,e_num_e_den_cum

    integer :: i,j,k
    real(rk) :: e_num,e_den

    k=1

    do i=1,size(search_for_up)
      if (e_num_walker(i)>1.e50_rk) then
        do j=k,size(search_in_up)
          if (search_for_up(i).eq.search_in_up(j).and.search_for_dn(i).eq.search_in_dn(j)) then

            if (hf_to_psit) then
              if (j==1) then
                if(i/=1) stop 'in linear_search_list_and_update i/=1 when j==1'
                e_num_walker(i)=psi_t_connected_e_loc_num(1)/psi_t_connected_e_loc_den(1)
                e_den_walker(i)=walk_wt(1)
              else
                e_num_walker(i)=psi_t_connected_e_loc_num(j)
                e_den_walker(i)=psi_t_connected_e_loc_den(j)
              endif
            else
              e_num_walker(i)=psi_t_connected_e_loc_num(j)
              e_den_walker(i)=psi_t_connected_e_loc_den(j)
            endif

            k=j+1
            exit

          elseif (abs(search_for_up(i))<abs(search_in_up(j)).or.(search_for_up(i)==search_in_up(j).and.abs(search_for_dn(i))<abs(search_in_dn(j)))) then

            e_num_walker(i) = 0._rk
            e_den_walker(i) = 0._rk
            k=j
            exit

          endif
        enddo ! end of linear search
        if (e_num_walker(i)>1.e50_rk) then ! This is needed if the last walker is not in psit_connected!
          e_num_walker(i) = 0._rk
          e_den_walker(i) = 0._rk
        endif
      endif ! use stored e_num,e_den

      e_num=e_num_walker(i)*walk_wt(i)
      e_den=e_den_walker(i)*walk_wt(i)

      if (e_num.ne.0._rk) then
        if (abs(e_den).lt.1.d-22)  e_den=abs(e_den) ! This is necessary because the sign of e_num depends on the sign of e_den, which can be either +0 or -0
        e_den_gen=e_den_gen+e_den
        e_den2_cum=e_den2_cum+e_den**2
        e_den_abs_cum=e_den_abs_cum+abs(e_den)
        e_num_gen=e_num_gen+e_num
        e_num2_cum=e_num2_cum+e_num**2
        e_num_abs_cum=e_num_abs_cum+e_num*sign(1._rk,e_den)
        e_num_e_den_cum=e_num_e_den_cum+e_num*e_den
!       write(6,'(''i,j,e_den,e_num,walk_wt(i)'',2i6,9es14.6)') i,j,e_den,e_num,walk_wt(i),e_num_walker(i),e_den_walker(i)
      endif

    enddo

  end subroutine linear_search_list_and_update

!====================================================================================================================================================================================================================
  subroutine linear_search_list(search_for_up,search_for_dn,search_in_up,search_in_dn,indices)
!====================================================================================================================================================================================================================

    ! Searches for all of the elements in search_for_up,dn (assumed to be sorted) in search_in_up,dn
    ! Returns indices, a vector of size(search_for_up), with
    ! elements equal to the indices in search_in_up,dn, if they
    ! are in that vector, and 0 otherwise
    ! Scales as size(search_for_up)+size(search_in_up)
    ! A Holmes, 5 Feb 2015

    implicit none

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: search_for_up(:),search_for_dn(:) ! List of dets (occupied ones) we are searching for (assumed to be in order)
    type(ik_vec),intent(in) :: search_in_up(:),search_in_dn(:) ! List of dets (with stored local energy pieces) we are searching in (assumed to be in order)
#else
    integer(ik),intent(in) :: search_for_up(:),search_for_dn(:) ! List of dets (occupied ones) we are searching for (assumed to be in order)
    integer(ik),intent(in) :: search_in_up(:),search_in_dn(:) ! List of dets (with stored local energy pieces) we are searching in (assumed to be in order)
#endif
    integer,intent(out) :: indices(:) ! Assumed to be size(search_for_up)

    integer :: i,j,k

    k=1
    indices(:) = 0

    do i=1,size(search_for_up)
      do j=k,size(search_in_up)
        if (search_for_up(i).eq.search_in_up(j).and.search_for_dn(i).eq.search_in_dn(j)) then
          indices(i) = j
          k=j+1
          exit
        elseif (abs(search_for_up(i))<abs(search_in_up(j)).or.(search_for_up(i)==search_in_up(j).and.abs(search_for_dn(i))<abs(search_in_dn(j)))) then
          k=j
          exit
        endif
      enddo
    enddo

  end subroutine linear_search_list

!====================================================================================================================================================================================================================
  subroutine linear_search(search_for_up,search_for_dn,search_in_up,search_in_dn,indices)
!====================================================================================================================================================================================================================

    ! Searches for all of the elements in search_for_up,dn (assumed to be sorted) in search_in_up,dn, and updates all the running sums when it finds them.
    ! Scales as size(search_for_up)+size(search_in_up)
    ! Not presently being used.
    ! A. Holmes, 27 Mar 2012

    implicit none

#ifdef NUM_ORBITALS_GT_127
    type(ik_vec),intent(in) :: search_for_up,search_for_dn ! List of dets we are searching for (assumed to be in order)
    type(ik_vec),intent(in) :: search_in_up(:),search_in_dn(:) ! List of dets we are searching in (assumed to be in order)
#else
    integer(ik),intent(in) :: search_for_up,search_for_dn ! List of dets we are searching for (assumed to be in order)
    integer(ik),intent(in) :: search_in_up(:),search_in_dn(:) ! List of dets we are searching in (assumed to be in order)
#endif
    integer,intent(out) :: indices ! Locations of the search_for elements in the search_in list; 0 if they are not in that list

    integer :: j,k

    k=1

!   do i=1,size(search_for_up)
      do j=k,size(search_in_up)
        if (search_for_up.eq.search_in_up(j).and.search_for_dn.eq.search_in_dn(j)) then
       !if (search_for_up(i).eq.search_in_up(j).and.search_for_dn(i).eq.search_in_dn(j)) then
          indices=j
          return
         !indices(i)=j
          k=j+1
          exit
        endif
      enddo
      indices=0
!       elseif (search_for_up<search_in_up(j).or.(search_for_up==search_in_up(j).and.search_for_dn(i)<search_in_dn(j))) then
!      !elseif (search_for_up(i)<search_in_up(j).or.(search_for_up(i)==search_in_up(j).and.search_for_dn(i)<search_in_dn(j))) then
!         indices=0
!        !indices(i)=0
!         k=j
!         exit
!       endif
!     enddo
!   enddo

  end subroutine linear_search

!!====================================================================================================================================================================================================================
!  function check_sorted_by_label(dets_up,dets_dn)
!!====================================================================================================================================================================================================================
!
!    ! Checks whether dets_up,dets_dn sorted by label.
!    ! A Holmes, 4 Sep 2012.
!
!    implicit none
!
!    type(ik_vec),intent(in) :: dets_up(:),dets_dn(:)
!    logical :: check_sorted_by_label
!
!    integer :: n
!    type(ik_vec),allocatable :: diff_up(:),diff_dn(:)
!
!    n = size(dets_up)
!    allocate(diff_up(n-1))
!
!    diff_up(1:n-1) = dets_up(1:n-1) - dets_up(2:n)
!
!    if (min(diff_up)<0) then
!      ! out of order
!      return
!    endif
!
!    diff_dn(1:n-1) = dets_dn(1:n-1) - dets_dn(2:n)
!
!
!  end function check_sorted_by_label
!

  subroutine from_upper_triangular_to_square(n,tri_nonzero_elems,tri_indices,tri_values,square_nonzero_elems,square_indices,square_values,num_tri_elems,num_square_elems)
    ! Convert a sparse matrix from upper triangular form to square matrix form
    ! A Holmes, 31 Jan 2014

    use types, only : i8b

    integer,intent(in) :: n
    integer(i8b),intent(in) :: tri_nonzero_elems(:),tri_indices(:)
    real(rk),intent(in) :: tri_values(:)
    integer(i8b),intent(out) :: square_nonzero_elems(:),square_indices(:)
    real(rk),intent(out) :: square_values(:)  ! Assume that indices, values are of size 2*size(tri_indices)-n
    integer(i8b),intent(in) :: num_tri_elems
    integer(i8b),intent(out) :: num_square_elems

    integer :: i,row
    integer :: filled(num_tri_elems)
    integer(i8b) :: cumulative(num_tri_elems)
    integer(i8b) :: tri_cumulative(num_tri_elems)

    num_square_elems = 0

    square_nonzero_elems(:) = tri_nonzero_elems(:) - 1
    do i=1, num_tri_elems
      square_nonzero_elems(tri_indices(i)) = square_nonzero_elems(tri_indices(i)) + 1
    enddo

    tri_cumulative(1) = 0
    do i=2,n
      tri_cumulative(i) = tri_cumulative(i-1) + tri_nonzero_elems(i-1)
    enddo

    cumulative(1) = 0
    do i=2,n
      cumulative(i) = cumulative(i-1) + square_nonzero_elems(i-1)
    enddo

    filled(:) = 0
    row = 1

    do row=1,n
      do i=tri_cumulative(row)+1,tri_cumulative(row)+tri_nonzero_elems(row)
        ! Fill in lower triangular part (including diagonal)
        filled(tri_indices(i)) = filled(tri_indices(i)) + 1
        square_indices(cumulative(tri_indices(i))+filled(tri_indices(i))) = row
        square_values(cumulative(tri_indices(i))+filled(tri_indices(i))) = tri_values(i)
        num_square_elems = num_square_elems + 1
        ! Fill in upper triangular part (not including diagonal)
        if (row/=tri_indices(i)) then
          filled(row) = filled(row) + 1
          square_indices(cumulative(row)+filled(row)) = tri_indices(i)
          square_values(cumulative(row)+filled(row)) = tri_values(i)
          num_square_elems = num_square_elems + 1
        endif
      enddo
    enddo

   !n_tot_nonzero = cumulative(n) + square_nonzero_elems(n)

  end subroutine from_upper_triangular_to_square


  subroutine from_square_to_upper_triangular(n,square_nonzero_elems,square_indices,square_values,tri_nonzero_elems,tri_indices,tri_values,num_square_elems,num_tri_elems)
    ! Convert a sparse matrix from square matrix form to upper triangular form
    ! A Holmes, 31 Jan 2014

    use types, only : i8b

    integer,intent(in) :: n
    integer(i8b),intent(in) :: square_nonzero_elems(:),square_indices(:)
    real(rk),intent(in) :: square_values(:)  ! Assume that indices, values are of size 2*size(tri_indices)-n
    integer(i8b),intent(out) :: tri_nonzero_elems(:),tri_indices(:)
    real(rk),intent(out) :: tri_values(:)
    integer,intent(in) :: num_square_elems
    integer,intent(out) :: num_tri_elems

    integer :: i,j,jold,write_tri,row
    integer :: cumulative(num_square_elems) !size(square_nonzero_elems))

    cumulative(1) = 0
    do i=2,n
      cumulative(i) = cumulative(i-1) + square_nonzero_elems(i-1)
    enddo

    ! j counts number of upper triangular elems on a row, jold is the value of j at the end of previous row
    j = 0
    jold = 0

    tri_nonzero_elems(n) = 1
    row = 1
    do row=1,n
      do i=cumulative(row)+1,cumulative(row)+square_nonzero_elems(row)
        if (i>1) then
          if (i-1==cumulative(row)) then
            tri_nonzero_elems(row-1) = j-jold
            if (write_tri==1) then
              write_tri = 0
              jold = j
            endif
          endif
        endif
        if (row==square_indices(i)) then
          write_tri = 1
        endif
        if (write_tri==1) then
          j = j+1
          tri_indices(j) = square_indices(i)
          tri_values(j) = square_values(i)
        endif
      enddo
    enddo
    num_tri_elems = j

  end subroutine from_square_to_upper_triangular


  subroutine from_square_to_band(n,square_nonzero_elems,square_indices,square_values,band_nonzero_elems,band_indices,band_values,num_square_elems,num_band_elems,index_min,index_max)
    ! Convert a sparse mabandx from square mabandx form to upper bandangular form
    ! A Holmes, 18 Mar 2014

    use types, only : i8b

    integer,intent(in) :: n
    integer(i8b),intent(in) :: square_nonzero_elems(:),square_indices(:)
    real(rk),intent(in) :: square_values(:)  ! Assume that indices, values are of size 2*size(band_indices)-n
    integer(i8b),intent(out) :: band_nonzero_elems(:),band_indices(:)
    real(rk),intent(out) :: band_values(:)
    integer(i8b),intent(in) :: num_square_elems
    integer(i8b),intent(out) :: num_band_elems
    integer(i8b),intent(in) :: index_min,index_max

    integer :: i,j,jold,row
    integer :: cumulative(num_square_elems)

    cumulative(1) = 0
    do i=2,n
      cumulative(i) = cumulative(i-1) + square_nonzero_elems(i-1)
    enddo

    ! j counts number of upper band elems on a row, jold is the value of j at the end of previous row
    j = 0
    jold = 0

    row = 1
    do row=1,n
      do i=cumulative(row)+1,cumulative(row)+square_nonzero_elems(row)
        if (i==cumulative(row)+1.and.row>1) then
          band_nonzero_elems(row-1) = j-jold
          jold = j
        endif
        if (square_indices(i)>=index_min.and.square_indices(i)<=index_max) then
          j = j+1
          band_indices(j) = square_indices(i)-index_min+1
          band_values(j) = square_values(i)
        endif
      enddo
    enddo
    band_nonzero_elems(n) = j-jold
    num_band_elems = j

  end subroutine from_square_to_band


  subroutine from_sparse_square_to_dense_square(n,square_nonzero,square_indices,square_values)
    ! Just prints out as a square matrix
    ! Not presently being used
    ! A Holmes, 17 Mar 2014

    integer,intent(in) :: n,square_indices(:),square_nonzero(:)
    real(rk),intent(in) :: square_values(:)

    integer :: i,j,ind
    real(rk),allocatable :: mat(:,:)

    allocate(mat(n,n))
    mat(:,:) = 0.0_rk
    ind = 0
    do i=1,n
      do j=1,square_nonzero(i)
        ind = ind + 1
        mat(i,square_indices(ind)) = square_values(ind)
        mat(square_indices(ind),i) = square_values(ind)
      enddo
    enddo

    write (6,*) "Matrix="; call flush(6)
    do i=1,n
      write (6,*) mat(i,1:n)
    enddo
    call flush(6)

  end subroutine from_sparse_square_to_dense_square

  subroutine from_sparse_tri_to_dense_square(n,tri_nonzero,tri_indices,tri_values)
    ! Just prints out as a square matrix
    ! Not presently being used
    ! A Holmes, 17 Mar 2014

    integer,intent(in) :: n,tri_indices(:),tri_nonzero(:)
    real(rk),intent(in) :: tri_values(:)

    integer :: i
    real(rk),allocatable :: mat(:,:)
    integer :: row

    allocate(mat(n,n))
    mat(:,:) = 0.0_rk
    row = 1
    do i=1,size(tri_indices)
      if (i>1) then
        if (tri_indices(i)<=tri_indices(i-1).or.(tri_nonzero(row)==1)) then
          row = row + 1
        endif
      endif
      mat(row,tri_indices(i)) = tri_values(i)
      mat(tri_indices(i),row) = tri_values(i)
    enddo

    write (6,*) "Matrix="; call flush(6)
    do i=1,n
      write (6,*) mat(i,1:n)
    enddo
    call flush(6)

  end subroutine from_sparse_tri_to_dense_square

  subroutine gen_hist(lo,hi,nbins,lbounds,bins)
    ! Generate histogram
    ! A Holmes, 21 Mar 2015

    use types, only : i8b

    real(rk),intent(in) :: lo,hi
    integer,intent(in) :: nbins
    real(rk),allocatable,intent(out) :: lbounds(:)
    integer(i8b),allocatable,intent(out) :: bins(:)
    integer :: i

    allocate(lbounds(nbins))
    allocate(bins(nbins))

    bins(:) = 0
    do i=1,nbins
      lbounds(i) = lo+(i-1)*(hi-lo)/(nbins-1)
    enddo

  end subroutine gen_hist


  subroutine add_to_hist(x,nbins,lbounds,bins)
    ! Add x to a histogram structure (nbins,lbounds,bins) already allocated by gen_hist
    ! A Holmes, 21 Mar 2015

    use types, only : i8b

    real(rk),intent(in) :: x
    integer,intent(in) :: nbins
    real(rk),intent(in) :: lbounds(:)
    integer(i8b),intent(inout) :: bins(:)
    integer :: ibin

    if (x<lbounds(1))  return

    ibin=min(nbins,int(1+(nbins-1)*(x-lbounds(1))/(lbounds(nbins)-lbounds(1)))) ! don't need the safeguards on next 2 lines
    !ibin=int(min(real(nbins,rk),1+(nbins-1)*(x-lbounds(1))/(lbounds(nbins)-lbounds(1))))
    !ibin=max(1,min(nbins,int(1+(nbins-1)*(x-lbounds(1))/(lbounds(nbins)-lbounds(1)))))
    if(ibin.le.0) then
      write(6,'(''nbins,lbounds(1),lbounds(nbins),x,ibin='',i12,3es12.4,i12)') nbins,lbounds(1),lbounds(nbins),x,ibin
      stop 'ibin.le.0 means x is NaN, but that should never happen'
    endif
    bins(ibin) = bins(ibin) + 1_i8b

  end subroutine add_to_hist


  subroutine get_occ_orbs_1det(det,occ)
  ! Get lists of which orbitals are occupied in given configuration
  ! A Holmes, 21 Mar 2015

#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: det
  type(ik_vec) :: tmp_det
#else
  integer(ik),intent(in) :: det
  integer(ik) :: tmp_det
#endif
 !integer,intent(in) :: norb
  integer,intent(out) :: occ(:)

  integer :: i_elec,i

    tmp_det = det
    i_elec = 0
    do while (tmp_det .ne. 0)
      i = trailz(tmp_det)+1
      i_elec = i_elec + 1
      occ(i_elec) = i
      tmp_det = ibclr(tmp_det,i-1)
    enddo

  end subroutine get_occ_orbs_1det

  subroutine get_occ_orbs_2dets(det_up,det_dn,occ_up,occ_dn)
  ! Get lists of which orbitals are occupied in given configuration
  ! A Holmes, 21 Mar 2015

#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: det_up,det_dn
  type(ik_vec) :: tmp_det
#else
  integer(ik),intent(in) :: det_up,det_dn
  integer(ik) :: tmp_det
#endif
 !integer,intent(in) :: norb
  integer,intent(out) :: occ_up(:),occ_dn(:)

  integer :: i_elec,i

    tmp_det = det_up
    i_elec = 0
    do while (tmp_det .ne. 0)
      i = trailz(tmp_det)+1
      i_elec = i_elec + 1
      occ_up(i_elec) = i
      tmp_det = ibclr(tmp_det,i-1)
    enddo

    tmp_det = det_dn
    i_elec = 0
    do while (tmp_det .ne. 0)
      i = trailz(tmp_det)+1
      i_elec = i_elec + 1
      occ_dn(i_elec) = i
      tmp_det = ibclr(tmp_det,i-1)
    enddo

  end subroutine get_occ_orbs_2dets


  subroutine get_occ_and_unocc_orbs(det_up,det_dn,norb,occ_up,occ_dn,unocc_up,unocc_dn)
  ! Get lists of which orbitals are occupied and unoccupied in given configuration
  ! A Holmes, 21 Mar 2015

#ifdef NUM_ORBITALS_GT_127
  type(ik_vec),intent(in) :: det_up,det_dn
#else
  integer(ik),intent(in) :: det_up,det_dn
#endif
  integer,intent(in) :: norb
  integer,intent(out) :: occ_up(:),occ_dn(:),unocc_up(:),unocc_dn(:)

  integer :: ne,nh,i

    ne = 0
    nh = 0
    do i=1,norb
      if (btest(det_up,i-1)) then
        ne = ne+1
        occ_up(ne) = i
      else
        nh = nh+1
        unocc_up(nh) = i
      endif
    enddo

    ne = 0
    nh = 0
    do i=1,norb
      if (btest(det_dn,i-1)) then
        ne = ne+1
        occ_dn(ne) = i
      else
        nh = nh+1
        unocc_dn(nh) = i
      endif
    enddo

  end subroutine get_occ_and_unocc_orbs


  subroutine setup_alias_rk(K,pdf,smaller,larger,J,q)
    ! Setup the alias method for efficient sampling of a discrete PDF
    ! Setup method takes O(K) time, but enables sampling in O(1) time
    ! Translated from Python code on https://hips.seas.harvard.edu/blog/2013/03/03/the-alias-method-efficient-sampling-with-many-discrete-outcomes/
    ! A Holmes, 24 Mar 2015
    ! Edited : A Holmes, 29 Sep 2015. Initialize J(i) = i, rather than
    ! J(i) = 1, so that if the random number exceeds a floating point
    ! approximation to P(i) = 1, i is returned (as it should be) rather
    ! than 1.

    integer,intent(in) :: K ! Number of discrete outcomes
    real(rk),intent(in) :: pdf(:) ! Normalized probability distribution
    integer,intent(inout) :: smaller(:),larger(:) ! Scratch variables; pass them in to avoid allocating them

    ! Output alias sampling quantities
    integer,intent(out) :: J(:)
    real(rk),intent(out) :: q(:)

    integer :: small,large
    integer :: n_s,n_l,i

      q(:) = 0._rk

      ! Sort data into outcomes with probability larger and smaller than 1/K

      smaller(:) = 0
      larger(:) = 0

      n_s = 0
      n_l = 0

      do i=1,K
        J(i) = i
        q(i) = K*pdf(i)
        if (q(i)<1.0) then
          n_s = n_s + 1
          smaller(n_s) = i
        else
          n_l = n_l + 1
          larger(n_l) = i
        endif
      enddo

      ! Loop through and create little binary mixtures that appropriately
      ! allocate the larger outcomes over the overall uniform mixture.

      do while (n_s > 0 .and. n_l > 0)
        small = smaller(n_s)
        large = larger(n_l)
        J(small) = large
        q(large) = q(large) + q(small) - 1._rk
        if (q(large) < 1._rk) then
          smaller(n_s) = large
          n_l = n_l - 1
        else
          n_s = n_s - 1
        endif
      enddo

  end subroutine setup_alias_rk


  subroutine setup_alias_real(K,pdf,smaller,larger,J,q)
    ! Setup the alias method for efficient sampling of a discrete PDF
    ! Setup method takes O(K) time, but enables sampling in O(1) time
    ! Translated from Python code on https://hips.seas.harvard.edu/blog/2013/03/03/the-alias-method-efficient-sampling-with-many-discrete-outcomes/
    ! A Holmes, 24 Mar 2015
    ! Edited : A Holmes, 29 Sep 2015. Initialize J(i) = i, rather than
    ! J(i) = 1, so that if the random number exceeds a floating point
    ! approximation to P(i) = 1, i is returned (as it should be) rather
    ! than 1.

    integer,intent(in) :: K ! Number of discrete outcomes
    real,intent(in) :: pdf(:) ! Normalized probability distribution
    integer,intent(inout) :: smaller(:),larger(:) ! Scratch variables; pass them in to avoid allocating them

    ! Output alias sampling quantities
    integer,intent(out) :: J(:)
    real,intent(out) :: q(:)

    integer :: small,large
    integer :: n_s,n_l,i

      q(:) = 0._rk

      ! Sort data into outcomes with probability larger and smaller than 1/K

      smaller(:) = 0
      larger(:) = 0

      n_s = 0
      n_l = 0

      do i=1,K
        J(i) = i
        q(i) = K*pdf(i)
        if (q(i)<1.0) then
          n_s = n_s + 1
          smaller(n_s) = i
        else
          n_l = n_l + 1
          larger(n_l) = i
        endif
      enddo

      ! Loop through and create little binary mixtures that appropriately
      ! allocate the larger outcomes over the overall uniform mixture.

      do while (n_s > 0 .and. n_l > 0)
        small = smaller(n_s)
        large = larger(n_l)
        J(small) = large
        q(large) = q(large) + q(small) - 1.0
        if (q(large) < 1._rk) then
          smaller(n_s) = large
          n_l = n_l - 1
        else
          n_s = n_s - 1
        endif
      enddo

  end subroutine setup_alias_real


  integer function sample_alias_rk(K,J,q)
    ! Use alias sampling to sample a PDF setup using setup_alias routine
    ! Samples in O(1) time
    ! Translated from Python code on https://hips.seas.harvard.edu/blog/2013/03/03/the-alias-method-efficient-sampling-with-many-discrete-outcomes/
    ! A Holmes, 24 Mar 2015

    use tools, only : random_int

    integer,intent(in) :: K ! Number of discrete outcomes

    ! Alias sampling quantities
    integer,intent(in) :: J(:)
    real(rk),intent(in) :: q(:)

    integer :: i
    real(rk) :: rannyu

      i = random_int(K)

      if (rannyu() < q(i)) then
        sample_alias_rk = i
      else
        sample_alias_rk = J(i)
      endif

  end function sample_alias_rk


  integer function sample_alias_real(K,J,q)
    ! Use alias sampling to sample a PDF setup using setup_alias routine
    ! Samples in O(1) time
    ! Translated from Python code on https://hips.seas.harvard.edu/blog/2013/03/03/the-alias-method-efficient-sampling-with-many-discrete-outcomes/
    ! A Holmes, 24 Mar 2015

    use tools, only : random_int

    integer,intent(in) :: K ! Number of discrete outcomes

    ! Alias sampling quantities
    integer,intent(in) :: J(:)
    real,intent(in) :: q(:)

    integer :: i
    real(rk) :: rannyu

      i = random_int(K)

      if (rannyu() < q(i)) then
        sample_alias_real = i
      else
        sample_alias_real = J(i)
      endif

  end function sample_alias_real


end module more_tools
