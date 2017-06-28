# This Makefile is to be used with gmake.  On Linux systems, gmake is
# the default make.

system := $(shell uname)
hostname := $(shell uname -n)
username := $(shell id -u -n)

# Libraries --------------------------------------------------------------

#LIBHOME=/home/cyrus/champ/lib
#LAPACK=-L$(LIBHOME)/lib2/lapack -llapack
#BLAS=-L$(LIBHOME)/lib2/blas -lblas
LIBS = -llapack -lblas #$(LAPACK) $(BLAS)
#LIBS =	-llapack-3

# For computer without MPI use gfortran or ifort (if present), and comment out -DMPI.
# For computer with MPI use mpif90
 FC=gfortran
#FC=ifort
#FC = mpif90
##FC = /data/alessandro/lib/bin/mpif90
#FC = /home/mjo98/mpich-install/bin/mpif90
#FC = /home/mjo98/mpich-install-gcc4.9.2/bin/mpif90
 FFLAGS   = -finit-local-zero -O3 -ffixed-line-length-132 -Wall #-DMAC
 F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall #-DMAC # -DNUM_ORBITALS_GT_127

ifeq ($(system),Darwin)
 FFLAGS   = -finit-local-zero -O3 -ffixed-line-length-132 -Wall -DMAC
 F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -DMAC # -DNUM_ORBITALS_GT_127
endif

ifeq ($(system),Linux)

# ifort flags
# -assume realloc_lhs is needed for automatic reallocation
#FFLAGS=-zero -r8 -pad -align all -fpp -ip -O3 -mtune=core2 -march=core2 -xT -extend_source -assume realloc_lhs
#F90FLAGS= -zero -r8 -pad -align all -O3 -fpp -ip -assume realloc_lhs
#F90FLAGS= -zero -r8 -pad -align all -O3 -fpp -pg -assume realloc_lhs
#F90FLAGS=-zero -r8 -pad -align all -O0 -fpp -ip -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG -assume realloc_lhs

# gfortran flags
# For debugging
#FFLAGS =   -finit-local-zero -O3 -ffixed-line-length-132 -Wall -g -fbounds-check -fbacktrace
#F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -g -fbounds-check -fbacktrace -DDEBUG # -DNUM_ORBITALS_GT_127
# For debugging and profiling
#FFLAGS =   -finit-local-zero -O3 -ffixed-line-length-132 -Wall -g -pg -fbounds-check -fbacktrace
#F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -g -pg -fbounds-check -fbacktrace # -DNUM_ORBITALS_GT_127
# For profiling with gprof or valgrind
# FFLAGS =   -finit-local-zero -O3 -ffixed-line-length-132 -Wall -g -pg
# F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -g -pg # -DNUM_ORBITALS_GT_127
# For production
 FFLAGS   = -finit-local-zero -O3 -ffixed-line-length-132 -Wall
 F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall # -DNUM_ORBITALS_GT_127

# LDFLAGS = -L/usr/local/tools/mpip/lib -lmpiP -lm#-static

ifeq ($(findstring cooley,$(hostname)), cooley)
FC = mpif90
# Cooley cluster at ANL

# INTEL COMPILERS - code needs changing to remove integer*16
# $HOME/.soft.cooley should have the following lines:
# +intel-composer-xe
# +mvapich2-intel
#
# FFLAGS = -zero -save -extend_source -w -r8 -O3 -pad
# F90FLAGS = -zero -save -extend_source -w -r8 -O3 -pad

# GNU COMPILERS
# $HOME/.soft.cooley should have the following 3 lines:
# +mvapich2
# +gcc-4.8.1
# @default
# Then type resoft

 FFLAGS   = -O3 -ffixed-line-length-132 -Wall
 F90FLAGS = -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall # -DNUM_ORBITALS_GT_127
endif

ifeq ($(findstring cab,$(hostname)), cab)

 FC = mpif90
# Cab cluster at LLNL
#For PARPACK, use this LIBS string and -DPARPACK
#LIBS = -llapack -lblas -L/g/g19/otten1/soft/arpack-ng -lparpack -larpack #$(LAPACK) $(BLAS)
#FFLAGS   = -finit-local-zero -O3 -ffixed-line-length-132 -Wall -DMPI
#F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -DMPI # -DNUM_ORBITALS_GT_127
 FFLAGS =   -finit-local-zero -O3 -ffixed-line-length-132 -Wall -g -fbounds-check -fbacktrace
 F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -g -fbounds-check -fbacktrace -DDEBUG  -DMPI  # -DNUM_ORBITALS_GT_127
endif

ifeq ($(findstring quartz,$(hostname)), quartz)

 FC = mpif90
# Cab cluster at LLNL
#For PARPACK, use this LIBS string and -DPARPACK
#LIBS = -llapack -lblas -L/g/g19/otten1/soft/arpack-ng -lparpack -larpack #$(LAPACK) $(BLAS)
#FFLAGS   = -finit-local-zero -O3 -ffixed-line-length-132 -Wall -DMPI
#F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -DMPI # -DNUM_ORBITALS_GT_127
 FFLAGS =   -finit-local-zero -O3 -ffixed-line-length-132 -Wall -g -fbacktrace -fbounds-check
 F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -g -fbacktrace -DDEBUG  -DMPI -fbounds-check #-DNUM_ORBITALS_GT_127
endif

# planck_compiler=g77
  planck_compiler=ifort

ifeq ($(hostname),planck)
# AMD Athlon cluster of Arias group
ifeq ($(planck_compiler),g77)
# to link g77 mpi libraries, make sure /usr/bin is ahead of /usr/local/bin
# FC = g77-3.2
  FC = g77-3.3
# FC = g77
  FFLAGS = -g -C -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check
# FFLAGS = -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check
# FFLAGS = -O3 -malign-double -ffixed-line-length-none -fno-automatic
# FFLAGS = -pg -O3 -malign-double -ffixed-line-length-none -fno-automatic
  LDFLAGS = -static
  FC_MPI = mpif77 -fc=g77
# FFLAGS_MPI = -O2 -malign-double -ffixed-line-length-none -fno-automatic -Wall -fbounds-check -DMPI
  FFLAGS_MPI = -O3 -malign-double -ffixed-line-length-none -fno-automatic -DMPI
  LDFLAGS_MPI = -static
endif

ifeq ($(planck_compiler),ifort)
# to link ifort mpi libraries, make sure /usr/local/bin is ahead of /usr/bin
# In the mpif77 line, the fc=ifc tells it to use ifc rather than the default g77
# fpp2 is the preprocessor
#{i|M|K|W} i=Pentium Pro, M-Pentium MMX, K=Pentium3, W=Pentium4,Xeon
# On Pentium4W results in 3% faster execution than not using it and 2% faster execution than using -axi, but uncertainties are equally big
# -O3 seems to be about the same a O2 for ifc (v7).  Have not tested ifort.
# -w turns off warnings
# -fpp2 for preprocessor (not needed for this code)
# Richard says to use -ip for interprocess (helps a lot with VASP)
# He says VASP did not run correctly with -fast or -xW
# He says to use :ib on mpirun to get infiniband network
# He found -O3WN -ipo -pad runs fastest, but there were several choices within the timing noise.
  FC = ifort
  FFLAGS = -g -CB -d4 -extend_source -w -r8 -fpp2 -O0
# FFLAGS = -extend_source -w -r8 -fpp2 -O2
# FFLAGS = -extend_source -w -r8 -O3N -ipo -pad
# FFLAGS = -extend_source -w -r8 -O3 -pad
# LDFLAGS = -Vaxlib -static
  LDFLAGS = -static
  FC_MPI = mpif77 -fc=ifort
# FFLAGS_MPI = -g -CB -extend_source -w -r8 -O2 -DMPI
# FFLAGS_MPI = -extend_source -w -r8 -O2 -DMPI
  FFLAGS_MPI = -extend_source -w -r8 -O3 -pad -DMPI
# LDFLAGS_MPI = -Vaxlib -static -DMPI
  LDFLAGS_MPI = -static -lmpi -DMPI
endif
endif

#qmc desktops at cornell
ifeq ($(findstring qmc,$(hostname)), qmc)
  FC=gfortran
# For debugging
 FFLAGS =   -finit-local-zero -O0 -ffixed-line-length-132 -Wall -g -fbounds-check -fbacktrace
 F90FLAGS = -finit-local-zero -O0 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -g -fbounds-check -fbacktrace
# For debugging and profiling
#FFLAGS =   -finit-local-zero -O3 -ffixed-line-length-132 -Wall -g -pg -fbounds-check -fbacktrace
#F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -g -pg -fbounds-check -fbacktrace # -DNUM_ORBITALS_GT_127
# For profiling with gprof or valgrind
# FFLAGS =   -finit-local-zero -O3 -ffixed-line-length-132 -Wall -g -pg
# F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -g -pg # -DNUM_ORBITALS_GT_127
# For production
#FFLAGS   = -finit-local-zero -O3 -ffixed-line-length-132 -Wall
#F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall # -DNUM_ORBITALS_GT_127
endif

#qmc1 at cornell
ifeq ($(hostname),qmc1)
  FC=ifort
  FFLAGS=-zero -r8 -pad -align all -fpp -ip -O3 -mtune=core2 -march=core2 -xT -extend_source
  FFLAGS_DEBUG=$(FFLAGS) -check all -check nooutput_conversion -debug extended -debug-parameters all -g -traceback -error_limit 1 -DDEBUG
endif

ifeq ($(hostname),qmc8)
  FC=gfortran
endif

ifeq ($(hostname), Evanger-linux)
  FC=gfortran
 # For TAU profiler.
 # FC=tau_f90.sh
 # F90FLAGS += -opt_verbos -optTauSelectFile="select.tau"
  # F90FLAGS += -DNUM_ORBITALS_GT_127
endif

# dft cluster
#ifeq ($(findstring login,$(hostname)), login)
ifeq ($(findstring dft,$(hostname)), dft)
# FC=gfortran
# FC=mpif90
# FC=/home/mjo98/mpich-install/bin/mpif90
  FC=/home/mjo98/mpich-install-gcc4.9.2/bin/mpif90
# For debugging
#FFLAGS =   -finit-local-zero -O0 -ffixed-line-length-132 -Wall -DMPI -g -fbounds-check -fbacktrace
#F90FLAGS = -finit-local-zero -O0 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -DMPI -g -fbounds-check -fbacktrace # -DNUM_ORBITALS_GT_127
# For debugging and profiling
#FFLAGS =   -finit-local-zero -O3 -ffixed-line-length-132 -Wall -DMPI -g -pg -fbounds-check -fbacktrace
#F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -DMPI -g -pg -fbounds-check -fbacktrace # -DNUM_ORBITALS_GT_127
# For profiling with gprof or valgrind
# FFLAGS =   -finit-local-zero -O3 -ffixed-line-length-132 -Wall -DMPI -g -pg
# F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -DMPI -g -pg -DNUM_ORBITALS_GT_127
# For production
 FFLAGS   = -finit-local-zero -O3 -ffixed-line-length-132 -Wall -DMPI
 F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -DMPI # -DNUM_ORBITALS_GT_127
endif

# Bridges computer at Pittsburgh
ifeq ($(findstring br,$(hostname)), br)
  FC=mpif90
# For debugging
#FFLAGS =   -finit-local-zero -O0 -ffixed-line-length-132 -Wall -DMPI -g -fbounds-check -fbacktrace
#F90FLAGS = -finit-local-zero -O0 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -DMPI -g -fbounds-check -fbacktrace # -DNUM_ORBITALS_GT_127
# For debugging and profiling
#FFLAGS =   -finit-local-zero -O3 -ffixed-line-length-132 -Wall -DMPI -g -pg -fbounds-check -fbacktrace
#F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -DMPI -g -pg -fbounds-check -fbacktrace # -DNUM_ORBITALS_GT_127
# For profiling with gprof or valgrind
# FFLAGS =   -finit-local-zero -O3 -ffixed-line-length-132 -Wall -DMPI -g -pg
# F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -DMPI -g -pg -DNUM_ORBITALS_GT_127
# For production
 FFLAGS   = -finit-local-zero -O3 -ffixed-line-length-132 -Wall -DMPI
 F90FLAGS = -finit-local-zero -O3 -ffree-line-length-none -cpp -x f95-cpp-input -Wall -DMPI # -DNUM_ORBITALS_GT_127
endif

ifeq ($(hostname),nanolab.cnf.cornell.edu)
  FC = ifort
# FFLAGS = -g -CB -d4 -zero -save -extend_source -w -r8 -fpp2 -O0
# FFLAGS = -zero -save -extend_source -w -r8 -fpp2 -O2
# FFLAGS = -zero -save -extend_source -w -r8 -O3 -ipo -pad
  FFLAGS = -zero -save -extend_source -w -r8 -O3 -pad
  F90FLAGS = -finit-local-zero -O3  -ffree-line-length-none -cpp -x f95-cpp-input -Wall -fbounds-check
# LDFLAGS = -Vaxlib -static
# LDFLAGS = -static -nothread
  LDFLAGS = -nothread
  FC_MPI = /usr/lam_intel/bin/mpif77
# FFLAGS_MPI = -g -CB -zero -save -extend_source -w -r8 -O0 -DMPI
# FFLAGS_MPI = -g -CB -zero -save -extend_source -w -r8 -O2 -DMPI
# FFLAGS_MPI = -zero -save -extend_source -w -r8 -O2 -DMPI
  FFLAGS_MPI = -zero -save -extend_source -w -r8 -O3 -pad -DMPI
# LDFLAGS_MPI = -Vaxlib -static -DMPI
# LDFLAGS_MPI = -static -nothread -lmpi -DMPI
  LDFLAGS_MPI = -nothread -lmpi -lpthread -DMPI
endif

#ifeq ($(username), cyrus)
#  F90FLAGS += -DNUM_ORBITALS_GT_127
#endif
ifeq ($(username), junhao)
  # F90FLAGS += -DNUM_ORBITALS_GT_127
  # F90FLAGS += -DSHARE_NODE
endif

endif

.SUFFIXES:
.SUFFIXES: .f90 .f95 .o .f .c

.c.o:
	$(CC) $(CC_FLAGS) -o $@ -c $<

.f.o:
	$(FC) $(FFLAGS) -o $@ -c $<
#	$(MFC) $(FFLAGS) -o $@ -c $<

.f90.o:
	$(FC) $(F90FLAGS) -o $@ -c $<
#	$(MFC) $(F90FLAGS) -o $@ -c $<

all: sqmc

# List the dependencies.  Frank generated these using makemake.pl
chemistry.o: commons/common_ham.o commons/common_imp.o commons/common_psi_t.o \
	commons/common_run.o commons/common_selected_ci.o generic_sort.o \
	mpi_routines.o more_tools.o overload.o tools.o types.o fhash_module.o
constants.o: types.o
do_walk.o: chemistry.o commons/common_ham.o commons/common_imp.o \
	commons/common_psi_t.o commons/common_run.o commons/common_walk.o \
	hamiltonian_mod.o heg.o hubbard.o more_tools.o overload.o \
	read_psi_trial.o semistoch.o tools.o types.o
generic_sort.o: overload.o types.o
hamiltonian_mod.o: chemistry.o commons/common_ham.o commons/common_imp.o \
	commons/common_psi_t.o commons/common_run.o \
	commons/common_selected_ci.o commons/common_walk.o heg.o hubbard.o \
	hci.o read_psi_trial.o semistoch.o tools.o types.o
heg.o: commons/common_ham.o commons/common_imp.o commons/common_psi_t.o \
	commons/common_run.o constants.o generic_sort.o more_tools.o \
	overload.o tools.o types.o
hubbard.o: commons/common_ham.o commons/common_imp.o commons/common_psi_t.o \
	commons/common_run.o commons/common_selected_ci.o constants.o \
	generic_sort.o more_tools.o overload.o tools.o types.o
hci.o: chemistry.o semistoch.o
more_tools.o: commons/common_ham.o commons/common_psi_t.o generic_sort.o \
	overload.o tools.o types.o
my_second.o: types.o
my_memory.o: tools.o
overload.o: types.o
read_psi_trial.o: commons/common_ham.o commons/common_psi_t.o \
	commons/common_run.o types.o
semistoch.o: chemistry.o commons/common_ham.o commons/common_imp.o \
	commons/common_psi_t.o commons/common_run.o \
	commons/common_selected_ci.o commons/common_walk.o generic_sort.o \
	heg.o hubbard.o more_tools.o overload.o tools.o types.o
sqmc_main.o: basic_tools.o do_walk.o
tools.o: commons/common_ham.o commons/common_run.o commons/common_walk.o overload.o types.o
mpi_routines.o: types.o commons/common_walk.o tools.o
test.o: overload.o types.o


commons/common_ham.o: types.o
commons/common_imp.o: types.o
commons/common_psi_t.o: types.o
commons/common_run.o: types.o
commons/common_selected_ci.o: types.o
commons/common_walk.o: types.o

SRCS =	hci.f90 test.f90 basic_tools.f90 commons/common_walk.o mpi_routines.f90 chemistry.f90 constants.f90 do_walk.f90 \
	generic_sort.f90 hamiltonian_mod.f90 heg.f90 hubbard.f90 \
	more_tools.f90 my_second.f90 my_memory.f90 overload.f90 read_psi_trial.f90 \
	semistoch.f90 sqmc_main.f90 tools.f90 types.f90 cdet.f90 det.f90 matinv.f90 \
	rannyu.f90 commons/common_ham.f90 commons/common_imp.f90 \
	commons/common_psi_t.f90 commons/common_run.f90 \
	commons/common_selected_ci.f90 commons/common_walk.f90 fhash_module.f90


OBJS =	hci.o basic_tools.o commons/common_walk.o rannyu.o mpi_routines.o chemistry.o constants.o do_walk.o generic_sort.o \
	hamiltonian_mod.o heg.o hubbard.o more_tools.o my_second.o my_memory.o overload.o \
	read_psi_trial.o semistoch.o sqmc_main.o tools.o types.o cdet.o det.o \
	matinv.o commons/common_ham.o commons/common_imp.o \
	commons/common_psi_t.o commons/common_run.o \
	commons/common_selected_ci.o fhash_module.o

sqmc: $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) -o $@ $(MPI)

merge_sort: merge_sort.o my_second.o
	$(FC) $(FFLAGS) $(LDFLAGS) merge_sort.o my_second.o $(LIBS) -o $@

test_sort: test_sort.o my_second.o
	$(FC) $(FFLAGS) $(LDFLAGS) test_sort.o my_second.o $(LIBS) -o $@

test: test.o overload.o types.o
	$(FC) $(FFLAGS) $(LDFLAGS) test.o overload.o types.o $(LIBS) -o $@

template: template.o types.o rannyu.o shell.o
	$(FC) $(FFLAGS) $(LDFLAGS) template.o types.o rannyu.o shell.o $(LIBS) -o $@

mem_avail: mem_avail.o
	$(FC) $(FFLAGS) $(LDFLAGS) mem_avail.o -o $@

alloc: alloc.o my_second.o
	$(FC) $(FFLAGS) $(LDFLAGS) alloc.o my_second.o -o $@

e2e_test: sqmc
	cd e2e_tests/heg && ../../sqmc -i i_det > o_det && python ../e2e_check.py o_det o_det_ref
	cd e2e_tests/heg && ../../sqmc -i i_st > o_st && python ../e2e_check.py o_st o_st_ref
	cd e2e_tests/c2 && ../../sqmc -i i_det > o_det && python ../e2e_check.py o_det o_det_ref
	cd e2e_tests/c2 && ../../sqmc -i i_st > o_st && python ../e2e_check.py o_st o_st_ref

clean:
	rm -f $(OBJS) *.mod

debug: FFLAGS += -g -fbounds-check -fbacktrace
debug: F90FLAGS += -g -fbounds-check -fbacktrace -DDEBUG
debug: all

profile: FFLAGS += -g -pg
profile: F90FLAGS += -g -pg
profile: all
