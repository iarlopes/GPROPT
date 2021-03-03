# Build LINKS (single-scale) or LINKS2 (coupled multi-scale) according to the targets
# specified:
# - `make`: single-scale, for debugging (default)
# - `make optm`: single-scale, optimised for performance
# - `make clean`: remove executables and object files
# - `make test_builds`: test which compilation types are successful (all available
#    compilers, develop and release, single and coupled multi-scale)
#
# To control the compiler used, set the variable `FC`:
# - `FC=ifort`: Intel Fortran (default)
# - `FC=gfortran`: GNU Fortran
# - `FC=nagfor`: NAG Fortran
# - `FC=pgfortran`: PGI Fortran
#
# Useful `make` options:
# `-B`: force rebuild of all targets
# `-j {N}`: parallel build with `N` threads; if `N` is absent, use all available threads
#
# Notes:
# 1) By default, this Makefile mirrors the source directory tree to store object files (in
#    `$OBJDIR`). The resulting executable is placed in directory `$BINDIR`.
#
# 2) By default, this Makefile serialises the output, which removes all colour that
#    compilers assign to error and warning messages. As a workaround, you can pipe the
#    output to the log colouriser `ccze` (install with `sudo apt install ccze`, on Ubuntu):
#    `make | ccze -A`
#
## History
# Name          | Date     | Comment
# ------------- | -------- | ---------------------------------------------------------------
# Igor Lopes    | Mar 2021 | Initial version (adapted from LINKS makefile)
#===========================================================================================
SHELL = /bin/bash
# Output colours
RED=$$(tput setaf 1)
YELLOW=$$(tput setaf 3)
GREEN=$$(tput setaf 2)
NORMAL=$$(tput sgr0)
# Find number of available processes
NPROC=$(shell grep -c ^processor /proc/cpuinfo)
# Serialise output (this is equivalent to setting --output-sync)
ifeq (0, $(words $(findstring $(MAKECMDGOALS), test_builds)))
  MAKEFLAGS += O
endif
#
# Compilation type: debug or release
# ---------------------------------------------------------------------
ifeq ($(findstring optm, $(MAKECMDGOALS)), optm)
  DEBUG = 0
  BUILD_TYPE += release
else
  DEBUG = 1
  BUILD_TYPE += develop
endif
#
# Executable file name and location
# ---------------------------------
BINDIR = ../bin
ALLBIN = $(BINDIR)/GPROPT_debug $(BINDIR)/GPROPT
ifeq ($(DEBUG),1)
  BIN = $(word 1, $(ALLBIN))
else ifeq ($(DEBUG),0)
  BIN = $(word 2, $(ALLBIN))
endif
#
# Directory where object (.o) and module (.mod) files are stored
# --------------------------------------------------------------
OBJROOT = ../obj
ifeq ($(DEBUG), 1)
  OBJDIR := $(addsuffix /debug, $(OBJROOT))
else
  OBJDIR := $(addsuffix /release, $(OBJROOT))
endif
#
# Fortran compiler setup
# ======================
FC = ifort
COMPILER = $(FC)

# Define compiler family and associated preprocessor symbol
ifeq ($(COMPILER),$(filter $(COMPILER),ifort))
  FCTYPE = intel
  PRE_SYM:=__INTEL_COMPILER
else ifeq ($(COMPILER),$(filter $(COMPILER),gfortran))
  FCTYPE = gnu
  PRE_SYM:=__GFORTRAN__
else ifeq ($(COMPILER),$(filter $(COMPILER),nagfor))
  FCTYPE = nag
  PRE_SYM:=NAGFOR
else ifeq ($(COMPILER),$(filter $(COMPILER),pgfortran pgf90 pgf95))
  FCTYPE = pgi
  PRE_SYM:=__PGI
endif
#
# Intel Fortran
# -------------
ifeq ($(FCTYPE),intel)
  ifeq ($(OS),Windows_NT) # Windows
    FFLAGS := /Qopenmp /heap-arrays:2096 /noerror-limit /module:$(OBJDIR)
    ifeq ($(DEBUG), 1) # Debug
      FFLAGS += /O0 /debug:full /traceback /check:all /gen-interfaces /warn:all \
                /Qtrapuv /Qfp-speculation=strict /fp:strict /Qprec /fpe=0 \
                /fpe-all=0 /Qftz /Qinit:snan,arrays,zero /nostandard-realloc-lhs \
								/check:noarg_temp_created
      # Other options that might be of interest:
                #/stand:f08 /gen-dep:$(OBJDIR)/deps.txt
    else # Release
      FFLAGS += /O3 /QxHost
    endif
  else # Other OSs
    FFLAGS := -qopenmp -heap-arrays 2096 -module $(OBJDIR) -nostandard-realloc-lhs # -stand f08
    # Extra flags for .f90 files only
    FFLAGS_F90 := -warn errors
    ifeq ($(DEBUG), 1) # Debug
      FFLAGS += -O0 -debug full -traceback -check all -gen-interfaces -warn all \
                -ftrapuv -fp-speculation=strict -mp1 -fp-model strict \
                -fpe=0 -fpe-all=0 -ftz -init=snan,arrays,zero -check noarg_temp_created
      # Other options that might be of interest:
                #-gen-dep=$(OBJDIR)/deps.txt
    else # Release
      # Using PGO (profiler-guided optimisation):
      # 1) Compile adding the flag `-prof-gen=threadsafe -prof-dir ../bin`
      # 2) Run some representative examples. The profiler will collect execution data to
      # optimise the code. These will be stored as .dyn/.dpi files in the directory
      # specified by the `-prof-dir` option above
      # 3) Recompile, now using `-prof-use -prof-dir ../bin` in place of `-prof-gen`
      FFLAGS += -O3 -xHost
      # Other options that might be of interest:
                #-fast -parallel -ipo
    endif
  endif
# # gfortran
# # --------
# else ifeq ($(FCTYPE),gnu)
#   FFLAGS := -fopenmp -ffree-line-length-none -fimplicit-none -J$(OBJDIR) -std=f2008 -pedantic-errors -fall-intrinsics
#   ifeq ($(DEBUG), 1) # Debug
#     FFLAGS += -O0 -g -ggdb -fbacktrace -fcheck=all -Wall -Wextra -fmax-errors=0 \
#               -finit-integer=-1000000 -finit-real=nan -finit-logical=true \
#               -finit-character=0 -ftrapv -fno-fast-math -fno-omit-frame-pointer \
#               -ffpe-trap=zero,overflow -fsanitize=undefined -fuse-ld=gold
#     # Other options that might be of interest:
#               #-ffpe-trap=invalid,zero,overflow -fsanitize=address
#   else # Release
#     # Using PGO (profiler-guided optimisation)
#     # 1) Compile adding the flag `-fprofile-generate`
#     # 2) Run some representative examples. The profiler will collect execution data to
#     # optimise the code. These will be stored as .gcda files in the same directory where
#     # the object files are created (by default)
#     # 3) Recompile, now using `-fprofile-use -fprofile-correction` in place of
#     # `-fprofile-generate`. Do not use `make clean` between both compilations, as this would
#     # delete the .gcda files! Instead, use the `-B` make option to force recompilation
#     FFLAGS += -O2 -march=native -mtune=native -flto=$(NPROC) -funroll-loops
#     # Other options that might be of interest:
#               #-Ofast -falign-functions=1 -falign-jumps=1 -falign-loops=1 -falign-labels=1
#               #-fassociative-math -freciprocal-math -fno-signed-zeros -fno-trapping-math
#               #-frename-registers
#   endif
# # NAG Fortran
# # -----------
# else ifeq ($(FCTYPE),nag)
#   FFLAGS := -kind=byte -openmp -mdir $(OBJDIR) -I$(OBJDIR) -f2008 \
#             -colour=error:black+green,warn:yellow,info:magenta -u -unsharedrts
#   ifeq ($(DEBUG), 1) # Debug
#     FFLAGS += -g -nan -ieee=nonstd -info -C=all
#   else # Release
#     FFLAGS += -O4 -Oblas -ieee=full
#   endif
# # PGI Fortran
# # -----------
# else ifeq ($(FCTYPE),pgi)
#   FFLAGS := -module $(OBJDIR) -I$(OBJDIR) -Mbackslash -mp
#   ifeq ($(DEBUG), 1) # Debug
#     FFLAGS += -g -O0 -Mbounds -traceback -Mchkstk -Mchkptr
#   else # Release
#     FFLAGS += -O3 -Munroll=c:1 -Mnoframe -Mlre -Mpre -Mflushz -Mcache_align -Mvect=simd
#   endif
else
  $(error Unsupported compiler!)
endif
#
# Add preprocessor passing
ifeq ($(FCTYPE),nag)
  FFLAGS += -fpp -I.
else
  FFLAGS += -cpp -I.
endif
#
# Pass debug build information
ifeq ($(DEBUG), 1)
	FFLAGS += -DDEBUG_BUILD
endif
#
# Linker flags and libraries
# ==========================
# In all cases, add the following
LDLIBS += -lpthread -lm -ldl
#
# Find source files and define associated object files
# ----------------------------------------------------
# Define a recursive wildcard function: (see http://stackoverflow.com/a/18258352)
# first argument - directory where files are to be searched recursively;
# second argument - list of extensions to be searched for
rwildcard=$(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2)$(filter $(subst *,%,$2),$d))
#
# Include all directories with source files
EXT = *.f90 *.F90 *.f
SRC := $(call rwildcard,,$(EXT))
#
# Object files: replace source extensions with .o, adding 'OBJDIR/' as prefix to path
OBJ = $(addprefix $(OBJDIR)/,$(addsuffix .o,$(basename $(SRC))))
#
# File creation/deletion on different OSs
# ---------------------------------------
# (fixPath converts Unix addresses (separator: '/') to Windows-style (separator: '\'))
ifeq ($(OS),Windows_NT)
  RM = del /Q
  MKDIR = mkdir
  fixPath = $(subst /,\,$1)
else
  RM = rm -f
  RMR = rm -rf
  MKDIR = mkdir -p
  fixPath = $1
endif
# Don't output compilation header with targets in list NOHEADER
NOHEADER := clean gen_dep upd_version test_builds
ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NOHEADER))))
  $(info +======================== Building GPROPT: $(BUILD_TYPE) ========================)
  ifdef IS_REPO
    $(info | Current branch: $(shell git rev-parse --abbrev-ref HEAD))
  endif
  $(info |       Compiler: $(FC))
  $(info |   Linker flags: $(LDFLAGS))
  $(info |    Linker libs: $(LDLIBS))
  $(info | Compiler flags: $(FFLAGS))
  $(info +=================================================================================)
endif
#
# Rules for .f90, .F90 (these go through the C++ pre-processor) and .f files
# --------------------------------------------------------------------------
.PHONY : all optm
all optm : $(BIN)
#	@true # This suppresses the "nothing to be done for ..." messages
#
$(OBJDIR)/%.o: %.f90
	@$(MKDIR) $(dir $@)
	@echo Compiling $< ...
	@$(FC) $(INC) $(FFLAGS) $(FFLAGS_F90) -c -o $@ $<
$(OBJDIR)/%.o: %.F90
	@$(MKDIR) $(dir $@)
	@echo Compiling $< ...
	@$(FC) $(INC) $(FFLAGS) $(FFLAGS_F90) $(CPPFLAGS) -c -o $@ $<
$(OBJDIR)/%.o: %.f
	@$(MKDIR) $(dir $@)
	@echo Compiling $< ...
	@$(FC) $(INC) $(FFLAGS) -c -o $@ $<
#
# Linking
# -------
#
$(BIN): $(OBJ)
	@$(MKDIR) $(dir $(BIN))
	@echo Linking $@ ...
	@$(FC) $(LDFLAGS) $(FFLAGS) -o $@ $(OBJ) $(LDLIBS)
	@echo "Done!"
	@echo "+=========================== Built LINKS: $(BUILD_TYPE) ========================"
	@if ! test -z "$(IS_REPO)"; then \
	  echo "| Current branch: $$(git rev-parse --abbrev-ref HEAD)"; \
	fi
	@echo "|       Compiler: $(FC)"
	@echo "|   Linker flags: $(LDFLAGS)"
	@echo "|    Linker libs: $(LDLIBS)"
	@echo "| Compiler flags: $(FFLAGS)"
	@echo "+================================================================================="
#
# Remove executables and all files in OBJROOT
# -------------------------------------------
.PHONY : clean
clean :
	@$(RMR) $(call fixPath,$(OBJROOT)/*)
	@$(RM)  $(call fixPath,$(ALLBIN))
	@$(RM)  .ac_fortran
	@echo "Cleanup complete!"
#
# Try building all versions of GPROPT with all available compilers
# ---------------------------------------------------------------
.PHONY : test_builds
# Find out all available compilers...
# ... single-scale build
HAS_GFORTRAN  := $(shell command -v gfortran 2> /dev/null)
HAS_IFORT     := $(shell command -v ifort 2> /dev/null)
HAS_PGFORTRAN := $(shell command -v pgfortran 2> /dev/null)
HAS_NAGFOR    := $(shell command -v nagfor 2> /dev/null)
#
# Function test_build tests a certain compilation type. Arguments:
#   1 - String containing compiler address (empty if not available)
#   2 - Make log filename (to be stored in directory $(BINDIR))
#   3 - Compiler name (to be echoed to stdout)
#   4 - Options passed to make for that type of compilation
test_build =                                                                     \
  if ! test -z "$(1)"; then                                                      \
    echo;                                                                        \
    MAKE_LOG="$(BINDIR)/$(2)";                                                   \
    $(MAKE) --no-print-directory clean >/dev/null 2>&1;                          \
    echo "$(3)...";                                                              \
    $(MAKE) --output-sync --no-print-directory $(4) -j > $$MAKE_LOG 2>&1;        \
    if [ $$? -eq 0 ]; then                                                       \
      echo "$(GREEN)Success!$(NORMAL) (log: file://$$(readlink -f $$MAKE_LOG))"; \
    else                                                                         \
      echo "$(RED)Failed$(NORMAL) (log: file://$$(readlink -f $$MAKE_LOG))";     \
    fi;                                                                          \
    if [ "$(3)" == "NAG Fortran" ]; then                                         \
      sed -i '/NAG Fortran Compiler/d' $$MAKE_LOG;                               \
    fi;                                                                          \
  fi
#
test_builds :
	@$(MKDIR) $(BINDIR)
	@$(RM)  $(BINDIR)/.make_*.log
	@echo "Test LINKS' builds"; echo "=================="
	@echo; echo; echo "Single-scale, develop"; echo "---------------------"
	@$(call test_build,$(HAS_GFORTRAN),.make_single_develop_gfortran.log,GNU Fortran,FC=gfortran)
	@$(call test_build,$(HAS_IFORT),.make_single_develop_ifort.log,Intel Fortran,FC=ifort)
	@$(call test_build,$(HAS_PGFORTRAN),.make_single_develop_pgfortran.log,PGI Fortran,FC=pgfortran)
	@$(call test_build,$(HAS_NAGFOR),.make_single_develop_nagfor.log,NAG Fortran,FC=nagfor)
	@echo; echo; echo "Single-scale, release"; echo "---------------------"
	@$(call test_build,$(HAS_GFORTRAN),.make_single_release_gfortran.log,GNU Fortran,optm FC=gfortran)
	@$(call test_build,$(HAS_IFORT),.make_single_release_ifort.log,Intel Fortran,optm FC=ifort)
	@$(call test_build,$(HAS_PGFORTRAN),.make_single_release_pgfortran.log,PGI Fortran,optm FC=pgfortran)
	@$(call test_build,$(HAS_NAGFOR),.make_single_release_nagfor.log,NAG Fortran,optm FC=nagfor)
