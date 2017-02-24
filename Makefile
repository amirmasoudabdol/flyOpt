# CODE VERSION ######################################################

VERSION = 10.0

# executables to make

# Build commands
# 2016, A.M.Abdol:
# 	Based on the METHOD variable provided to this file, the build order will change
# 	and different binary file will be build, fly_ss or fly_ess.

ssRule=ssa
essRule=essa
ssFolder=ss
essFolder=ess

ifeq ($(METHOD), -DSS)
	execFile = fly_ss
	methodRule=$(ssRule)
	methodFolder=$(ssFolder)
	methodINCLUDES=-I../ss
else ifeq ($(METHOD), -DESS)
	execFile = fly_ess
	methodRule=$(essRule)
	methodFolder=$(essFolder)
	methodINCLUDES=-I../ess
endif

# FLAGS FOR -v FOR ALL EXECUTABLES ##################################
# this passes user and host name, compiler and version to the com-
# pile so it can be printed in the -v message of each executable

USRFLAG  = -DUSR=\"$(USER)\"
HOSTFLAG = -DMACHINE=\"$(HOST)\"
COMPFLAG = -DCOMPILER=\"$(CC)\"
FLAGFLAG = -DFLAGS=\"optimized\"
VERSFLAG = -DVERS=\"$(VERSION)\"

VFLAGS = $(USRFLAG) $(HOSTFLAG) $(COMPFLAG) $(FLAGFLAG) $(VERSFLAG) 

# find out about which architecture we're on and set compiler 
# accordingly

SSFLAGS = -DSS
ESSFLAGS = -DESS

OSTYPE=linux-gnu
#ifeq ($(OSTYPE),linux-gnu)

# This should recognize most linux machines by looking for the string "linux" in $OSTYPE 
# The gcc version specified manually due to the miss configuration of my Mac
ifneq (,$(findstring linux,$(OSTYPE)))
	CC = gcc
	DEBUGFLAGS = $(DEBUGFLAGS)
	PROFILEFLAGS = $(PROFILEFLAGS)
	FLYEXECS = unfold printscore scramble $(execFile)
	SUNDIALS = /usr/local
endif

ifeq ($(OSTYPE),osf1)
	CC = cc
endif

#echo $(CC)

# find the compiler and set corresponding flags

ifeq ($(CC),cc)
	CCFLAGS = -std1 -fast -DALPHA_DU -DNDEBUG 
	DEBUGFLAGS = -std1 -O0 -g
	PROFILEFLAGS = -O2 -g1 
	LIBS = -lm -ldxml -lgsl -lgslcblas
	FLIBS = $(LIBS)
	KCC = /bin/kcc
	KFLAGS = -ckapargs=' -scalaropt=3 -optimize=5 -roundoff=3 -arl=4 '
# uncomment 2 lines below if you don't want kcc
#	KCC = $(CC)
#	KFLAGS = $(CCFLAGS)
endif

ifeq ($(CC),icc)
# lucas flags for profiling and debugging 2-6-03
# 	CCFLAGS = -O3 -DNDEBUG
        CCFLAGS = -O3 -xW -tpp7 -ipo 
        PRECFLAGS = -mp -prec_div -pc80 
#       DEBUGFLAGS = -g  -inline_debug_info -O1  -xW -tpp7
        DEBUGFLAGS = -g  -inline_debug_info -O0
#	DEBUGFLAGS = -g
#       PROFILEFLAGS = -prof_dir profile_data -prof_gen -O2  -xW -tpp7
	PROFILEFLAGS = -p -qp -O2 -xW -tpp7
#       USEPROFFLAGS = -prof_use  -prof_dir profile_data  -O3 -xW -tpp7 -ipo -opt_report
	LIBS = -limf -lgsl -lgslcblas
r/local	LIBS = -lm
	FLIBS = -limf -lgsl -lgslcblas -static
	KCC = $(CC)
	KFLAGS = $(CCFLAGS)
        export ICC = "yes"
endif

ifeq ($(CC),gcc)
  	CCFLAGS = -Wall -m64 -O2 -std=gnu99 -DHAVE_SSE2 $(METHOD) 
   	PROFILEFLAGS = -g -pg -O2 -DHAVE_SSE2
	LIBS = -lm -lgsl -lgslcblas -lsundials_cvode -lsundials_nvecserial -L$(SUNDIALS)/lib
	FLIBS = -lm -lgsl -lgslcblas -lsundials_cvode -lsundials_nvecserial -L$(SUNDIALS)/lib
	KCC = $(CC)
	KFLAGS = $(CCFLAGS)
endif

# debugging?

ifdef DEBUG
	CCFLAGS = $(DEBUGFLAGS)
	FLAGFLAG = -DFLAGS=\"debugging\"
else
	DEBUG = "Off"
endif

ifdef PROFILE
	CCFLAGS = $(PROFILEFLAGS)
	FLAGFLAG = -DFLAGS=\"profiling\"
	KCC = $(CC)
	KFLAGS =
else
	PROFILE = "Off"
endif

ifdef USEPROFILE
        CCFLAGS = $(USEPROFFLAGS)
endif

# export all variables that Makefiles in subdirs need
# 2012 july 25, I needed to add the location of my sundials libraries - A. Crombach

export INCLUDES = -I. -I../utils -I/usr/local/include -I$(SUNDIALS)/include $(methodINCLUDES) -I/usr/include/malloc
export CFLAGS = -std=gnu99 $(CCFLAGS) $(INCLUDES)
export VFLAGS
export CC
export KCC
export KFLAGS
export LIBS
export FLIBS
export FLYEXECS

export SSFLAGS
export ESSFLAGS
export SAFLAGS
export METHOD
export execFile

#define targets

fly: utl 
	cd fly && $(MAKE)

deps: 
	cd fly && $(MAKE) -f basic.mk Makefile && chmod +w Makefile

utl: $(methodRule)
	cd utils && make

$(methodRule):
	cd $(methodFolder) && make

clean:
	rm -f core* *.o *.il
	rm -f */core* */*.o */*.il
	rm -f fly/unfold fly/printscore fly/scramble
	rm -f fly/fly_ss fly/fly_ess

veryclean:
	rm -f core* *.o *.il
	rm -f */core* */*.o */*.il */*.slog */*.pout */*.uout
	rm -f fly/unfold fly/printscore fly/scramble
	rm -f fly/fly_ss fly/fly_ess
	rm -f utils/gen_deviates
	rm -f fly/Makefile
	rm -f ss/*.o
	rm -f ess/*.o

cleanoutput:
	rm -f input/*.state
	rm -f *.csv
	rm -f *.out

help:
	@echo "make: this is the Makefile for fly code"
	@echo "      always 'make deps' first after a 'make veryclean'"
	@echo ""
	@echo "      the following targets are available:"
	@echo "      utl:       make object files in the utils directory only"
	@echo "      fly:       compile the fly code (which is in 'fly')"
	@echo "      clean:     gets rid of cores and object files"
	@echo "      veryclean: gets rid of executables and dependencies too"
	@echo ""
	@echo "      your current settings are:"   
	@echo "      compiler:  $(CC)"
	@echo "      flags:     $(CFLAGS)"
	@echo "      debugging: $(DEBUG)"
	@echo "      profiling: $(PROFILE)"
	@echo "      os type:   $(OSTYPE)"
	@echo ""


