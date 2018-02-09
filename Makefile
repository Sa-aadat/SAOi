# SAOi:
# SAOi: Sat May 02 09:04:51 SAST 2009, Albert Groenwold, Stellenbosch
# SAOi: Mon 27 Jul 2017 14:37:50 SAST, Sa-aadat Parker, UCT 
# SAOi: Makefile
# SAOi:
#

#
.PHONY: clean proper
#
#
# Select the new Gnu fortran95 compiler ##############################
#   F77 = gfortran-4.1
#
# Select the new Gnu fortran95 compiler ##############################
#   F77 = gfortran-4.2
#
# Select the new Gnu fortran95 compiler ##############################
#   F77 = gfortran-4.3
#
# Select the new Gnu fortran95 compiler ##############################
    F77 = gfortran
#
# Select the Sun Studio Compiler (after: module add sunstudio12) #####
#   F77 = f95  -f77=input,output  # requires F77 compliant I/O
#
# Select the old intel compiler ######################################
#   F77 = ifc
#
# Select the new intel compiler on older platforms ###################
#   F77 = ifort -no-ipo
#
# Select the new intel compiler ######################################
#   F77 = ifort
#
# Select the default new g95 compiler ################################
#     F77 = g95
#
# Select the new g95 compiler on 32 bit machines######################
#   F77 = g95-32
#
# Select the new g95 compiler on 64 bit machines; 32 bit integers ####
#   F77 = g95-64-32
#
# Select the new g95 compiler on 64 bit machines; 64 bit integers ####
#   F77 = g95-64-64
#
# Specify some optional flags ########################################
#   FFLAGS = -fast    # full optimization for ifort on 64 bit machines
#
# Specify some other optional flags ##################################
#     FFLAGS = -O2
#
# Specify some other optional flags ##################################
#   FFLAGS = -O3 
#
# Specify some other optional flags ##################################
#   FFLAGS = -O -pg
#
# Specify some more optional flags ###################################
#   FFLAGS = -O -mcmodel=medium
#
# Or specify debugging flag ##########################################
#   FFLAGS = -g
#
# Or specify debugging flag and check everything with ifort ##########
#   FFLAGS = -g -check all -warn all
#
# Or specify debugging flag and check everything with gfortran #######
    FFLAGS = -ggdb -Wall -mcmodel=medium -m64 -Wconversion -Wextra
    #-frange-check -Warray-bounds -Waliasing -Wsurprising -fbounds-check -Wline-truncation -Wconversion -Wextra -Wno-unused -fbounds-check -frange-check
#
#   -ffpe-trap=invalid,zero,overflow,underflow,inexact,denormal-Warray-temporaries
#
# Or specify debugging flag and check everything with g95 ############
#   FFLAGS = -g  -Wall -Wextra -Wno-unused -fbounds-check  -frange-check
# 
# Specify C compiler
    GC = gcc
#    
# Specify flags for the C compiler
    CFLAGS = -g -Wall -Warray-bounds -fbounds-check -Wconversion -Wno-unused  -Wno-unused-function -mcmodel=medium -m64


COMFLAG =

ARCH = $(shell uname -m)

OS = $(shell uname)

MYHOME = $(shell ~/)

#
.SUFFIXES : .o .f90 .f .f03

#
#  Include files for algorithms

include Repo/tmp/Makefile.nag

include Repo/tmp/Makefile.cplex

include Repo/tmp/Makefile.galahad


#
#  Object files for the SAOi algorithm

files1 =         approximations.f Gradients.f split.f           \
                 SAOi_dense.f lbfgsb24.f Functions.f timer.f    \
                 Initialize.f dlbfgsb24f.f Utils_Dense.f Falk.f \
                 diaHess.f defaults.f uncstr.f conserve.f       \
                 trustfilter.f NAGdums.f qcauchy.f formats.f    \
                 enforce.f dgalQP.f diaHessUser.f Special.f     \
                 testNaN.f cplexqp.f dlbfgsb24u.f smmp1.f       \
                 SAOi_sparse.f Utils_Sparse.f dgalQPs.f
                 
files2 =         dlsqpd.f90 dqpad.f90 dqpbd.f90 dqpcd.f90

files3 =         SAOi.f 

files4 =         dcplex.c

files5 =         Fem.f03 FemInit.f03 geom.f03 main.f03          \
                 Fem_g.f03 Fem_gg.f03

object_files =   approximations.o Gradients.o split.o           \
                 SAOi_dense.o lbfgsb24.o Functions.o timer.o    \
                 Initialize.o dlbfgsb24f.o Utils_Dense.o Falk.o \
                 diaHess.o  defaults.o uncstr.o conserve.o      \
                 trustfilter.o NAGdums.o qcauchy.o formats.o    \
                 enforce.o dgalQP.o diaHessUser.o Special.o     \
                 testNaN.o cplexqp.o dlbfgsb24u.o dlsqpd.o      \
                 dqpad.o dqpbd.o dqpcd.o smmp1.o                \
                 SAOi_sparse.o Utils_Sparse.o dgalQPs.o         \
                 SAOi.o  Fem.o FemInit.o geom.o dcplex.o        \
                 main.o Fem_g.o Fem_gg.o
                 
.c.o:
	   $(GC) -c $(CFLAGS) $*.c

.f.o:
	   $(F77) -c $(FFLAGS) $*.f

.f90.o:
	   $(F77) -c $(FFLAGS) $*.f90
	   
.f03.o:
	   $(F77) -c $(FFLAGS) $*.f03

SAOi: $(object_files)
	   $(F77) $(FFLAGS) $(cplexflag) $(COMFLAG) -fopenmp -o SAOi $(object_files) $(cplexlib) $(naglib) $(gallib64d) -lblas -llapack

nag:
	cp Repo/nag/NAGdums.f NAGdums.f
	cp Repo/nag/Makefile.nag Repo/tmp/Makefile.nag
	make

no-nag:
	cp Repo/bare/NAGdums.f NAGdums.f
	cp Repo/bare/Makefile.nag Repo/tmp/Makefile.nag
	make

cplex:
	cp Repo/cplex/dcplex.c dcplex.c
	cp Repo/cplex/cplexqp.f cplexqp.f
	cp Repo/cplex/Makefile.cplex Repo/tmp/Makefile.cplex
	make

no-cplex:
	cp Repo/bare/dcplex.c dcplex.c
	cp Repo/bare/cplexqp.f cplexqp.f
	cp Repo/bare/Makefile.cplex Repo/tmp/Makefile.cplex
	make

galahad64:
	cp Repo/galahad/dlsqpd.f90 dlsqpd.f90
	cp Repo/galahad/dqpad.f90 dqpad.f90
	cp Repo/galahad/dqpbd.f90 dqpbd.f90
	cp Repo/galahad/dqpcd.f90 dqpcd.f90
	cp $(galmod64d)/galahad_lsqp_double.mod .
	cp $(galmod64d)/galahad_qpa_double.mod .
	cp $(galmod64d)/galahad_qpb_double.mod .
	cp $(galmod64d)/galahad_qpc_double.mod .
	cp Repo/galahad/Makefile.galahad Repo/tmp/Makefile.galahad
	make

no-galahad:
	cp Repo/bare/dlsqpd.f90 dlsqpd.f90
	cp Repo/bare/dqpad.f90 dqpad.f90
	cp Repo/bare/dqpbd.f90 dqpbd.f90
	cp Repo/bare/dqpcd.f90 dqpcd.f90
	rm -f galahad_lsqp_double.mod
	rm -f galahad_qpa_double.mod
	rm -f galahad_qpb_double.mod
	rm -f galahad_qpc_double.mod
	cp Repo/bare/Makefile.galahad Repo/tmp/Makefile.galahad
	make

bare:
	cp Repo/bare/enforce.f enforce.f
	cp Repo/bare/dlsqpd.f90 dlsqpd.f90
	cp Repo/bare/dqpad.f90 dqpad.f90
	cp Repo/bare/dqpbd.f90 dqpbd.f90
	cp Repo/bare/dqpcd.f90 dqpcd.f90
	rm -f galahad_lsqp_double.mod
	rm -f galahad_qpa_double.mod
	rm -f galahad_qpb_double.mod
	rm -f galahad_qpc_double.mod
	cp Repo/bare/Makefile.galahad Repo/tmp/Makefile.galahad
	cp Repo/bare/NAGdums.f NAGdums.f
	cp Repo/bare/Makefile.nag Repo/tmp/Makefile.nag
	cp Repo/bare/dcplex.c dcplex.c
	cp Repo/bare/cplexqp.f cplexqp.f
	cp Repo/bare/Makefile.cplex Repo/tmp/Makefile.cplex
# 	cp Examples/unimodal/default/*  .
	make

full:
	cp Repo/cplex/dcplex.c dcplex.c
	cp Repo/cplex/cplexqp.f cplexqp.f
	cp Repo/cplex/Makefile.cplex Repo/tmp/Makefile.cplex
	cp Repo/galahad/dlsqpd.f90 dlsqpd.f90
	cp Repo/galahad/dqpad.f90 dqpad.f90
	cp Repo/galahad/dqpbd.f90 dqpbd.f90
	cp Repo/galahad/dqpcd.f90 dqpcd.f90
	cp $(galmod64d)/galahad_lsqp_double.mod .
	cp $(galmod64d)/galahad_qpa_double.mod .
	cp $(galmod64d)/galahad_qpb_double.mod .
	cp $(galmod64d)/galahad_qpc_double.mod .
	cp Repo/galahad/Makefile.galahad Repo/tmp/Makefile.galahad
	make

test:
	cp Repo/bare/enforce.f .
	make galahad
	./test.sh

galahad:
   ifeq ($(OS),Linux)
        ifeq ($(ARCH),x86_64)
	       make galahad64
        endif
   else
	@echo '  '
	@echo ' OS is not Linux 64 ... We cannot help at this stage, sorry. '
	@echo '  '
   endif

approximations.o : approximations.f ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
SAOi.o           : SAOi.f           ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
Gradients.o      : Gradients.f      ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
split.o          : split.f          ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
SAOi_dense.o     : SAOi_dense.f     ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
lbfgsb24.o       : lbfgsb24.f       ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
dlbfgsb24f.o     : dlbfgsb24f.f     ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
dlbfgsb24u.o     : dlbfgsb24u.f     ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
Functions.o      : Functions.f      ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
Initialize.o     : Initialize.f     ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
Utils_Dense.o    : Utils_Dense.f    ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
Falk.o           : Falk.f           ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
diaHess.o        : diaHess.f        ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
diaHessUser.o    : diaHessUser.f    ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
defaults.o       : defaults.f       ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
uncstr.o         : uncstr.f         ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
conserve.o       : conserve.f       ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
trustfilter.o    : trustfilter.f    ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
qcauchy.o        : qcauchy.f        ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
timer.o          : timer.f                                                  Makefile
NAGdums.o        : NAGdums.f                                                Makefile
formats.o        : formats.f                                                Makefile
testNaN.o        : testNaN.f                                                Makefile
enforce.o        : enforce.f        ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
dgalQP.o         : dgalQP.f         ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
dqpad.o          : dqpad.f90        ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
dqpbd.o          : dqpbd.f90        ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
dqpcd.o          : dqpcd.f90        ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
Wolfe.o          : Wolfe.f          ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
Special.o        : Special.f        ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile # setsea.h
cplexqp.o        : cplexqp.f        ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
dqpbox.o         : dqpbox.f         ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
dcplex.o         : dcplex.c                                                 Makefile
smmp1.o          : smmp1.f          ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
SAOi_sparse.o    : SAOi_sparse.f    ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
Utils_Sparse.o   : Utils_Sparse.f   ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
dgalQPs.o        : dgalQPs.f        ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
FemInit.o        : FemInit.f03      ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
Fem.o            : Fem.f03          ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
main.o           : main.f03         ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
geom.o           : geom.f03         ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
Fem_g.o          : Fem_g.f03        ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile
Fem_gg.o         : Fem_gg.f03       ctrl.h ctrl_get.inc ctrl_set.inc size.h Makefile

clean :                # also cleans the test bed
	rm -f *.AXM
	rm -f *.DSP
	rm -f *.ESK
	rm -f *.GRF
	rm -f *.HIS
	rm -f *.L?
	rm -f *.NON
	rm -f *.OUT
	rm -f *.TVU
	rm -f *.DSC
	rm -f *.COR
	rm -f *.GEN
	rm -f *.GHI
	rm -f *.ID
	rm -f *.KE
	rm -f *.LM
	rm -f *.SAP
	rm -f *.TE
	rm -f *.LE
	rm -f *.RSX
	rm -f fort.*
	rm -f core
	rm -f edsap
	rm -f eigen
	rm -f *.gz
	rm -f *~
	rm -f *.OPT
	rm -f *.out* out* ttt
	rm -f DIRlog.dat
	rm -f DIRmissed.dat
	rm -f DIRresult.dat
	rm -f *.summary
	rm -f junk*
	rm -f TowerDervs
	rm -f TowerFuncs
	rm -f Dump_approx_stuff
	rm -f DesVars
	rm -f Elem*.dat
	rm -f Prob*.dat
	rm -f *.bak
	rm -f temp*
	rm -f Temp*
	rm -f mbb*
	rm -f ptwist*
	rm -f fplt*
	rm -f *.sub.po*
	rm -f *.sub.o*
	rm -f halfcyl*
	rm -f iterate.dat
	rm -f CUTEr.output
	rm -f AUTOMAT.d
	rm -f ELFUN.f
	rm -f GROUP.f
	rm -f EXTER.f
	rm -f makemycuter
	rm -f installmycuter
	rm -f OUTSDIF.d
	rm -f RANGE.f
	rm -f siflist*
	rm -f *.spc
	rm -f *.SPC
	rm -f saoimin
	rm -f *.o
	rm -f callgrind.*
	rm -f dltopo8
	rm -f ptopo8
	rm -f *.inf
	rm -f *.lds
	rm -f *.str
	rm -f Jacob*.ind
	rm -f sctopo* topo tower *.inp *.out *.dat
	rm -f DesCons *.2.f stamp*  .*swp* *.f1 sdsaoi* saoi* *.lp 
# 	cp Repo/bare/enforce.f enforce.f
# 	cp Repo/bare/dlsqpd.f90 dlsqpd.f90
# 	cp Repo/bare/dqpad.f90 dqpad.f90
# 	cp Repo/bare/dqpbd.f90 dqpbd.f90
# 	cp Repo/bare/dqpcd.f90 dqpcd.f90
	rm -f galahad_lsqp_double.mod
	rm -f galahad_qpa_double.mod
	rm -f galahad_qpb_double.mod
	rm -f galahad_qpc_double.mod
# 	cp Repo/bare/Makefile.galahad Repo/tmp/Makefile.galahad
# 	cp Repo/bare/NAGdums.f NAGdums.f
# 	cp Repo/bare/Makefile.nag Repo/tmp/Makefile.nag
# 	cp Repo/bare/dcplex.c dcplex.c
# 	cp Repo/bare/cplexqp.f cplexqp.f
# 	cp Repo/bare/Makefile.cplex Repo/tmp/Makefile.cplex
	rm -f *.case *.geo *.msh *.matid *.ndbnd *.ndlds TopOpt.*

proper : clean
	rm -f saoid.a
	rm -f *.mod
	rm -f t
	rm -f SAOi
