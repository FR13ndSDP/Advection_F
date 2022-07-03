AMREX_HOME ?= ../..

DEBUG = FALSE
DIM = 3
COMP = gcc
PRECISION = DOUBLE
USE_MPI = FALSE
USE_OMP = FALSE
USE_CUDA = FALSE

BL_NO_FORT = TRUE

#programs = fsnapshot
programs = test

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

multiple_executables = $(addsuffix .$(machineSuffix).ex, $(programs))
default: $(multiple_executables)

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

clean::
	$(SILENT) $(RM) $(multiple_executables)
