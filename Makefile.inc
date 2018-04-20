# Common definitions

FC = gfortran
LD = gfortran

# Put whatever is needed to compile against umfpack
LIBS = -L/usr/lib/x86_64-linux-gnu -lumfpack

OPTIM = -cpp -pipe -mtune=native -march=native -O2
NORMAL = -cpp -mtune=native -march=native -O2 -fexternal-blas
DEBUG = -cpp -std=f2008 -pedantic -Wconversion -Wall \
	-Wcharacter-truncation -Wunderflow -g -fbounds-check -fbacktrace \
	-fimplicit-none -fdump-core -ffpe-trap=invalid,zero
# -ffpe-trap=invalid,zero,underflow,denormal : denormal, unerflow causes problem
# with iso_c_bindings ???
# Don't track overflow error with lapack, because it creates overflows and work
# with them.
DEBUG_OPTIM = $(DEBUG) -Wextra -Warray-temporaries -ffree-line-length-0 \
	-fcheck=all -finit-real=nan

FFLAGS = $(DEBUG)
LDFLAGS = $(FFLAGS)

AR	= ar
ARFLAGS	=
RANLIB	= ranlib

RM	= rm

ECHO	= echo

# .SILENT :

# Regles implicites
.SUFFIXES:
.SUFFIXES: .o .f90 .mod

%.o : %.mod
