EXES =  spin_flip_2d
MODULES =   constants.mod	    \
            c_functions.mod     \
			field_constants.mod \
            hdf5_utilities.mod  \
            user_input.mod      \
			utilities.mod       \
			rfield_class.mod	\
			cfield_class.mod	\
			eigen_vectors.mod   \
            spin_flip_procedures.mod

MODULE_OBJECTS = $(MODULES:.mod=.f90.o)
ROOT_BUILD_DIR := ${CURDIR}

# include Makefile.inc.cray if on a Cray
# otherwise include Makefile.inc.mkl
ifeq ($(strip $(BC_HOST)),mustang)
$(info Using Mustang build config)
include ./Make_include/Makefile.inc.mustang
endif

ifeq ($(strip $(BC_HOST)),gaffney)
$(info Using Gaffney build config)
include ./Make_include/Makefile.inc.gaffney
LDFLAGS = -lmpi
endif

ifeq ($(strip $(BC_HOST)),centennial)
$(info Using Centennial build config)
include ./Make_include/Makefile.inc.centennial
LDFLAGS = -lmpi
endif

ifeq ($(strip $(BC_HOST)),onyx)
$(info Using Onyx build config)
include ./Make_include/Makefile.inc.onyx
$(info Using Cray ftn)
endif

ifeq ($(strip $(BC_HOST)),narwhal)
$(info Using Narwhal build config)
include ./Make_include/Makefile.inc.narwhal
$(info Using Cray ftn)
endif

.PHONY: all
all: ${EXES} ${TESTS}

.PHONY: exes
exes: ${EXES}

.PHONY: tests
tests:  $(MODULE_OBJECTS)
	cd TEST && $(MAKE)

c_functions.mod: called_from_fortran.c.o

spin_flip_2d: main.f90.o ${MODULE_OBJECTS} called_from_fortran.c.o
	${MPIF90} -o $@ $^ ${LDFLAGS} ${SCALAPACK} ${HDF5_LIBS}
	chmod 755 $@

main.f90.o: main.f90 ${MODULE_OBJECTS}
	${MPIF90} -c -o $@ $< ${FFLAGS} -I ${HDF5_INC}

%.c.o: %.c
	$(CC) -c -o $@ $< ${CFLAGS}

%.f90.o: %.f90
	${MPIF90} -c -o $@ $< ${FFLAGS} -I ${HDF5_INC}

%.mod: %.f90.o

.PHONY: debug
debug:	FFLAGS = -g -Wall
debug:	${EXES}

.PHONY: clean
clean:
	rm -f *.mod *.MOD *.f90.* *.c.o ${EXES}
	cd TEST && $(MAKE) clean

%.f90:


.PRECIOUS: %.mod
