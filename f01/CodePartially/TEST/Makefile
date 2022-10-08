TESTS = test_user_input     \
        test_user_input_mpi \
        test_mkdir          \
        test_rfield		    \
        test_cfield		    \
        test_eigen_vectors  \
        test_ev             \
        test_functions
		# test_lapack		    \
		test_scalapack	    \
		test_scala_module   \
		test_scatter_gather 

COMMON_MODULES = ../user_input.f90.o ../field_constants.f90.o ../constants.f90.o ../utilities.f90.o
FIELD_MODULES  = ../rfield_class.f90.o ../cfield_class.f90.o

TEST_FFLAGS = $(FFLAGS) -I.. -cpp -g

.PHONY:
all: ${TESTS}

test_user_input:    test_user_input.f90.o \
                    ${COMMON_MODULES}
	${MPIF90} -o $@ $^ ${LDFLAGS}

test_user_input_mpi:    test_user_input_mpi.f90.o \
                    ${COMMON_MODULES}
	${MPIF90} -o $@ $^ ${LDFLAGS}

test_ev:    test_ev.f90.o           \
            ../eigen_vectors.f90.o  \
            ../utilities.f90.o      \
            ${COMMON_MODULES}       \
            ${FIELD_MODULES}
	${MPIF90} -o $@ $^ ${LDFLAGS}

test_functions:     test_functions.f90.o            \
                    ../spin_flip_procedures.f90.o   \
                    ../eigen_vectors.f90.o          \
                    ${COMMON_MODULES}               \
                    ${FIELD_MODULES}
	${MPIF90} -o $@ $^ ${LDFLAGS}

test_cfield:    test_cfield.f90.o   \
                ${COMMON_MODULES}   \
                ${FIELD_MODULES}
	${MPIF90} -o $@ $^ ${LDFLAGS}

test_rfield:    test_rfield.f90.o       \
                ../rfield_class.f90.o   \
                 ${COMMON_MODULES}
	${MPIF90} -o $@ $^ ${LDFLAGS}

test_eigen_vectors: test_eigen_vectors.f90.o    \
                    ../eigen_vectors.f90.o      \
                    ../utilities.f90.o          \
                    ${COMMON_MODULES}
	${MPIF90} -o $@ $^ ${LDFLAGS} ${SCALAPACK}

test_lapack: test_lapack.f90.o
	${MPIF90} -o $@ $< ${LAPACK} ${LDFLAGS}

test_scalapack: test_scalapack.f90.o
	${MPIF90} -o $@ $< ${SCALAPACK} ${LDFLAGS}

test_scala_module: ${MODULE_OBJECTS} random_module.f90.o
	${MPIF90} -o $@ $^ ${SCALAPACK} ${LDFLAGS}

test_scatter_gather: ${MODULE_OBJECTS} test_scatter_gather.f90.o
	${MPIF90} -o $@ $^ ${SCALAPACK} ${LDFLAGS}

test_mkdir: ../c_functions.f90.o ../called_from_fortran.c.o test_mkdir.f90.o
	${F90} -o $@ $^ -I.. ${LDFLAGS}

%.f90.o: %.f90 ${MODULES}
	${MPIF90} -c -o $@ $< ${TEST_FFLAGS} 

%.f90:

.PRECIOUS: %.mod

.PHONY: clean
clean:
	rm -f *.mod *.f90.o ${TESTS}

