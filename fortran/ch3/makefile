include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
FFLAGS += -g

poisson: poisson.o
	-${FLINKER} -o poisson poisson.o  ${PETSC_LIB}
	${RM} poisson.o

distclean:
	@rm -f *~ poisson *tmp

