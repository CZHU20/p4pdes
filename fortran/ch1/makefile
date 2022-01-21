include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
FFLAGS +=

e: e.o
	-${FLINKER} -o e e.o  ${PETSC_LIB}
	${RM} e.o

distclean:
	@rm -f *~ e *tmp

