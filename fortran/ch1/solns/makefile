include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
FFLAGS +=

expx: expx.o
	-${FLINKER} -o expx expx.o  ${PETSC_LIB}
	${RM} e.o

distclean:
	@rm -f *~ expx *tmp

