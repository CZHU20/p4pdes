include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
FFLAGS += -g

vecmatksp: vecmatksp.o
	-${FLINKER} -o vecmatksp vecmatksp.o  ${PETSC_LIB}
	${RM} vecmatksp.o

sparsemat: sparsemat.o
	-${FLINKER} -o sparsemat sparsemat.o  ${PETSC_LIB}
	${RM} sparsemat.o

tri: tri.o
	-${FLINKER} -o tri tri.o  ${PETSC_LIB}
	${RM} tri.o

distclean:
	@rm -f *~ vecmatksp sparsemat *tmp

