include Makefile.inc

OBJS 	=   prec.o sort.o matrices.o list.o mat_list.o mat_csr.o \
	umfpack_calls.o mat_csc.o fe_mesh_quentin.o p1.o fonction.o
MODS 	=   prec.mod sort.mod matrices.mod list.mod mat_list.mod mat_csr.mod \
	umfpack_calls.mod mat_csc.mod fe_mesh_quentin.mod p1.mod fonction.mod

EXE	= laplace2 test_csc
EXE_OBJS = laplace2.o test_csc.o

all: $(OBJS) $(EXE_OBJS)
	-for m in $(EXE); do (\
		$(LD) $(LDFLAGS) -o $$m.exe $$m.o $(OBJS) $(LIBS) ); done

clean :
	-$(RM) -f $(OBJS) $(EXE_OBJS)

distclean : clean
	-$(RM) -f $(MODS) $(EXE)
	-for m in $(EXE); do (\
		$(RM) -f $$m.exe ); done

doc : force_look
	doxygen
	cd doc/latex; $(MAKE) $(MFLAGS)

force_look :
	true

# Les dépendances... faites à la main
list.o 		: prec.o
mat_csc.o	: prec.o matrices.o list.o mat_list.o umfpack_calls.o sort.o
mat_csr.o	: prec.o matrices.o list.o mat_list.o
mat_list.o	: prec.o matrices.o
matrices.o	: prec.o
sort.o		: prec.o
umfpack_calls.o : prec.o
fe_mesh_quentin.o	: vtkCellType.inc prec.o
p1.o  		: fe_mesh_quentin.o prec.o mat_list.o fonction.o mat_csc.o mat_csr.o matrices.o

laplace2.o 	: prec.o fe_mesh_quentin.o mat_list.o mat_csr.o mat_csc.o p1.o

.f90.o:
	$(FC) $(FFLAGS) -c $<
.f90.mod:
	$(FC) $(FFLAGS) -c $<
