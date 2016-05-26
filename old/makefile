# Start by cleaning the list of suffixes
.SUFFIXES:

# We want to specify suffix rules for .f90 .f .o .x files. 
.SUFFIXES: .f90 .f .o .x

# Fortran compiler
FC = ifort 
#FC = f95

# Optimization
OFLAG=-O3

.f90.o:; $(FC) -c -o $@ $(OFLAG) $*.f90
.f.o:; $(FC) -c -o $@ $(OFLAG) $*.f
.o.x:; $(FC) -o $@ $*.o

EXE = addvacn.x ave-mesh.x band-procar-p1.x dos-procar.x neig.x relax-new.x \
      bmfit.x output-poscar.x

execution: $(EXE)

# The target bmfit.x depends on the dependence bmfit.o, mina.o and exerfc.o.
bmfit.x: bmfit.o mina.o exerfc.o
	$(FC) -o bmfit.x bmfit.o mina.o exerfc.o 

# Type 
# make -f MakefileName install
# to install *.x in the directory ~/bin (Please make sure that the directory ~/bin exists).
install:
	mv *.x ~/bin

# Type 
# make -f MakefileName clean
# to remove *.x *.o in the current directory.
clean:
	rm -f *.x *.o
