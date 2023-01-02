FC = mpif90
OPT = -ffpe-trap=invalid,zero,overflow -fbounds-check -g3 -O0 -fstack-protector-all -finit-real=snan -fbacktrace
EXE = exec

$(EXE): parameters.o exchange.o functions.o charge.o matrix.o solver.o main.o
	mpif90 -o $@ $^
%.o: %.f90
	mpif90 -c $<

clean:
	rm -f *.o *.mod $(EXE)
