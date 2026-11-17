# using ifort and MKL
FC = ifort -r8 -O
LIBS = -mkl=sequential

%.o: %.f90
	$(FC) -c $<


OBJGET := basis.o gradient.o pes_shell.o test.o


test.x : $(OBJGET)
	$(FC) -o $@ $^ $(LIBS)

clean:
	rm -f *.o *.mod *.x	
