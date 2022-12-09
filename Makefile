# Make file to run the project
FC=gfortran
FFLAGS= -fcheck=all -ffpe-trap=invalid,zero,overflow -O0 -fbacktrace -g
SRC=module.f90 convection.f90
OBJ=${SRC:.f90=.o}

%.o : %.f90
	$(FC) $(FFLACS) -o $@ -c $<

run: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

exec: run
	run

plot: exec
	python Visualization.py

vis:
	python Visualization.py

clean:
	@del *.avi *.exe *.mod *.o results\*.txt results\*.csv figures_cosine\*.png figures_random\*.png