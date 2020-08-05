fc = gfortran
fflags += -O3 -ffree-line-length-none
parallelflags += -fopenmp
debugflags += -g -fcheck=all -Wall -ffree-line-length-none -fbacktrace

#modules must be in appropriate order
modules = src/modules/declarations.f90 src/modules/nodes.f90 src/modules/rbf.f90 src/modules/input.f90 src/modules/fluid_vars.f90 src/modules/BCs.f90 src/modules/output.f90 src/modules/MINPACK_lmder_lmdif_lmstr.f
exec = src/RBF-LM_solver.f90
parallel_exec = src/RBF-LM_solver_parallel.f90

LM_solver : $(modules) $(exec)
	$(fc) $(fflags) $(modules) $(exec) -o LM_solver
	make clean
	
parallel_LM_solver : $(modules) $(exec)
	$(fc) $(fflags) $(parallelflags) $(modules) $(parallel_exec) -o LM_solver
	make clean

debug : $(modules) $(exec)
	$(fc) $(debugflags) $(modules) $(exec) -o LM_solver

clean : 
	rm *.mod

delete :
	rm LM_solver

