LAPACK_DIR = /workspace/tools/lapack
F90 = gfortran
FTOOL = params.f90 tool.f90 
FFLAGS = -L$(LAPACK_DIR)/lib -llapack -lblas

eigen: $(FTOOL) naive.f90 
	$(F90)  $(FTOOL) naive.f90 $(FFLAGS) -o eigen

clean:
	rm -f eigen params.mod subprog.mod
