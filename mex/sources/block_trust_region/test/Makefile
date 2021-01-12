FC = gfortran
FCFLAGS = -g -Wall

all: dulmage_mendelsohn_test trust_region_test

dulmage_mendelsohn_test: dulmage_mendelsohn.o dulmage_mendelsohn_test.o
	$(FC) $(FCFLAGS) $^ -o $@

dulmage_mendelsohn_test.o: dulmage_mendelsohn_test.f08 dulmage_mendelsohn.mod
	$(FC) $(FCFLAGS) -c $< -o $@

dulmage_mendelsohn.mod: dulmage_mendelsohn.o

dulmage_mendelsohn.o: ../dulmage_mendelsohn.f08
	$(FC) $(FCFLAGS) -c $< -o $@

trust_region_test: trust_region.o trust_region_test.o
	$(FC) $(FCFLAGS) $^ -o $@ -llapack -lblas

trust_region_test.o: trust_region_test.f08 trust_region.mod
	$(FC) $(FCFLAGS) -c $< -o $@

trust_region.mod: trust_region.o

trust_region.o: ../trust_region.f08 lapack.mod
	$(FC) $(FCFLAGS) -c $< -o $@

blas_lapack.o: ../../blas_lapack.F08
	$(FC) $(FCFLAGS) -c $< -o $@

lapack.mod: blas_lapack.o

clean:
	rm -f *.o *.mod dulmage_mendelsohn_test trust_region_test
