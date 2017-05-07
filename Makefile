# your FORTRAN compiler
FCOMP=gfortran
F90=$(FCOMP)

# set here link reference needed for BLAS & LAPACK
BLASLAPACK=-L/usr/lib -llapack -lblas  

# compilation options
opt = -Ofast

# JSONFORTRAN lib
jsonfortrandir=lib/jsonfortran
jsonfortran=-L$(jsonfortrandir) -ljsonfortran
jsonfortranlist=$(jsonfortrandir)/json_module.o
# AGMG lib
agmgdir=lib/agmg
listagmg=$(agmgdir)/dagmg.o $(agmgdir)/dagmg_mumps.o
#
lib=$(BLASLAPACK) $(jsonfortran)
export FCOMP
export F90
export opt
externallist=$(listagmg)
locallist=types.o common_func.o error_handler.o iovtk.o structure_generator.o fvm_iem_matrix.o load_configuration.o process_results.o



diffusion_model: $(locallist) main.o
	cd $(agmgdir);make dseq
	$(FCOMP) $(opt) -o diffusion_model $(locallist) main.o $(externallist) $(lib) 

test_simple: $(locallist) test_simple.o
	cd $(agmgdir);make dseq
	$(FCOMP) $(opt) -o test_simple $(locallist) test_simple.o $(externallist) $(lib) 

test_simple.o: test_simple.f95
	$(FCOMP) $(opt) -c test_simple.f95
main.o: main.f95
	$(FCOMP) $(opt) -c main.f95

load_configuration.o: load_configuration.f95
	$(FCOMP) $(opt) -c load_configuration.f95 -I$(jsonfortrandir)
	
iovtk.o: iovtk.f95
	$(FCOMP) $(opt) -c iovtk.f95

structure_generator.o: structure_generator.f95
	$(FCOMP) $(opt) -c structure_generator.f95

fvm_iem_matrix.o: fvm_iem_matrix.f95
	$(FCOMP) $(opt) -c fvm_iem_matrix.f95

process_results.o: process_results.f95
	$(FCOMP) $(opt) -c process_results.f95
	
error_handler.o: error_handler.f95
	$(FCOMP) $(opt) -c error_handler.f95

common_func.o: common_func.f95
	$(FCOMP) $(opt) -c common_func.f95

types.o: types.f95
	$(FCOMP) $(opt) -c types.f95
	
clean:
	rm -f *.o *.mod
	rm -f diffusion_model
	rm -f test_simple
