all: modules main
modules:
	mkdir -p obj mod
	gfortran -c mersene.f90 -o obj/mersene.o -Jmod
	gfortran -c utils.f90 -o obj/utils.o -Jmod
	gfortran -c incidence_structure.f90 -o obj/incidence_structure.o -Jmod
	gfortran -c randgen.f -o obj/randgen.o -Jmod
main: modules
	mkdir -p bin/
	mkdir -p mod/
	gfortran obj/*.o bibd_ca.f90 -o bin/bibd_ca -Imod/
clean:	
	rm -rf obj/ bin/ mod/
