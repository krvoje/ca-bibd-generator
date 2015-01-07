all: normal
modules:
	gfortran -fsyntax-only mersene.f90
	gfortran -fsyntax-only incidence_structure.f90
normal: modules
	gfortran incidence_structure.f90 randgen.f bibd_ca.f90 mersene.f90 -o bibd_ca
inline: modules
	gfortran -finline-functions randgen.f bibd_ca.f90 mersene.f90 -o bibd_ca_inline
clean:	
	rm bibd bibd_inline mtmod.mod incidence_structure.mod
