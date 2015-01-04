all: normal inline
module:
	gfortran -fsyntax-only mersene.f90
normal: module
	gfortran randgen.f bibd.f90 mersene.f90 -o bibd
inline: module
	gfortran -finline-functions randgen.f bibd.f90 mersene.f90 -o bibd_inline
clean:	
	rm bibd bibd_inline mtmod.mod
