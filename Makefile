TARGET=target
BIN=$(TARGET)/bin
OBJ=$(TARGET)/obj
MOD=$(TARGET)/mod
FORTRAN_SRC=src/main/fortran
all: modules main
modules:
	mkdir -p $(OBJ) $(MOD)
	gfortran -c $(FORTRAN_SRC)/mersene.f90 -o $(OBJ)/mersene.o -J$(MOD)
	gfortran -c $(FORTRAN_SRC)/utils.f90 -o $(OBJ)/utils.o -J$(MOD)
	gfortran -c $(FORTRAN_SRC)/incidence_structure.f90 -o $(OBJ)/incidence_structure.o -J$(MOD)
	gfortran -c $(FORTRAN_SRC)/randgen.f -o $(OBJ)/randgen.o -J$(MOD)
main: modules
	mkdir -p $(BIN)/
	mkdir -p $(MOD)/
	gfortran $(OBJ)/*.o $(FORTRAN_SRC)/bibd_ca.f90 -o $(BIN)/bibd_ca -I$(MOD)/
clean:	
	rm -rf $(TARGET)
