TARGET=target
BIN=$(TARGET)
OBJ=$(TARGET)/obj
MOD=$(TARGET)/mod
SRC=src
all: modules main
modules:
	mkdir -p $(OBJ) $(MOD)
	gfortran -c $(SRC)/mersene.f90 -o $(OBJ)/mersene.o -J$(MOD)
	gfortran -c $(SRC)/utils.f90 -o $(OBJ)/utils.o -J$(MOD)
	gfortran -c $(SRC)/incidence_structure.f90 -o $(OBJ)/incidence_structure.o -J$(MOD)
	gfortran -c $(SRC)/randgen.f -o $(OBJ)/randgen.o -J$(MOD)
main: modules
	mkdir -p $(BIN)/
	mkdir -p $(MOD)/
	gfortran $(OBJ)/*.o $(SRC)/bibd_ca.f90 -o $(BIN)/bibd_ca -I$(MOD)/
clean:	
	rm -rf $(TARGET)
