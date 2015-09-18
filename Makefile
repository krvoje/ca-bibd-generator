TARGET=target
BIN=$(TARGET)
OBJ=$(TARGET)/obj
MOD=$(TARGET)/mod
SRC=src
SWITCHES=-O3
all: modules main
modules:
	mkdir -p $(OBJ) $(MOD)
	gfortran -c $(SRC)/mersene.f90 -o $(OBJ)/mersene.o -J$(MOD) $(SWITCHES)
	gfortran -c $(SRC)/utils.f90 -o $(OBJ)/utils.o -J$(MOD) $(SWITCHES)
	gfortran -c $(SRC)/incidence_structure.f90 -o $(OBJ)/incidence_structure.o -J$(MOD) $(SWITCHES)
	gfortran -c $(SRC)/n_opt.f90 -o $(OBJ)/n_opt.o -J$(MOD) $(SWITCHES)
	gfortran -c $(SRC)/randgen.f -o $(OBJ)/randgen.o -J$(MOD) $(SWITCHES)
	gfortran -c $(SRC)/tabu_list.f90 -o $(OBJ)/tabu_list.o -J$(MOD) $(SWITCHES)
main: modules
	mkdir -p $(BIN)/
	mkdir -p $(MOD)/
	gfortran $(OBJ)/*.o $(SRC)/bibd_ca.f90 -o $(BIN)/bibd_ca -I$(MOD)/ $(SWITCHES)
clean:	
	rm -rf $(TARGET)
