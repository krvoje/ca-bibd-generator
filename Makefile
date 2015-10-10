TARGET=target
TARGET_JAVA=target/java
BIN=$(TARGET)
OBJ=$(TARGET)/obj
MOD=$(TARGET)/mod
SRC_FORTRAN=src/fortran
SRC_JAVA=src/java
SWITCHES=-O3
all: java
fortran: modules main
java:
	mkdir -p $(TARGET_JAVA)
	javac $(SRC_JAVA)/org/krvoje/bibd/*.java -d $(TARGET_JAVA)
	echo "Main-Class: org.krvoje.bibd.RandomCABIBD" > $(TARGET_JAVA)/MANIFEST.MF
	jar cvfm $(BIN)/bibd_ca.jar	$(TARGET_JAVA)/MANIFEST.MF -C $(TARGET_JAVA)/ .
modules:
	mkdir -p $(OBJ) $(MOD)
	gfortran -c $(SRC_FORTRAN)/mersene.f90 -o $(OBJ)/mersene.o -J$(MOD) $(SWITCHES)
	gfortran -c $(SRC_FORTRAN)/utils.f90 -o $(OBJ)/utils.o -J$(MOD) $(SWITCHES)
	gfortran -c $(SRC_FORTRAN)/incidence_structure.f90 -o $(OBJ)/incidence_structure.o -J$(MOD) $(SWITCHES)
	gfortran -c $(SRC_FORTRAN)/n_opt.f90 -o $(OBJ)/n_opt.o -J$(MOD) $(SWITCHES)
	gfortran -c $(SRC_FORTRAN)/randgen.f -o $(OBJ)/randgen.o -J$(MOD) $(SWITCHES)
main: modules
	mkdir -p $(BIN)/
	mkdir -p $(MOD)/
	gfortran $(OBJ)/*.o $(SRC_FORTRAN)/bibd_ca.f90 -o $(BIN)/bibd_ca -I$(MOD)/ $(SWITCHES)
clean:	
	rm -rf $(TARGET)
