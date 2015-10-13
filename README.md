# Balanced incomplete block design incidence matrix generator

An old university project. 

Tries to generate a balanced incomplete block design (http://mathworld.wolfram.com/BlockDesign.html) incidence matrix, using a cellular automaton grid. Works nice for small matrices :)

Java and Fortran.

The Fortran version builds with gfortran 4.7.2., and currently hosts an older version of the algorithm.

Usage:

    ./target/bibd_ca v k lambda

For instance to generate a couple of examples:

    make    
    java -jar ./target/bibd_ca.jar 7 3 1
    java -jar ./target/bibd_ca.jar 7 4 2
    java -jar ./target/bibd_ca.jar 9 3 1
    java -jar ./target/bibd_ca.jar 13 3 1
    java -jar ./target/bibd_ca.jar 15 3 1
