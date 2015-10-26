# Balanced Incomplete Block Design Incidence Matrix Generator

An implementation of an [evolutionary algorithm](https://en.wikipedia.org/wiki/Evolutionary_algorithm) that tries to generate a *Balanced Incomplete Block Design* (see [BIBD](http://mathworld.wolfram.com/BlockDesign.html) ) incidence matrix. Currently works nice for small matrices :)

Requires JDK 1.6 to build. To generate a couple of examples:

    make    
    java -jar ./target/bibd_ca.jar 7 3 1
    java -jar ./target/bibd_ca.jar 7 4 2
    java -jar ./target/bibd_ca.jar 9 3 1
    java -jar ./target/bibd_ca.jar 13 3 1
    java -jar ./target/bibd_ca.jar 15 3 1
    java -jar ./target/bibd_ca.jar 19 3 1
    java -jar ./target/bibd_ca.jar 21 3 1

The project also hosts an older version of the algorithm written in Fortran (builds with gfortran 4.7.2.)
