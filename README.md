# Balanced incomplete block design incidence matrix generator

An old university project. 

Tries to generate a balanced incomplete block design (http://mathworld.wolfram.com/BlockDesign.html) incidence matrix. Uses a cellular automaton with a help of an n-optimisation algorithm. Works nice for small matrices.

Builds with gfortran 4.7.2.

Usage:

    ./target/bin/bibd_ca v k lambda [optimisation-steps = 2]

For instance to generate a Fano plane matrix:

    make
    # To use a 2-opt improvement on each iteration
    ./target/bin/bibd_ca 7 3 1
    # For a completely random generation
    ./target/bin/bibd_ca 7 3 1 0
