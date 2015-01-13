# Balanced incomplete block design incidence matrix generator

An old university project. 

Tries to generate a balanced incomplete block design (http://mathworld.wolfram.com/BlockDesign.html) incidence matrix. Uses a cellular automaton with a help of an n-optimisation algorithm. Works nice for small matrices.

Builds with gfortran 4.7.2.

Usage:
    ./bin/bibd_ca v k lambda optimisation-steps

For instance to generate a Fano plane matrix:

    make
    # To use a 2-opt improvement on each iteration
    ./bin/bibd_ca 7 3 1 2
    # For a completely random generation
    ./bin/bibd_ca 7 3 1 0
