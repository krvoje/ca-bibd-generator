# Balanced incomplete block design incidence matrix generator

An old university project. 

Tries to generate a balanced incomplete block design (http://mathworld.wolfram.com/BlockDesign.html) incidence matrix. Uses a cellular automaton with a help of an n-optimisation algorithm. Works nice for small matrices.

Builds with gfortran 4.7.2.

Usage:

    ./target/bibd_ca v k lambda [optimisation-steps = 2]

For instance to generate a Fano plane matrix, its complement, and a Steiner triplet of order 9:

    make
    # To use a 2-opt improvement on each iteration
    ./target/bibd_ca 7 3 1
    ./target/bibd_ca 7 4 2
    ./target/bibd_ca 9 3 1

    # For a completely random generation
    ./target/bibd_ca 7 3 1 0
    ./target/bibd_ca 7 4 2 0
    ./target/bibd_ca 9 3 1 0
