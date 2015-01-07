# A cellular automaton BIBD incidence matrix generator

An old university project. Tries to generate a BIBD (http://mathworld.wolfram.com/BlockDesign.html) incidence matrix, using a cellular automaton whose cells are the cells of the incidence matrix.

For instance to generate a Fano plane matrix:

    make
    ./bibd_ca 7 3 1 2
