# Balanced Incomplete Block Design Incidence Matrix Generator

An implementation of an [evolutionary algorithm](https://en.wikipedia.org/wiki/Evolutionary_algorithm) that tries to generate a [Balanced Incomplete Block Design](http://mathworld.wolfram.com/BlockDesign.html) [incidence matrix](https://en.wikipedia.org/wiki/Incidence_matrix). Currently works nice for small matrices :)

To build do:
   
    sbt: assembly

To se a couple of examples on generation:

    sbt "runMain org.krvoje.bibd.Program 7 3 1"
    sbt "runMain org.krvoje.bibd.Program 7 4 2"
    sbt "runMain org.krvoje.bibd.Program 9 3 1"
    sbt "runMain org.krvoje.bibd.Program 13 3 1"
    sbt "runMain org.krvoje.bibd.Program 15 3 1"
    sbt "runMain org.krvoje.bibd.Program 19 3 1"
    sbt "runMain org.krvoje.bibd.Program 21 3 1"
    
Depending on my current fiddling with the algorithm, the larger cases may or may not work.

The project also hosts an older version of the algorithm written in Fortran (builds with gfortran 4.7.2.), but this is outdated.

## So what's this, and what does it do?

This an old university project of mine, which I got hooked on after rediscovering it in an obscure folder of my hard drive. The code tries to employ an evolutionary algorithm to search for an incidence structure that satisfies specific conditions. The search space is represented by column permutations of an incidence matrix, and the population is given by the incidences themselves (the cells of the matrix).

## An evolutionary algorithm

Evolutionary algorithms aim to solve problems by mimicking natural selection. Thus the jargon borrows from biology terms like: *population*, *generation*, *fitness* and *mutation*. An evolutionary algorithm (very roughly) does the following steps:
    
    generate initial *population*
    until STOP_CRITERIA is met
        breed the next generation by some *mutation* procedure
        compute *fitness* of each individual
        based on *fitness* determine which individuals survive

They are a [metaheuristic](https://en.wikipedia.org/wiki/Metaheuristic) (i.e. *guessing*) technique, generally useful for solving problems with a large solution space, where solving them using more exact techniques is infeasible.

## BIBD (Balanced Incomplete Block Design)

A *Balanced Incomplete Block Design* is a particular type of [*incidence structure*](https://en.wikipedia.org/wiki/Incidence_structure) (e.g. graph) that satisfies specific conditions. 

An *incidence structure* consists of *vertices* and *blocks* (e.g. points and lines). *Blocks* contain *vertices*. Two *vertices* are said to be *incident* if a *block* contains them.

In addition, *BIBD* satisfies the properties:
- each *block* contains *k* *vertices*
- each *vertex* is contained in *r* blocks
- each 2 *vertices* are contained together exactly in *lambda* points

The parameters satisfy the equations:
    
    v r = b k
    lambda (v - 1) = r (k - 1)

(*v* is for *vertices*, *b* is for *blocks)

## The actual algorithm

Incidence structures lend themselves well to matrix representation. If we use a (*v* x *b*) matrix where *rows* are *vertices*, and *cols* are *blocks*, we codify incidences (or their absence) using zeros and ones.


                                               b
                                               l
                                               o
                                               c
                                               k
                                               |
                                          0 1 |0| 1 0 0 1
                                          0 0 |1| 1 1 0 0
                                          1 0 |0| 1 0 1 0
                                          0 0 |0| 0 1 1 1
                               vertex ->  1 0 |1| 0 0 0 1
                                          0 1 |1| 0 0 1 0
                                          1 1 |0| 0 1 0 0

The idea of the algorithm is to use a kind of a [cellular automaton](https://en.wikipedia.org/wiki/Cellular_automaton) grid, where the cells are the fields of the incidence matrix.

The algorithm does the following basic steps:
    
    generate a random incidence matrix
    while not BIBD:
        pick two cells in a row
        compute fitness
        if the cells are eligible for change swap them
