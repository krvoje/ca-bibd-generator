MODULE rnd
    
    CONTAINS

        REAL FUNCTION rand(x)
            INTEGER x
            REAL rendom

            CALL init_random_seed()
            CALL random_number(rendom)

            rand=rendom*x
            RETURN
        END FUNCTION

        REAL FUNCTION randint(x)
            INTEGER x
            REAL rendom

            CALL init_random_seed()
            CALL random_number(rendom)

            randint=int(rendom*x)
            RETURN
        END FUNCTION

        SUBROUTINE init_random_seed()
                    INTEGER :: i, n, clock
                    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

                    CALL RANDOM_SEED(size = n)
                    ALLOCATE(seed(n))

                    CALL SYSTEM_CLOCK(COUNT=clock)

                    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
                    CALL RANDOM_SEED(PUT = seed)

                    DEALLOCATE(seed)
        END subroutine

END MODULE
