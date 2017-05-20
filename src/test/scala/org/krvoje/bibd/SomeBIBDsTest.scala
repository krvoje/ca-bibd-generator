package org.krvoje.bibd

import utest._

object SomeBIBDsTest extends TestSuite {

    def tests = TestSuite {
        'SteinerTriplets {
            MutationCABIBD(7, 3, 1).findBIBD
            MutationCABIBD(9, 3, 1).findBIBD

            MutationCABIBD(13, 3, 1).findBIBD
            MutationCABIBD(15, 3, 1).findBIBD

            MutationCABIBD(19, 3, 1).findBIBD
            MutationCABIBD(21, 3, 1).findBIBD

            MutationCABIBD(25, 3, 1).findBIBD
            MutationCABIBD(27, 3, 1).findBIBD

            MutationCABIBD(31, 3, 1).findBIBD
            MutationCABIBD(33, 3, 1).findBIBD
        }

        'gacabibd {
            GACABIBD(7, 3, 1).findBIBD
            GACABIBD(9, 3, 1).findBIBD

            GACABIBD(13, 3, 1).findBIBD
            GACABIBD(15, 3, 1).findBIBD

            GACABIBD(19, 3, 1).findBIBD
            GACABIBD(21, 3, 1).findBIBD

            GACABIBD(25, 3, 1).findBIBD
            GACABIBD(27, 3, 1).findBIBD

            GACABIBD(31, 3, 1).findBIBD
            GACABIBD(33, 3, 1).findBIBD
        }
    }

}
