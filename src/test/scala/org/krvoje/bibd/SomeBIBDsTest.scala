package org.krvoje.bibd

import org.scalatest.FeatureSpec

class SomeBIBDsTest extends FeatureSpec {

    implicit val lastChange = LastChange(now)

    feature("Steiner triplets") {
        mutation(7, 3, 1)
        mutation(9, 3, 1)
        mutation(13, 3, 1)
        mutation(15, 3, 1)
        mutation(19, 3, 1)
        mutation(21, 3, 1)
        mutation(25, 3, 1)
        mutation(27, 3, 1)
        mutation(31, 3, 1)
        mutation(33, 3, 1)
        ga(7, 3, 1)
        ga(9, 3, 1)
        ga(13, 3, 1)
        ga(15, 3, 1)
        ga(19, 3, 1)
        ga(21, 3, 1)
        ga(25, 3, 1)
        ga(27, 3, 1)
        ga(31, 3, 1)
        ga(33, 3, 1)
    }

    private def mutation(v: Int, k: Int, lambda: Int) =
        scenario(s"(Mutation) Steiner triplet 2-($v, $k, $lambda)") {MutationCABIBD(v,k,lambda).findBIBD}
    private def ga(v: Int, k: Int, lambda: Int) =
        scenario(s"(GA) Steiner triplet 2-($v, $k, $lambda)") {GACABIBD(v,k,lambda).findBIBD}
}
