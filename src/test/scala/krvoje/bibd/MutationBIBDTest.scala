package krvoje.bibd

import org.scalatest.FeatureSpec

class MutationBIBDTest extends FeatureSpec {

    //implicit val lastChange = ReferenceFrame(now)

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
    }

    private def mutation(v: Int, k: Int, lambda: Int): Unit =
        scenario(s"(Mutation) Steiner triplet 2-($v, $k, $lambda)") {MutationCABIBD(v,k,lambda).findBIBD}
}
