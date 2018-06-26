package org.krvoje.bibd

import org.scalatest.FeatureSpec

class GABIBDTest extends FeatureSpec {

    //implicit val lastChange = ReferenceFrame(now)

    feature("Steiner triplets") {
        ga(7, 3, 1)
        ga(9, 3, 1)
        ga(13, 3, 1)
        ga(15, 3, 1)
        //ga(19, 3, 1)
        //ga(21, 3, 1)
        //ga(25, 3, 1)
        //ga(27, 3, 1)
        //ga(31, 3, 1)
        //ga(33, 3, 1)
    }

    private def ga(v: Int, k: Int, lambda: Int) =
        scenario(s"(GA) Steiner triplet 2-($v, $k, $lambda)") {GACABIBD(v,k,lambda).findBIBD}
}
