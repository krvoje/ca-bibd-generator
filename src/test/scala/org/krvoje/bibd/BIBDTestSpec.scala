package org.krvoje.bibd

import org.specs2.mutable.Specification

class BIBDTestSpec extends Specification {

  "Steiner Triplets" >> {
    "They should halt" >> {
      RandomCABIBD.main(Array("7", "3", "1"))
      RandomCABIBD.main(Array("9", "3", "1"))

      RandomCABIBD.main(Array("13", "3", "1"))
      RandomCABIBD.main(Array("15", "3", "1"))

      RandomCABIBD.main(Array("19", "3", "1"))
      RandomCABIBD.main(Array("21", "3", "1"))

      RandomCABIBD.main(Array("25", "3", "1"))
      RandomCABIBD.main(Array("27", "3", "1"))

      RandomCABIBD.main(Array("31", "3", "1"))
      RandomCABIBD.main(Array("33", "3", "1"))
      true === true
    }
  }

}
