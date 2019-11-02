package krvoje.bibd

import org.scalatest.FeatureSpec

class MutationBIBDTest extends FeatureSpec {
  feature("Steiner triplets") {
    steiner(7)
    steiner(9)
    steiner(13)
    steiner(15)
    steiner(19)
  }
  feature("Hadamard plane") {
    hadamard(1)
    hadamard(2)
    hadamard(3)
  }
  feature("Affine Plane") {
    affine(2)
    affine(3)
    affine(4)
  }
  feature("Projective plane") {
    projective(1)
    projective(2)
    projective(3)
  }

  private def affine(n: Int): Unit = mutation("Affine plane")(n * n, n, 1)
  private def hadamard(n: Int): Unit = mutation("Hadamard design")(4*n+3, 2*n+1, n)
  private def projective(n: Int): Unit = mutation("Projective plane")(n*n + n + 1, n + 1, 1)
  private def steiner(n: Int): Unit = mutation("Steiner triplet")(n, 3, 1)

  private def mutation(name: String)(v: Int, k: Int, lambda: Int): Unit =
    scenario(s"(Mutation) $name 2-($v, $k, $lambda)") {MutationCABIBD(v,k,lambda).findBIBD}
}
