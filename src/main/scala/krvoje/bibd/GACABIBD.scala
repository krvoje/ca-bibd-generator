package krvoje.bibd

import krvoje.bibd.util.Implicits._

import scala.collection.mutable.ListBuffer

case class GACABIBD(vertices: Int,
  blocksPerVertex: Int,
  lambda: Int,
  logMe: Boolean = true)(implicit val lastChange: ReferenceFrame = ReferenceFrame(now)) {

  var generations = 0

  var population = ListBuffer.empty[MutationCABIBD]

  val proto = MutationCABIBD(vertices = vertices,
    blocksPerVertex = blocksPerVertex,
    lambda = lambda,
    logMe = false)

  val MAX_POP_SIZE = proto.is.v

  val rf = ReferenceFrame(now)
  for(i <- 0 until MAX_POP_SIZE)
    population += MutationCABIBD(vertices = vertices,
      blocksPerVertex = blocksPerVertex,
      lambda = lambda,
      logMe = false)(ReferenceFrame(now))

  log("")
  log(s"(v=${population(0).is.v}, k=${population(0).is.k}, λ=$lambda, b=${population(0).is.b}, r=${population(0).is.r})")

  def findBIBD: IncidenceStructure = {
    var found: Option[IncidenceStructure] = None
    while (found.isEmpty) {
      //assert(population.size == maxPopulationSize)
      population.foreach {mutant =>
        if (checkBIBD(mutant.is)) {
          found = mutant.findBIBD
        }
        mutant.mutate(isCheckForStaleness = true)
      }
      if(found.nonEmpty) return found.get
    }
    throw new RuntimeException("This part of the code should be unreachable")
  }

  def kill(mutant: MutationCABIBD) = {
    population -= mutant
  }
  def birth(ancestors: Seq[MutationCABIBD]) = {
    population += MutationCABIBD(vertices, blocksPerVertex, lambda, false)(ReferenceFrame(now))
  }

  private def checkBIBD(is: IncidenceStructure): Boolean = {
    if(is.isBIBD) {
      //System.out.print(CLEAR_SCREEN)
      log(s"Generations: $generations. Population: ${population.size}")
      log(s"An incidence matrix for (${is.v}, ${is.k}, ${is.lambda}) found!")
      log(is.toString)

      true
    } else {
      false
    }
  }

  private def log(str: String) = if(logMe) println(this.getClass.getSimpleName + ":" + str)
}
