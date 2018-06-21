package org.krvoje.bibd

import scala.collection.mutable.ListBuffer

case class GACABIBD(vertices: Int,
  blocksPerVertex: Int,
  lambda: Int,
  logMe: Boolean = true)(implicit val lastChange: ReferenceFrame = ReferenceFrame(now)) {

  var generations = 0

  var population = ListBuffer.empty[MutationCABIBD]

  val rf = ReferenceFrame(now)
  for(i <- 0 until 10)
    population += MutationCABIBD(vertices = vertices,
      blocksPerVertex = blocksPerVertex,
      lambda = lambda,
      logMe = false)(ReferenceFrame(now))

  val proto = population(0)

  log("")
  log(s"(v=${population(0).is.v}, k=${population(0).is.k}, λ=$lambda, b=${population(0).is.b}, r=${population(0).is.r})")

  def findBIBD: IncidenceStructure = {
    while (true) {
      //assert(population.size == maxPopulationSize)
      population.foreach {mutant =>
        if (checkBIBD(mutant.is)) return mutant.findBIBD
        mutant.mutate(checkForStaleness = true)
      }

      population = population.map {
        mutant =>
          if(mutant.stale) {
            MutationCABIBD(vertices, blocksPerVertex, lambda, false)(ReferenceFrame(now))
          } else {
            mutant
          }
      }

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
