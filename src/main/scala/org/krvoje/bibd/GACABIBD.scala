package org.krvoje.bibd

import scala.collection.mutable.ListBuffer
import scala.util.Random

case class GACABIBD(vertices: Int, blocksPerVertex: Int, lambda: Int, logMe: Boolean = true)(implicit val lastChange: LastChange) {

  var generations = 0

  val population = ListBuffer.empty[MutationCABIBD]

  population += MutationCABIBD(vertices, blocksPerVertex, lambda, false)(LastChange(now))
  val proto = population(0)

  val populationSizeSto = Stochasticity(min = 1, max = proto.is.v * proto.is.b * proto.is.r * proto.is.k * proto.is.lambda, increment = 1)
  val heuristicSto = Stochasticity(min = population(0).is.b, max = population(0).is.maxHeuristicDistance, increment = 1)

  log("")
  log(s"(v=${population(0).is.v}, k=${population(0).is.k}, λ=$lambda, b=${population(0).is.b}, r=${population(0).is.r})")

  def findBIBD: IncidenceStructure = {
    while (true) {
      //assert(population.size == maxPopulationSize)
      for (mut <- population) {
        if (checkBIBD(mut.is)) return mut.findBIBD
        if (Random.nextInt(mut.is.maxHeuristicDistance) < heuristicSto.value()) {
          lastChange.timestamp = now
          population -= mut
          if(population.size < populationSizeSto.value()) population += MutationCABIBD(vertices, blocksPerVertex, lambda, false)(LastChange(now))
        } else {
          mut.mutate(checkForStaleness = false)
        }
      }

      generations += 1
    }
    throw new RuntimeException("This part of the code should be unreachable")
  }

  private def checkBIBD(is: IncidenceStructure): Boolean = {
    if(is.isBIBD) {
      //System.out.print(CLEAR_SCREEN)
      log("Generations: " + generations)
      log(s"An incidence matrix for (${is.v}, ${is.k}, ${is.lambda}) found!")
      log(is.toString)

      true
    } else {
      false
    }
  }

  private def log(str: String) = if(logMe) println(this.getClass.getSimpleName + ":" + str)
}
