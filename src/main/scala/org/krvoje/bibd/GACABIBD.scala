package org.krvoje.bibd

import scala.collection.mutable.ListBuffer
import scala.util.Random

case class GACABIBD(vertices: Int, blocksPerVertex: Int, lambda: Int) {

  var generations = 0
  val populationSize = vertices * blocksPerVertex

  val population = ListBuffer.empty[MatrixIncidenceStructure]
  for(i <- 0 until populationSize)
    population += MatrixIncidenceStructure(vertices, blocksPerVertex, lambda)

  implicit val lastChange = Moment(System.currentTimeMillis())
  val heuristicSto = Stochasticity(min = 1, max = population(0).maxHeuristicDistance, increment = 1)

  println()
  println(s"(v=${population(0).v}, k=${population(0).k}, λ=$lambda, b=${population(0).b}, r=${population(0).r})")

  def findBIBD: IncidenceStructure = {
    while (true) {
      assert(population.size == populationSize)
      for (is <- population) {
        if (checkBIBD(is)) return is
        if (Random.nextInt(is.maxHeuristicDistance) < heuristicSto.value()) {
          lastChange.when = System.currentTimeMillis()
          population -= is
          population += MatrixIncidenceStructure(vertices, blocksPerVertex, lambda)
        }
      }

      generations += 1
    }
    throw new RuntimeException("This part of the code should be unreachable")
  }

  private def checkBIBD(is: IncidenceStructure): Boolean = {
    if(is.isBIBD) {
      //System.out.print(CLEAR_SCREEN)
      System.out.println("Generations: " + generations)
      System.out.println(s"An incidence matrix for (${is.v}, ${is.k}, ${is.lambda}) found!")
      System.out.print(is)

      true
    } else {
      false
    }
  }
}
