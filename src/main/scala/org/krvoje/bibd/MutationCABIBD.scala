package org.krvoje.bibd

import scala.util.Random

case class MutationCABIBD(vertices: Integer, blocksPerVertex: Integer, lambda: Integer) {

  val rnd = Random

  var maxChangeFactor: BigInt = 0
  var minChangeFactor: BigInt = 0
  var staleCellsCount: BigInt = 0
  var generations: BigInt = 0
  var iterations: BigInt = 0
  var unchanged: BigInt = 0
  implicit val lastNonStaleChange = Moment(System.currentTimeMillis())
  var lastChange = System.currentTimeMillis()

  val is = MatrixIncidenceStructure(vertices, blocksPerVertex, lambda)
  val changeFactor = Array.ofDim[Int](is.v, is.b)

  println()
  println(s"(v=${is.v}, k=${is.k}, λ=$lambda, b=${is.b}, r=${is.r})")

  val unchangedTreshold = Stochasticity(min = 1, max = is.v * is.b * is.r * is.k * is.lambda, increment = 1)
  val stalenessTreshold = Stochasticity(min = 1, max = is.v * is.b, increment = 1)

  calculateChangeFactors

  // If this never halts, the program does not halt
  def findBIBD: IncidenceStructure = {
    var activeRow = 0
    var dormantRow = 0
    var cfa = 0
    var cfd = 0
    var rcf = 0
    var doWeChange = false

    var col = -1
    while(true)
    {
      col = (col + 1) % is.b
      //col = rnd.nextInt(is.b)

      activeRow = randomActiveInCol(is, col)
      dormantRow = randomDormantInCol(is, col)

      if(is.isBIBD) {
        //System.out.print(CLEAR_SCREEN)
        System.out.println("Generations: " + generations)
        System.out.println("Iterations: " + iterations)
        System.out.println(s"An incidence matrix for (${is.v}, ${is.k}, ${is.lambda}) found!")
        System.out.print(is)

        return this.is; // Sparta
      }

      this.iterations += 1

      cfa = changeFactor(activeRow)(col)
      cfd = changeFactor(dormantRow)(col)
      rcf = rnd.nextInt(max(1,maxChangeFactor).intValue())

      doWeChange = (rcf < cfa && rcf < cfd) || stale

      if (doWeChange) {
        is.flip(activeRow, col)
        is.flip(dormantRow, col)
        calculateChangeFactors
        unchanged = 0
        generations += 1
        if(!stale) lastNonStaleChange.when = System.currentTimeMillis()
        lastChange = System.currentTimeMillis()
      }
      else {
        unchanged += 1
      }
    }

    // Never reached
    throw new RuntimeException("This part of the code should be unreachable")
  }

  def calculateChangeFactors: Unit = {
    maxChangeFactor = 0
    staleCellsCount = 0
    forIndex(0, is.v) { row =>
      forIndex(0, is.b) { col =>
        changeFactor(row)(col) = calculateChangeFactor(row, col)
        maxChangeFactor = max(changeFactor(row)(col), maxChangeFactor)
        minChangeFactor = min(changeFactor(row)(col), minChangeFactor)
        if(changeFactor(row)(col) <= 0) staleCellsCount += 1
      }
    }
  }

  def calculateChangeFactor(row: Int, col: Int): Int = {
    var changeFactor = 0
    var delta = 0

    // Max increment of: (v-1) * r - lambda
    forIndex(0, is.v) { otherRow =>
      if(otherRow != row) {
        delta = math.abs(is.lambda - is.rowIntersection(row, otherRow))
        if (is.lambda < is.rowIntersection(row, otherRow)) {
          if (is.active(row, col) && is.active(otherRow, col))
            changeFactor += delta
          else if (is.active(row, col) && is.dormant(otherRow, col))
            changeFactor -= delta
        }
        if (is.lambda > is.rowIntersection(row, otherRow)) {
          if (is.dormant(row, col) && is.active(otherRow, col))
            changeFactor += delta
          else if (is.dormant(row, col) && is.dormant(otherRow, col))
            changeFactor -= delta
        }
      }
    }

    // Max increment of: (b - k)
    delta = is.k - is.sumInCol(col)
    if(is.active(row,col))
      changeFactor -= delta
    else if(is.dormant(row,col))
      changeFactor += delta

    // Max increment of: (v - r)
    delta = is.r - is.sumInRow(row)
    if(is.active(row,col))
      changeFactor -= delta
    else if(is.dormant(row,col))
      changeFactor += delta

    return changeFactor
  }

  def randomActiveInCol(is: IncidenceStructure, col: Int): Int = {
    var row = rnd.nextInt(is.v)
    while(!is.active(row, col))
      row = rnd.nextInt(is.v)
    return row
  }

  def randomDormantInCol(is: IncidenceStructure, col: Int): Int = {
    var row = rnd.nextInt(is.v)
    while(!is.dormant(row, col))
      row = rnd.nextInt(is.v)
    return row
  }

  private def stale: Boolean = {
    assert(staleCellsCount <= is.v * is.b)
    unchanged > unchangedTreshold.value() &&
      staleCellsCount > stalenessTreshold.value()
  }
}