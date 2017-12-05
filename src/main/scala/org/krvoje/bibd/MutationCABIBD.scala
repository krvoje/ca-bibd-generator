package org.krvoje.bibd

import scala.util.Random

case class MutationCABIBD(vertices: Integer, blocksPerVertex: Integer, lambda: Integer, logMe: Boolean = true)(implicit val lastChange: LastChange) {

  val rnd = Random

  var maxChangeFactor: BigInt = 0
  var minChangeFactor: BigInt = 0
  var staleCellsCount: BigInt = 0
  var generations: BigInt = 0
  var iterations: BigInt = 0
  var unchanged: BigInt = 0

  val is = MatrixIncidenceStructure(vertices, blocksPerVertex, lambda)
  val changeFactor = Array.ofDim[Int](is.v, is.b)

  log("")
  log(s"(v=${is.v}, k=${is.k}, λ=$lambda, b=${is.b}, r=${is.r})")


  val inc = Stochasticity(1, 500)
  val blocksInc = Stochasticity(1, is.b, 1)
  val blocksSto = Stochasticity(1, is.v * is.b, increment = blocksInc.copy().value)
  val unchangedTreshold = Stochasticity(1, is.v * is.b * is.r * is.k * is.lambda, increment = blocksInc.copy().value, stochasticity = blocksSto.copy())
  val stalenessTreshold = Stochasticity(1, is.v * is.b, increment = blocksInc.copy().value, stochasticity = inc.copy())

  // If this never halts, the program does not halt
  def findBIBD: IncidenceStructure = {
    calculateChangeFactors

    var col = -1
    while(true)
    {
      if (is.isBIBD) {
        log("Generations: " + generations)
        log("Iterations: " + iterations)
        log(s"An incidence matrix for (${is.v}, ${is.k}, ${is.lambda}) found!")
        log(is.toString)

        return this.is; // Sparta
      }

      col = (col + 1) % is.b
      mutate(col)
      //col = rnd.nextInt(is.b)
    }

    // Never reached
    throw new RuntimeException("This part of the code should be unreachable")
  }

  def mutate(checkForStaleness: Boolean): Unit = {
    for(col <- 0 until is.b) mutate(col, checkForStaleness)
  }

  def mutate(col: Int, checkForStaleness: Boolean = true): Unit = {
    val activeRow = randomActiveInCol(is, col)
    val dormantRow = randomDormantInCol(is, col)

    this.iterations += 1

    val cfa = changeFactor(activeRow)(col)
    val cfd = changeFactor(dormantRow)(col)
    val rcf = rnd.nextInt(max(1, maxChangeFactor).intValue())

    val doWeChange = (rcf < cfa && rcf < cfd) || checkForStaleness && stale

    if (doWeChange) {
      is.flip(activeRow, col)
      is.flip(dormantRow, col)
      calculateChangeFactors
      unchanged = 0
      generations += 1
      if (!stale) lastChange.timestamp = now
    }
    else {
      unchanged += 1
    }
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

  def stale: Boolean = {
    assert(staleCellsCount <= is.v * is.b)
    unchanged > unchangedTreshold.value() &&
      staleCellsCount > stalenessTreshold.value()
  }

  private def log(str: String) = if(logMe) println(this.getClass.getSimpleName + ":" + str)
}