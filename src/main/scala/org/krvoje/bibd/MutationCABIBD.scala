package org.krvoje.bibd

import scala.util.Random

case class MutationCABIBD(
  vertices: Integer,
  blocksPerVertex: Integer,
  lambda: Integer,
  logMe: Boolean = true
)(implicit val rf: ReferenceFrame = ReferenceFrame(now)) {

  val is = MatrixIncidenceStructure(vertices, blocksPerVertex, lambda)

  val maxPossibleChangeFactor: BigInt = (is.v-1) * is.r - lambda + (is.b-is.k) + (is.v-is.r)
  var maxChangeFactor: BigInt = 0
  var minChangeFactor: BigInt = 0
  var staleCellsCount: BigInt = 0

  var generations: BigInt = 0
  val changeFactor = Array.ofDim[Int](is.v, is.b)

  log("")
  log(s"(v=${is.v}, k=${is.k}, λ=$lambda, b=${is.b}, r=${is.r})")


  /** How long we wait between two changes (iterations) */
  val changeInterval = Stochasticity(
    min = 1
    , max = is.v * is.b //* is.r * is.k * is.lambda
    , changeInterval = is.b
  )

  /** How much of iterations unchanged means a matrix is stale */
  val unchangedTreshold = Stochasticity(
    min = is.b
    , max = is.v * is.b * is.v//is.r * is.k * is.lambda
    , changeInterval = changeInterval.copy()
  )

  /** How many stale cells means that a matrix is stale (stale means that they are unable to change)*/
  val stalenessTreshold = Stochasticity(
    min = is.b
    , max = is.v * is.b
    , changeInterval = changeInterval.copy()
  )

  // If this never halts, the program does not halt
  def findBIBD: IncidenceStructure = {
    calculateChangeFactors

    var col = -1
    while(true)
    {
      if (is.isBIBD) {
        log("Generations: " + rf.generation)
        log("Iterations: " + rf.currentIteration)
        log("Stale resets: " + rf.staleResets)
        log(s"An incidence matrix for (${is.v}, ${is.k}, ${is.lambda}) found!")
        log(is.toString)

        return this.is; // Sparta
      }

      col = (col + 1) % is.b
      mutate(col)
    }

    // Never reached
    throw new RuntimeException("This part of the code should be unreachable")
  }

  def mutate(checkForStaleness: Boolean): Unit = {
    for(col <- 0 until is.b) mutate(col, checkForStaleness)
  }

  def mutate(col: Int, checkForStaleness: Boolean = true): Unit = {
    calculateChangeFactors
    val activeRow = randomActiveInCol(is, col)
    val dormantRow = randomDormantInCol(is, col)

    val cfa = changeFactor(activeRow)(col)
    val cfd = changeFactor(dormantRow)(col)
    val rcf = Random.nextInt(max(1, maxChangeFactor).intValue())

    val ripeForChange = (rcf < cfa && rcf < cfd)
    val isStale = checkForStaleness && stale
    val doWeChange = ripeForChange || isStale

    if (doWeChange) {
      is.flip(activeRow, col)
      is.flip(dormantRow, col)
      calculateChangeFactors
      rf.changed(!ripeForChange && isStale)
    }
    rf.iterate
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
    var row = Random.nextInt(is.v)
    while(!is.active(row, col)) {
      row = Random.nextInt(is.v)
    }
    row
  }

  def randomDormantInCol(is: IncidenceStructure, col: Int): Int = {
    var row = Random.nextInt(is.v)
    while(!is.dormant(row, col)){
      row = Random.nextInt(is.v)
    }
    row
  }

  def stale: Boolean = {
    rf.unchanged >= unchangedTreshold && staleCellsCount > stalenessTreshold
  }

  private def log(str: String) = if(logMe) println(this.getClass.getSimpleName + ":" + str)
}