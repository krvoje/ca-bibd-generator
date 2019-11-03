package krvoje.bibd

import java.io.{BufferedWriter, File, FileWriter}

import krvoje.bibd
import Implicits._

import scala.util.Random

case class MutationCABIBD(
  vertices: Integer,
  blocksPerVertex: Integer,
  lambda: Integer,
  isLogged: Boolean = true,
  isLoggedInFile: Boolean = false,
)(implicit val rf: ReferenceFrame = ReferenceFrame(now())) {

  val is = new IncidenceStructure(vertices, blocksPerVertex, lambda)
  val globalSB = new StringBuilder()

  val maxPossibleChangeFactor: Int = (is.v-1) * is.r - lambda + (is.b-is.k) + (is.v-is.r)
  var maxChangeFactor: Int = 0
  var minChangeFactor: Int = 0

  var generations: Int = 0
  val changeFactor: Array[Array[Int]] = Array.ofDim[Int](is.v, is.b)

  log("")
  log(s"(v=${is.v}, k=${is.k}, λ=$lambda, b=${is.b}, r=${is.r})")

  val unchangedSto = Stochasticity(
    min = is.v * is.b * 2,
    max = is.v * is.b * 3,
    increment = Stochasticity(1, is.v * is.b),
    changeThreshold = Stochasticity(1, is.b),
  )

  // If this never halts, the program does not halt
  def findBIBD: IncidenceStructure = {
    calculateChangeFactors()

    var col = -1
    while(true)
    {
      if (is.isBIBD) {
        log("Generations: " + rf.generation)
        log("Iterations: " + rf.currentIteration)
        log("Stale resets: " + rf.staleResets)
        log(s"An incidence matrix for 2-(${is.v}, ${is.k}, ${is.lambda}) found!")
        log(unchangedSto.value().toString)
        log(is.toString)

        if (isLoggedInFile) writeToFile(s"2-(${is.v}, ${is.k}, ${is.lambda}).tmp")

        return this.is; // Sparta
      }

      col = (col + 1) % is.b
      mutate(col)
    }

    // Never reached
    throw new RuntimeException("This part of the code should be unreachable")
  }

  def mutate(): Unit = {
    for(col <- 0 until is.b) mutate(col)
  }

  def mutate(col: Int): Unit = {
    calculateChangeFactors()
    val activeRow = randomActiveInCol(is, col)
    val dormantRow = randomDormantInCol(is, col)

    val cfa = changeFactor(activeRow)(col)
    val cfd = changeFactor(dormantRow)(col)
    val rcf = Random.nextInt(max(1, maxChangeFactor+1).intValue)

    val ripeForChange = (rcf < cfa && rcf < cfd)
    val isStale = isMatrixStale
    val doWeChange = ripeForChange || isStale

    if (doWeChange) {
      is.flip(activeRow, col)
      is.flip(dormantRow, col)
      calculateChangeFactors()
      rf.changed(!ripeForChange && isStale)
    }
    rf.iterate()
  }

  def calculateChangeFactors(): Unit = {
    maxChangeFactor = 0
    minChangeFactor = 0
    forIndex(0, is.v) { row =>
      forIndex(0, is.b) { col =>
        changeFactor(row)(col) = calculateChangeFactor(row, col)
        maxChangeFactor = max(changeFactor(row)(col), maxChangeFactor)
        minChangeFactor = min(changeFactor(row)(col), minChangeFactor)
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

  def isMatrixStale: Boolean = {
    rf.unchanged > unchangedSto
  }

  private def log(str: String): String = if(isLogged) {
    val msg: String = this.getClass.getSimpleName + ":" + str
    println(msg)
    globalSB.append(s"$msg\n")
    msg
  } else {
    ""
  }

  private def writeToFile(filename: String): Unit = {
    val file = new File(s"target/$filename")
    val bw = new BufferedWriter(new FileWriter(file))
    bw.write(globalSB.toString)
    bw.close()
  }
}
