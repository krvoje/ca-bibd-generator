package org.krvoje.bibd

import scala.util.Random

object RandomCABIBD extends App {

  val rnd = Random

  var maxChangeFactor =0
  var minChangeFactor =0
  var generations = 0
  var iterations = 0
  var animate = false

  val CLEAR_SCREEN = "\u001B[2J"

  if(args.length != 3 && args.length != 4) {
    System.out.println("Usage: ")
    System.out.println("java -jar ca-bibd.jar v k lambda [animate]")
    System.exit(0)
  }

  val vertices = Integer.parseInt(args(0))
  val blocksPerVertex = Integer.parseInt(args(1))
  val lambda = Integer.parseInt(args(2))
  if(args.length == 4) {
    if(args(3) == "animate") {
      animate = true
    }
  }

  val is = new MatrixIncidenceStructure(vertices, blocksPerVertex, lambda)
  val changeFactor = Array.ofDim[Int](is.v, is.b)
  randomize(is)
  getBibd // If this never halts, the program does not halt

  //System.out.print(CLEAR_SCREEN)
  System.out.println("Generations: " + generations)
  System.out.println("Iterations: " + iterations)
  System.out.println(s"An incidence matrix for (${is.v}, ${is.k}, ${is.lambda}) found!")
  System.out.print(is)

  def randomize(is: IncidenceStructure) {
    forIndex(0, is.k) { i =>
      forIndex(0, is.b) { col =>
        var row = rnd.nextInt(is.v)
        while(is.active(row, col))
          row = rnd.nextInt(is.v)
        is.setIncidence(row, col, true)
      }
    }
    calculateChangeFactors
  }

  def getBibd: IncidenceStructure = {
    var unchanged = 0
    var lastChange = System.currentTimeMillis()

    var activeRow = 0
    var dormantRow = 0
    var cfa = 0
    var cfd = 0
    var rcf = 0
    var doWeChange = false
    var stale = false

    var col = -1
    while(true)
    {
      col = (col + 1) % is.b

      activeRow = randomActiveInCol(is, col)
      dormantRow = randomDormantInCol(is, col)

      if(is.isBIBD) {
        return this.is; // Sparta
      }

      this.iterations += 1

      cfa = changeFactor(activeRow)(col)
      cfd = changeFactor(dormantRow)(col)
      rcf = rnd.nextInt(math max(1,maxChangeFactor))

      stale = unchanged > is.v * is.b * is.r * is.k

      doWeChange = (rcf < cfa && rcf < cfd) ||
        stale ||
        maxChangeFactor <= 0

      if (doWeChange) {
        is.flip(activeRow, col)
        is.flip(dormantRow, col)
        calculateChangeFactors
        //System.out.println("Unchanged: " + unchanged + ", maxChangeFactor: " + maxChangeFactor)
        unchanged = 0
        lastChange = System.currentTimeMillis()
        generations += 1
        if(animate) {
          sleep(50)
          System.out.print(CLEAR_SCREEN)
          System.out.println(is)
        }
      }
      else {
        unchanged += 1
      }
    }
    return null // Never reached
  }

  def calculateChangeFactors: Unit = {
    maxChangeFactor = 0
    forIndex(0, is.v) { row =>
      forIndex(0, is.b) { col =>
        changeFactor(row)(col) = calculateChangeFactor(row, col)
        maxChangeFactor = math.max(changeFactor(row)(col), maxChangeFactor)
        minChangeFactor = math.min(changeFactor(row)(col), minChangeFactor)
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

  def sleep(millis: Long): Unit = {
    Thread.sleep(millis)
  }
}
