package org.krvoje.bibd

import scala.util.Random

case class MatrixIncidenceStructure(val v: Int, val k: Int, val lambda: Int, val randomizeMe: Boolean = true) extends IncidenceStructure{

  val (r: Int, b: Int) = computeParams(v, k, lambda)

  private val _incidences = Array.ofDim[Int](v, b)
  private val _rowIntersection = Array.ofDim[Int](v, v)
  private val _colIntersection = Array.ofDim[Int](b, b)
  private val _sumInRow = Array.ofDim[Int](v)
  private val _sumInCol = Array.ofDim[Int](b)
  private var _sumTotal = 0
  private val _sumIdeal = v*r
  private var _heuristicDistance = 0
  val maxHeuristicDistance = b * (v - k) + (v * v - v) * (r - lambda)

  updateCache()
  if(randomizeMe) randomize

  override def rowIntersection(row1: Int, row2: Int): Int = _rowIntersection(row1)(row2)

  override def colIntersection(col1: Int, col2: Int): Int = _colIntersection(col1)(col2)

  override def sumInRow(row: Int): Int = _sumInRow(row)

  override def sumInCol(col: Int): Int = _sumInCol(col)

  override def dormant(row: Int, col: Int): Boolean = incidences(row, col) == 0

  override def active(row: Int, col: Int): Boolean = incidences(row, col) != 0

  override def setIncidence(row: Int, col: Int, alive: Boolean): Unit = {
    if(!(alive && active(row, col))) flip(row, col)
  }

  override def incidences(row: Int, col: Int): Int = _incidences(row)(col)

  override def isBIBD: Boolean = {
    if(_sumTotal != _sumIdeal) {
      return false
    }

    forIndex(0, v) { row =>
      if(sumInRow(row) != r) return false
    }

    forIndex(0, b) { col =>
      if(sumInCol(col) != k) return false
    }

    forIndex(0, v - 1) { row =>
      forIndex(row + 1, v) { otherRow =>
        if(_rowIntersection(row)(otherRow) != lambda) return false
      }
    }

    return true
  }

  override def flip(row: Int, col: Int): Unit = {
    val newVal = if(active(row, col)) 0 else 1
    _incidences(row)(col) = newVal
    val increment = if(newVal == 0) -1 else 1

    _sumTotal += increment
    _sumInRow(row) += increment
    _sumInCol(col) += increment
    _rowIntersection(row)(row) += increment
    _colIntersection(col)(col) += increment

    forIndex(0, v) { otherRow =>
      if (_incidences(otherRow)(col) == 1) {
        _rowIntersection(otherRow)(row) += increment
        _rowIntersection(row)(otherRow) += increment
      }
    }

    forIndex(0, b) {otherCol =>
      if (_incidences(row)(otherCol) == 1) {
        _colIntersection(otherCol)(col) += increment
        _colIntersection(col)(otherCol) += increment
      }
    }
    calculateHeuristicDistance
  }

  override def heuristicDistance: Int = _heuristicDistance


  private def computeParams(v: Int, k: Int, lambda: Int) = {
    val rDouble = lambda.toDouble * (v-1) / (k-1)
    val bDouble = rDouble * v / k

    if(!rDouble.isValidInt || !bDouble.isValidInt)
      throw new RuntimeException(s"Invalid BIBD parameters. ($v, $k, $lambda)")

    (rDouble.toInt, bDouble.toInt)
  }

  private def updateCache(): Unit = {
    _sumTotal = 0

    forIndex(0, v) { row =>
      _sumInRow(row) = 0
      forIndex(0, b) { col =>
        _sumInRow(row) += _incidences(row)(col)
      }
      _sumTotal += _sumInRow(row)
    }

    forIndex(0, b) { col =>
      _sumInCol(col) = 0
      forIndex(0, v) { row =>
        _sumInCol(col) += _incidences(row)(col)
      }
    }

    forIndex(0, v) { row1 =>
      forIndex(0, v) { row2 =>
        _rowIntersection(row1)(row2) = 0
        forIndex(0, b) { col =>
          _rowIntersection(row1)(row2) += _incidences(row1)(col) * _incidences(row2)(col)
        }
      }
    }

    calculateHeuristicDistance()
  }

  private def calculateHeuristicDistance(): Unit = {
    _heuristicDistance = 0
    forIndex(0, v) { row =>
      _heuristicDistance += math.abs(sumInRow(row) - r)
    }
    forIndex(0, b) { col =>
      _heuristicDistance += math.abs(sumInCol(col) - k)
    }
    forIndex(0, v) { row1 =>
      forIndex(0, v) { row2 =>
        if(row1 != row2) {
          _heuristicDistance += math.abs(_rowIntersection(row1)(row2) - lambda)
        }
      }
    }
  }

  def randomize = {
    forIndex(0, k) { i =>
      forIndex(0, b) { col =>
        var row = Random.nextInt(v)
        while(active(row, col))
          row = Random.nextInt(v)
        setIncidence(row, col, true)
      }
    }
  }

  override def toString: String = {
    val sb = StringBuilder.newBuilder
    forIndex(0, v) {
      row => sb.append("\n")
        forIndex(0, b) {
          col => sb.append(_incidences(row)(col))
        }
    }
    sb.append("\n")
    sb.toString()
  }
}
